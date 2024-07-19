"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import time
import traceback
from multiprocessing import Queue

from warpdemux.read_until import ReadUntilClient

from warpdemux.file_proc import load_model
from warpdemux.models.dtw_svm import confidence_margin
from warpdemux.sig_proc import extract_adapter, mad_winsor, normalize, segment_signal
from warpdemux.live_balancing.balancer import BarcodeBalancers
from warpdemux.live_balancing.config_parser import BalancingConfig, ModelConfig
from warpdemux.live_balancing.utils import ReadObject, Result

# TODO make worker functions into classes


def segmentation_worker(
    input_queue: Queue,
    output_queue: Queue,
    config: ModelConfig,
) -> None:
    while True:
        read_obj = input_queue.get()
        start_time = time.time()

        if read_obj is None:  # Use None as a signal to stop.
            break

        try:
            adapter_sig = extract_adapter(
                read_obj.data_arr,
                adapter_start=0,
                adapter_end=read_obj.polya_start,
                extract_padding=config.spc.sig_extract.padding,
            )
            # required
            norm_adapter_sig = mad_winsor(
                adapter_sig,
                outlier_thresh=config.spc.sig_norm.outlier_thresh,
                window_size=config.spc.sig_norm.winsor_window,
            )

            # optional
            norm_adapter_sig = normalize(
                norm_adapter_sig,
                config.spc.sig_extract.normalization,
            )

            segment_avgs, _ = segment_signal(
                norm_adapter_sig,
                num_events=config.spc.segmentation.num_events,
                min_obs_per_base=config.spc.segmentation.min_obs_per_base,
                running_stat_width=config.spc.segmentation.running_stat_width,
            )
            segment_avgs = segment_avgs.reshape(1, -1)

            if not segment_avgs.size:
                # TODO
                print("no segments")

                continue

            norm_segment_avgs = normalize(
                segment_avgs, method=config.spc.segmentation.normalization
            )

            read_obj.data_arr = norm_segment_avgs[
                :, -config.spc.segmentation.barcode_num_events :
            ]  # barcode_fpt

            read_obj.time_per_step.append(time.time() - start_time)

            output_queue.put(read_obj)
        except:
            print("Exception in segmentation worker")
            traceback.print_exc()
            return


def classification_worker(
    input_queue: Queue, output_queue: Queue, config: ModelConfig
) -> None:
    try:
        model = load_model(config.model_name)
    except:
        print("Exception in loading model")
        traceback.print_exc()
        return

    while True:
        try:
            read_obj = input_queue.get()  # type:ReadObject
            if read_obj is None:  # Use None as a signal to stop.
                break

            start_time = time.time()

            y_pred, y_prob = model.predict(
                read_obj.data_arr,
                nproc=1,
            )

            read_obj.data_arr = y_prob.reshape(1, -1)
            read_obj.time_per_step.append(time.time() - start_time)
            read_obj.is_outlier = y_pred == -1

            output_queue.put(read_obj)

        except:
            print("Exception in classification worker")
            traceback.print_exc()
            return


def balance_worker(
    input_queue: Queue,
    output_queue: Queue,
    config: BalancingConfig,
    client: ReadUntilClient,
    balancers: BarcodeBalancers,
) -> None:
    while True:
        try:
            read_obj = input_queue.get()  # type:ReadObject
            start_time = time.time()

            if read_obj is None:  # Use None as a signal to stop.
                break

            label = None
            decision = None

            balancer = balancers.get_balancer(read_obj.channel)

            balancer_name = balancer.name

            if read_obj.is_outlier:
                label = "outlier"
                if balancer.balance_type == "reject_all":
                    decision = "reject"
                else:
                    decision = "retain"
                reason = "outlier_detected"

            elif confidence_margin(read_obj.data_arr) < config.pred_conf_threshold:
                # not sent to balancer, so won't show up in per barcode reporting
                if balancer.balance_type == "reject_all":
                    decision = "reject"
                else:
                    decision = "retain"
                label = "unclassified"
                reason = "below_conf_pred"

            else:  # can't be outlier and has enough confidence
                barcode_idx = read_obj.data_arr.argmax()
                label = str(barcode_idx)
                decision, reason = balancers.balance(
                    barcode_idx, read_obj.read_id, read_obj.channel
                )
                if decision:
                    decision = "retain"
                else:
                    decision = "reject"

            if decision == "reject":
                if (
                    read_obj.chunk_length - read_obj.polya_start
                ) > config.max_signal_after_polya:
                    # chunk has accumulated too much data after polyA
                    # don't reject it to avoid pores getting blocked
                    decision = "retain"
                    reason = "too_late_to_reject"
                    # when this happens a lot, the reporting balance will be quite off
                    # TODO: add a counter for this / report it

                else:
                    reject_duration = balancers.get_balancer(
                        read_obj.channel
                    ).reject_duration
                    if reject_duration is None:
                        reject_duration = config.reject_duration

                    client.unblock_read(
                        read_obj.channel,
                        read_obj.read_number,
                        duration=reject_duration,
                    )

            read_obj.time_per_step.append(time.time() - start_time)
            output_queue.put(
                Result(
                    read_obj,
                    label,
                    decision,
                    reason,
                    balancer.name,
                    balancer.balance_type,
                )
            )

        except:
            print("Exception in balancing worker")
            traceback.print_exc()
            return
