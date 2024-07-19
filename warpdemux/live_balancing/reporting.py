"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import os
import time
import traceback
from typing import Union

import numpy as np
import pandas as pd

from warpdemux.live_balancing.balancer import BarcodeBalancer, BarcodeBalancers
from warpdemux.live_balancing.utils import Result

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


class ProcessedCounters:
    def __init__(self, num_bcs: int):
        self.AF = 0
        self.AU = 0
        self.AC = 0
        self.AN = 0
        self.RF = 0
        self.RU = 0
        self.RC = 0
        self.RN = 0
        self.AC_per_bc = np.zeros(num_bcs)
        self.RC_per_bc = np.zeros(num_bcs)
        self.total = 0

    def _update_one(self, result: Result):
        self.total += 1
        if result.decision == "reject":
            if result.barcode_label == "unclassified":
                self.RU += 1
            elif result.barcode_label == "outlier":
                self.RN += 1
            elif result.barcode_label == "failed_read":
                self.RF += 1
            else:
                self.RC += 1
                self.RC_per_bc[int(result.barcode_label)] += 1
        else:
            if result.barcode_label == "failed_read":
                self.AF += 1
            elif result.barcode_label == "unclassified":
                self.AU += 1
            elif result.barcode_label == "outlier":
                self.AN += 1
            else:
                self.AC += 1
                self.AC_per_bc[int(result.barcode_label)] += 1

    def _update_many(self, counter_obj: "ProcessedCounters"):
        self.total += counter_obj.total
        self.AC += counter_obj.AC
        self.AF += counter_obj.AF
        self.AU += counter_obj.AU
        self.AN += counter_obj.AN
        self.RC += counter_obj.RC
        self.RF += counter_obj.RF
        self.RU += counter_obj.RU
        self.RN += counter_obj.RN
        for i in range(self.AC_per_bc.size):
            self.AC_per_bc[i] += counter_obj.AC_per_bc[i]
        for i in range(self.RC_per_bc.size):
            self.RC_per_bc[i] += counter_obj.RC_per_bc[i]

    def update(self, result: Union[Result, "ProcessedCounters"]):
        if isinstance(result, Result):
            self._update_one(result)
        elif isinstance(result, ProcessedCounters):
            self._update_many(result)


class ProcessedCountersContainer:
    def __init__(self, balancers: BarcodeBalancers):
        self._balancers = balancers
        num_bcs = balancers.config_list[0].num_bcs
        self._counters = {}
        for name in balancers.names:
            self._counters[name] = ProcessedCounters(num_bcs)

    def get_counters(self, name):
        return self._counters[name]

    def get_counters_dict(self):
        return self._counters

    def get_counters_df(self):
        return pd.DataFrame([c.__dict__ for c in self._counters.values()])

    def update(self, result):
        if isinstance(result, Result):
            # get name of balancer based on channel
            self._counters[result.balancer_name].update(result)
        elif isinstance(result, ProcessedCountersContainer):
            for name, counter in result.get_counters_dict().items():
                self._counters[name].update(counter)

    @property
    def names(self):
        return list(self._counters.keys())


def report_worker(input_queue, config, skip_stats, balancers: BarcodeBalancers):
    fpath = os.path.join(config.save_path, f"barcode_balancing_{config.runid}.csv")
    print_width = 9

    try:
        stats = ProcessedCountersContainer(balancers)
        stop_condition = False

    except:
        print("Exception in report worker")
        traceback.print_exc()
        return

    def _process_queue():
        res = []
        _stats = ProcessedCountersContainer(balancers)
        stop_condition = False

        if input_queue.empty():
            return None, None, False

        while True:
            try:
                item = input_queue.get(block=False)  # type: Result
            except:
                break  # queue is now empty

            if item is None:  # process is being closed
                stop_condition = True
                break

            res.append(item.to_dict())

            if item.balance_type == "failed_read":
                continue  # don't count failed reads

            _stats.update(item)

        if not len(res):  # can happen if the only item in the queue was None
            return None, None, stop_condition

        return pd.DataFrame(res), _stats, stop_condition

    def report_skip_stat_balance():
        # TODO: split in positive and negative missed_obs vals
        m = 0
        s = 0
        if len(skip_stats["missed_obs"]) > 2:
            m = np.mean(skip_stats["missed_obs"])
            s = np.std(skip_stats["missed_obs"])

        print(f"Missed observations (mean, std): {m:.1f}, {s:.1f}")
        print(f"Skipped reads due to missed start: {skip_stats['missed_reads']}")
        print(
            f"Skipped reads due to not finding poly A: {skip_stats['too_long_reads']}"
        )
        print(
            "Skipped reads due to not matching expected signal:"
            f" {skip_stats['not_real_read']}"
        )
        print("")

    def report_current_balance(balancer: "BarcodeBalancer"):
        print("--------------------------------------")
        print(f"Balancer: {balancer.name} ({balancer.balance_type})")
        print("--------------------------------------")

        t = time.strftime("%H:%M:%S")

        # main count dict contains the decisions we made (per adapter)
        counters = stats.get_counters(balancer.name)

        # TODO: use 1-indexed true barcode class, rather than index
        labels = " ".join(
            [f"{i:02d}".ljust(print_width) for i in range(counters.AC_per_bc.size)]
        )

        print("")
        print(f"[{config.runid}] --- {t} ")
        print("")

        print("Processed reads")
        print(f"Accepted Classified: {counters.AC} / {counters.total}")
        print(f"Accepted Unclassified: {counters.AU} / {counters.total}")
        print(f"Accepted Failed: {counters.AF} / {counters.total}")
        print(f"Accepted Noise: {counters.AN} / {counters.total}")
        print(f"Rejected Classified: {counters.RC} / {counters.total}")
        print(f"Rejected Unclassified: {counters.RU} / {counters.total}")
        print(f"Rejected Failed: {counters.RF} / {counters.total}")
        print(f"Rejected Noise: {counters.RN} / {counters.total}")
        print("")

        counts = " ".join([f"{int(c)}".ljust(print_width) for c in counters.AC_per_bc])
        frac = " ".join(
            [
                f"{c/counters.AC if counters.AC else 0:.2f}".ljust(print_width)
                for c in counters.AC_per_bc
            ]
        )

        print("Processed reads; #Accepted per BC")
        print("Fraction represents accepted processed reads per bc / processed reads")
        print(f"{labels}")
        print(f"{counts}")
        print(f"{frac}")
        print("")

        counts = " ".join([f"{int(c)}".ljust(print_width) for c in counters.RC_per_bc])
        frac = " ".join(
            [
                f"{c/counters.RC if counters.RC else 0:.2f}".ljust(print_width)
                for c in counters.RC_per_bc
            ]
        )

        print("Processed reads; #Rejected per BC")
        print("Fraction represents rejected processed reads per bc / processed reads")

        print(f"{labels}")
        print(f"{counts}")
        print(f"{frac}")
        print("")

        if balancer.balance_type == "read_count":
            matched_read_counts = balancer.pod5_read_sets.matched_read_count_per_bc
            all_matched = np.sum(matched_read_counts)

            counts = " ".join(
                [f"{int(c)}".ljust(print_width) for c in matched_read_counts]
            )
            frac = " ".join(
                [
                    f"{c/all_matched if all_matched else 0:.2f}".ljust(print_width)
                    for c in matched_read_counts
                ]
            )

            print(f"Detected pod5 reads per BC")
            print(
                "Fraction represents classified pod5 reads per bc / classified reads in"
                " pod5"
            )

            print(f"{labels}")
            print(f"{counts}")
            print(f"{frac}")
            print("")

            unmatched_classif_counts = (
                balancer.pod5_read_sets.unmatched_seen_read_count_per_bc
            )

            counts = " ".join(
                [f"{int(c)}".ljust(print_width) for c in unmatched_classif_counts]
            )
            frac = " ".join(
                [
                    f"{unmatched/(matched+unmatched) if (matched+unmatched) else 0:.2f}".ljust(
                        print_width
                    )
                    for matched, unmatched in zip(
                        matched_read_counts, unmatched_classif_counts
                    )
                ]
            )

            print(f"Detected adapters without pod5 per BC")
            print(
                "Fraction represents classified reads without pod5 entry per bc / all"
                " classified reads"
            )

            print(f"{labels}")
            print(f"{counts}")
            print(f"{frac}")
            print("")

            print(
                "Estimated average read length (kbases): "
                f"{balancer.pod5_read_sets.matched_estimated_read_length:.2f}",
            )
            print("Estimated average read length per BC (kbases)")
            print(f"{labels}")
            print(
                " ".join(
                    [
                        f"{c:.2f}".ljust(print_width)
                        for c in balancer.pod5_read_sets.matched_estimated_read_length_per_bc
                    ]
                )
            )
            print("")

        elif balancer.balance_type == "base_normalization":
            matched_base_counts = balancer.pod5_read_sets.matched_kbases_per_bc
            all_matched = np.sum(matched_base_counts)

            # TODO: use 1-indexed true barcode class, rather than index
            labels = " ".join(
                [f"{i:02d}".ljust(print_width) for i in range(matched_base_counts.size)]
            )

            counts = " ".join(
                [f"{int(c)}".ljust(print_width) for c in matched_base_counts]
            )
            frac = " ".join(
                [
                    f"{c/all_matched if all_matched else 0:.2f}".ljust(print_width)
                    for c in matched_base_counts
                ]
            )

            print(f"Detected pod5 kbases per BC")
            print(
                "Fraction represents classified pod5 kbases per bc / classified kbases"
                " in pod5"
            )

            print(f"{labels}")
            print(f"{counts}")
            print(f"{frac}")
            print("")

            print(
                "Estimated average read length (kbases): "
                f"{balancer.pod5_read_sets.matched_estimated_read_length:.2f}",
            )
            print("Estimated average read length per BC (kbases)")
            print(f"{labels}")
            print(
                " ".join(
                    [
                        f"{c:.2f}".ljust(print_width)
                        for c in balancer.pod5_read_sets.matched_estimated_read_length_per_bc
                    ]
                )
            )
            print("")

        else:
            pass

        if balancer.blacklist.sum():
            # TODO: use 1-indexed true barcode class, rather than index
            indices = np.argwhere(balancer.blacklist).ravel()
            bl = [f"bc{str(i).zfill(2)}" for i in indices]
            print(f"Blacklist: [{', '.join(bl)}]")
            print("")

        if balancer.ignorelist.sum():
            # TODO: use 1-indexed true barcode class, rather than index
            indices = np.argwhere(balancer.ignorelist).ravel()
            il = [f"bc{str(i).zfill(2)}" for i in indices]
            print(f"Ignorelist: [{', '.join(il)}]")
            print("")

    def current_balance_stats(balancer: "BarcodeBalancer"):
        counters = stats.get_counters(balancer.name)

        balancer_stats = dict(
            time=time.strftime("%H:%M:%S"),
            balancer_name=balancer.name,
            balance_type=balancer.balance_type,
            total_processed=counters.total,
            AC=counters.AC,
            AU=counters.AU,
            AF=counters.AF,
            AN=counters.AN,
            RC=counters.RC,
            RU=counters.RU,
            RF=counters.RF,
            RN=counters.RN,
            **{f"AC{i:02d}": int(c) for i, c in enumerate(counters.AC_per_bc)},
            **{
                f"AC{i:02d}_frac": c / counters.AC if counters.AC else 0
                for i, c in enumerate(counters.AC_per_bc)
            },
            **{f"RC{i:02d}": int(c) for i, c in enumerate(counters.RC_per_bc)},
            **{
                f"RC{i:02d}_frac": c / counters.RC if counters.RC else 0
                for i, c in enumerate(counters.RC_per_bc)
            },
        )

        if balancer.balance_type == "read_count":
            matched_read_counts = balancer.pod5_read_sets.matched_read_count_per_bc
            all_matched = np.sum(matched_read_counts)

            unmatched_classif_counts = (
                balancer.pod5_read_sets.unmatched_seen_read_count_per_bc
            )

            balancer_stats.update(
                dict(
                    **{
                        f"pod5{i:02d}": int(c)
                        for i, c in enumerate(matched_read_counts)
                    },
                    **{
                        f"pod5{i:02d}_frac": c / all_matched if all_matched else 0
                        for i, c in enumerate(matched_read_counts)
                    },
                    **{
                        f"no_pod5{i:02d}": int(c)
                        for i, c in enumerate(unmatched_classif_counts)
                    },
                    **{
                        f"no_pod5{i:02d}_frac": c1 / (c1 + c2) if (c1 + c2) else 0
                        for i, (c1, c2) in enumerate(
                            zip(matched_read_counts, unmatched_classif_counts)
                        )
                    },
                    est_avg_kb=balancer.pod5_read_sets.matched_estimated_read_length,
                    **{
                        f"est_avg_kb{i:02d}": v
                        for i, v in enumerate(
                            balancer.pod5_read_sets.matched_estimated_read_length_per_bc
                        )
                    },
                )
            )
        elif balancer.balance_type == "base_normalization":
            matched_base_counts = balancer.pod5_read_sets.matched_kbases_per_bc
            all_matched = np.sum(matched_base_counts)

            balancer_stats.update(
                dict(
                    **{f"pod5_kb{i:02d}": v for i, v in enumerate(matched_base_counts)},
                    **{
                        f"pod5_kb{i:02d}_frac": c / all_matched if all_matched else 0
                        for i, c in enumerate(matched_base_counts)
                    },
                    est_avg_kb=balancer.pod5_read_sets.matched_estimated_read_length,
                    **{
                        f"est_avg_kb{i:02d}": v
                        for i, v in enumerate(
                            balancer.pod5_read_sets.matched_estimated_read_length_per_bc
                        )
                    },
                )
            )

        if balancer.blacklist.sum():
            balancer_stats.update(
                **{
                    f"blacklist{i:02d}": int(v)
                    for i, v in enumerate(balancer.blacklist)
                }
            )
        else:
            balancer_stats.update(
                {f"blacklist{i:02d}": 0 for i in range(balancer.num_bcs)}
            )

        if balancer.ignorelist.sum():
            balancer_stats.update(
                **{
                    f"ignorelist{i:02d}": int(v)
                    for i, v in enumerate(balancer.ignorelist)
                }
            )
        else:
            balancer_stats.update(
                {f"ignorelist{i:02d}": 0 for i in range(balancer.num_bcs)}
            )
        return balancer_stats

    def save_current_balance(balancer: "BarcodeBalancer"):
        _fpath = os.path.join(
            config.save_path,
            f"balancing_stats__{balancer.name}__{config.runid}.csv",
        )

        pd.DataFrame().from_dict(
            current_balance_stats(balancer), orient="index"
        ).T.to_csv(
            _fpath,
            mode="a",
            index=False,
            header=not os.path.exists(_fpath),
        )

    def save_classified_read_ids_without_pod5(balancer: "BarcodeBalancer"):
        if not balancer.is_watcher:
            return

        _fpath = os.path.join(
            config.save_path,
            f"classified_adapters_without_pod5_entry__{balancer.name}__{config.runid}.csv",
        )

        balancer.pod5_read_sets.unmatched_seen_df.to_csv(_fpath, index=False)

    def report_time_stats():
        df = pd.read_csv(
            os.path.join(config.save_path, f"barcode_balancing_{config.runid}.csv"),
            index_col=False,
        )

        print("")
        m = df.time_polya_detect.mean()
        s = df.time_polya_detect.std()
        print(f"PolyA detection time: {m:.4f} +/- {s:.4f}")

        m = df.time_segmentation.mean()
        s = df.time_segmentation.std()
        print(f"Segmentation time: {m:.4f} +/- {s:.4f}")

        m = df.time_classification.mean()
        s = df.time_classification.std()
        print(f"Classification time: {m:.4f} +/- {s:.4f}")

        m = df.time_balancing.mean()
        s = df.time_balancing.std()
        print(f"Balancing time: {m:.4f} +/- {s:.4f}")

        m = df.total_compute_time.mean()
        s = df.total_compute_time.std()
        print(f"Total compute time: {m:.4f} +/- {s:.4f}")

        m = df.total_time.mean()
        s = df.total_time.std()
        print(f"Total time: {m:.4f} +/- {s:.4f}")
        print("")

    while not stop_condition:
        time.sleep(config.save_every_sec)

        try:
            df, _stats, stop = _process_queue()

            if df is not None:
                # TODO: split results csv into multiple files as size increases?
                df.round(6).to_csv(
                    fpath,
                    mode="a",
                    index=False,
                    header=not os.path.exists(fpath),
                )

                stats.update(_stats)
                for name in stats.names:
                    balancer = balancers.get_balancer(name)
                    report_current_balance(balancer)
                    save_current_balance(balancer)
                    save_classified_read_ids_without_pod5(balancer)

                print("--------------------------------------")
                print("")

                report_skip_stat_balance()

        except:
            print("Exception in report worker")
            traceback.print_exc()
            return

        stop_condition = stop

    print("--------------------------------------")
    print("Closing watcher, final stats.")
    print("--------------------------------------")

    report_time_stats()
