"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import argparse
import time
from concurrent.futures import ThreadPoolExecutor

from warpdemux import read_until
from warpdemux.live_balancing.session import Session


def setup_read_until_client(mk_port: int) -> read_until.ReadUntilClient:
    return read_until.ReadUntilClient(
        mk_port=mk_port,
        one_chunk=False,
        cache_type=read_until.AccumulatingCache(size=5120),
        calibrated_signal=True,
        # critical: the following two settings only work with wdx branch of read_until. Otherwise, disable filter_strands
        filter_strands=True,
        prefilter_classes={"adapter"},
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mk_port",
        type=int,
        help="port to connect to MinKNOW on",
    )
    parser.add_argument(
        "config_path",
        type=str,
        help="path to the config file",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    read_until_client = setup_read_until_client(args.mk_port)

    analysis = Session(
        read_until_client,
        args.config_path,
    )

    with ThreadPoolExecutor(1) as executor:
        executor.submit(analysis.run, batch_size=args.batch_size)

    time.sleep(args.time_to_run)
    analysis.stop()


if __name__ == "__main__":
    main()
