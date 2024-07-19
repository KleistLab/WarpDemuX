"""
WarpDemuX

Copyright (c) 2023 by Wiep K. van der Toorn
Contact: w.vandertoorn@fu-berlin.de

"""

import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from warpdemux.live_balancing.session import Session

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config_file",
    type=str,
    default=None,
)


class DummyRead:
    def __init__(self, raw_data, number, id):
        self.raw_data = raw_data
        self.number = number
        self.id = id

        self.start_sample = 0
        self.chunk_start_sample = 0
        self.chunk_length = 6000


class DummyClient:
    def __init__(
        self,
        max_c: int = 800,
    ):
        testdir = os.path.join(os.path.dirname(__file__), "test")
        self.data = np.load(os.path.join(testdir, "test_signals.npy"))
        self.c = 0
        self.max_c = max_c  # 800

        self.is_running = True
        self.queue_length = self.max_c
        self.signal_dtype = np.float32

        read_ids = np.concatenate(
            [
                np.loadtxt(os.path.join(testdir, f"read_ids{i}.txt"), dtype=str)
                for i in range(3)
            ]
        )
        read_ids = read_ids[: int(0.5 * self.max_c)]
        read_ids = np.insert(
            read_ids,
            np.random.randint(0, read_ids.size, self.max_c - read_ids.size),
            [f"{i:03d}" for i in np.arange(self.max_c - read_ids.size)],
        )

        self.data = [
            (np.random.randint(1, 127), DummyRead(self.data[i][:6000], i, read_ids[i]))
            for i in range(self.max_c)
        ]

    def get_read_chunks(
        self,
        batch_size=1,
        last=False,
        min_chunk_length=3000,
        *args,
        **kwargs,
    ):
        for i in range(self.c, min(self.c + batch_size, self.max_c)):
            chunk = self.data[i]
            time.sleep(0.01)
            self.queue_length -= 1
            yield chunk
        self.c += batch_size

        if self.c >= self.max_c:
            self.is_running = False

    def stop_receiving_read(*args, **kwargs):
        pass

    def unblock_read(*args, **kwargs):
        pass

    def run(self):
        pass

    def reset(self):
        pass


class DummySession(Session):
    def _check_client(self):
        return True


def debug_test():
    config_path = os.path.join(
        os.path.dirname(__file__),
        "test",
        "config_only_read_count.toml",
    )

    read_until_client = DummyClient(800)
    session = DummySession(read_until_client, config_path)
    with ThreadPoolExecutor(1) as executor:
        future = executor.submit(session.run, batch_size=1)
    future.result()  # waits for the run to finish

    config_path = os.path.join(
        os.path.dirname(__file__),
        "test",
        "config_only_adapter_count.toml",
    )

    read_until_client = DummyClient(800)
    session = DummySession(read_until_client, config_path)
    with ThreadPoolExecutor(1) as executor:
        future = executor.submit(session.run, batch_size=1)
    future.result()  # waits for the run to finish

    config_path = os.path.join(
        os.path.dirname(__file__),
        "test",
        "config_only_base_normalization.toml",
    )

    read_until_client = DummyClient(800)
    session = DummySession(read_until_client, config_path)
    with ThreadPoolExecutor(1) as executor:
        future = executor.submit(session.run, batch_size=1)
    future.result()  # waits for the run to finish

    config_path = os.path.join(
        os.path.dirname(__file__),
        "test",
        "config_only_none.toml",
    )

    read_until_client = DummyClient(800)
    session = DummySession(read_until_client, config_path)
    with ThreadPoolExecutor(1) as executor:
        future = executor.submit(session.run, batch_size=1)
    future.result()  # waits for the run to finish

    config_path = os.path.join(
        os.path.dirname(__file__),
        "test",
        "config_only_reject_all.toml",
    )

    read_until_client = DummyClient(800)
    session = DummySession(read_until_client, config_path)
    with ThreadPoolExecutor(1) as executor:
        future = executor.submit(session.run, batch_size=1)
    future.result()  # waits for the run to finish

    config_path = os.path.join(
        os.path.dirname(__file__),
        "test",
        "config.toml",
    )

    read_until_client = DummyClient(800)
    session = DummySession(read_until_client, config_path)
    with ThreadPoolExecutor(1) as executor:
        future = executor.submit(session.run, batch_size=1)
    future.result()  # waits for the run to finish


if __name__ == "__main__":
    args = parser.parse_args()

    if args.config_file is None:
        debug_test()

    else:
        read_until_client = DummyClient(800)
        session = DummySession(read_until_client, args.config_file)
        with ThreadPoolExecutor(1) as executor:
            future = executor.submit(session.run, batch_size=1)
        future.result()
