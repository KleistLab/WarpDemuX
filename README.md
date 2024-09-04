# WarpDemuX

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![DOI:10.1101/2024.07.22.604276](http://img.shields.io/badge/DOI-10.1101/2024.07.22.604276-blue.svg)](https://doi.org/10.1101/2024.07.22.604276)
[![DOI:10.21203/rs.3.rs-4783223/v1](http://img.shields.io/badge/10.21203/rs.3.rs-4783223/v1-blue.svg)](https://doi.org/10.21203/rs.3.rs-4783223/v1)

We introduce WarpDemuX, an ultra-fast and highly accurate adapter-barcoding and demultiplexing approach. WarpDemuX enhances speed and accuracy by fast processing of the raw nanopore signal, the use of a light-weight machine-learning algorithm and the design of optimized barcode sets. Additionally, integrating WarpDemuX into sequencing control software enables real-time enrichment of target molecules through barcode-specific adaptive sampling.

## Installation

```{bash}
# Clone this repository with the adapt submodule
git clone --recursive https://github.com/KleistLab/WarpDemuX.git [path/to/store/WarpDemuX]

# Create a new conda environment using the environment.yml file
mamba env create -n WDX -f [path/to/store/WarpDemuX]/environment.yml
conda activate WDX
```

WarpDemuX depends on ADAPTed, our tool for adapter and poly(A) tail detection. ADAPTed is included as a submodule for now. To make sure WDX can access this module correctly, you need to install WDX in editable mode. This can be done by using pip's '-e' flag:

```
pip install -e [path/to/store/WarpDemuX]

# or, for barcode-specific adaptive sampling
cd [path/to/store/WarpDemuX]
pip install -e '.[live-demux]'
```

Alternatively, if you don't wish to install in editable mode, you can install the warpdemux and adapted packages separately:

```
cd [path/to/store/WarpDemuX]
pip install .

cd warpdemux/adapted
pip install .
```

## Usage

```{bash}
conda activate WDX
warpdemux demux --input INPUT [INPUT ...] --output OUTPUT --model_name MODEL_NAME
```

Where INPUT are file(s) or directory(s) to be demultiplexed (pod5 files) and OUTPUT is the path to where the run output folder should be created. MODEL_NAME is the name of the model to use for demultiplexing. See the [Models](#models) section for a list of available models.

For further instructions and options, run `warpdemux --help` and `warpdemux demux --help`.

### File formats

WarpDemuX only supports the pod5 file format. Please convert your fast5/slow5/blow5 files to pod5 format before running WarpDemuX.
  
## Models

### Available models

- WDX-DPC_rna002_v0_4_0
- WDX4_rna002_v0_4_0
- WDX6_rna002_v0_4_0
- WDX8_rna002_v0_4_0
- WDX10_rna002_v0_4_0
- WDX12_rna002_v0_4_2

### WDX-DPC_rna002_v0_4_0

- For 4 samples, using the original DeePlexiCon barcodes.

### WDX4_rna002_v0_4_0

- For 4 samples. Uses barcodes WDX_bc04, WDX_bc05, WDX_bc06, and WDX_bc08.

### WDX4_rna002_v0_4_0

- For 6 samples. Uses barcodes WDX_bc01, WDX_bc03, WDX_bc05, WDX_bc06, WDX_bc07, and WDX_bc11.

### WDX8_rna002_v0_4_0

- For 8 samples. Uses barcodes WDX_bc01, WDX_bc03, WDX_bc05, WDX_bc06, WDX_bc07, WDX_bc09, WDX_bc11 and WDX_bc12.

### WDX10_rna002_v0_4_0

- For 10 samples. Uses barcodes WDX_bc01, WDX_bc02, WDX_bc03, WDX_bc05, WDX_bc06, WDX_bc07, WDX_bc09, WDX_bc10, WDX_bc11, and WDX_bc12.

### WDX12_rna002_v0_4_2

- For 12 samples. Uses all 12 WDX barcodes.

## Tips

If the length of the polyA tail is important to you, keep in mind that you need to preload a sufficient amount of the signal into memory. This will slow down the adapter detection code. Alternatively, you can re-analyse reads for which `polya_truncated` is True with higher `max_obs_adapter` settings.

## Live balancing

### Adapter classification

WarpDemuX classifies raw adapter signals to adapter labels and assigns a meta-label to each signal. The meta-label can be one of the following:

- 'classified': the signal is confidently classified as one of the training adapters
- 'unclassified': the signal is not confidently classified as one of the training adapters
- 'outlier': the signal is considered an outlier with respect to the training signal distribution
- 'failed': the signal is does not satisfy the requirements to be considered an adapter signal

Based on the classification and the active balancing strategy, WarpDemuX decides whether to accept or reject the read.

### Balancing strategies

WarpDemuX currently supports 5 balancing strategies:

- `none`: no balancing is done, all reads (classified, unclassified, outliers and failed) are accepted
- `reject_all`: no balancing is done, all reads (classified, unclassified, outliers and failed) are rejected
- `adapter_count`: barcodes are balanced based on the number of classified adapter signals (default). Outliers and failed reads are not part of the adapter balance and thus automatically accepted.
- `base_normalization`: barcodes are balanced based on the number of sequenced kbases. Outliers and failed reads are not part of the adapter balance and thus automatically accepted.
- `read_count`: barcodes are balanced based on the number of sequenced reads. Outliers and failed reads are not part of the adapter balance and thus automatically accepted.

In addition, WarpDemuX supports blacklisting barcodes. Barcodes on the blacklist are automatically rejected and not considered when computing the balance.

### config.toml

The live balancing is configured using a config file. Several example config files can be found in `warpdemux/live_balancing/test`.

#### `[model]` section

- `model_name` should be one of the available models, see [Available models](#available-models).

#### `[flowcell]` section

- `flowcell_type` should be one of `flongle`, `minion`.
- `min_channel` and `max_channel` specify the range of channels to use. Default is 1 and max channel on the flowcell.

#### `[processing]` section

- `nproc_segmentation` specifies the number of processes to use for segmentation. Default is 2.
- `nproc_classification` specifies the number of processes to use for classification. Default is 4.

#### `[acquisition]` section

- `max_missed_start_offset` specifies the maximum number of missed offsets to allow. Default is 1.
- `max_chunk_size` is the maximum number of observations to watch a read for and try to detect the poly(A) boundary.
    If the read has grown larger than `max_chunk_size`, no polyA detection is attempted and the read will not
    be further processed. Default value is 25000.
- `min_chunk_size` is the minimum number of observations to watch a read for before submitting it to the processing queue.
    Smaller read chunks are left to accumulate signal before attempting processing. Default value is 5000.
- `min_adapter_length` is the minimum length of the adapter. Used in the detection the adapter-poly(A) boundary.
    Defaults to `min_chunk_size`.

#### `[balancing]` section

- `pred_conf_threshold` determines the confidence threshold for the classifier. If the confidence score is below this threshold, the read is considered to be 'unclassified'. See the section on balancing strategies to see how unclassified reads are handled. Default is 0.2.
- `reject_duration` is the duration in seconds to reverse the channel charge when rejecting a read. Default is 0.1.

#### `[[balancers]]` sections

You can use multiple balancers tracking different channels to compare and execute different adaptive sampling strategies.
Each balancer should have its own `[[balancers]]` section in the config file. Note that when using only one balancer,
the section should still be named `[[balancers]]`. When using multiple balancers (= when you have multiple `[[balancers]]` sections),
`channel_frac` or `channel_num` should be specified for each balancer. These parameters are mutually exclusive.

Note that, after all `[[balancers]]` sections and respective `channel_frac` or `channel_num` are processed,
any remaining channels are assigned to a balancer with `balance_type = 'none'`.
This means that an extra balancer may get created, or that an existing balancer of balance type
`none` may get more channels than specified by `channel_frac` or `channel_num`.
This is done to make sure that all channels are assigned to a balancer.

- `channel_frac` is the fraction of the flowcell channels that is assigned to the balancer.
    The number of channels is calculated as `int(channel_frac * channels_on_flowcell)`.
    Pay attention that this can lead to the number of channels being rounded down,
    it may be more intuitive to use `channel_num` instead.
    The channels are per balancer will not be in a consecutive range, rather they are chosen randomly to limit the impact of local flowcell irregularities and air bubbles.

- `channel_num` is the number of flowcell channels that is assigned to the balancer. The channels are chosen randomly.

- `balance_type` can be one of the following:
  - `none`: no balancing is done
  - `reject_all`: all barcodes are rejected
  - `adapter_count`: barcodes are balanced based on the number of classified adapter sequences (default)
  - `base_normalization`: barcodes are balanced based on the number of sequenced kbases
  - `read_count`: barcodes are balanced based on the number of sequenced reads

- `min_stat` is the minimum number of 'unit' to see of an bc, before balancing is considered.
    The meaning of 'unit' depends on the balancing type. unit=counts for all balancing types,
    except for type 'base_normalization', for which it is kbases. Default is 100.

- `reject_duration` is the balancer specific duration in seconds to reverse the channel charge when rejecting a read.
    If not specified, this falls back to the balancing.reject_duration value, see `[balancing]` section.

- `pod5_watch_dir` is required for balance types`read_count` and `base_normalization` and specifies the pod5 output
    directory of the run (`pod5_watch_dir = "/path/to/pod5_dir"`). The pod5 files written to this directory during the run are
    monitored to calculate the number of sequenced reads/kbases per barcode. Path should be a relative path to root of the WDX repo,
    or an absolute path.

- `pod5_check_interval` specifies the interval at which the pod5 files are checked, e.g.`pod5_check_interval = 1` (default 1 second)

- `watch_barcodeXX` is a boolean flag that can be set to false (`watch_barcodeXX = false`) to place a barcode on the
    ignore list. Reads for barcodes on the ignore list are automatically accepted and not considered when computing the balance.
    Setting barcodes on the ignore list is usefull when certain barcodes were not used the experiment that you're sequencing.
    Default is true for all barcodes.
    Barcodes are 0-indexed and should be 02d formatted, i.e.: watch_barcode00 - watch_barcode11.

- `blacklist_barcodeXX` is a boolean flag that can be set to true (`blacklist_barcodeXX = true`) to place a barcode on the
    blacklist. Reads for barcodes on the blacklist are automatically rejected and not considered when computing the balance.
    Default is false for all barcodes.
    Barcodes are 0-indexed and should be 02d formatted, i.e.: blacklist_barcode00 - blacklist_barcode11.

- `max_barcodeXX` can be used to add an upper bound for a barcode (in 'unit' as explained above) using`max_barcodeXX = XYZ`
    Setting a max will put the barcode on the blacklist when the bound is reached.
    The default value is `inf` for all barcodes.
    Barcodes are 0-indexed and should be 02d formatted, i.e.: max_barcode00 - max_barcode11.

- `watch_for_missing` is a boolean that specifies whether to keep an eye out for experimental errors:
    if almost no reads (`< min_stat`) were seen for a barcode after `wait_to_see` seconds
    the barcode is added to the ignore list such that the balance does not stagnate.
    Default is true.

- `wait_to_see` specifies when (how many seconds into the run) to decide if any barcodes are missing (see `watch_for_missing`).
    Make sure that this value is larger than the time it takes to see the first `min_stat` reads/kbases for a barcode,
    else all barcodes may end up on the ignore list and that is not what you want.
    Given that the pore scans take 4 minutes, we advise to set this to at least 900 seconds (default) for `min_stat=100`.

#### `[reporting]` section

- `save_path` specifies the directory to save the balancing log and results to. Default is `./results`.
- `save_every_sec` specifies the interval in seconds at which the balancing log and results are saved. Default is 30.

#### Notes

- Note that for balance type `reject_all`  the blacklist and upper bounds are ignored (ignore list can be used)
- Note that for balance type `none` the ignore list and upper bounds are ignored (blacklist can be used)

### Testing

```{bash}
conda activate WDX
cd [path/to/store/WarpDemuX]

python -m warpdemux.live_balancing.dummy

python -m warpdemux.live_balancing.dummy --config_file warpdemux/live_balancing/test/config_only_adapter_count.toml
```

## Licensing Information

### WarpDemuX segmentation and read_until modules

The WarpDemuX segmentation module is based on code from tombo (<https://github.com/nanoporetech/tombo>). The WarpDemuX read_until module is based on code from read_until_api (<https://github.com/nanoporetech/read_until_api>), version 3.4.1. These modules have retained their original licenses. The full text of the licenses can be found in the `segmentation` and `read_until` directories.

### Project License

This project, with the exeption of the `segmentation` and `read_until` module, is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0). You can view the full text of the license at the following link:
<https://creativecommons.org/licenses/by-nc/4.0/legalcode>

- The code in `warpdemux/segmentation` is licensed under the terms of the Mozilla Public License 2.0 (MPL 2.0).
- The code in `warpdemux/read_until` is licensed under the terms of the Mozilla Public License 2.0 (MPL 2.0).

### Dependencies Licenses

- **MPL 2.0**: `pod5`, `vbz-h5py-plugin`, `read-until`, `minknow_api`
  - Licensed under the Mozilla Public License 2.0. Full license text available at [MPL 2.0 License](https://www.mozilla.org/en-US/MPL/2.0/).
- **Apache 2.0**: `dtaidistance`
  - Licensed under the Apache License 2.0. Full license text available at [Apache 2.0 License](https://www.apache.org/licenses/LICENSE-2.0).
- **BSD 3-Clause**: `scikit-learn`, `scipy`, `read_until`, `minknow_api`
  - Licensed under the BSD 3-Clause License. Full license text available at [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause).
- **MIT**: `attrs`, `toml`, `tqdm`
  - Licensed under the MIT License. Full license text available at [MIT License](https://opensource.org/licenses/MIT).

Please ensure compliance with each license's terms and conditions.

### Patent Information

An international priority patent application was filed jointly by RKI, HZI and FU Berlin on April 26, 2024, at the European Patent Office (EPO) under number PCT/EP2024/061629.
