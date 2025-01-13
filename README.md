# WarpDemuX

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Preprint: BioRXiV](https://img.shields.io/badge/BioRXiV-doi:10.1101/2024.07.22.604276-blue)](https://doi.org/10.1101/2024.07.22.604276)
[![Preprint: RSquare](https://img.shields.io/badge/Research_Square-doi:10.21203/rs.3.rs--4783223/v1-blue)](https://doi.org/10.21203/rs.3.rs-4783223/v1)
[![Zenodo](https://img.shields.io/badge/Zonodo-doi:10.5281/zenodo.14154555-green)](https://doi.org/10.5281/zenodo.14154555)

WarpDemuX is an ultra-fast and high-accuracy adapter-barcoding and demultiplexing tool for Nanopore direct RNA sequencing. It enhances speed and accuracy through:

- Fast processing of raw signals
- Light-weight machine learning algorithms
- Optimized barcode sets
- Real-time enrichment capabilities through barcode-specific adaptive sampling

Currently supports both SQK-RNA002 and SQK-RNA004 chemistries.


## Quick Start

Installation takes approximately 10 minutes:

```{bash}
# Clone this repository with the adapt submodule
git clone --recursive https://github.com/KleistLab/WarpDemuX.git [path/to/store/WarpDemuX]

# Create a new conda environment using the environment.yml file
# We advise to use mamba instead of conda for speed
mamba env create -n WDX -f [path/to/store/WarpDemuX]/environment.yml
mamba activate WDX

# install in editable mode
pip install -e [path/to/store/WarpDemuX]

# For barcode-specific adaptive sampling
cd [path/to/store/WarpDemuX]
pip install -e '.[live-demux]'

```

WarpDemuX depends on ADAPTed, our tool for adapter and poly(A) tail detection. ADAPTed is included as a submodule for now. To make sure WDX can access this module correctly, you need to install WDX in editable mode. See above.

If you don't wish to install in editable mode, you can install the warpdemux and adapted packages separately:

```
cd [path/to/store/WarpDemuX]
pip install .

cd warpdemux/adapted
pip install .
```

## Basic Usage

```{bash}
conda activate WDX
warpdemux demux -i INPUT [INPUT ...] -o OUTPUT -m MODEL_NAME -j NCORES
```

Where:
- `INPUT`: Pod5 file(s) or directory(s) to demultiplex
- `OUTPUT`: Path for run output folder
- `MODEL_NAME`: Model to use (see [Models](#models) section)
- `NCORES`: Number of cores for parallel processing

Make sure to set `NCORES`, especially when running on a cluster. If not specified, WarpDemuX might incorrectly assume all cores of the node are available to you, which affects the memory allocation and may lead to memory issues.

For further instructions and options, run `warpdemux --help` and `warpdemux demux --help`.

### Advanced Usage Options
#### Fingerprinting

To preprocess your data to obtain the fingerprint for each read, without running barcode prediction, use:

```{bash}
warpdemux prep -i INPUT [INPUT ...] -o OUTPUT -m MODEL_NAME -j NCORES ...
```

Use `--export /path/to/valid/config.toml` to specify a custom config file for the preprocessing step, e.g. when you are developing a new model.

#### Handling Short Reads in RNA004


The CNN adapter detection method used for RNA004 has limitations when processing short reads where the adapter signal comprises approximately half of the total signal length. This is due to interference with signal normalization. There are two approaches to handle this:


1. **Enable LLR Fallback (During Initial Run, default)**
   - By default, WarpDemuX enables automatic fallback to the more sensitive LLR method for short reads
   - Note: This can significantly increase runtime if your data contains many adapter-only reads
   - You can disable this behavior by setting the `--export cnn_boundaries.fallback_to_llr_short_reads=false` runtime argument
   
2. **Process Failed Reads Separately**
   - Run WarpDemuX with default or disabled LLR fallback settings first
   - Then use the retry functionality to process failed reads (see below)
   - This may be preferential when you know upfront that you have many short reads in your data

#### Continuing Interrupted Runs
Use the `continue` subcommand to resume an interrupted run:

```{bash}
warpdemux continue /path/to/previous/run/output/WDX[n_barcodes]_[chemistry]_[version]_[UUID]
```

To modify performance parameters before continuing:
1. Locate the `command.json` file in the previous run's output folder
2. Edit performance parameters (e.g., batch size, number of processes)
3. Important: Do not modify processing parameters as this may cause failures


#### Recovering Failed Reads

For RNA004 runs, you can attempt to recover failed reads using a more sensitive detection method:

```{bash}
warpdemux retry /path/to/previous/run/output/WDX[n_barcodes]_[chemistry]_[version]_[UUID]
```

Recovered reads are added to the main output folder. Note that retrying also improves yield for non-short reads.

## Barcodes

WarpDemuX uses custom barcode sequences embedded within the RTA (Reverse Transcription Adapter) during library preparation. 

### Available Barcodes

| Barcode ID | Barcode sequence    |
|------------|-------------------|
| Barcode 1  | TTTTTACTGCCAGTGACT |
| Barcode 2  | AGGGGAGAGAGCCCCCCC |
| Barcode 3  | CACGTCATTTTCCACGTC |
| Barcode 4  | GGAGGCCAGGCGGACCGA |
| Barcode 5  | ACGGACCTTTTGACTTAA |
| Barcode 6  | TATTGCATACTGCGCCGC |
| Barcode 7  | CCACGGAGGGAGGATTGG |
| Barcode 8  | TTACCGGCAGTGACGGAC |
| Barcode 9  | CGAGATTGCATCCCCCCC |
| Barcode 10 | TACCACCTGCCGGCGGCC |
| Barcode 11 | GCCCGCCGGGGGAGAAGC |
| Barcode 12 | TTTTTTTTACCGGCAGTT |

For detailed information about barcode design, optimization, and performance characteristics, please refer to our manuscript.


## Models

### Available models


WarpDemuX RNA004 models are optimized for reads with poly(A) tails. For datasets with many short poly(A) tails, you have two options:

1. Enable LLR fallback with `--export cnn_boundaries.fallback_to_llr_short_reads=true`
2. Switch to LLR detection entirely (slower but better handles short poly(A) tails)

Note: CNN detection is only available for RNA004, while LLR detection is the standard method for RNA002.

| Model Name | Chemistry | Samples | Barcodes Used |
|------------|-----------|----------|---------------|
| WDX-DPC_rna002_v0_4_4 | RNA002 | 4 | Original DeePlexiCon barcodes |
| WDX4_rna002_v0_4_4 | RNA002 | 4 | WDX_bc04, WDX_bc05, WDX_bc06, WDX_bc08 |
| WDX6_rna002_v0_4_4 | RNA002 | 6 | WDX_bc01, WDX_bc03, WDX_bc05, WDX_bc06, WDX_bc07, WDX_bc11 |
| WDX8_rna002_v0_4_4 | RNA002 | 8 | WDX_bc01, WDX_bc03, WDX_bc05, WDX_bc06, WDX_bc07, WDX_bc09, WDX_bc11, WDX_bc12 |
| WDX10_rna002_v0_4_4 | RNA002 | 10 | WDX_bc01, WDX_bc02, WDX_bc03, WDX_bc05, WDX_bc06, WDX_bc07, WDX_bc09, WDX_bc10, WDX_bc11, WDX_bc12 |
| WDX12_rna002_v0_4_4 | RNA002 | 12 | All 12 WDX barcodes |
| WDX4_rna004_v0_4_4 | RNA004 | 4 | WDX_bc03, WDX_bc04, WDX_bc05, WDX_bc07 |

**Coming soon:**
- RNA004 tRNA models
- RNA004 non-polyA models
- RNA004 larger barcode sets

## Target Accuracy Modes


WarpDemuX features a flexible accuracy control system that lets you optimize the balance between prediction accuracy and data yield. Each model includes barcode-specific calibrated thresholds that account for unique characteristics and error patterns.

### Available Target Accuracies

- 95.0%
- 97.5%
- 99.0%
- 99.5%
- 99.9%

### How It Works

The confidence thresholds are determined through a barcode-specific calibration procedure:


1. Baseline error rates are established for each barcode using calibration datasets
2. Model- and barcode-specific confidence thresholds are computed for each target accuracy level
3. Predictions are filtered using these thresholds to maintain accuracy targets

### Using Target Accuracy Modes

Currently, you can apply target accuracy filtering post-prediction:

1. Group results by predicted barcode
2. Apply the desired confidence threshold from `target_accuracy_thresholds/*.csv`
3. Filter predictions that meet or exceed the threshold

*Note: Direct integration of accuracy modes into the main workflow is planned for a future release.*

### Choosing the Right Target

Consider these factors when selecting a target accuracy:

- **Higher Accuracy (99.0-99.9%)**
  - Fewer but more reliable assignments
  - Best for applications requiring high precision
  - Recommended for critical analyses

- **Lower Accuracy (95.0-97.5%)**
  - Larger dataset with more assignments
  - Suitable for exploratory analyses
  - Better for applications that can tolerate some noise

Choose based on your experimental requirements and the relative importance of precision versus dataset size.

## Performance

### Runtime
- Scales linearly with number of cores
- RNA004: ~2-3 minutes per 100,000 reads (8 cores, standard laptop)

### Memory Requirements
- Recommended: 2GB RAM per core (with default minibatch size of 1000)
- Example: 16GB RAM for 8 cores

## Expected input and output 

### Input
WarpDemuX exclusively supports the pod5 file format. You must convert any fast5/slow5/blow5 files to pod5 format before processing.

### Output


WarpDemuX creates an output directory named `WDX[n_barcodes]_[chemistry]_[version]_[UUID]` with the following structure:

```
WDX[n_barcodes]_[chemistry]_[version]_[UUID]/
├── failed_reads/ # Contains statistics for reads where adapter detection failed
├── predictions/ # Contains barcode predictions for successful reads
└── detected_boundaries/ # Optional: Created when --save_boundaries true
└── fingerprints/ # Optional: Created when --save_fpts true
```

### Predictions Output
The `predictions/` directory contains `barcode_predictions_[INDEX].csv` files, with each file containing `batch_size_output` reads (configurable via command line). For each successfully processed read, the following information is recorded:

| Column Name | Description |
|------------|-------------|
| read_id | Pod5 file read ID (**Note**: May differ from basecalling BAM file read ID, see [split reads](#split-reads)) |
| predicted_barcode | Predicted barcode (`-1` indicates noise classification) |
| confidence_score | Prediction confidence score |
| p01-p12 | Individual probability scores for each barcode |
| p-1 | Probability score for noise class |

Note: Available probability columns (p01-p12) vary by model.


### Failed Reads Output
The `failed_reads/` directory contains `failed_reads_[INDEX].csv` files with detailed signal statistics for reads where adapter detection failed. Key metrics include:

**Signal Information:**
- read_id: Pod5 file read ID
- signal_len: Total signal length
- preloaded: Length of analyzed preloaded signal

**Adapter Statistics:**
- adapter_start/end/len: Position and length
- adapter_mean/std/med/mad: Signal statistics (pA)

**PolyA Tail Statistics:**
- polya_start/end/len: Position and length
- polya_mean/std/med/mad: Signal statistics (pA)
- polya_candidates: Alternative end positions

**RNA Transcript Statistics:**
- rna_preloaded_start/len: Position and length (post-polya)
- rna_preloaded_mean/std/med/mad: Signal statistics (pA)

**Detection Method Results:**
- llr_*: LLR method detection results
- cnn_*: CNN method detection results
- mvs_*: MVS method detection results

**Quality Metrics:**
- adapter_med_dt/mad_dt: Dwell time statistics. Median/median absolute deviation for the dwell time per translocation event in the adapter.
- fail_reason: Specific cause of detection failure

### Boundaries
When `--save_boundaries true` is set, successfully detected boundaries are saved to `detected_boundaries_[INDEX].csv` files in the `detected_boundaries/` directory. These files contain the same columns as failed reads output, except for the `fail_reason` column.

## Basecalling and split reads

With newer versions of Dorado (the basecaller), single reads from the pod5 file may be split into one or more reads during basecalling. See [this issue](https://github.com/nanoporetech/dorado/issues/848) for more details. When this happens:

1. Each child read is assigned a new read ID in the basecalling output
2. The original pod5 read ID is stored in the `pi:Z` field of the BAM file

### Linking predictions to basecalled reads
Currently, users need to manually link WarpDemuX predictions with basecalled split read outputs using the `pi:Z` field in the BAM file. We will provide tooling to automatically handle this linking given a BAM file and an output directory soon.

## PolyA tail length estimation

WarpDemuX calculates various signal statistics for each read during adapter detection validation. These statistics include an estimate of the poly(A) tail length. By default, these statistics are not saved, but can be enabled using the `--save_boundaries true` command line option.

### Accuracy considerations

The poly(A) tail length estimation has some important limitations to note:

- The default settings are optimized for demultiplexing speed rather than accurate poly(A) measurement
- Longer poly(A) tails may be underestimated
- Results should be considered approximate estimates rather than precise measurements

### Improving accuracy

If accurate poly(A) tail length estimation is important for your analysis, you can adjust the settings:

1. Increase the amount of signal preloaded for analysis by setting a higher `max_obs_trace` value:
```{bash}
warpdemux demux [...] --export core.max_obs_trace=XYZ
```
2. Always enable boundary saving to capture the measurements:
```{bash}
warpdemux demux [...] --save_boundaries true
```
Note that increasing `max_obs_trace` will result in slower processing times, so consider the tradeoff between accuracy and speed based on your specific needs.

## Barcode-based adaptive sampling (Live Balancing)

Live balancing has been executed as a proof-of-concept for RNA002. 

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

## Testing

WarpDemuX has been tested on Ubuntu 20.04.6 LTS and Ubuntu 22.04.4 LTS, as well as WSL 2 on Windows 11.

### Demux

```{bash}
warpdemux demux -i [path/to/WarpDemuX]/test_data/demux/4000_rna004.pod5 -j 8
```

### Live balancing

```{bash}
conda activate WDX
cd [path/to/WarpDemuX]

python -m warpdemux.live_balancing.dummy --config_file warpdemux/../test_data/live_balancing/config_only_read_count.toml
```

## Troubleshooting


### Compilation Errors

When you first run WarpDemuX, it automatically compiles the Cython code. If you encounter a compilation error like:


Try these solutions in order:

1. **Load GCC Module** (if using a cluster/HPC system):
   ```bash
   module avail              # List available modules
   module load GCC          # Load the GCC module
   warpdemux demux [...]    # Run WarpDemuX
   ```

2. **Create Clean Environment** (if module loading doesn't work):
   ```bash
   # Create new environment with minimal dependencies
   mamba create -n WDX2 python=3.8
   mamba activate WDX2
   
   # Install WarpDemuX and ADAPTed via pip
   pip install [path/to/store/WarpDemuX]
   pip install [path/to/store/WarpDemuX]/warpdemux/adapted
   
   # Run WarpDemuX
   warpdemux demux [...]
   ```

If you continue to experience issues, please open a GitHub issue with your system details and the full error message.

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

## How to Cite

If WarpDemuX has been helpful for your research, please cite our work:

> Demultiplexing and barcode-specific adaptive sampling for nanopore direct RNA sequencing  
> van der Toorn W, Bohn P, Liu-Wei W, Olguin-Nava M, Smyth RP, von Kleist M  
> *bioRxiv* 2024.07.22.604276; doi: [https://doi.org/10.1101/2024.07.22.604276](https://doi.org/10.1101/2024.07.22.604276)

You can also cite our Research Square preprint:

> Demultiplexing and barcode-specific adaptive sampling for nanopore direct RNA sequencing  
> Max von Kleist, Wiep van der Toorn, Wang Liu-Wei et al.  
> *Research Square* PREPRINT (Version 1), 02 September 2024  
> doi: [https://doi.org/10.21203/rs.3.rs-4783223/v1](https://doi.org/10.21203/rs.3.rs-4783223/v1)
