## Unreleased

### Changed

- command line argument `--input`/`-i` is now required for the `demux` and `prep` subcommands.
- command line argument `--model_name`/`-m` is now required for the `demux` and `prep` subcommands.
- The progress bar now updates every second, rather than every 10 seconds.
- You can now change certain parameters at the command line when running `continue`. Available parameters are: `num_proc`, `batch_size_output`, `minibatch_size` and `model_name`. The same holds for the `predict` subcommand, see below.

### Added

- The `predict` subcommand has been added. It can be used to predict barcodes for previously preprocessed reads (subcommand `prep`).

### Removed

- Removed the `read_id_csv_colname` argument. When providing an `read_id_csv` file, the file has to contain a `read_id` column.
- The `retry` subcommand has been removed. `adapted` now handles retrying failed reads during the initial run. You can turn automatic retrying on by setting the `fallback_to_llr` parameter in the `cnn_boundaries` or `start_peak_detection` sections of the config file.

## [v0.4.5] - 2024-12-09

### Changed

- WarpDemuX now depends on ADAPTed v0.2.4. This includes a new detection workflow that uses the inital peak in the adapter signal as baseline for the RNA start detection.
- Date and time of the run have been added to the output directory name.

### Added

- `--export` flag to pass particular configuration variables or a full config file.
- `prep` subcommand to preprocess reads for WarpDemuX, this includes the detection of the boundaries and the fingerprinting of the barcode sequences.

### Fixed

- Missing `bottleneck` dependency for editable WarpDemuX installation.

## [v0.4.4] - 2024-11-13

### Changed

- WarpDemuX now depends on ADAPTed v0.2.3. This includes a new detection workflow that uses a CNN to predict the boundaries of the adapter and polyA signals. The CNN method replaces the previous llr workflow as the primary detection method and provides faster detection. The CNN workflow comes with a depency on torch.
- Barcode fingerprints are not saved by default. Use `--save_fpts` to save them.
- Boundaries are not saved by default. Use `--save_boundaries` to save them. In `continue` mode for previous runs in which `save_boundaries` was `False`, previously failed reads will be reanalysed together with the remaining reads. This means that the percentage of succesfully detected adapters will be artificially deflated.
- `batch_size` was renamed to `batch_size_output`, it now specifies the number of reads per output file rather than the number of reads per thread.
- The file processing module has been improved: no itermediate files are created anymore, barcode predictions are streamed directly to output files. Nearly all functions in this module have been rewritten. A seperate thread handles the preloading of the reads into memory (io operations in big batch sizes, enqueuing in minibatches), the queue/worker combinations explicitely takes care to enqueue at most `num_proc` minibatches to limit memory usage.
- Added short names for several parser arguments: `-i` for `--input`, `-o` for `--output`, `-b` for `--batch_size_output`, `-m` for `--model_name`.
- `config.batch.bidx_passed/failed` is now `config.batch.batch_idx_pass/fail`.
- Reads are no longer indexed before running WarpDemuX. This saves time and memory, but means that the progress bar now only shows the number of reads processed rather than the percentage of the total.
- Results directories now include the model name (WDX4-12)and a random UUID.
- the `pod5` executable installed with WarpDemuX is used to calculate the total number of reads to process parallel to the demux process. At the start of demultipexing the progress bar does not have a 'total' value. This is updated once the `pod5` process has finished.
- The parser flag for setting the `num_proc` was changed to `-j`.
- The file processing module has been rewritten to use a separate thread for each type of worker, and to enqueue minibatches of reads rather than individual reads. This reduces memory usage and allows for more efficient I/O.
- Signal normalization has been simplified. Outliers (spikes) are now clipped to a MAD-based estimate, rather than locally imputed.
- Major refactoring of the `file_proc` module.


### Removed

- The subcommands `fpts` and `resegment` have been removed. WarpDemuX now only supports demux from pod5 (`demux`) or continuing from an incomplete previous WarpDemuX run (`continue`).
- Several parser arguments were removed: `minibatch_size`, `DEBUG`.
- The `SigNormConfig` section of the `SigProcConfig` has been removed. Normalization is handled by ADAPTed.

### Added

- Retry functionality. In addition to continuing from a previous run, WarpDemuX can now retry processing failed reads from a previous run (`warpdemux retry`). For this, a slightly slower adapter/polyA detection workflow is used (see ADAPTed `combined_detect_llr2`). Automatic retrying will be integrated into the main detection workflow in a future release.
- RNA004 support (WDX4), with barcode-specific target accuracy confidence score filtering (`target_accuracy_thresholds` folder).
- Test data for live balancing dummy test and demultiplexing has been added in the `test_data` folder.
- Extensive updates to the documentation.

### Fixed

- WarpDemuX now correctly loads the chemistry-specific config files. Previously, whole sections of the chemistry-specific (ADAPTed) config file were overwritten by the WarpDemuX config file.
- Errors in live balancing dummy test have been fixed.


## [v0.4.3] - 2024-09-11

### Added

- Added a `--chemistry` flag to the parser. This can be used to specify the latest default config TOML file.
- command.json file is now saved in the output folder and contains the command used to run WarpDemuX.
- Logging: process outputs are now logged to the `warpdemux.log` file and to stdout.
- WarpDemuX now supports continuing from a previous (incomplete) run using the `continue <continue_from_path>` subcommand.

### Fixed

- setup.py: added a `config_files` folder to the package data.
- signal preload size is now explicitly updated based on the `max_obs_trace` value after initializing the config object.

### Changed

- SigProgConfig loader functions now load the config in a two-step approach. First, they load the version&chemistry specific adapted config file, then they update the config with the chemistry-specific or user-provided config file. This allows to remove redundant parameters from the warpdemux config files.
`
- Config names now include the translocation speed of the chemistry, ie: `rna002_70bps@v0.4.3`.
- The created output subdirectory path now contains the version of warpdemux used.
- Now depends on adapted v0.2.2.
- Event segmentation now relies on true peak detection in the t-scores, rather than sort-and-select. The `num_events` highest peaks in the calculated t-statistics, as detected with a minimal distance of `min_obs` are returned as changepoints.
- `save_dwell_times` in parser is False by default (True before).
- The output directory is now named after the version of WarpDemuX and a random UUID rather than the current date and time.

### Removed

- The `--create_subdir` argument is removed. The output directory is now always created in the specified output folder.
