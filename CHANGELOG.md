
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
