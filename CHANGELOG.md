
## [v0.4.3]

### Added

- Added a `--chemistry` flag to the parser. This can be used to specify the latest default config TOML file.

### Fixed

- setup.py: added a `config_files` folder to the package data.
- signal preload size is now explicitly updated based on the `max_obs_trace` value after initializing the config object.

### Changed

- SigProgConfig loader functions now load the config in a two-step approach. First, they load the version&chemistry specific adapted config file, then they update the config with the chemistry-specific or user-provided config file. This allows to remove redundant parameters from the warpdemux config files.
`
- Config names now include the translocation speed of the chemistry, ie: `rna002_70bps@v0.4.3`.
- The created output subdirectory path now contains the version of warpdemux used.
- Now depends on adapted 029a1a513995c3ee878f74dc594a15735d977e68
