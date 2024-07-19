# Dev Setup Install Instructions

1. install miniconda
2. Add conda-forge and bioconda channels: <https://bioconda.github.io/>
3. create and activate a python 3.8 conda environment
4. Install additional packages: `conda/mamba install jupyterlab plotly python-kaleido=0.1.0 pandas numpy bottleneck joblib scikit-learn=1.3.1 cython=0.29.36 dtaidistance tqdm toml`
5. `pip install pod5 minknow_api`
    (or on windows `pip install pod5 minknow_api`)
6. Clone WarpDemuX github `git clone https://github.com/wvandertoorn/WarpDemuX`

## Other notes

### Replaying Flongle runs on MinION simulated devices

Sadly MinKNOW does not support the addition of a Flongle simulated device.
However, for dRNA sequencing of IVT RNA Flongles are the best choice, and
we have performed most of our runs to generate training data and validate
on FLO-FLG001 flow cells.

Replaying of Flongle recordings itself is not an issue - the issue lies in
the requirement to have a physical MinION Mk1b + Flongle adapter plugged
into the PC. Attempting to replay a Flongle file on a MinION simulated
device results in an error when simulation.py attempts to set the calibration
data, because the batch file only contains data for channels 1-126, while
the set_calibration function of the gRPC interface expects all channels of
the MinION device to have values. Notably, when no calibration data can be
found in the replay file, MinKNOW by itself sets them all to (0,2048). Thus,
we decided to manually set the calibration of channels 127+ that we don't have
data for anyway by adding the following code to
MinKNOW/ont-python/Lib/site-packages/bream4/utility/simulation.py (line 33):

```python
if len(calibration) < len(device.get_channel_list()):
    device.logger.info("Replaying a Flongle data on MinION device. Adding empty calibrations.")
    for channel in device.get_channel_list():
        if channel not in calibration:
            calibration[channel] = (0, 2048)
```

Next is to add a simulated MinION position. To do this activate the conda
environment, then run
`python -m minknow_api.examples.manage_simulated_devices --add MS00000`
to add a simulated MinION device.

Next you can click on Start Run, select the simulated device, select the
appropriate flow cell (this will determine which sequencing protocols you see) -
select the correct sequencing protocol and add your bulk fast5 file
to replay. Then you can start the run (we recommend recording the optional
event and classification data, but raw data does not need to be re-recorded).

Note that after the run starts there is one more thing to take care off:
The timing of pore scans between Flongle and MinION flow cells is not synced,
as such pores will abberantly be disabled by MinKNOW at the start of the run.
To reactivate them, a manual pore scan at minute 4-5 is required.

### minknow_api Handshake failed

If when trying to add a simulated device the error "handshake failed"
appears, the grpcio installation appears to be faulty. Remove all versions
of grpcio from the conda environment (also via pip), remove minknow_api,
then reinstall minknow_api. This will also install a compatible version
of grpcio.

### Adding/Editing custom sequencing protocols

Custom sequencing protocols can be added by placing them in the
<MinKNOW/conf/package/sequencing> directory. After adding/editing a protocol,
in MinKNOW press Start, then the 3 dots on the right and "reload scripts".

## Live Adapative Sampling Checklist

0. Have poly-adenylated RNA

1. Copy the proper protocol and channel_states file with shorter break_read_seconds
 time over to MinKNOW/conf/package/sequencing/

2. Open MinKNOW, click Start, press three buttons and reload scripts (just to be sure)

3. Configure the sequencing, selecting the new protocol. Then turn off reserve reads,
set pore scan interval to 4h (or more).

4. In the Output screen, select the correct output folder, select pod5 read output,
set number of reads to file to 50 (or 100). Enable Bulk output, including raw data.

5. Start the run, wait ~3-4 min until pore scan is finished.

6. Add pod5 directory (it will only be generated after pore scan is finished) to the
correct warpdemux .toml file.

7. Start WarpDemuX analysis
