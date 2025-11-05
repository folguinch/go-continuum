# GoContinuum

[![DOI](https://zenodo.org/badge/177511811.svg)](https://zenodo.org/badge/latestdoi/177511811)
[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/pylint-dev/pylint)

This program uses an asymmetric sigma clip algorithm to find line-free channels
in the spectrum of the maximum of a data cube. It then uses these channels to
produce a continuum subtracted ms and quality assurance continuum images. 

## Installation

Using `pip`:
```bash
pip install git+https://github.com/folguinch/goco-helpers
pip install git+https://github.com/folguinch/go-continuum
```

### Dependencies

Dependencies are installed by the steps above. These include the following
python packages:
* `numpy`
* `scipy`
* `matplotlib`
* `astropy`
* `casadata`
* `casatasks`

For parallel cleaning, a standalone version of [`CASA`](https://casa.nrao.edu/)
is needed, and the following path needs to be defined:
```bash
export MPICASA="/path/to/bin/mpicasa -n {0} /path/to/casa/bin/casa"
```

### Developer installation

First install [`goco-helpers`](https://github.com/folguinch/goco-helpers),
which provides helper function and scripts. Then install `go-continuum`
following the preferred method through
[`poetry`](https://python-poetry.org/):
```bash
git clone git@github.com:folguinch/go-continuum.git
cd go-continuum
poetry install
```

It can also be installed through `pip`:
```bash
git clone git@github.com:folguinch/go-continuum.git
cd go-continuum
pip install -e .
```

## Usage

The current version can be run as:
```bash
python -m go_continuum.goco config_file.cfg
```

Unlike previous versions, the directory structure of the outputs is generated
by `goco`, while the input structure is determined by the user and passed
along through the configuration file.

### Generating a `config` file

The configuration file tells the program where the data is located and allows
us to customize the different steps in the program. An example with the
available sections and options can be found in
[`goco-helpers`](https://github.com/folguinch/goco-helpers). They also
provide a script to generate a configuration file from the available data.
Fist copy one of the examples in `goco-helpers`, and modify it to generate a 
template (see the 
[documentation of `goco-helpers`](https://folguinch.github.io/goco-helpers/goco-helpers/goco_helpers.html)).
Then generate the configuration file:
```bash
python -m goco_helpers.config_generator template.cfg uvdata1.ms uvdata2.ms ...
```
The script will generate configuration files per science targets in the provided
measurement sets, identifying the name name of the sources and other parameters
(e.g. a guess of the cell/pixel and optimal image sizes).

## Algorithm implementation

At this point self-calibration has not been implemented within the code. We
recommend applying self-calibration tables before running `goco`.

After setting-up the environment, `goco` will compute dirty cubes per SPW to
find a representative spectrum. Dirty images can be cropped to in order to
save memory. Previously calculated dirty cubes can be 
provided in the configuration file under the `dirty` section with options
per SPW `image_name_spwX` where `X` is the SPW number.

To find a representative spectrum and obtain line-free channels:
- [x] Search for the maximum value in each SPW (neglecting the 10 channels at 
the edges)
- [x] Combine maximum values, reject the farthest point from the others and
take the average among the other values to define the *representative source*.
- [x] Extract the spectra at the *representative source*.
- [x] Find the continuum channels from the spectrum of the
*representative source* (see [`AFOLI`]() documentation for default values):
    1. Generates a masked spectrum with a mask containing: the edges of the
    spectrum (`extremes` parameter), requested flagged channels (`flagchans`),
    and invalid values (`invalid_values`).
    2. Runs sigma-clip in the spectrum to filter out line emission/absorption
    (`sigma`, `censtat` and `niter` parameters).
    3. It dilates the mask by requested amount (`dilate`).
    4. Removes masked section below a minimum width (`min_width`).
    5. Masks small gaps between masked sections (`min_gap`).
    6. Saves statistics to a table if requested.
    7. Retrieve the masking levels by changing the number of iterations for
    sigma-clipping
    8. Finds and saves the masked channel ranges in CASA format.
    9. Plots the spectrum and step statistics.
    10. If requested, it writes channel range files for different levels
    factors of the real continuum from sigma-clip.

For continuum then:
- [x] Split the visibilities and average channels:
    1. Split using the channels in the mask (files with `cont_avg`).
    2. Split without any masking (files with `cont_all`).
- [x] If more than 1 EB, concatenate the visibilities.
- [x] Make quality control images using automatic PB cleaning for both
  visibility sets above.
- [x] Clean using auto-masking.

Cleaning steps use the `tclean` parameter `nsigma` to find the stopping
threshold (`nsigma` option in `continuum` section). The `auto-masking`
parameters are determined from the recommended values based on the 75%
baseline set by the `b75` option.

For lines then:
- [x] Run `uvcontsub` (in CASA) for each EB using the mask to subtract the 
    continuum.
- [x] Apply calibration table if needed.
- [x] If more than 1 EB, concatenate the visibilities.
- [ ] Clean the contsub data.

