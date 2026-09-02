# ISOTOPIA

ISOTOPIA is a software package for the prediction of radio-isotope production. It uses cross sections from the IAEA medical isotope data library and the TENDL nuclear data library and calculates radioactive yields as a function of irradiation time and other irradiation characteristics.

## Documentation and reference

A description of the code and its options can be found in the [ISOTOPIA User Manual (pdf)](https://github.com/arjankoning1/isotopia/blob/main/doc/isotopia.pdf).

The reference to be used for ISOTOPIA is:

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155, 1 (2019).

## Installation

### Prerequisites

The following are the prerequisites for compiling and using ISOTOPIA:

- GNU make
- a recent Fortran compiler, such as GNU Fortran (gfortran)
- the `isotopia.libs/` cross-section database for calculations and sample cases
- git, only when ISOTOPIA is downloaded using `git clone`

### Downloads

ISOTOPIA and its cross-section library can be downloaded in one of the following ways.

#### 1. Frozen version ISOTOPIA-2.2 (December 2025)

The frozen ISOTOPIA-2.2 distribution is available from the [TALYS page](https://nds.iaea.org/talys/). It can be retrieved by clicking on the download links or with

```bash
curl -LO https://nds.iaea.org/talys/codes/isotopia.tar
tar zxf isotopia.tar

curl -LO https://nds.iaea.org/talys/codes/isotopia.libs.tar
tar zxf isotopia.libs.tar
```

Keep the code and cross-section library as sibling directories:

```text
parent_directory/
├── isotopia/
└── isotopia.libs/
```

This version is fixed and will not change.

#### 2. Latest beta version without git

Users who do not have git can download snapshots of the current `main` branches directly from GitHub:

```bash
curl -L \
  -o isotopia-main.tar.gz \
  https://github.com/arjankoning1/isotopia/archive/refs/heads/main.tar.gz

tar zxf isotopia-main.tar.gz
mv isotopia-main isotopia

curl -L \
  -o isotopia.libs-main.tar.gz \
  https://github.com/arjankoning1/isotopia.libs/archive/refs/heads/main.tar.gz

tar zxf isotopia.libs-main.tar.gz
mv isotopia.libs-main isotopia.libs
```

Keep both directories in the same parent directory:

```text
parent_directory/
├── isotopia/
└── isotopia.libs/
```

The downloaded snapshots contain the latest versions of the `main` branches at the time of download. To obtain newer versions later, download the snapshots again.

#### 3. Latest beta version using git

Users with git can clone both repositories with

```bash
git clone https://github.com/arjankoning1/isotopia.git
git clone https://github.com/arjankoning1/isotopia.libs.git
```

Again, keep both repositories in the same parent directory.

The advantage of this method is that the local installations can subsequently be updated with

```bash
cd isotopia
git pull --ff-only

cd ../isotopia.libs
git pull --ff-only
```

By default, ISOTOPIA derives `isotopia.libs/` from the parent directory of `ISOTOPIA_DIR`. If the library is stored elsewhere, use the existing `crosspath` input keyword:

```text
crosspath /path/to/isotopia.libs/
```

A separate `ISOTOPIA_LIBS` environment variable is therefore not required.

### Installation instructions

#### 1. Frozen version ISOTOPIA-2.2

For the frozen tar distribution:

```bash
cd isotopia
./install_isotopia.bash
```

An alternative is:

```bash
cd isotopia/source
make
```

The frozen distribution retains its own installation scripts and settings.

#### 2. Latest beta version

The installation procedure is identical whether the latest beta version was obtained as a GitHub tar snapshot or using `git clone`.

From the `isotopia/` directory, run

```bash
./install_isotopia.bash
```

which automatically executes the `Makefile` in `isotopia/source`.

An alternative is:

```bash
cd isotopia/source
make
```

The executable is installed as

```text
isotopia/bin/isotopia
```

For the latest beta version, the default compiler is `gfortran`. When `gfortran` is used and no `FFLAGS` are supplied, the Makefile uses:

```text
-w -O3 -ffp-contract=off
```

For other compilers, no default compiler flags are imposed.

Compiler and compilation options can be passed through `install_isotopia.bash`, for example:

```bash
# GNU Fortran
./install_isotopia.bash FC=gfortran FFLAGS="-O3 -ffp-contract=off"

# Intel Fortran
./install_isotopia.bash FC=ifx FFLAGS="-O3"
```

Set:

```bash
export ISOTOPIA_DIR="/Users/koning/isotopia"
export PATH="$ISOTOPIA_DIR/bin:$PATH"
export ISOTOPIA_USER="Your Name"
```

These lines can be added to `~/.zshrc` or `~/.profile`.

If setting `ISOTOPIA_DIR` is not possible, edit `code_dir` in `source/machine.f90` and rebuild ISOTOPIA.

## The ISOTOPIA package

The `isotopia/` directory contains:

- `README.md` this README file
- `LICENSE` the License file
- `install_isotopia.bash` the installation script
- `source/` the Fortran source code of ISOTOPIA and the Makefile
- `bin/` the executable after successful installation
- `files/` abundance and decay data files needed for the calculations
- `doc/` the tutorial in PDF format
- `samples/` the input and output files of the sample cases and the `verify` script used to run the sample cases

The separate sibling `isotopia.libs/` repository contains the cross-section database.

In total, about 20 GB of free disk space is required for ISOTOPIA and its data.

## Sample cases

The sample cases assume:

```text
parent_directory/
├── isotopia/
└── isotopia.libs/
```

To run the sample cases:

```bash
cd samples
./verify
```

From the top-level ISOTOPIA directory, the same test can be started with:

```bash
make -C source check
```

`make check` sets `ISOTOPIA_DIR` automatically and uses the sibling `isotopia.libs/` database.

ISOTOPIA runs as:

```bash
isotopia < isotopia.inp > isotopia.out
```

assuming that `isotopia/bin` has been added to `PATH`.

## License and Copyright

This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
