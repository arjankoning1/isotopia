# ISOTOPIA

ISOTOPIA is a software package for the prediction of radio-isotope production. It uses cross sections from the IAEA medical isotope data library and the TENDL nuclear data library and calculates radioactive yields as a function of irradiation time and other irradiation characteristics.

## Documentation and reference

A description of the code and its options can be found in the [ISOTOPIA User Manual (pdf)](https://github.com/arjankoning1/isotopia/blob/main/doc/isotopia.pdf).

The reference to be used for ISOTOPIA is:

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155, 1 (2019).

## Installation

### Prerequisites

The following are the prerequisites for compiling and using ISOTOPIA:

- git (only if the package is downloaded via GitHub)
- GNU make
- a recent Fortran compiler, such as GNU Fortran (gfortran)
- the `isotopia.libs/` cross-section database for calculations and sample cases

### Downloads

#### 1. Download the tar files (frozen version ISOTOPIA-2.2)

This is available at the the [TALYS page](https://nds.iaea.org/talys/), and can be downloaded by clicking on the download link or
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

#### 2. Using git (latest beta version)

```bash
git clone https://github.com/arjankoning1/isotopia.git
git clone https://github.com/arjankoning1/isotopia.libs.git
```

Again, keep both repositories in the same parent directory.

By default, ISOTOPIA derives `isotopia.libs/` from the parent directory of `ISOTOPIA_DIR`. If the library is stored elsewhere, use the existing `crosspath` input keyword:

```text
crosspath /path/to/isotopia.libs/
```

A separate `ISOTOPIA_LIBS` environment variable is therefore not required.

### Installation instructions

#### 1. For the tar file (frozen version ISOTOPIA-2.2)

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

#### 2. For the git version (latest beta version)

```bash
cd isotopia
./install_isotopia.bash
```

which automatically executes the `Makefile` in `isotopia/source`.

An alternative is:

```bash
cd isotopia/source
make
```

For the git version, the default compiler is `gfortran`. When `gfortran` is used and no `FFLAGS` are supplied, the Makefile uses:

```text
-w -O3 -ffp-contract=off
```

For other compilers, no default compiler flags are imposed.

Compiler and compilation options can be passed through `install_isotopia.bash`, for example:

```bash
./install_isotopia.bash FC=gfortran FFLAGS="-O3 -ffp-contract=off"
./install_isotopia.bash FC=ifx FFLAGS="-O3"
```

The executable is installed as `isotopia/bin/isotopia`.

Set:

```bash
export ISOTOPIA_DIR="/Users/koning/isotopia"
export PATH="$ISOTOPIA_DIR/bin:$PATH"
export ISOTOPIA_USER="Your Name"
```

These lines can be added to `~/.zshrc` or `~/.profile`.

If setting `ISOTOPIA_DIR` is not possible, edit `code_dir` in `source/machine.f90` and rebuild ISOTOPIA.

For the modern git version, `code_build.bash` and `path_change.bash` are no longer required and can be removed after adopting the new installer, Makefile and `machine.f90`.

## The ISOTOPIA package

The `isotopia/` directory contains:

+ `README.md` is this README file
+ `LICENSE` is the License file
+ `install_isotopia.bash` is  the installation script
+ `source/` contains the Fortran source code of ISOTOPIA and the Makefile
+ `bin/` contains the executable after successful installation
+ `files/` contains abundance and decay data files needed for the calculations
+ `doc/` contains the tutorial in pdf format
+ `samples/` contains the input and output files of the sample cases, and the *verify* script for the user to run the sample cases

The separate sibling `isotopia.libs/` repository contains the cross-section database. In total, about 20 GB of free disk space is required for ISOTOPIA and its data.

## Sample cases

The sample cases assume:

```text
parent_directory/
├── isotopia/
└── isotopia.libs/
```

Run:

```bash
cd samples
./verify
```

or from the top-level directory:

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
