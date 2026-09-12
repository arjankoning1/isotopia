# Photon benchmark from IAEA-TECDOC-2051

`medical_photo.pdf` does not contain a set of tabulated physical yields like
TECDOC-1211. Its sole directly numerical production-yield benchmark is Table 2:
the end-of-irradiation Mo-99 inventory from twenty 1 mm x 12 mm diameter,
98%-enriched Mo-100 disks irradiated for 6.5 days by a 42 MeV, 8 kW electron
beam. The 0 mm Ta result is 1174.1 GBq.

The input here approximates that stack as a homogeneous 2 cm Mo-100 slab with
the same 1.13097 cm2 frontal area (about 23.1 g), sets the electron current to
8 kW / 42 MeV = 0.19047619 mA, and requests the activity at EOI. Its `fgamma
0.3` is ISOTOPIA's generic electron-to-photon conversion assumption.

Run it with:

```bash
ISOTOPIA_DIR=/Users/koning/isotopia \
ISOTOPIA_BIN=/Users/koning/isotopia/bin/isotopia \
./run_cases.bash
```

This is deliberately a **model comparison**, not a validation target. The
TECDOC value is a FLUKA calculation for a disk stack and electron transport;
ISOTOPIA uses a thick-target Kramers bremsstrahlung spectrum and does not model
the tantalum converter, inter-disk transport, or detailed target cooling. Thus
the 1 mm and 2 mm Ta results in Table 2 (1224.7 and 1229.1 GBq) cannot be
represented by separate ISOTOPIA inputs without an externally supplied photon
spectrum.

The report's other yields are narrative, geometry-specific estimates (for
example the Ra-225 and Xe-123 discussions), rather than enough controlled input
data for a numerical benchmark. Their production routes can still be explored
with ISOTOPIA samples, but they are not added as misleading ratios here.
