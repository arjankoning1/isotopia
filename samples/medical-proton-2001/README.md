# Proton yield benchmarks from IAEA-TECDOC-1211

These are ISOTOPIA inputs for every proton route in Chapter 5 of
`medical_proton.pdf` for which the report tabulates a physical yield.  Each
case uses the highest endpoint in its corresponding table and `Eback 1.e-6`.
This represents the report's integral yield: the beam slows from the listed
endpoint to below the reaction threshold.  It is not a chosen production
window.

The table's physical yield is in GBq/C.  Compare it with ISOTOPIA's
`Initial production rate [GBq/(mA.h)]`, divided by 3.6, since 1 mAh = 3.6 C.
`cases.tsv` records the report table, product, endpoint, and (where
unambiguous) the tabulated endpoint value.  Run all cases without leaving
activation files in the source directories with:

```bash
ISOTOPIA_DIR=/Users/koning/isotopia \
ISOTOPIA_BIN=/Users/koning/isotopia/bin/isotopia \
./run_cases.bash
```

The comparison is a data-library benchmark, not an exact regression test:
the 2001 values use the report's recommended cross sections, Ziegler stopping
powers, and then-current decay data, whereas ISOTOPIA uses its current medical
and TENDL libraries.  Do not compare the reported `A1` or `A2` columns to the
initial-yield line.  Those are the 1-hour and saturation activities at 1 uA.

Table 5.1.13b provides two incompatible fitted evaluations, so its reference
cell is intentionally labelled rather than given a single value. The report
discusses `Rb(p,xn)Sr-82` but provides no recommended yield table, so it is
not a benchmark case here.
