# Reactor benchmarks from IAEA-TECDOC-1340

`medical_reactor.pdf` is an operational production manual rather than a
tabulated evaluated-yield database. This suite therefore includes only recipes
that state a target mass, thermal neutron flux, irradiation time, and measured
or expected activity: natural KBr to Br-82, natural chromium to Cr-51, and
Au-197 to Au-198.

The KBr input converts the reported 93 mg KBr capsule to 62.3 mg elemental
bromine. The targets are treated as thin (`selfshield n`) because the report
does not specify the active powder geometry needed to calculate attenuation.
The final activity extracted by the runner includes the stated cooling time.

Run the suite with:

```bash
ISOTOPIA_DIR=/Users/koning/isotopia \
ISOTOPIA_BIN=/Users/koning/isotopia/bin/isotopia \
./run_cases.bash
```

ISOTOPIA uses its built-in HFR neutron spectrum, with `fluxtotal` set equal to
the reported thermal-flux value. This is an approximation: the report supplies
only a thermal-flux scalar, not the full neutron spectrum, so `fluxtotal` is not
strictly the same physical quantity. Epithermal contributions, target geometry,
and reactor-position corrections are also unavailable. Ratios are consequently
useful operational comparisons, not strict nuclear-data regressions. The Cr-51
report value is a lower bound (`>30 GBq`), so its reported ratio should be read
as an upper bound.

Other entries in the manual either give a yield range, a figure only, a
fission-product process, or omit enough irradiation/source detail that a single
ISOTOPIA ratio would be misleading.
