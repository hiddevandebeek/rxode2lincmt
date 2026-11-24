# Carry pairs a linCmt() model can carry at once

With an rxode2 that keeps the carry in its solve buffer the carry is
sized to each model, so there is no cap and this is
`.Machine$integer.max`; with an older rxode2 it is the fixed array's 8.
nlmixr2est keeps a model's carry within it.

## Usage

``` r
.linCmtCarryMaxPairs()
```

## Value

integer
