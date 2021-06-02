# Tests

The implementation of the AMIS algorithm is tested in two ways:

- `test_AMIS.R`: A fixed seed test that verifies that the output of the `amis`
  function remains consistent with previous implementions. This test
  only performs a couple of iterations of the AMIS algorithm and is
  relatively short to run. It is most useful during development.
- `test_ecfd.R` A test that computes the empirical cumulant distribution function
  (ECDF) from the sampled weighted parameters and compares it to the
  ECDF computed from the reference prevalence data. This test is
  longer to run but provide good assurance that the package is running
  the AMIS correctly.
  
Both tests are integration tests, in the sense that they test the
top-level function `amis`, as opposed to its building blocks. As of
June, 2nd 2021, current funding does not allow for the implementation
of unit tests. 

## Prerequisites

To run the tests, you must have the NTD Modelling Consortium trachoma
model installed.  See [installing the NTD trachoma
model](#installing-the-ntd-trachoma-model).

## Fixed seed regression test

Test `test_AMIS.R` runs two iteration of the AMIS algorithm for a
single pixel group (see [pixel groups](#pixel-groups)) and compares the resulting
parameters and weights to trusted data. The random generator seed is
fixed to a known value (default `1`) for the numbers to be
reproducible.

The test can be run from the top-level package directory (the
directory that contains `DESCRIPTION`) as follows

```R
R --no-save < tests/test_AMIS.R
```

The reference data is generated from a trusted implementation of AMIS
used as a starting point for the development of this package. For more
information on how to generate the reference data, see [Generating the
test reference data](#generating-the-test-reference-data) below.


### Generating the test reference data

The fixed seed test `test_AMIS.R` reads in a geostatistical map in
`test_data/prevalence_map.csv`. This map can be re-generated using
`scripts/generate_test_map.csv`.

The expected results of the AMIS algorithm can be generated using
`scripts/AMIS_five_iterations_orig.R`. This script performs five
iterations of the AMIS algorithm and records the result of each
iteration is a separate file `test_data/param_iteration_?.csv` where
`?` is the iteration number.

## Empirical CDF based test

Script `test_ecdf.R` runs the AMIS algorithm for the NTD trachoma
model for a single pixel group. The script compares the Empirical
Cumulant Distribution Function (ECDF) for both the weighted infection
prevalence sampled through the AMIS algorithm and the initial
prevalence data on which the weights calculation is based.

The test can be run from the top-level package directory (the
directory that contains `DESCRIPTION`) as follows

```shell
R --vanilla < test_ecdf.R [plot_dir]
```

It generates one plot per map pixel in the pixel group.  If argument
`plot_dir` is unspecified, plots are written in directory
`tests/ecdf_plots`.

### Expected output

![typical expected ecdf plot](expected_ecdf_plot.png)

# Installing the NTD trachoma model

Navigate to the package top-level (the directory that contains `DESCRIPTION`), _e.g._

```shell
cd /path/to/trachomAMIS/
```
Create a python 3 virtual environment and activate it:

```shell
python3 -m venv .venv
source .venv/bin/activate
```

Install the NTD trachoma model:

```shell
pip install tests/ntd-model-trachoma/
```
