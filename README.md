# Description

This package provides an implementation of the Adaptive Multiple
Importance Sampling algorithm, as described in

_Integrating Geostatistical Maps And Transmission Models Using Adaptive Multiple Importance Sampling_
Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez, T. Deirdre Hollingsworth, Simon E.F. Spencer
medRxiv 2020.08.03.20146241; doi: https://doi.org/10.1101/2020.08.03.20146241 

Currently, the implementation is limited to 1D models.

# Installation

Make sure you have the package [devtools](https://devtools.r-lib.org/)
installed. Then

```R
devtools::install_github("OxfordRSE/trachomAMIS")
```

# Usage

The package exports a single function `amis` that takes a
geostatistical map, a model and returns sampled parameters and their
associated weights.

```R
param_and_weights <- trachomAMIS::amis(geo_map, model, amis_params)
```

- `geo_map`: A matrix representing the geostatistical map, with one
  row per pixel.
- `model`: A function implementing the model. It can be anything, as
  long as it conforms to a specific interface. See defining a model
  function.
- `amis_params`: A list containing the parameters for the AMIS algorithm.
  - `nsamples`: The number of sample parameters to draw at each iteration.
  - `target_ess`: The target effective sample size.
  - `T`: The maximum number of iterations.
  - `delta`: The Randon-Nikodym parameter.
  
## Defining a model function

The `amis` function expects its argument `model_func` to be a function with the
following interface

```
observables <- model_func(seeds, parameters)
```

- `parameters`: A matrix of parameter values (`double`)
- `seeds`: A vector of seeds (`integer`)

Function `model_func` is expected to run the model for each pair
(`seed`, `parameter`) and return the corresponding values for the
observable (_e.g._ infection prevalence).

## Wrapping a model

If your model does not comply with the above interface, you can
provide a wrapper function. The following example illustrates how it
can be done in the case where the function implementing the model
reads and writes values in files and expects several additional
parameters to be specified:

```R
wrapped_model <- function(seeds, parameters) {
	# write input on disk
	input_file <- "beta_values.csv"
    write_model_input(seeds, parameters, input_file)

	# run model
    run_model(input_file, mda_file, output_file, infect_output,
               SaveOutput = F, OutSimFilePath = NULL,
               InSimFilePath = NULL)

	# read results and return obsvervable
    res <- read.csv(output_file)
    return(100 * res[, dim(res)[2]])
  }
```
