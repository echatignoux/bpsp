# BP-SP Model Evaluation

This repository contains the source code and documentation required to
evaluate and fit the BP-SP model.

## Repository Structure

### 1. `vignettes/`
Contains detailed guides and examples for using the model:
* **`eval_bpsp.html`**: A comprehensive example demonstrating how to fit the BP-SP model to a simulated dataset.
* **`simulate-data.html`**: Documentation detailing the methodology used to generate the simulated datasets.

### 2. `inst/stan/`
The core Stan implementations:
* **`bp_sp.stan`**: The main Stan program.
* **`sp_fun.stan`**: A collection of Stan functions used within the
  main program, including the BP likelihood and spatial analysis
  functions (adapted from
  [ConnorDonegan/Stan-IAR](https://github.com/ConnorDonegan/Stan-IAR/)). 

### 3. `R/`
* **`format_data.R`**: An R utility function designed to preprocess
  and format raw data for compatibility with Stan. 

### 4. `data/`
Includes demonstration data:
* **`dt_cu.rda`**: Baseline values (incidence rates, Se, and PPV by
  age) derived from Corpus Uteri (CU) cancer estimates in France
  (2013–2019) using hospitalization proxies. 
* **`shap_geo.rda`**: A geographical shapefile (`sf` format) covering
  a subset of 704 French postal codes. 
* **`dt_sim.rda`**: A simulated dataset generated according to the
  steps described in the `simulate-data` vignette. 
