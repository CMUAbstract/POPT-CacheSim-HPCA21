# Cache Simulators for P-OPT

Repo for **P-OPT: Practical Optimal Cache Replacement for Graph Analytics** (HPCA 2021)

This repo contains cache simulators for the following cache replacement policies:
* LRU
* DRRIP
* P-OPT
* Transpose-driven Belady's OPT replacement policy (T-OPT)


## Repo Organization

* `simulators`   : Pin-based cache simulators for the different replacement policies
* `applications` : Annotated graph applications designed to work with the simulators
* `scripts`      : Helper scripts to launch experiments (and plot results)

## Requirements

The cache simulators are pintools designed to work with Pin-2.14 (Download script provided in `scripts/download_pin.py`)

The pintools require an older version of g++ (v4.9) 

All the applications and scripts were tested on Debian Stretch with g++-6.3.0 and python (v3.5.3)

## Usage Instructions

The following scripts should be run from the `scripts` directory:

* **Step 1**: Run `download_pin.py` (to install Pin) followed by `download_and_build_graphs.py` (to build all inputs)  
* **Step 2**: Run `run_cache_sims.py` to start the cache simulations
* **Step 3**: Run `plot_llcmiss_red.py` to plot LLC miss reduction from different policies 

## Adding More Policies/Applications

`simulators/README.md` contains information on defining new LLC replacement policies and the `applications/README.md`
directory provides information for creating new applications to evaluate the policies on. 

## Feedback 

Please create an [issue](https://github.com/CMUAbstract/POPT-CacheSim-HPCA21/issues) for reporting bugs or providing feedback.

