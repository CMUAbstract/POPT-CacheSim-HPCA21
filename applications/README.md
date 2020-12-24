# (Instrumented) Applications

Applications whose memory access pattern is analyzed via the 
cache simulator. 

## Adding simulator hooks

All simulator-related hooks begin with the prefix `PIN_`. 

To add a new simulator hook, look into _function replacement_ in the 
simulator README. The hook should be defined in **both** the application
and the simulator

The `popt` and `opt-ideal` versions of applications require more information
exchange between the application and the simulator. To compare the different 
hooks required across different versions, compare `pr.cc` across different
directories.

### Hooks for different types of graph applications

The repo currently includes different versions of `pr.cc` which is representative of 
Sparse Matrix Dense Vector (SpM-DV) multiplication.

We will be including more applications that are examples of popular graph computation patterns:
* `pr_delta.cc`: An application using Push-Pull direction switching where the pull phase is 
                 simulated (representative of SpM-SpV multiplication). The `popt` and `opt-ideal`
                 versions need to compute rereferences for two irregularly-accessed data structures 
                 (_srcData_ and _frontier_)
* `cc_sv.cc`   : A Push version of SpM-DV. The `popt` and `opt-ideal` versions will use the CSC as
                 the transpose that is used for computing rereferences
