**A Comprehensive and Generalized Wood Combustion Model**

This repository contains the UDF source code for a mesoscale wood combustion model. The wood combustion model developed by the authors (see publications listed below) incorporates pyrolysis, char oxidation, and gas phase reactions in a coupled manner. This UDF code provides pyrolysis and char oxidation source terms in the Navier-Stokes equations. A methane combustion mechanism should be used along with this code to simulate gas phase reactions. The code also alters the wood density and thermal conductivity. This code can be used in any CFD simulations wherein wood degradation / combustion is involved. This model has been developed for the combustion of a mixed hardwood species. The pyrolysis, char oxidation reaction rates and energy source terms for a different wood species (eg. a softwood species) need to be experimentally obtained before using the current wood combustion model. 

Theoretical model development and validation data can be found in the publications listed below:

1. Banagiri S., Khadakkar I., Parameswaran M., Meadows J., Lattimer B. Y., A mesoscale CFD model to simulate wood combustion, XI International Seminar on Fire and Explosive Hazards.

2. Banagiri S., Parameswaran M., Khadakkar I., Meadows J., Lattimer B. Y., A reduced wood pyrolysis mechanism for evaluating solid and gas phase parameters, Fuel 381 (2025) 133416, https://doi.org/10.1016/j.fuel.2024.133416.

3. Banagiri S, Meadows J, Lattimer BY. A computational fluid dynamics model to estimate local quantities in firebrand char oxidation. Journal of Fire Sciences. 2023;41(6):241-268. doi:10.1177/07349041231195847.

**Usage Instructions:**

In order to use this code, you will need access to Ansys Fluent and Visual Studio for compiling the C++ code. 
1. You will need to mesh a geometry with a wood (solid zone) and a fluid zone. Conjugate heat transfer needs to be implemented in order to use the code.
2. If you are using the UDF source code for the first time, then you will need to compile and build a UDF library in the same folder where you have saved your source code and Fluent case file. To do this, launch the Native Tools command prompt for Visual Studio (VS) and change directory (cd) to the path where you saved your source code and Fluent case file. Now, launch Ansys Fluent from the command prompt in the **same directory as the Fluent case file and the UDF source code**.  
3. You will need to hook the UDF source terms in your Fluent case setup. This [link](https://www.afs.enea.it/project/neptunius/docs/fluent/html/udf/node140.htm) provides a detailed explanation of how to hook UDFs to your Fluent case.
4. The current UDF source code has different UDFs for setting solid properties (DEFINE_PROPERTY), source terms in the solid and fluid domains (DEFINE_SOURCE), initializing the calculation (DEFINE_INIT), functions to displays errors / warnings (DEFINE_EXECUTE_ON_LOADING), and functions executing at the end of time step (DEFINE_EXECUTE_AT_END). These need to be hooked separately. The detailed procedure for hooking these UDFs is also provided in this [link](https://www.afs.enea.it/project/neptunius/docs/fluent/html/udf/node140.htm).
5. Be sure to check the solid zone and fluid ID's and solid and fluid interface ID's in your case setup. These need to match with the ID's declared in the global variables.
6. Initialize and run the calculation as usual.
7. The execution workflow is as follows: After hooking the UDFs, the DEFINE_EXECUTE_ON_LOADING macro gets executed. This macro displays custom warning/error messages in the Fluent GUI. The DEFINE_INIT macro gets executed once you initialize the Fluent case file. Before running the calculations, Fluent executes the DEFINE_SOURCE macro. This macro supplies the pyrolysis and char oxidation species and energy source terms. These source terms are supplied to the governing equations and impact the flowfield, reaction rates, and heat transfer. As you run the calculations, Fluent performs inner iterations to solve the N-S equations plus additional transport, species, and radiative transport equations. After Fluent performs inner iterations, the DEFINE_PROPERTY macro gets executed. This macro sets the wood density and thermal conductivity as a function of a progress variable. After executing the DEFINE_PROPERTY macro, Fluent performs outer iterations to solve the governing equations. After performing outer iterations, Fluent executes the DEFINE_EXECUTE_AT_END macro. This macro does not directly alter the governing equations; however, this macro consists of several user defined memories (C_UDMIs) which are used in the DEFINE_PROPERTY and DEFINE_SOURCE macros.  
    
**Common Pitfalls/Tips for executing the UDFs:**

1. This source code consists of many user-defined memory allocations (C_UDMIs). These C_UDMIs are used in DEFINE_SOURCE and DEFINE_PROPERTY macros. Be sure to allocate enough user defined memory **prior to initializing the case**. Failure to allocate enough user defined memory will result in SIGSEV segmentation errors. This [link](https://www.afs.enea.it/project/neptunius/docs/fluent/html/udf/node103.htm) describes how to allocate user defined memory.
2. Flow property functions (C_R(c,t) and C_T(c,t) for cell density and temperature) are not available to Fluent before the case is initialized. Therefore, you cannot use these functions in the DEFINE_EXECUTE_ON_LOADING macro.
3. If you want to use DEFINE_SPECIFIC_HEAT to alter the wood specific heat, Fluent does not provide the capability to use / alter cell properties (i.e., C_T(c,t), C_R(c,t) etc.) within this macro.
4. User defined memories on a boundary face (F_UDMIs) cannot be plotted as report plots or visualized using contour plots. Only C_UDMIs may be post-processed using report or contours plots.
5. If you are calculating flow gradients or reconstruction gradients using UDFs, then you must enable temporary retainment of memory using solve/set/advanced/retain-temporary-solver-mem yes in your Fluent journal file.
6. Be carefuly while altering or using global variables. These variables are **not** erased when a simulation is paused , interrupted, or re-initialized. The only way to reset the global variables is to unhook and hook the UDFs.
7. The current UDF is parallelized. It is a good practice to parallize UDF macros (even if you are running them in serial). This [link](https://www.afs.enea.it/project/neptunius/docs/fluent/html/udf/node212.htm) provides detailed instructions on how to do this.

**Citation:**

In order to refer to the wood combustion model, please cite the following articles:

```
@article{char_oxidation_model,
author = {Shrikar Banagiri and Joseph Meadows and Brian Y Lattimer},
title ={A computational fluid dynamics model to estimate local quantities in firebrand char oxidation},

journal = {Journal of Fire Sciences},
volume = {41},
number = {6},
pages = {241-268},
year = {2023},
doi = {10.1177/07349041231195847},

URL = { 
    
        https://doi.org/10.1177/07349041231195847
    
    

},
eprint = { 
    
        https://doi.org/10.1177/07349041231195847
    
}
}

@conference{mesoscale_model,
author = {Shrikar Banagiri and Ishanee Khadakkar and Manjunath Parameswaran and Joseph Meadows and Brian Y Lattimer},
title ={A mesoscale CFD model to simulate wood combustion},
year = {2025},
month = {2025-06-15}
publisher = {XI International Seminar on Fire and Explosive Hazards. Rome, Italy},
language = {en}
}

```

While referring to the pyrolysis mechanism in the wood combustion model, please cite the following article:

```
@article{BANAGIRI2025_pyro_mech,
title = {A reduced wood pyrolysis mechanism for evaluating solid and gas phase parameters},
journal = {Fuel},
volume = {381},
pages = {133416},
year = {2025},
issn = {0016-2361},
doi = {https://doi.org/10.1016/j.fuel.2024.133416},
url = {https://www.sciencedirect.com/science/article/pii/S0016236124025651},
author = {Shrikar Banagiri and Manjunath Parameswaran and Ishanee Khadakkar and Joseph Meadows and Brian Y. Lattimer},
keywords = {Wood pyrolysis, Reduced kinetics, Thermogravimetry, Thermal degradation, Combustion},
abstract = {A reduced wood pyrolysis mechanism comprising of three parallel, first-order reactions was developed to model solid phase heat flow and gas phase combustion parameters. To this end, simultaneous thermal analyzer (STA) studies on a mixed hardwood sawdust sample were conducted for different heating rates (5, 10, and 20 K/min). A nonlinear least-squares model-fitting algorithm was used to extract the Arrhenius kinetic parameters from these experiments. Gaseous species mass flow rates at 1 K/min and 5 K/min were measured by a gas chromatograph (GC) coupled with the STA. The concept of volatile fractions was introduced to correlate the species mass flow rates to the pseudo-component (hemicellulose, cellulose, and lignin) mass loss rates. Sensible specific heats of individual pseudo-components and char were found by using their mass fractions and a progress variable for the degradation. Using the kinetic parameters, gas phase mass fractions, and sensible specific heats, the heats of decomposition of the pseudo-components were evaluated. These parameters formed the inputs of a solid phase heat flow model. This solid phase model was validated against differential scanning calorimetry (DSC) data. The gas phase composition generated by the reduced mechanism was validated against micro-combustion calorimeter (MCC) experiments at different heating rates. The proposed mechanism can be used in wood combustion studies wherein the solid phase heat flow and gas phase combustion processes are coupled. Furthermore, this methodology can be extended to the degradation of any composite material.}
}
```

**Acknowledgment:**

This work was supported by the United States Department of Energy (DOE) through the Office of Energy Efficiency and Renewable Energy (EERE), Project Number: DE-EE0009770. This project is a collaborative effort between the University of Alabama and Virginia Tech. 
