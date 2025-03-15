**A Comprehensive and Generalized Wood Combustion Model**

This repository contains the UDF source code for the mesoscale wood combustion model. This UDF code provides pyrolysis and char oxidation source terms in the Navier-Stokes equations. The code also alters the wood density and thermal conductivity.

Additional information can be found in the publications listed below:


1. Banagiri S., Khadakkar I., Parameswaran M., Meadows J., Lattimer B. Y., A mesoscale CFD model to simulate wood combustion, XI International Seminar on Fire and Explosive Hazards.

2. Banagiri S., Parameswaran M., Khadakkar I., Meadows J., Lattimer B. Y., A reduced wood pyrolysis mechanism for evaluating solid and gas phase parameters, Fuel 381 (2025) 133416, https://doi.org/10.1016/j.fuel.2024.133416.


**Usage Instructions:**

In order to use this code, you will need access to Ansys Fluent. 
1. You will need to mesh a geometry with a wood (solid zone) and a fluid zone. Conjugate heat transfer needs to be implemented in order to use the code.
2. You will need to hook the UDF source terms in your Fluent case setup. This [link](https://www.afs.enea.it/project/neptunius/docs/fluent/html/udf/node140.htm) provides a detailed 
