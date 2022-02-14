# DigitalEnergyTestbed Modelica Models

This repository contains the [Modelica](https://modelica.org/) models developed as part of research project [DigitalEnergyTestbed](https://energieforschung.at/projekt/offene-testumgebung-zur-evaluierung-von-digitalisierungsloesungen-fuer-integrierte-strom-waermenetze/)


## Package SubstationTeststand

Folder [`substation-teststand`](./substation-teststand) contains package *SubstationTeststand*, which models a test stand for district heating substations.
The package requires the [Modelica IBPSA Library v3.0.0](https://github.com/ibpsa/modelica-ibpsa/releases/tag/v3.0.0) to be loaded.


## Package TestbedExample

Folder [`testbed-example`](./testbed-example) contains package *TestbedExample*, which implements test cases that demonstrate the use of a district heating substation test stand as part of a larger testbed.
These models are used for the [setup of a virtual testbed implementation](https://github.com/AIT-IES/detb-lablink-example).
The package requires the [Modelica IBPSA Library v3.0.0](https://github.com/ibpsa/modelica-ibpsa/releases/tag/v3.0.0), the [Modelica DisHeatLib v1.0.0](https://github.com/AIT-IES/DisHeatLib/releases/tag/v1.0.0) and package *SubstationTeststand* to be loaded.
