# DyeWorkflow

The workflow aims to calculate dye's optical and redox properties. Screening through the predicted properties can result in a number of dyes that fulfill the criteria of the specific system.

The workflow consists of four parts:
1. Preotimisaton
2. Optimisation
3. Calculation of solution-phase Gibbs free energies for neutral and oxidized dye
4. Calculation of excitation energies and oscillator strengths

**Preotimization** 

Assuming that the starting structures are far from the potential energy surface (PES) minimum, the workflow preoptimizes the starting structures using the density functional tight-binding method.

**Optimisation**

Optimization follows preoptimization. From the same starting point (preoptimized structure), workflow optimizes the geometry of neutral and oxidized dye. We optimize the dye's geometry in a vacuum with GFN1-xTB. Successfully optimized geometries are exported to as .xyz files.

**Calculation of solution-phase Gibbs free energies**

The solution-phase Gibbs free energy for both oxidation states consists of two steps. 
First is the ADF calculation with COSMO to create the dye's .coskf file of the solute (the dye). The second step is COSMO-RS calculation that uses .coskf of a solvent (available in ADFCRS database for common solvents) and previously calculated .coskf of the solute.

**Calculation of excitation energies and oscillator strengths**

Excitation energies are calculated using the sTD-DFT method on the neutral geometries.


The files required before starting the workflow are the starting molecular coordinates (.xyz) and the COSMO-RS file (.coskf). 
