# MATLAB-FEA-of-plate-with-circular-hole
# 2D FEM Analysis of a Plate with Circular Hole

## Overview

This project presents a **2D Finite Element Method (FEM)** analysis of a tensile-loaded plate containing a circular hole. The FEM solver is developed in **MATLAB** using **CST (Constant Strain Triangle)** elements.

The main objective is to study the **displacement, strain, stress distribution, and stress concentration** around the circular hole and to validate the numerical results using analytical elasticity and **Abaqus FEA**.

## Methodology

The analysis follows the standard FEM procedure:

1. Create the 2D plate geometry with a central circular hole.
2. Generate a triangular mesh using **CST elements**.
3. Define material properties such as Young's modulus and Poisson's ratio.
4. Apply fixed boundary conditions and tensile loading.
5. Assemble the element stiffness matrices into the global stiffness matrix.
6. Solve the FEM system:

   `K × U = F`

7. Calculate element strain and stress from the obtained nodal displacements.
8. Calculate **von Mises stress** and identify the maximum stress near the hole.
9. Calculate the **stress concentration factor (Kt)**.
10. Compare FEM results with the **Kirsch analytical solution**.
11. Cross-verify the MATLAB FEM results with **Abaqus FEA**.

## FEM Formulation

The project uses **3-node CST triangular elements** under **plane-stress conditions**.

The main FEM equations are:

`K × U = F`

`ε = B × U`

`σ = D × ε`

where:

- `K` = Global stiffness matrix
- `U` = Nodal displacement vector
- `F` = Applied force vector
- `B` = Strain-displacement matrix
- `D` = Material/constitutive matrix
- `ε` = Strain
- `σ` = Stress

## Stress Concentration

The stress concentration factor is calculated as:

`Kt = σmax / σnominal`

For an infinite plate with a circular hole under uniaxial tension, the **Kirsch analytical solution** gives:

`Kt = 3`

The FEM result is compared with this theoretical value to evaluate the accuracy of the numerical solution.

## Validation

The MATLAB FEM results are validated through:

- **Kirsch analytical elasticity solution**
- **Abaqus FEA**
- **Mesh refinement/convergence study**

The comparison focuses mainly on **stress distribution, maximum stress, stress concentration factor, and displacement**.

## Software & Technologies

- **MATLAB**
- **Finite Element Method (FEM)**
- **CST (Constant Strain Triangle)**
- **Abaqus**
- **Plane Stress**
- **Elasticity**
- **Kirsch Solution**
