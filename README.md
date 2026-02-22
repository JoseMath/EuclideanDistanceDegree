# EuclideanDistanceDegree

A Macaulay2 package for computing Euclidean distance degrees and solving nearest point problems.

This repository contains a Macaulay2 package implementing symbolic and numerical methods for computing Euclidean Distance Degrees (ED degrees) of algebraic varieties and for solving associated nearest point problems.

## Overview

Given a variety X ⊂ CC^n or PP^{n−1} and a data point u, the Euclidean Distance Degree counts the number of complex critical points of the squared distance function ||x − u||^2 restricted to X.

This package:
- Builds the Lagrange multiplier system for the ED problem.
- Provides symbolic and numerical methods to compute ED‑degrees.
- Offers homotopy continuation tools (via Bertini) for solving the nearest point problem.

## Repository Structure

### Core package
- EuclideanDistanceDegree.m2 — main package file exporting core functions.

### Building equations
- LagrangeEquations.m2 — constructs ED critical equations.
- conormal.m2 — conormal/dual constructions.

### Symbolic computation algorithms
- EDD_Determinantal.m2 - using the minors of the augmented Jacobian matrix to find the ideal of critical points
- EDD_LeftKernel.m2 — using the left kernel of the augmented Jacobian matrix to find an ideal that defines critical points along with a left kernel vector
- EDD_ProjectiveLeftKernel.m2 — similiar to EDD_LeftKernel but includes additional homogenization steps for a homogeneous system of equations.

### Numerical homotopy methods
- EDD_Numerical.m2 — numerical pipeline 
- EDD_ParameterHomotopy.m2 — parameter homotopy for weighted and unweighted problems.

### Examples & demos
- EDD_demo.m2
- ExamplesSymbolicEDDegree.m2
- ExamplesNumericEDDegree.m2
- ExampleEDDProjectiveHomotopy.m2 — detailed demo of projective homotopies.

## Installation & Setup

In Macaulay2:
```
restart
needsPackage "EuclideanDistanceDegree"
-- Optional:
-- needsPackage "Bertini"
```

Load examples:
```
load "ExamplesSymbolicEDDegree.m2"
load "ExampleEDDProjectiveHomotopy.m2"
```

## Usage Overview

1. Define a model.
2. Define the ideal of critical points using provided constructors.
3. Compute ED-degree symbolically, or
4. Solve nearest‑point problems numerically via parameter homotopy continuation.



## Comments on how to run Bertini (v1.7) on macOS (Used in the numerical methods)

This comment explains how to bypass the macOS Gatekeeper warning and open **bertini** if you trust the source.

### Download and allow the App from System Settings

1. After downloading Bertini from its homepage https://bertini.nd.edu/download.html
2. Attempt to open **bertini** by double-clicking it.
3. macOS will display a warning:
   > "Apple cannot verify that this app is free from malware."
4. Open **System Settings**.
5. Navigate to **Privacy & Security**.
6. Scroll down until you see the message:
   > "bertini was blocked because it is not from an identified developer."
7. Click **Allow Anyway**.
8. Try opening **bertini** again.
9. A similar warning will appear, but now with an **Open** button.
10. Click **Open** to launch the app.

