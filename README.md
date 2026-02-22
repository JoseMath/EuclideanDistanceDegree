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


# How to run Bertini (v1.7) on macOS (Sequoia 15.6.1)

The third-party software Bertini is called in our numerical continuation methods. This comment explains how to bypass the macOS Gatekeeper warning and open **bertini** if you trust the source.

## Download and allow the App from System Settings

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


## Troubleshooting: unable to run Bertini in the terminal

Bertini is usually run through the terminal and relies on input text files. Our Macaulay2 package, when using numerical methods, streamlines the creation of the input files and automates calls to Bertini. In this section we explain some our personal best practices and other ways to troubleshoot common issues.  

### 1. If macOS still blocks Bertini

You may need to remove Apple's quarantine flag:

```bash
xattr -rd com.apple.quarantine /Applications/BertiniApple_v1.7
```

Then attempt to run it again.

---
### 2. Fixing `zsh: command not found: bertini` on macOS

This guide explains how to resolve the error:

```
zsh: command not found: bertini
```

This happens when macOS cannot find the **Bertini** executable in your system `PATH`.

---

#### 2a. Find where Bertini is installed

You must locate the `bertini` executable. Common locations include:

- `/Applications/BertiniApple_v1.7/bertini`
- `/usr/local/bin/bertini`
- `~/Applications/BertiniApple_v1.7/bertini`
- A folder where you unzipped it (often `~/Downloads`)

Try searching with:

```bash
find /Applications -name bertini 2>/dev/null
```

or:

```bash
find ~ -name bertini 2>/dev/null
```

---

#### 2b. Test running Bertini using the full path

Once you find it, try running it directly. Example:

```bash
/Applications/BertiniApple_v1.7/bertini 
```

or 

```bash
/Applications/BertiniApple_v1.7/bertini /Applications/BertiniApple_v1.7/examples/zero_dim/basic_example/input 
```

If this works, Bertini is installed correctly; your shell just doesn't know where it is yet. The former will likely return something like this 
```
ERROR: 'input' does not exist!!!

Bertini will now exit due to this error.
```
because the input file does not exist in the same directory as Bertini.

---

#### 2c. Add Bertini to your PATH

To run Bertini without specifying the path, you must add the directory containing `bertini` to your `PATH`.

Example path:

```
/Applications/BertiniApple_v1.7
```

Edit your shell configuration:

```bash
open -e ~/.zshrc
```

Add this line at the bottom (update folder as needed):

```bash
export PATH="/Applications/BertiniApple_v1.7:$PATH"
```

Save and reload configuration:

```bash
source ~/.zshrc
```

Test:

```bash
bertini
```

---
## How to tell Macaulay2 how to find Bertini?

There is a distinction between path and PATH that can sometimes lead to confusion and cause errors when Macaulay2 is trying to find Bertini. 

---

### What Is the System PATH?

The **PATH** is an environment variable used by the operating system's shell (e.g., `zsh`, `bash`) to determine where to look for executable programs.

When you run a command such as:

```bash
bertini input
```

the shell searches for `bertini` in each directory listed in your PATH.

For example, 
```
/usr/local/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin
```

In Macaulay2 you get see the PATH by using the command

```
getenv "PATH"
```

If `bertini` is not in the PATH seen by Macaulay2, you get errors like:
  ```
  zsh: command not found: bertini
  ```

---

### path versus PATH 
Macaulay2 has a separate mechanism for locating its own **.m2 packages**. These are files written in the Macaulay2 language.
This path list is *not* used to locate external programs. The path for packages usually includes directories like:

```
./
Library/Application Support/Macaulay2/code/
Library/Application Support/Macaulay2/local/share/Macaulay2/
/opt/homebrew/share/Macaulay2/
```

These directories contain:
- Package source files (`.m2`)
- Documentation metadata
- Compiled Macaulay2 libraries

But, to run Bertini we need to search for system executables. 

### Specifying the correct path to Bertini
You can make `bertini` discoverable by Macaulay2 using the Configuration option when loading the package.

```
needsPackage("Bertini", Configuration=>{"BERTINIexecutable"=>"/Applications/BertiniApple_v1.7/bertini"})
```
or putting it in a directory in the system PATH.


---




