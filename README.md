# FVM-High-Older-DOLINC-Scheme

This repo provides code for our paper [Order-lifted data inversion/retrieval method of neighbor cells to implement general high-order schemes in unstructured-mesh-based finite-volume solution framework](https://doi.org/10.1016/j.jcp.2025.113735) ([arXiv version](https://arxiv.org/abs/2404.08952)), implemented in the OpenFOAM framework.
* Authors: Hao Guo \[[Google Scholar](https://scholar.google.com/citations?hl=zh-CN&user=ZOhz0b8AAAAJ)\], Boxing Hu, Peixue Jiang, Xiaofeng Ma, Yinhai Zhu

<p align="center">
  <img align="center" width="400" src="/docs/RTI4096.png">
</p>
<p align="center" > Density contours of Rayleigh–Taylor instability problem using high-order DOLINC scheme. </p>

<p align="center">
  <img align="center" width="400" src="/docs/R2D4096.png">
</p>
<p align="center"> Density contours of two-dimensional Riemann problem using high-order DOLINC scheme. </p>


## Abstract:

This study introduces an order-lifted inversion/retrieval method for implementing high-order schemes within the framework of an unstructured-mesh-based finite-volume method. This method defines a special representation called the data order-lifted inversion of neighbor cells (DOLINC) differential, which transforms the degrees of freedom of wide templates into differentials of various orders stored in local grid cells. Furthermore, to retrieve the original far-field information without bias during the reconstruction/interpolation of face values, the corresponding accurate inversion formulas are derived based on the defined DOLINC differentials. The order-lifted inversion method can be applied to multi-dimensional polyhedral-mesh solvers by considering the influence of grid non-uniformity on high-order schemes. It seamlessly accommodates multi-process parallel computing for high-order methods without requiring special consideration for the boundary interface. This method not only enhances the numerical accuracy of second-order finite-volume methods, but also demonstrates a significant computational-speed advantage over similar methods. A series of benchmark cases, including the linear advection, Burgers, and Euler equations, are comprehensively validated to assess the practical performance of the method. The results indicate that the unstructured-mesh high-order schemes implemented based on this method achieve theoretical accuracy in practical computations and substantially reduce computational costs compared with methods that increase grid resolution.


## Code structure:

* `src/` contains the source code of several DOLINC schemes and a related solver.
  * `DOLINCSchemes/` contains the code of all DOLINCSchemes to compile the basic library `libDOLINCSchemes.so`.
    * `schemes/` contains the code of fixed stencil schemes.
    * `doliSchemes/` contains the code of non-oscillatory  schemes and related abstract classes.
  * `EulerFoam/` contains the code of a density-based compressible-flow solver with high-order time integration.
* `examples/` contains a set of examples of using the high-order schemes.
  * `Riemann2D/` is a classical version of two-dimensional Riemann Problems.
  * `RTInstability/` is the benchmark of Rayleigh–Taylor instability.


## Compilation:

Run `Allwmake` to compile both the library of DOLINC schemes and the Euler equation solver.


## Usage:

For main DOLINC library `libDOLINCSchemes.so`, below schemes are available.

* 3rd-order upwind stencil scheme (USR3): superior alternative to `linearUpwind` scheme.
* 4th-order central stencil scheme (CSR4): superior alternative to `linear` scheme.
* ENO schemes based on 2nd-order (UENO2, first-order DOLINC differentials) and 3rd-order stencils (UENO3, second-order DOLINC differentials): superior alternative to TVD-type schemes.

The solver EulerFoam is similar to rhoCentrolFoam but with additional functions. Two types of time scheme is integrated:
* 1st-order forward Euler scheme (ForwardEuler): same as rhoCentralFoam.
* 3rd-order accurate TVD-type Runge–Kutta schemes (TVDRungeKutta3): better accuracy and robustness.

See files and scripts in example cases for more instruction.


### Guide for your own case:

Add below line to `system/controlDict` for linking to `.so`:

```cpp
libs            (DOLINCSchemes);
```

If using a density-based solver, add below lines to `system/fvSchemes`:

```cpp
// Only for demonstration
interpolationSchemes
{
    scheme          UENO3 grad(U) 1.2;  // 1.0 of original ENO scheme
    schemeV         UENO3 grad(U) 1.2;

    reconstruct(rho) $scheme;
    reconstruct(U)  $schemeV;
    reconstruct(T)  $scheme;
}
```

If using the EulerFoam, the time advance scheme should also be set:

```cpp
advanceScheme       ForwardEuler;  // ForwardEuler TVDRungeKutta3
```

If using a pressure-based solver, add below lines to `system/fvSchemes`:

```cpp
// Only for demonstration
divSchemes
{
    schemeV         UENO3 grad(U) 1.2;

    div(phi,U)      Gauss $schemeV;
    div(phi,T)      Gauss USR3 grad(T);
    div(phi,e)      Gauss CSR4 grad(T);
}
```


## Examples:

So far, we have updated two examples as demonstrations of the use of high-order DOLINC schemes for density-based solvers and pressure based solvers.

### 2D Riemann Problem

This is a classical case for high-order scheme resolution test.

Adjust the settings in `Allrun` script and run the script.

```bash
cd examples/Riemann2D
./Allclean ; ./Allrun >& log &
```

### Rayleigh--Taylor instability

This case is performed for demonstration of pressure-based solvers.

Adjust the settings in `Allrun` script and run the script.

```bash
cd examples/RTInstability
./Allclean ; ./Allrun >& log &
```

## Requirements:

OpenFOAM v2106


## Citation

If you find this repo useful in your research, please consider citing our paper: [Order-lifted data inversion/retrieval method of neighbor cells to implement general high-order schemes in unstructured-mesh-based finite-volume solution framework](https://doi.org/10.1016/j.jcp.2025.113735).

```
@article{Guo2025,
    title={Order-lifted data inversion/retrieval method of neighbor cells to implement general high-order schemes in unstructured-mesh-based finite-volume solution framework},
    author={Guo, Hao and Hu, Boxing and Jiang, Peixue and Ma, Xiaofeng and Zhu, Yinhai},
    year={2025},
    journal={Journal of Computational Physics},
    volume={524},
    pages={113735},
    DOI={10.1016/j.jcp.2025.113735},
}
```


## License

**FVM-High-Older-DOLINC-Scheme** is published under the **GNU GPL Version 3** license, which can be found in the LICENSE file.


## Problems

If you find any bugs in the code or have trouble in compiling/running **FVM-High-Older-DOLINC-Scheme** in your machine, you are very welcome to [open an issue](https://github.com/Fracturist/FVM-High-Older-DOLINC-Scheme/issues) in this repository.


## Future plan

More schemes and examples/tutorials are in progress.

