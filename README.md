# PrimitiveSolver
PrimitiveSolver is a library that provides a simple interface for constructing and using general equations of state (EOSes) in GRMHD simulations. It requires the [NumTools](https://github.com/jfields7/num-tools) library to work. The code itself is based on Wolfgang Kastaun's [RePrimAnd](https://github.com/wokast/RePrimAnd) library, particularly its algorithm for recovering the primitive variables of the GRMHD equations (see <https://arxiv.org/abs/2005.01821>). It is intended to interact with [GR-Athena++](https://arxiv.org/abs/2101.08289), but its dependencies are limited to the `AthenaArray` data structure and the names of a few data types.

To compile PrimitiveSolver, run `make`.

To install PrimitiveSolver, run `make install`. It may require sudo privileges. The default install location is in `/usr/local/`. To change this, update `INSTALL_DIR` in the Makefile. Header files will be installed in `include/PrimitiveSolver`, and the library `libPrimitiveSolver.a` will be installed in `lib`.

To uninstall PrimitiveSolver, run `make uninstall`, which also may require sudo privileges.

To run unit tests, run `make test`.

# How to use PrimitiveSolver
## Setting up an Equation of State
PrimitiveSolver utilizes a policy-based design. An equation of state consists of two parts: an `EOSPolicy` and an `ErrorPolicy`. The former provides methods for calculating various thermodynamic quantities, and the latter provides methods for handling various errors that might occur during a simulation, such as what to do if the density falls below an artificial atmosphere or the primitive recovery routine fails. These are wrapped up in a single template class, EOS:

```c++
template<typename EOSPolicy, typename ErrorPolicy> class EOS;
```

For example, to create an ideal gas that simply resets unphysical variables to the atmosphere, one would simply write

```c++
EOS<IdealGas, ResetFloor> eos;
```

In order to align with most EOS tables for neutron stars, an `EOS` object revolves around three thermodynamic variables: number density in the rest frame (`n`), temperature (`T`), and particle density fractions (such as the fraction of electrons, quarks, etc.) (`Y`). The rationale for this is that most EOS tables calculate energy density as functions of `n`, `T`, and `Y`, so any other quantity requires an expensive root find. Other thermodynamic quantities can be calculated using the following functions:
```c++
Real GetEnergy(Real n, Real T, Real *Y) // total fluid energy density
Real GetPressure(Real n, Real T, Real *Y) // fluid pressure
Real GetEntropy(Real n, Real T, Real *Y) // entropy per baryon
Real GetEnthalpy(Real n, Real T, Real *Y) // enthalpy per baryon
Real GetSoundSpeed(Real n, Real T, Real *Y) // speed of sound in the fluid
Real GetSpecificEnergy(Real n, Real T, Real *Y) // specific energy
```

**Note that the specific entropy and enthalpy returned by `EOS` are per baryon, not per mass.** 

Because the Valencia formulation typically used in GRMHD is based on pressure instead of temperature, two utility functions are provided to calculate the temperature:
```c++
Real GetTemperatureFromE(Real n, Real e, Real *Y) // temperature from energy density
Real GetTemperatureFromP(Real n, Real p, Real *Y) // temperature from pressure
```

Additionally, `EOS` also provides functions to access to two member variables:
```c++
Real GetNSpecies() // Number of particle species
Real GetBaryonMass() // Mass of a baryon as used by the EOS
```
Both the number of particle fractions and baryon mass are fixed by the EOS and cannot be changed at runtime in any of the provided EOSes (and user-implemented EOSes should generally maintain this policy).
