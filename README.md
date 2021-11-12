# PrimitiveSolver
PrimitiveSolver is a library that provides a simple interface for constructing and using general equations of state (EOSes) in GRMHD simulations. It adapts the root-finding methods from the [NumTools](https://github.com/jfields7/num-tools) library but no longer requires it as an explicit dependency. The code itself is based on Wolfgang Kastaun's [RePrimAnd](https://github.com/wokast/RePrimAnd) library, particularly its algorithm for recovering the primitive variables of the GRMHD equations (see <https://arxiv.org/abs/2005.01821>). It is intended to interact with [GR-Athena++](https://arxiv.org/abs/2101.08289), but it has no strict dependencies on the code.

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
Real GetEnergy(Real n, Real T, Real *Y); // total fluid energy density
Real GetPressure(Real n, Real T, Real *Y); // fluid pressure
Real GetEntropy(Real n, Real T, Real *Y); // entropy per baryon
Real GetEnthalpy(Real n, Real T, Real *Y); // enthalpy per baryon
Real GetSoundSpeed(Real n, Real T, Real *Y); // speed of sound in the fluid
Real GetSpecificEnergy(Real n, Real T, Real *Y); // specific energy
```

**Note that the specific entropy and enthalpy returned by `EOS` are per baryon, not per mass.** 

Because the Valencia formulation typically used in GRMHD is based on pressure instead of temperature, two utility functions are provided to calculate the temperature:
```c++
Real GetTemperatureFromE(Real n, Real e, Real *Y); // temperature from energy density
Real GetTemperatureFromP(Real n, Real p, Real *Y); // temperature from pressure
```

Additionally, `EOS` also provides functions to access to two member variables:
```c++
Real GetNSpecies(); // Number of particle species
Real GetBaryonMass(); // Mass of a baryon as used by the EOS
```
Both the number of particle fractions and baryon mass are frequently fixed by the specific EOS, so `EOS` does not provide functionality to change them. Specific EOSes may provide this functionality, but it should be neither expected nor required in the general case.

There is one more function provided by `EOS` that is used specifically by the primitive solver:
```c++
Real GetMinimumEnthalpy();
```
This function, which should be a constant-time operation, returns the global minimum of the enthalpy per baryon.

## Using the Primitive Solver
The primitive solver is embedded in a class called `PrimitiveSolver`. Because it must store a pointer to the `EOS` object, it is also a template class:
```c++
template<typename EOSPolicy, typename ErrorPolicy> class PrimitiveSolver;
```

The constructor takes a single argument to an `EOS` pointer, so an instantiation would look something like this:
```c++
PrimitiveSolver<IdealGas, ResetFloor> ps(&eos);
```

Because all the heavy lifting should be done by the `EOS` object itself, `PrimitiveSolver` only has a few member functions:
```c++
bool PrimToCon(Real prim[NPRIM], Real cons[NCONS],
               Real bu[NMAG], Real g3d[NSPMETRIC],
               Real g3u[NSPMETRIC]);
bool ConToPrim(Real prim[NPRIM], Real cons[NCONS],
               Real bu[NMAG], Real g3d[NSPMETRIC],
               Real g3u[NSPMETRIC]);
```
The `PrimToCon()` function converts a set of primitive variables `prim = (rho, v^i, P, T, Y)` with a magnetic field `bu = B^i` into conserved variables `cons = (D, S_i, tau, D_Y)`. The `ConToPrim()` function does the inverse. It takes a set of conserved variables and a magnetic field and inverts them to get the primitive variables. In addition to `prim`, `cons`, and `bu`, both functions also require the spatial metric, `g3d`, and the inverse spatial metric, `g3u`. There are a couple notes here:
1. Note that none of the variables are densitized. Because there are a variety of differing opinions on when to densitize and undensitize, how to interpolate GR quantities to align with MHD quantities, etc., the code simply leaves it to the user to handle these things.
2. `NPRIM` is *not* the same size as `NCONS`. Because `NPRIM` also stores the temperature, it has one extra variable.

Additionally, `PrimitiveSolver` has two functions for accessing constant member variables:
```c++
EOS<EOSPolicy, ErrorPolicy> *const GetEOS() const; // Get a pointer to the EOS used by this PrimitiveSolver.
const int GetNSpecies() const; // Get the number of particle fractions this PrimitiveSolver expects.
```

# Adding a New Equation of State
The `EOS` class implements all equations of state via policy classes. It expects the policy to define the following (protected) methods:
```c++
Real TemperatureFromE(Real n, Real e, Real *Y);  // Temperature from energy density
Real TemperatureFromP(Real n, Real p, Real *Y);  // Temperature from pressure
Real Energy(Real n, Real T, Real *Y);            // Energy density from temperature
Real Pressure(Real n, Real T, Real *Y);          // Pressure from temperature
Real Entropy(Real n, Real T, Real *Y);           // Entropy per baryon  (NOT per mass!)
Real Enthalpy(Real n, Real T, Real *Y);          // Enthalpy per baryon (NOT per mass!)
Real SoundSpeed(Real n, Real T, Real *Y);        // Sound speed from temperature
Real SpecificEnergy(Real n, Real T, Real *Y);    // Specific energy from temperature
Real MinimumEnthalpy();                          // Global minimum enthalpy per baryon (NOT per mass!)
```
It also expects the following member (protected) variables:
```c++
int n_species;  // Number of particle species, should not exceed MAX_SPECIES
Real mb;        // Baryon mass, should be fixed by the EOS
Real max_n;     // Maximum number density, should be fixed by the EOS
Real min_n;     // Minimum number density, should be fixed by the EOS
Real max_T;     // Maximum temperature, should be fixed by the EOS
Real min_T;     // Minimum temperature, should be fixed by the EOS
```

The protected variables are all provided by the `EOSPolicyInterface` class, so the easiest way to make a new EOS policy is to inherit from it:
```c++
class NewEOSPolicy : public EOSPolicyInterface {
  protected:
    /// Implement all the methods above
    /// ...
  public:
    /// Any additional methods specific to the EOS, such as an adiabatic index, should be made available here.
};
```
