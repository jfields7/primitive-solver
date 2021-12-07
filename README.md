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

# Unit Tests
It is strongly recommended that any new `EOSPolicy` or `ErrorPolicy` classes have unit tests to go with them. There is a basic unit testing framework included in the repository under the `tests` folder. There are folders for three kinds of tests: `error` includes tests for the `ErrorPolicy`, `eos` contains tests for a specific `EOSPolicy`, and `primitive` contains tests for running `PrimitiveSolver` with a specific `EOSPolicy`.

## Running and Adding New Unit Tests
Tests can be run from the repository root by running `make test`. This runs a separate Makefile inside `tests` that compiles each test and runs it. Inside `tests\Makefile`, each set of tests has its own executable (again, the framework is very basic). Suppose we wrote a new `EOSPolicy` for the Synge relativistic ideal gas. Then we could add the following recipe:
```make
eos_synge : $(OBJ_FILES)
	$(call make_test,eos,test_synge)
```
The `make_test` function takes two arguments: a directory (relative to `tests`) and a source file. In this case, we're saying that we have a test at `./eos/test_synge.cpp`. From there, `make_test` will automatically compile the test with all the correct libraries and header files.

The next step is to make sure that the recipe is actually applied correctly. So, we need to append the test to the `tests` recipe inside the Makefile:
```make
tests: dirs ... eos_synge
```

Lastly, we need to make sure that the test is executed when we type `make test` in the repository root. There is one more recipe, `run`, in the test Makefile for this. Append the new test to run as follows, making sure to add the trailing `; \` to the line before:
```make
run: dirs tests
	...
	...; \
	./test_synge
```

Now when we compile with `make test`, the new tests will automatically run.

## Unit Testing Framework
The unit testing framework can be included in test source code via the `testing.hpp` header file. It provides a few auxiliary functions:
```c++
void PrintGreen(const std::string& text);
void PrintRed(const std::string& text);
void PrintBold(const std::string& text);
void PrintError(double expected, double actual);
double GetError(double expected, double actual);
```
The functions do exactly what their names suggest: `PrintGreen()` prints a string of text to the terminal in bright green, `PrintRed()` prints in red, and `PrintBold()` prints bold text. `PrintError()` prints out what the expected value was, what the actual value is, and what the relative error between the two is. Finally, `GetError()` performs a relative error calculation.

The heart of the unit testing framework is the `UnitTests` class. There are three functions of interest:
```c++
UnitTests(const std::string& name);

template<class ... Types>
void RunTest(bool (*test)(Types...), std::string name, Types ... args);

void PrintSummary();
```

The `UnitTests()` constructor takes a string which will be used as the name for this test suite. `PrintSummary()` should be called at the end in order to tell the user how many tests passed and which tests, if any, failed.

The `RunTest()` function takes a little bit of explanation. This is a variadic template function, which makes it possible for `UnitTests` to run any test without needing to know the number of arguments beforehand. The function accepts the following parameters: a pointer to a test function with an arbitrary number of arguments that returns a `bool`, a name for the test, and then any arguments needed by the test function. Suppose we want to run the test `bool TestTemperatureFromEnergy(Real n, Real e, Real *Y)` from `testfunctions.hpp`. Then we would type the following:
```c++
tester.RunTest(&TestTemperatureFromEnergy, "Temperature from Energy Test", n, e, Y);
```
We have three arguments after the name of the test because `TestTemperatureFromEnergy()` needs three arguments.

Alternatively, if we built a test that takes no arguments but confirmed that an object's constructor works as expected, we would do the following instead:
```c++
tester.RunTest(&TestConstructor, "Construction Test");
```
Note that because `TestConstructor()` takes no arguments, we call `RunTest()` with only two arguments.

There is one noticeable limitation on test functions constructed this way: variadic templates don't work very well with references, so all parameters should be passed by value or by pointer.

## Writing Unit Tests
Good unit tests should test every non-trivial component of the code. For an `EOSPolicy`, this roughly translates to every EOS function, such as `TemperatureFromE`, `TemperatureFromP`, and so forth. Good tests will also test a variety of inputs, including both typical and extreme inputs, boundary cases, and outright wrong inputs. For example, the `PiecewisePolytrope` tests include densities at every specified polytrope, densities for an ideal gas component below the first polytrope, and continuity tests.

In order to make this easier, there are some utility functions provided in `testfunctions.hpp` and `primitive_utility.hpp`:
```c++
// testfunctions.hpp
bool TestTemperatureFromEnergy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos,
  Real n, Real T, Real *Y, const Real tol);

bool TestTemperatureFromPressure(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos,
  Real n, Real T, Real *Y, const Real tol);

bool TestEnthalpy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos,
  Real n, Real T, Real *Y, const Real tol);

bool TestSpecificEnergy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos,
  Real n, Real T, Real *Y, const Real tol);

bool TestConToPrim(Primitive::PrimitiveSolver<EOSPolicy, ErrorPolicy>* ps, 
  Real prim[NPRIM], Real cons[NCONS], Real bu[NMAG], 
  Real gd[NMETRIC], Real gu[NMETRIC], const Real tol);
```

```c++
// primitive_utility.hpp

// Functions for modifying the metric.
void MinkowskiMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]);

void SchwarzschildMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]);

void ScrewballMinkowskiMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]);

// Functions for modifying the velocity.
void ZeroVelocity(Real prim[NPRIM]);

void StrongVelocity(Real prim[NPRIM]);

// Functions for modifying the magnetic field.
void ZeroField(Real bu[NMAG]);

void StrongField(Real bu[NMAG]);

// Functions for non-zero particle fractions.
void ParticleFractions(Real prim[NPRIM], int s);
```

Generally speaking, all unit tests should consist of a single source file. Additional source files, if needed, can be put in `tests/src/` and `tests/include/`. Supposing we had an `EOSPolicy` called `SyngeGas`, we could write some tests for it as follows:
```c++
#include <eos.hpp>
#include <ps_types.hpp>
#include <synge_gas.hpp>
#include <do_nothing.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

using namespace Primitive;

int main(int argc, char* argv[]) {
  UnitTests tester{"Synge Gas EOS"};

  EOS<SyngeGas, DoNothing> eos;
  const Real tol = 1e-12; // error tolerance for floating-point arithmetic.
  Real n = 1.345e2; // number density
  Real T = 4.985e3; // temperature
  Real *Y = nullptr; // empty, SyngeGas doesn't depend on particle fractions.

  // Check that the temperature is self-consistent with the energy.
  tester.RunTest(&TestTemperatureFromEnergy<SyngeGas, DoNothing>,
                 "Temperature from Energy Test",
                 &eos, n, T, Y, tol);

  // Check that the temperature is self-consistent with the pressure.
  tester.RunTest(&TestTemperatureFromPressure<SyngeGas, DoNothing>,
                 "Temperature from Pressure Test",
                 &eos, n, T, Y, tol);

  tester.PrintSummary();

  return 0;
}

```

*Note*: It is _strongly_ recommended that any `EOSPolicy` class have tests both for the EOS itself and for its interactions with `PrimitiveSolver`.
