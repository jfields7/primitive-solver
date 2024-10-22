# Point Debugger
This is a tool to test inversions for specific conserved states. The debugger can be
compiled from the main directory of PrimitiveSolver using `make build_debugger`, then
run with `./point-debugger -i <input.par>`. A typical input file looks something like
the following:

```
# EOS policy parameters
[EOS]
eos_policy = IdealGas # Can also be PiecewisePolytrope or EOSCompOSE
code_units = Nuclear # Set to GeometricSolar if mb != 1, such as in EOSCompOSE

## IdealGas parameters
gamma = 2

## PiecewisePolytrope parameters
npieces = 2
density1 = 1
density2 = 10
gamma1 = 2.0
gamma2 = 1.8
mb = 1
P0 = 100
rho_min = 0
gamma_thermal = 1.5

## EOSCompOSE parameters
table = SFHo.h5


# Error policy parameters
[Error]
error_policy = ResetFloor # Other option is DoNothing
dfloor = 1e-15   # rho floor in code units
tfloor = 1e-1    # temperature floor
dthreshold = 1.0 # density threshold
vmax   = 0.999   # Maximum velocity
bsqmax = 1e6     # Maximum magnetization
y1floor = 0.0    # atmosphere for first particle species
#y<n>floor = 0.0 # floor for nth particle species, if present

# PrimitiveSolver parameters
[PrimitiveSolver]
iterations = 30    # number of root-solver iterations
tol        = 1e-15 # tolerance for root solver

# Physical variable state
[State]
D   = 1.0
Sx  = 0.0
Sy  = 0.0
Sz  = 0.0
tau = 1.0
Dy1 = 0.5
#Dy<n> = 0.0

Bx  = 0.01
By  = 0.
Bz  = 0. 

gxx = 1.0
gxy = 0.0
gxz = 0.0
gyy = 1.0
gyz = 0.0
gzz = 1.0
```
