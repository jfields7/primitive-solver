# PrimitiveSolver
A library describing a general EOS solver for GRMHD. It requires the (NumTools)[https://github.com/jfields7/num-tools] library to work.

To compile PrimitiveSolver, run `make`.

To install PrimitiveSolver, run `make install`. It may require sudo privileges. The default install location is in `/usr/local/`. To change this, update `INSTALL_DIR` in the Makefile. Header files will be installed in `include/PrimitiveSolver`, and the library `libPrimitiveSolver.a` will be installed in `lib`.

To uninstall PrimitiveSolver, run `make uninstall`, which also may require sudo privileges.

To run unit tests, run `make test`.
