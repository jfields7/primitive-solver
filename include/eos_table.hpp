//! \file eos_table.hpp
//  \brief Defines EOSTable, which stores information from a tabulated
//         equation of state.

#include <string>

#include <eos_units.hpp>

class EOSTable {
  public:
    enum TableVariables {
      ETN = 0,    //! number density [cgs]
      ETT = 1,    //! temperature [cgs]
      ETNVARS = 2
    };

    EOSTable(std::string fname);
};
