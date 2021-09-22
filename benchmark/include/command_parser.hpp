#ifndef COMMAND_PARSER_HPP
#define COMMAND_PARSER_HPP
//! \file command_parser.hpp
//  \brief An OS-independent parser for command line arguments.

#include <string>
#include <vector>
#include <map>
#include <utility>

#include <ps_types.hpp>

class CommandParser {
  public:
    enum class CommandType {
      INT,
      REAL,
      STR,
      BOOL
    };

    enum class ErrorCode {
      SUCCESS,
      UNKNOWN,
      FORMAT,
      MISSING,
      NOPARAM,
      BADPARAM
    };

    struct Error {
      ErrorCode code;
      int position;
      std::string arg;
    };

  private:
    struct Argument {
      CommandType type; // Type for the argument
      bool set; // If the argument has been initialized or not.
      std::string val; // The current value for this argument.
      std::string name; // Argument name
      bool has_short; // If the argument has a short form or not.
    };

    std::map<std::string, Argument> arguments;

    // Helper function for adding arguments to the map.
    void AddArgument(std::string name, CommandType type, bool required, 
                     std::string def_val, bool has_short);

    std::string GetValue(std::string name, CommandType type);

    std::string GetTypeName(CommandType type);

    /// Chop off the hyphens, extending the argument if a short
    /// form is provided.
    std::string FormatArgument(char argv[]);

    // Find the long form for a variable from the short form.
    std::string FindLongForm(char arg);

  public:
    CommandParser() = default;
    ~CommandParser() = default;

    // Add an argument to the argument map with the indicated type.
    void AddString(std::string name, bool required, std::string def_val = "", bool has_short = false);
    void AddReal(std::string name, bool required, Real def_val = 0.0, bool has_short = false);
    void AddInteger(std::string name, bool required, int def_val = 0, bool has_short = false);
    void AddBoolean(std::string name, bool has_short = true);

    // Get the current value from the argument map for the specified
    // argument.
    std::string GetString(std::string name);
    bool GetBoolean(std::string name);
    Real GetReal(std::string name);
    int GetInteger(std::string name);

    // Taking the console input, parse all the arguments and update
    // the argument map accordingly.
    Error ParseArguments(int argc, char* argv[]);

    // Get the error name
    static std::string GetErrorCodeName(ErrorCode code);
};

#endif
