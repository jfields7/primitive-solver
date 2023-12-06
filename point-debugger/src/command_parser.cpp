//! \file command_parser.cpp
//  \brief Implementation file for CommandParser

#include <command_parser.hpp>
#include <stdexcept>
#include <sstream>
#include <locale>

void CommandParser::AddArgument(std::string name, CommandType type,
                                bool required, std::string def_val,
                                bool has_short) {
  // Make sure the argument is alphanumeric.
  for (unsigned int i = 0; i < name.length(); i++) {
    if (!std::isalnum(name[i])) {
      std::stringstream ss;
      ss << "Argument " << name << " is not alphanumeric.";
      throw std::invalid_argument(ss.str());
    }
  }

  // If the argument already exists, throw an exception.
  if (arguments.find(name) != arguments.end()) {
    std::stringstream ss;
    ss << "Argument " << name << " already exists.";
    throw std::invalid_argument(ss.str());
  }

  // If the argument has a short form and there is already another
  // short form variable, we throw an exception.
  if (has_short) {
    std::string long_form = FindLongForm(name[0]);
    if (long_form.compare("") != 0) {
      std::stringstream ss;
      ss << "Argument " << name << " shares the same short form as argument "
         << long_form << ".";
      throw std::invalid_argument(ss.str());
    }
  }

  // If the argument does not exist, we can add it to the argument map.
  Argument arg{type, !required, def_val, name, has_short};
  arguments.emplace(name, arg);
}

void CommandParser::AddString(std::string name, bool required, std::string def_val,
    bool has_short) {
  AddArgument(name, CommandType::STR, required, def_val, has_short);
}

void CommandParser::AddReal(std::string name, bool required, Real def_val,
    bool has_short) {
  AddArgument(name, CommandType::REAL, required, std::to_string(def_val), has_short);
}

void CommandParser::AddInteger(std::string name, bool required, int def_val,
    bool has_short) {
  AddArgument(name, CommandType::INT, required, std::to_string(def_val), has_short);
}

void CommandParser::AddBoolean(std::string name, bool has_short) {
  // Note that all booleans have a default value of false.
  AddArgument(name, CommandType::BOOL, false, std::string("f"), has_short);
}

std::string CommandParser::GetTypeName(CommandParser::CommandType type) {
  switch(type) {
    case CommandType::INT:
      return "integer";
      break;
    case CommandType::REAL:
      return "real";
      break;
    case CommandType::STR:
      return "string";
      break;
    case CommandType::BOOL:
      return "boolean";
      break;
  }

  return "UNDEFINED";
}

std::string CommandParser::GetValue(std::string name, CommandParser::CommandType type) {
  // If the argument does not exist, throw an exception.
  if (arguments.find(name) == arguments.end()) {
    std::stringstream ss;
    ss << "Argument " << name << " does not exist.";
    throw std::invalid_argument(ss.str());
  }

  Argument &arg = arguments[name];

  if (arg.type != type) {
    std::stringstream ss;
    ss << "Argument " << name << " is not of type " << GetTypeName(type) << ".";
    throw std::invalid_argument(ss.str());
  }

  if (arg.set == false) {
    std::stringstream ss;
    ss << "Argument " << name << " has not been initialized yet!";
    throw std::domain_error(ss.str());
  }
  
  return arg.val;
}

std::string CommandParser::GetString(std::string name) {
  return GetValue(name, CommandType::STR);
}

bool CommandParser::GetBoolean(std::string name) {
  std::string val = GetValue(name, CommandType::BOOL);
  if (val.compare("f") == 0) {
    return false;
  }
  if (val.compare("t") == 0) {
    return true;
  }

  return false;
}

int CommandParser::GetInteger(std::string name) {
  std::string val = GetValue(name, CommandType::INT);
  return std::stoi(val);
}

Real CommandParser::GetReal(std::string name) {
  std::string val = GetValue(name, CommandType::REAL);
  return std::stod(val);
}

std::string CommandParser::FindLongForm(char arg) {
  for (auto it : arguments) {
    if (arg == it.second.name[0]) {
      if (it.second.has_short) {
        return it.second.name;
      }
    }
  }

  return "";
}

std::string CommandParser::FormatArgument(char arg[]) {
  std::string argument{arg};

  // If the argument is shorter than two characters,
  // it can't possibly be a real argument.
  int length = argument.length();
  if (length < 2) {
    std::stringstream ss;
    ss << "Argument " << argument << " is too short.";
    throw std::invalid_argument(ss.str());
  }

  // If the first two characters are hyphens, the
  // argument must be in long form, and we can go ahead
  // and chop off the hyphens and return the argument.
  if (argument.compare(0, 2, "--") == 0) {
    return argument.substr(2, std::string::npos);
  }

  // If we only have one hyphen, then we need to search
  // for a long-form equivalent.
  if (argument.compare(0, 1, "-") == 0) {
    if (argument.length() != 2) {
      std::stringstream ss;
      ss << "Argument " << argument << " is too long to be short-form.";
      throw std::invalid_argument(ss.str());
    }

    std::string long_form = FindLongForm(argument[1]);
    if (long_form.compare("") == 0) {
      std::stringstream ss;
      ss << "Could not expand argument " << argument << ".";
      throw std::invalid_argument(ss.str());
    }

    return long_form;
  }

  std::stringstream ss;
  ss << "Invalid argument format: " << argument << ".";
  throw std::invalid_argument(ss.str());

  return "";
}

CommandParser::Error CommandParser::ParseArguments(int argc, char* argv[]) {
  // Loop over the arguments, skipping the first,
  // since it's just the name of the program.
  for (int i = 1; i < argc; i++) {
    std::string arg;
    // Try to format the argument. If that fails,
    // we return a format error.
    try {
      arg = FormatArgument(argv[i]);
    }
    catch (std::invalid_argument& e){
      Error error{ErrorCode::FORMAT, i, argv[i]};
      return error;
    }
    
    // Try to find the argument. If that fails,
    // return an unknown argument error.
    if (arguments.find(arg) == arguments.end()) {
      Error error{ErrorCode::UNKNOWN, i, arg};
      return error;
    }

    // Now we can read in the argument.
    Argument& arg_obj = arguments[arg];

    // If the argument is a boolean, then we just need
    // to set it and we're done.
    if (arg_obj.type == CommandType::BOOL) {
      arg_obj.set = true;
      arg_obj.val = "t";
    }
    else {
      // Otherwise, we need to read in the next argument,
      // too.
      if ( i == argc-1) {
        // This is the case that there is no following argument.
        Error error{ErrorCode::NOPARAM, i, arg};
        return error;
      }
      else {
        // Load the next argument as the parameter for this
        // argument.
        std::string parameter = std::string(argv[i+1]);
        // Make sure that it's a valid argument, of course.
        try {
          switch(arg_obj.type) {
            case CommandType::INT: {
              // We try to convert to an integer and check if
              // it throws an exception.
              std::stoi(parameter);

              break;
            }
            case CommandType::STR:
              break;
            case CommandType::REAL: {
              // We try to convert to a Real and check if
              // it throws an exception.
              std::stod(parameter);
              break;
            }
            default:
              break;
          }
        }
        catch (std::invalid_argument& e) {
          Error error{ErrorCode::BADPARAM, i, arg};
          return error;
        }
        // If we've made it this far, the parameter should be good.
        arg_obj.set = true;
        arg_obj.val = parameter;

        // We need to increment the index manually to skip the parameter.
        i++;
      }
    }
  }

  // Check that all required arguments have been set.
  for (auto it : arguments) {
    if (it.second.set == false) {
      Error error{ErrorCode::MISSING, 0, it.second.name};
      return error;
    }
  }

  Error error{ErrorCode::SUCCESS, 0, ""};
  return error;
}

std::string CommandParser::GetErrorCodeName(ErrorCode code) {
  switch(code) {
    case ErrorCode::SUCCESS:
      return std::string("Success");
      break;
    case ErrorCode::UNKNOWN:
      return std::string("Unknown Argument");
      break;
    case ErrorCode::FORMAT:
      return std::string("Format Error");
      break;
    case ErrorCode::MISSING:
      return std::string("Missing Argument");
      break;
    case ErrorCode::NOPARAM:
      return std::string("Missing Parameter");
      break;
    case ErrorCode::BADPARAM:
      return std::string("Invalid Parameter");
      break;
  }

  return std::string("Missing Code");
}
