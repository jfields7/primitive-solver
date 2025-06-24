#include <paramreader.hpp>
#include <fstream>
#include <iostream>
#include <ctype.h>

using namespace std;

// ParamReader {{{
ParamReader::ParamReader(){

}
// }}}

// ~ParamReader{{{
ParamReader::~ParamReader(){

}
// }}}

// readFile {{{
ParamReader::ParamResult ParamReader::readFile(string fname){
  // Try to load the file.
  ifstream file;
  try{
    file.open(fname.c_str(),ifstream::in);
  }
  catch(ifstream::failure &e){
    return BAD_FILENAME;
  }
  if(!file.is_open()){
    return BAD_FILENAME;
  }

  string section;

  // We need to loop until we hit the end of the file.
  while(!file.eof()){
    // Try to read a line. Ignore initial whitespace.
    char line[256];
    while(file.peek() == ' '){
      file.get();
    }
    // Grab the rest of the line.
    file.getline(line, 256);
    // Strings are easier to work with, so we convert the line to a string.
    string str = string(line);
    if(str.length() == 0){
      continue;
    }

    // Treat # and ; as comments.
    size_t pos = str.find('#');
    if(pos != string::npos){
      str = str.substr(0,pos);
    }
    pos = str.find(';');
    if(pos != string::npos){
      str = str.substr(0,pos);
    }
    if(str.length() == 0){
      continue;
    }

    // Chop whitespace off the end.
    while(str.length() > 0 && str.back() == ' '){
      str = str.substr(0,str.length()-1);
    }
    if(str.length() == 0){
      continue;
    }

    // Now we can check if this is a section.
    if(str.front() == '['){
      if(str.back() != ']'){
        return SYNTAX_ERROR;
      }
      // Extract only the name.
      str = str.substr(1,str.length()-2);
      // Make sure this section doesn't already exist. If so, return an error.
      if(data.find(str) != data.end()){
        return MULTIPLE_DEFINITIONS;
      }
      // Make a new map for the section.
      data[str]=map<string, string>();
      section = str;

    }
    else{
      // If it's not in a section, there's a naked variable, which is an error.
      if(data.size() == 0){
        return UNSECTIONED_PARAMETER;
      }

      ParamResult result = parseParameter(str, section);
      if(result != SUCCESS){
        return result;
      }
    }
  }

  file.close();

  return SUCCESS;
}
// }}}

// parseParameter {{{
ParamReader::ParamResult ParamReader::parseParameter(std::string &str, std::string &section){
  // There shouldn't be any whitespace at the beginning or the end. We can
  // go ahead and just look for where the equals sign is.
  size_t pos = str.find('=');
  if(pos == string::npos || pos == 0 || pos == str.length()-1){
    return SYNTAX_ERROR;
  }

  // Subdivide the string into a parameter and a value.
  string parameter = str.substr(0,pos);
  string value = str.substr(pos+1,string::npos);

  // Trim whitespace off the end of the parameter and off the beginning of the value.
  while(parameter.back() == ' '){
    parameter.pop_back();
  }
  while(value.front() == ' '){
    value.erase(0,1);
  }

  // Make sure the only characters left in the parameter are alphanumeric or have underscores.
  for(size_t i = 0; i < parameter.length(); i++){
    if(!isalnum(parameter[i]) && (parameter[i]!='_')){
      return INVALID_PARAMETER;
    }
  }
  // Make sure the only characters left in the value are also alphanumeric. Allow a '{' at the front
  // and a '}' at the back for arrays. Other permitted characters include '_', '.', '-', ',', and 
  // whitespace.
  for(size_t i = 0; i < value.length(); i++){
    if(!isalnum(value[i]) && (value[i]!='_') && (value[i] != '.') && (value[i] != '-') && (value[i] != ' ')
       && (value[i] != ',')){
      if((i != 0 && value[i]!='{') || (i!=value.length()-1 && value[i]!='}')){
        return INVALID_VALUE;
      }
    }
  }

  // If we've made it this far, we can add the parameter and its value to the 
  // data map. Multiple definitions are not an issue; if the parameter already
  // exists, its last value is overwritten.
  data[section][parameter] = value;

  return SUCCESS;
}
// }}}

// stringToVector {{{
ParamReader::ParamResult ParamReader::stringToVector(std::vector<std::string>& vec, std::string &str){
  // Make sure we have a legal string first.
  if(str.size() == 0){
    return SYNTAX_ERROR;
  }
  // Check that the first character is a '{' character. If so,
  // pop it off.
  if(str.front() != '{'){
    return SYNTAX_ERROR;
  }
  str.erase(0,1);
  // Check that the last character is a '}' character. If so,
  // pop it off.
  if(str.back() != '}'){
    return SYNTAX_ERROR;
  }
  str.pop_back();
  // Sort through all the arguments.
  while(str.size() > 0){
    size_t pos = str.find(',');
    std::string arg = str.substr(0, pos);
    str.erase(0, pos);
    str.erase(0, 1);
    // We should only have legal strings entering at this point because of the checks that occur
    // in parseParameter.
    vec.push_back(arg);
  }

  return SUCCESS;
}
// }}}

// hasSection {{{
bool ParamReader::hasSection(string section){
  return (data.find(section) != data.end());
}
// }}}

// hasParameter {{{
bool ParamReader::hasParameter(string section, string parameter){
  if(!hasSection(section)){
    return false;
  }
  return (data[section].find(parameter) != data[section].end());
}
// }}}

// readAsString {{{
string ParamReader::readAsString(string section, string parameter){
  if(!hasParameter(section,parameter))
  {
    cout << "Warning! " << section << " : " << parameter << " does not exist!\n";
    return string("");
  }
  else{
    return data[section][parameter];
  }
}
// }}}

// readAsDouble {{{
double ParamReader::readAsDouble(string section, string parameter){
  if(!hasParameter(section,parameter)){
    cout << "Warning! " << section << " : " << parameter << " does not exist!\n";
    return 0.0;
  }
  double result;
  try{
    result = stod(data[section][parameter]);
  }
  catch(...){
    cout << "Warning! " << section << " : " << parameter << " could not be read as a double!\n";
    return 0.0;
  }

  return result;
}
// }}}

// readAsInt {{{
int ParamReader::readAsInt(string section, string parameter){
  if(!hasParameter(section,parameter)){
    cout << "Warning! " << section << " : " << parameter << " does not exist!\n";
    return 0;
  }
  int result;
  try{
    result = stoi(data[section][parameter]);
  }
  catch(...){
    cout << "Warning! " << section << " : " << parameter << " could not be read as an integer!\n";
    return 0;
  }

  return result;
}
// }}}

// clearData {{{
void ParamReader::clearData(){
  data.clear();
}
// }}}
