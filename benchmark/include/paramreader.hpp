#ifndef PARAM_READER_HPP
#define PARAM_READER_HPP

#include <string>
#include <map>
#include <vector>

/******************************************************************************
 *
 * Class Name: ParamReader
 * Author: Jacob Fields
 * Date Modified: 29 April 2020
 *
 * Description: A very simple class for reading in parameter files. It is
 *              designed to be 100% self-contained and dependent only on the
 *              C++ STL.
 *
 *****************************************************************************/

class ParamReader{
  public:
    /**
     * An enum indicating the outcome of an operation.
     */
    typedef enum ParamResult_t {
      SUCCESS,
      BAD_FILENAME,
      SYNTAX_ERROR,
      BAD_TYPE,
      MULTIPLE_DEFINITIONS,
      UNSECTIONED_PARAMETER,
      INVALID_PARAMETER,
      INVALID_VALUE,
    } ParamResult;
  private:
    /**
     * Data is stored in a map of maps. The first map is the sections contained
     * in the file, and the second map contains all the parameters within each
     * section. Data is stored in std::vectors to facilitate reading arrays.
     */
    std::map<std::string, std::map<std::string,std::string>> data;

    /**
     * Parse a parameter from a string. If it's valid, this string is
     * automatically added to the specified section in the data map. If it's
     * invalid, nothing is added and an error is returned.
     */
    ParamResult parseParameter(std::string &str, std::string &section);

    /**
     * Convert a string into an array. It must contain exactly one '{' and '}'
     * each, which must respectively be the first and last characters.
     */
    ParamResult stringToVector(std::vector<std::string>& vec, std::string &str);
  public:
    /**
     * Constructor
     */
    ParamReader();
    /**
     * Destructor
     */
    ~ParamReader();

    /**
     * Read a file off the disk and load it into memory.
     * @param fname - A string containing the filename
     * @return - A ParamResult indicating success or why it failed.
     */
    ParamResult readFile(std::string fname);

    /**
     * Check if a particular section exists.
     * @param section - The section in question.
     * @return - True if the section was found, false otherwise.
     */
    bool hasSection(std::string section);

    /**
     * Check if a particular parameter exists in a section.
     * @param section - The section containing the parameter.
     * @param parameter - The parameter in question.
     * @return - True if the parameter was found, false otherwise.
     */
    bool hasParameter(std::string section, std::string parameter);

    /**
     * Read a parameter value as a string.
     * @param section - The section containing the parameter.
     * @param parameter - The parameter in question.
     * @return - The value of the parameter if it exists. Otherwise,
     *           it returns the string "NULL".
     */
    std::string readAsString(std::string section, std::string parameter);

    /**
     * Read a parameter value as a double.
     * @param section - The section containing the parameter.
     * @param parameter - The parameter in question.
     * @return - The value of the parameter if it exists. Otherwise,
     *           it returns 0.0.
     */
    double readAsDouble(std::string section, std::string parameter);

    /**
     * Read a parameter value as an integer.
     * @param section - The section containing the parameter.
     * @param parameter - The parameter in question.
     * @return - The value of the parameter if it exists. Otherwise,
     *           it returns 0.
     */
    int readAsInt(std::string section, std::string parameter);

    /**
     * Read a parameter value as an array of strings.
     * @param section - The section containing the parameter.
     * @param parameter - The parameter in question.
     *
     */

    /**
     * Clear all the data that's been read in.
     */
    void clearData();
};
#endif
