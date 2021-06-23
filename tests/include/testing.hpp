#ifndef TESTING_H
#define TESTING_H
//! \file testing.h
//  \brief A header file for a basic testing framework.
#include <string>
#include <vector>

/// Utility functions used by the testing framework.
void PrintGreen(const std::string& text);
void PrintRed(const std::string& text);

/// A class for handling unit tests.
class UnitTests {
  private:
    std::string suite_name;
  public:
    int test_count;
    int tests_passed;
    std::vector<std::string> tests_failed;
    UnitTests(std::string name);

    void RunTest(bool (*test)(), std::string name);

    void PrintSummary();
};

#endif
