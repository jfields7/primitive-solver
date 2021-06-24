//! \file testing.cpp
//  \brief Implementation of the testing framework.

#include <testing.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

void PrintGreen(const std::string& text) {
  std::cout << "\033[1;32m" << text << "\033[0m";
}

void PrintRed(const std::string& text) {
  std::cout << "\033[1;31m" << text << "\033[0m";
}

void PrintBold(const std::string& text) {
  std::cout << "\033[1m" << text << "\033[0m";
}

void PrintError(double expected, double actual) {
  double error = GetError(expected, actual);
  std::cout << "  Expected: " << expected << "\n";
  std::cout << "  Actual: " << actual << "\n";
  std::cout << "  Error: " << error << "\n";
}

double GetError(double expected, double actual) {
  return fabs((actual - expected)/actual);
}

UnitTests::UnitTests(std::string name){
  test_count = 0;
  tests_passed = 0;
  tests_failed = std::vector<std::string>();
  suite_name = name;

  std::stringstream str;
  str << "Starting unit test suite: " << suite_name << "\n\n";
  PrintBold(str.str());
}

void UnitTests::PrintSummary() {
  std::stringstream str;
  str << "Test suite summary: " << suite_name << "\n";
  PrintBold(str.str());
  std::cout << tests_passed << "/" << test_count << " tests passed.\n";
  if (tests_passed == test_count) {
    PrintGreen("All tests passed!\n");
  }
  else {
    PrintBold("Tests failed:\n");
    for (std::string test : tests_failed) {
      std::stringstream ss;
      ss << test << "\n";
      PrintRed(ss.str());
    }
  }
  std::cout << "\n";
}
