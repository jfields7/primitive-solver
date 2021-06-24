#ifndef TESTING_H
#define TESTING_H
//! \file testing.h
//  \brief A header file for a basic testing framework.
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <functional>

/// Printing functions used by the testing framework.
void PrintGreen(const std::string& text);
void PrintRed(const std::string& text);
void PrintBold(const std::string& text);
void PrintError(double expected, double actual);
double GetError(double expected, double actual);


/// A class for handling unit tests.
class UnitTests {
  private:
    std::string suite_name;
  public:
    int test_count;
    int tests_passed;
    std::vector<std::string> tests_failed;
    UnitTests(std::string name);

    template<class T, class ... Types>
    //void RunTest(bool (*test)(Types...), std::string name, Types ... args) {
    void RunTest(std::function<bool(T, Types...)> test, std::string name, T eos, Types ... args) {
      test_count++;
      std::cout << "Running test: " << name << "\n";
      std::stringstream str;
      if (test(eos, args...)) {
        tests_passed++;
        str << name << " passed!\n";
        PrintGreen(str.str());
      }
      else {
        tests_failed.push_back(name);
        str << name << " failed!\n";
        PrintRed(str.str());
      }
      std::cout << "\n";
    }
    template<class ... Types>
    void RunTest(std::function<bool(Types...)> test, std::string name, Types ... args) {
      test_count++;
      std::cout << "Running test: " << name << "\n";
      std::stringstream str;
      if (test(args...)) {
        tests_passed++;
        str << name << " passed!\n";
        PrintGreen(str.str());
      }
      else {
        tests_failed.push_back(name);
        str << name << " failed!\n";
        PrintRed(str.str());
      }
      std::cout << "\n";
    }

    void PrintSummary();
};

#endif
