#pragma once
#include "common_FEBio.h"
#include <string>
#include <pybind11/pybind11.h>
#include <fstream>

template <typename ReturnType, typename... Args>
ReturnType pyfunctioncall(FEModel* fem, const std::string& logs_string_firstline, std::ofstream& outFile, pybind11::module_& pymodule, const std::string& function_name, Args&&... args)
{
    try {
        return pymodule.attr(function_name.c_str())(std::forward<Args>(args)...).template cast<ReturnType>();
    }
    catch (const std::exception& e) {
        std::string logs_string = logs_string_firstline + "\n" + e.what();
        write_to_log_2(fem, logs_string, outFile);
        throw;
    }
    catch (...) {
        std::string logs_string = logs_string_firstline + "\n" + "Unknown exception caught";
        write_to_log_2(fem, logs_string, outFile);
        throw;
    }
}

template <typename... Args>
void pyfunctioncall(FEModel* fem, const std::string& logs_string_firstline, std::ofstream& outFile, pybind11::module_& pymodule, const std::string& function_name, Args&&... args)
{
    try {
        pymodule.attr(function_name.c_str())(std::forward<Args>(args)...);
    }
    catch (const std::exception& e) {
        std::string logs_string = logs_string_firstline + "\n" + e.what();
        write_to_log_2(fem, logs_string, outFile);
        throw;
    }
    catch (...) {
        std::string logs_string = logs_string_firstline + "\n" + "Unknown exception caught";
        write_to_log_2(fem, logs_string, outFile);
        throw;
    }
}