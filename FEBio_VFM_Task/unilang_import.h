#pragma once
#include <cstdint>

// if in windows Define DLLIMPORT as __declspec(dllimport), else Define it as nothing
#ifdef _WIN32
#define DLLIMPORT __declspec(dllimport)
#else
#define DLLIMPORT
#endif

extern "C" DLLIMPORT int unilang_load_file(char* filename);

extern "C" DLLIMPORT int unilang_regist_int_value(const char* objectname, int* object);
extern "C" DLLIMPORT int unilang_regist_double_value(const char* objectname, double* object);
extern "C" DLLIMPORT int unilang_regist_bool_value(const char* objectname, bool* object);
extern "C" DLLIMPORT int unilang_regist_pointer_value(const char* objectname, void** object);
extern "C" DLLIMPORT int unilang_regist_uintptr_t_value(const char* objectname, std::uintptr_t* object);

extern "C" DLLIMPORT int unilang_regist_function_var(const char* objectname, void** object, const char* result_type, const int parameters_number, const char** parameters_type);