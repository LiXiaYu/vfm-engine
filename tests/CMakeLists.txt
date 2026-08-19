add_executable(vfm_jacobian_policy_tests
    unit/jacobian_policy_tests.cpp)

target_include_directories(vfm_jacobian_policy_tests
    PRIVATE
        "${PROJECT_SOURCE_DIR}/FEBio_VFM_Task/include")

target_compile_features(vfm_jacobian_policy_tests PRIVATE cxx_std_20)

add_test(NAME jacobian_policy COMMAND vfm_jacobian_policy_tests)
