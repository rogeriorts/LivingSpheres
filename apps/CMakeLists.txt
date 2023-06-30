add_executable(app app.cpp ../include/LivingSPheres/simulator.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(app PUBLIC OpenMP::OpenMP_CXX)
endif()