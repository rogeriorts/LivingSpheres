# Testing library
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.10)
FetchContent_MakeAvailable(catch2)
# Adds Catch2::Catch2

# Tests need to be added as executables first
add_executable(
  test_cell_movements_due_external_forces 
  test_cell_movements_due_external_forces.cpp
  ../include/LivingSPheres/simulator.cpp
  )
add_executable(
  test_grid 
  test_grid.cpp 
  ../include/LivingSPheres/simulator.cpp
)

# I'm using C++17 in the test
#target_compile_features(testlib PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(test_cell_movements_due_external_forces PRIVATE  Catch2::Catch2)
target_link_libraries(test_grid PRIVATE  Catch2::Catch2)



# If you register a test, then ctest and make test will run it. You can also run
# examples and check the output, as well.
add_test(NAME test_cell_movements_due_external_forces COMMAND $<TARGET_FILE:test_cell_movements_due_external_forces>) # Command can be a target
add_test(NAME test_grid COMMAND $<TARGET_FILE:test_grid>) # Command can be a target