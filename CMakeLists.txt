cmake_minimum_required(VERSION 3.10.2)

project(LivingSpheres)

include(FetchContent)

add_subdirectory(apps)

enable_testing()

add_subdirectory(tests)

