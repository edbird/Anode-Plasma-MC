CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(AnodePlasmaMC)

ADD_EXECUTABLE(AnodePlasmaMC main.cpp)

# add the vector 3 library folder
INCLUDE_DIRECTORIES("~/Vector3")
