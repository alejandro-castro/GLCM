# CMakeLists.txt

# Older versions of CMake are likely to work just fine but, since
# I don't know where to cut off I just use the version I'm using
cmake_minimum_required(VERSION "3.17")

# name of this example project
project(glcm)

# Tell compiler to use C++ 14 features which is needed because
# Clang version is often behind in the XCode installation
set(CMAKE_CXX_STANDARD 14)

# configure the necessary common CMake environment variables
# needed to include and link the OpenCV program into this
# demo project, namely OpenCV_INCLUDE_DIRS and OpenCV_LIBS
find_package( OpenCV REQUIRED )

# tell the build to include the headers from OpenCV
include_directories( ${OpenCV_INCLUDE_DIRS} )
# Includes and libraries to use manage mat files
include_directories("/Applications/MATLAB_R2019b.app/extern/include")
link_directories("/Applications/MATLAB_R2019b.app/bin/maci64" "/Applications/MATLAB_R2019b.app/extern/lib/maci64")

# Includes user libraries
include_directories("${CMAKE_SOURCE_DIR}/Libraries")
# specify the executable target to be built
add_executable(glcm main.cpp ${CMAKE_SOURCE_DIR}/Libraries/glcm.cpp ${CMAKE_SOURCE_DIR}/Libraries/CVMatlabInterface.cpp ${CMAKE_SOURCE_DIR}/Libraries/LinkedList.cpp)

# tell it to link the executable target against OpenCV
target_link_libraries(glcm ${OpenCV_LIBS} mat mx) #I didn't need the mex library, maybe I will in the future
