cmake_minimum_required(VERSION 2.8)
project(main)

set(CMAKE_BUILD_TYPE  "Debug")
find_package( OpenCV REQUIRED)
## 为使用 sophus，需要使用find_package命令找到它
#find_package(Sophus REQUIRED)
#include_directories(${Sophus_INCLUDE_DIRS})

# Eigen
#find_package(Eigen REQUIRED)
#include_directories(${Eigen_INDLUDE_DIRS})
include_directories("/usr/include/eigen3")
add_executable(main   polygon.cpp)

#set(OpenCV_DIR "/home/whf/opencv/build/share/OpenCV")
include_directories(
        ${OpenCV_INCLUDE_DIRS}
)

# 链接OpenCV库
target_link_libraries(main  ${OpenCV_LIBS})
#target_link_libraries(main ${Sophus_LIBRARIES})
