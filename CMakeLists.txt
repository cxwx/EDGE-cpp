cmake_minimum_required(VERSION 3.10)
project(EDGE_CPP VERSION 1.0)

# C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# 查找 GAMERA（通过 Homebrew 安装）
find_path(GAMERA_INCLUDE_DIR gamera/Utils.h
    HINTS /opt/homebrew/opt/gamera/include /usr/local/opt/gamera/include)

find_library(GAMERA_LIBRARY gamera
    HINTS /opt/homebrew/opt/gamera/lib /usr/local/opt/gamera/lib)

# 其他依赖
find_package(Eigen3 REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP)

# 包含目录
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${GAMERA_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${GSL_INCLUDE_DIRS}
)

# 源文件（排除独立程序）
file(GLOB_RECURSE SOURCES
    "src/utils/*.cpp"
    "src/core/*.cpp"
    "src/numerics/*.cpp"
)

# 主程序源文件
set(MAIN_SRC "src/main.cpp")

# 可执行文件
add_executable(edge_cpp ${MAIN_SRC} ${SOURCES})

# 链接库
target_link_libraries(edge_cpp
    ${GAMERA_LIBRARY}
    ${GSL_LIBRARIES}
    ${GSL_CBLAS_LIBRARIES}
    ${OpenMP_CXX_LIBRARIES}
)

# 编译选项（OpenMP 支持需要特殊处理）
if(OpenMP_FOUND)
    target_compile_options(edge_cpp PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang>:-Xpreprocessor -fopenmp>
        $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fopenmp>
    )
    target_link_options(edge_cpp PRIVATE
        $<$<LINK_LANG_AND_ID:CXX,AppleClang,Clang>:-Xpreprocessor -fopenmp>
        $<$<LINK_LANG_AND_ID:CXX,GNU>:-fopenmp>
    )
endif()

# 测试可执行文件
add_executable(test_edge tests/test_basic.cpp)
target_link_libraries(test_edge
    ${GAMERA_LIBRARY}
    ${GSL_LIBRARIES}
    ${GSL_CBLAS_LIBRARIES}
)

# 坐标测试可执行文件
add_executable(test_coordinates tests/test_coordinates.cpp
                                       src/utils/EdgeUtils.cpp
                                       src/core/Diffusion.cpp
                                       src/core/Spectrum.cpp
                                       src/core/Coordinates.cpp
                                       src/numerics/Interpolation.cpp)
target_link_libraries(test_coordinates
    ${GAMERA_LIBRARY}
    ${GSL_LIBRARIES}
    ${GSL_CBLAS_LIBRARIES}
)

# 完整EDGE可执行文件
add_executable(edge_complete src/edge_complete.cpp
                             src/utils/EdgeUtils.cpp
                             src/core/Diffusion.cpp
                             src/core/Spectrum.cpp
                             src/core/Coordinates.cpp
                             src/core/MultiSource.cpp
                             src/numerics/Interpolation.cpp)
target_link_libraries(edge_complete
    ${GAMERA_LIBRARY}
    ${GSL_LIBRARIES}
    ${GSL_CBLAS_LIBRARIES}
)

# 简化演示程序可执行文件
add_executable(edge_simple src/edge_simple.cpp
                           src/utils/EdgeUtils.cpp
                           src/core/Diffusion.cpp
                           src/core/Spectrum.cpp
                           src/core/Coordinates.cpp
                           src/core/MultiSource.cpp
                           src/numerics/Interpolation.cpp)
target_link_libraries(edge_simple
    ${GAMERA_LIBRARY}
    ${GSL_LIBRARIES}
    ${GSL_CBLAS_LIBRARIES}
)

# 安装规则
install(TARGETS edge_cpp edge_simple DESTINATION bin)

# 打印配置信息
message(STATUS "GAMERA include: ${GAMERA_INCLUDE_DIR}")
message(STATUS "GAMERA library: ${GAMERA_LIBRARY}")
message(STATUS "Eigen3 include: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "OpenMP found: ${OpenMP_FOUND}")
