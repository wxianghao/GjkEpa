add_executable(testgjk testgjk.cpp)
target_link_libraries(testgjk PRIVATE gjkepa Eigen3::Eigen)
target_include_directories(testgjk PRIVATE ${PROJECT_SOURCE_DIR}/include)