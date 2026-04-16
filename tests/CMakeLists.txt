add_executable(testgjk testgjk.cpp)
target_link_libraries(testgjk PRIVATE gjkepa Eigen3::Eigen)
target_include_directories(testgjk PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(testepa testepa.cpp)
target_link_libraries(testepa PRIVATE gjkepa Eigen3::Eigen)
target_include_directories(testepa PRIVATE ${PROJECT_SOURCE_DIR}/include)