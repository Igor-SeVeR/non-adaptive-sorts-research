add_executable(${PROJECT_NAME} experiment.cpp generator.cpp generator.h)

target_include_directories(${PROJECT_NAME} PUBLIC ../adaptive_sorts ../nonadaptive_sorts)
target_link_directories(${PROJECT_NAME} PRIVATE ../adaptive_sorts ../nonadaptive_sorts)
target_link_libraries(${PROJECT_NAME} PRIVATE adaptive_sorts nonadaptive_sorts)