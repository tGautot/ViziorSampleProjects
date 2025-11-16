set( IMGUI_DIR  ../external/imgui )

set( IMGUI_SOURCES 
  ${IMGUI_DIR}/imgui.cpp 
  ${IMGUI_DIR}/imgui_demo.cpp 
  ${IMGUI_DIR}/imgui_draw.cpp 
  ${IMGUI_DIR}/imgui_tables.cpp 
  ${IMGUI_DIR}/imgui_widgets.cpp
  ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp 
  ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp)

add_executable(circle_path main.cpp ${IMGUI_SOURCES})

target_link_libraries(circle_path Vizior::Vizior)
target_include_directories(circle_path PRIVATE ${IMGUI_DIR} ${IMGUI_DIR}/backends)

