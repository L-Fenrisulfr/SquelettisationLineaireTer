add_executable(main MACOSX_BUNDLE main.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(main ${VTK_LIBRARIES})
else()
  target_link_libraries(main vtkHybrid vtkWidgets)
endif()

