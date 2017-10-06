add_executable(${PROJECT_NAME}
  df_main.f  
  df_rnd.f  
#  JetOutput-0.7.f  
  jet_output_hdf5.for
  hydro_3d_hdf5.for
  oscarnxtev.f  
  parameter_inp.f  
  rad_fnc.f  
  rad_inp.f  
  qhat.f
  )


target_link_libraries(${PROJECT_NAME} ${HDF5_FORTRAN_LIB})

# install compiled binary
install (TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
