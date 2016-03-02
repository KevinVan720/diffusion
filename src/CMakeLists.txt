add_executable(${PROJECT_NAME}
  df_main.f  
  df_rnd.f  
  JetOutput-0.7.f  
  oscarnxtev.f  
  parameter_inp.f  
  rad_fnc.f  
  rad_inp.f  
  qhat.f
  )

# install compiled binary
install (TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
