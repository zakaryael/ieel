#
# Include and configure the simulation programs
#
#
set(program_APPS
  test_buckling
  )

foreach (program ${program_APPS})
    install_ieel_application(${program} bin)
endforeach ()

