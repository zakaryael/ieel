#
# Include and configure the simulation programs
#
#
set(program_APPS
  test_buckling
  test_swim
  test_buckling2D
  test_swim2D
  swim2D_cellflow
  learn2D_noflow
  learn2D_cellflow
  swim2D_learned_cellflow
  )

foreach (program ${program_APPS})
    install_ieel_application(${program} bin)
endforeach ()

