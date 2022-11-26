#
# Include and configure the simulation programs
#
#
set(program_APPS
  test_buckling
  test_swim
  test_buckling2D
  learn2D_cellflow
  learn2D_reduced
  policy_gradient
  traj_snaps
  swim
  swimmer_traj
  markovianization
  bootstrapping_cellflow
  semi_differential_sarsa
  pg_with_et
  )

foreach (program ${program_APPS})
    install_ieel_application(${program} bin)
endforeach ()

