Module SolutePartrace

Contains
  !****************************************************************
  ! transformation from concentration [micromol / cm^3] to osmotic head [cm]
  ! Fac = factor for transformation
  Subroutine Salinity
    USE typedef
    USE DoussanMat, ONLY: PHs_osmotic
    USE SolData, ONLY: SoilSoluteConcentration
    IMPLICIT NONE
    Real(dp), Parameter :: Fac=-50.0

    Phs_osmotic = Fac * SoilSoluteConcentration

  End Subroutine Salinity
  !********************************************************************************
End Module SolutePartrace
