--
-- -----------------------------------------------------------------------------
-- Office PVTN (PVTN) Data Section Version 2013.1 Jun 14 2013
-- -----------------------------------------------------------------------------
--
-- File: 11_PVT.INC
-- Created on: Feb-18-2016 at: 11:52:00
--
-- *****************************************************************************
-- *                                  WARNING                                  *
-- *                THIS FILE HAS BEEN AUTOMATICALLY GENERATED.                *
-- *          ANY ATTEMPT TO EDIT MANUALLY MAY RESULT IN INVALID DATA.         *
-- *****************************************************************************
--
-- OFFICE-PVTN-HEADER-DATA
--
-- -----------------------------------------------------------------------------
-- Office PVTN Keywords
-- -----------------------------------------------------------------------------
--
ECHO
---------------------
-- SURFACE DENSITIES OF RESERVOIR FLUIDS
--
--        OIL  WATER    GAS
DENSITY
-- 
-- Fluid Densities at Surface Conditions
-- 
        49.1     63.02   0.06054
/
 
----------------------
-- PVT PROPERTIES OF DRY GAS (NO VAPOURISED OIL)
-- WE WOULD USE PVTG TO SPECIFY THE PROPERTIES OF WET GAS
--
--   PGAS   BGAS   VISGAS
PVDG
-- 
-- Dry Gas PVT Properties (No Vapourised Oil)
-- 
        14.7   166.666     0.008
       264.7    12.093    0.0096
       514.7     6.274    0.0112
      1014.7     3.197     0.014
      2014.7     1.614    0.0189
      2514.7     1.294    0.0208
      3014.7      1.08    0.0228
      4014.7     0.811    0.0268
      5014.7     0.649    0.0309
      9014.7     0.386     0.047
/
 
-- PVT PROPERTIES OF LIVE OIL (WITH DISSOLVED GAS)
-- WE WOULD USE PVDO TO SPECIFY THE PROPERTIES OF DEAD OIL
--
-- FOR EACH VALUE OF RS THE SATURATION PRESSURE, FVF AND VISCOSITY
-- ARE SPECIFIED. FOR RS=1.27 AND 1.618, THE FVF AND VISCOSITY OF
-- UNDERSATURATED OIL ARE DEFINED AS A FUNCTION OF PRESSURE. DATA
-- FOR UNDERSATURATED OIL MAY BE SUPPLIED FOR ANY RS, BUT MUST BE
-- SUPPLIED FOR THE HIGHEST RS (1.618).
--
--   RS      POIL  FVFO  VISO
PVTO
-- 
-- Live Oil PVT Properties (Dissolved Gas)
-- 
          0.14       800       1.2      1.09 /
          0.18      1000     1.225      1.01 /
          0.22      1200      1.24      0.96 /
          0.26      1400      1.28       0.9 /
           0.3      1600     1.285      0.85 /
          0.35      1800      1.31      0.77 /
           0.4      2000     1.335      0.74 /
          0.45      2200      1.36       0.7 /
           0.5      2400     1.385      0.66 /
          0.55      2600      1.41      0.62 /
           0.6      2800      1.43      0.58 /
          0.65      3000      1.46      0.54 /
          0.71      3200     1.485      0.51 /
          0.84      3600      1.55      0.45 /
          0.98      4000      1.62       0.4 /
          1.13      4400     1.675      0.37 /
          1.31      4800      1.73      0.35 /
         1.405      5000      1.76     0.335
                    5200      1.75      0.34
                    5600      1.74      0.36
                    6000      1.73      0.37
                    6400      1.71      0.39
                    6800       1.7       0.4 /
/
 
ECHO
------------------- PVT -----------------------
-- ROCK COMPRESSIBILITY
--
--    REF. PRES   COMPRESSIBILITY
ROCK
-- 
-- Rock Properties
-- 
        14.7    3e-006                                        
/
 
