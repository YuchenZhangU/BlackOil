--
-- -----------------------------------------------------------------------------
-- Office SCAL (SCAL) Data Section Version 2013.1 Jun 14 2013
-- -----------------------------------------------------------------------------
--
-- File: 11_SCAL.INC
-- Created on: Feb-18-2016 at: 11:52:00
--
-- *****************************************************************************
-- *                                  WARNING                                  *
-- *                THIS FILE HAS BEEN AUTOMATICALLY GENERATED.                *
-- *          ANY ATTEMPT TO EDIT MANUALLY MAY RESULT IN INVALID DATA.         *
-- *****************************************************************************
--
-- OFFICE-SCAL-HEADER-DATA
--
-- -----------------------------------------------------------------------------
-- Office SCAL Keywords
-- -----------------------------------------------------------------------------
--
ECHO
-- Report Levels for Grid Section Data
-- RPTGRID
-- 'DX'
-- 'DY'
-- 'DZ'
-- 'PERMX'
-- 'PERMY'
-- /
--Gas RELATIVE PERMEABILITY AND CAPILLARY PRESSURE ARE TABULATED AS
-- A FUNCTION OFGas SATURATION.
-----------------
-- SIMILARLY FOR GAS
--
--  SGAS   KRG   PCOG
SGFN
-- 
-- Gas Saturation Functions
-- 
             0           0           0
           0.1         0.1           0
           0.2         0.2           0
           0.3         0.3           0
           0.4         0.4           0
           0.5         0.5           0
           0.6         0.6           0
           0.7         0.7           0
           0.8         0.8           0
/
 
----------------
-- OIL RELATIVE PERMEABILITY IS TABULATED AGAINST OIL SATURATION
-- FOR OIL-WATER AND OIL-GAS-CONNATE WATER CASES
--
--  SOIL   KROG
SOF2
-- 
-- Oil Saturation Functions
-- 
             0           0
           0.2           0
           0.3      0.0109
           0.4     0.04375
           0.5   0.0984375
           0.6       0.175
           0.7   0.2734375
           0.8     0.39375
           0.9   0.5359375
             1         0.7
/
 
