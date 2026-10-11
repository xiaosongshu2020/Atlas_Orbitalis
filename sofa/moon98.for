      SUBROUTINE iau_MOON98 ( DATE1, DATE2, PV )
*+
*  - - - - - - - - - - -
*   i a u _ M O O N 9 8
*  - - - - - - - - - - -
*
*  Approximate geocentric position and velocity of the Moon.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  n.b. Not IAU-endorsed and without canonical status.
*
*  Given:
*     DATE1    d       TT date part A (Notes 1,4)
*     DATE2    d       TT date part B (Notes 1,4)
*
*  Returned:
*     PV     d(3,2)    Moon p,v, GCRS (au, au/d, Note 5)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways, among
*     others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*     The limited accuracy of the present algorithm is such that any
*     of the methods is satisfactory.
*
*  2) This routine is a full implementation of the algorithm
*     published by Meeus (see reference) except that the light-time
*     correction to the Moon's mean longitude has been omitted.
*
*  3) Comparisons with ELP/MPP02 over the interval 1950-2100 gave RMS
*     errors of 2.9 arcsec in geocentric direction, 6.1 km in position
*     and 36 mm/s in velocity.  The worst case errors were 18.3 arcsec
*     in geocentric direction, 31.7 km in position and 172 mm/s in
*     velocity.
*
*  4) The original algorithm is expressed in terms of "dynamical time",
*     which can either be TDB or TT without any significant change in
*     accuracy.  UT cannot be used without incurring significant errors
*     (30 arcsec in the present era) due to the Moon's 0.5 arcsec/sec
*     movement.
*
*  5) The result is with respect to the GCRS (the same as J2000.0 mean
*     equator and equinox to within 23 mas).
*
*  6) Velocity is obtained by a complete analytical differentiation
*     of the Meeus model.
*
*  7) The Meeus algorithm generates position and velocity in mean
*     ecliptic coordinates of date, which the present routine then
*     rotates into GCRS.  Because the ecliptic system is precessing,
*     there is a coupling between this spin (about 1.4 degrees per
*     century) and the Moon position that produces a small velocity
*     contribution.  In the present routine this effect is neglected
*     as it corresponds to a maximum difference of less than 3 mm/s
*     and increases the RMS error by only 0.4%.
*
*  References:
*
*     Meeus, J., Astronomical Algorithms, 2nd edition, Willmann-Bell,
*     1998, p337.
*
*     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G. & Laskar, J., Astron.Astrophys., 1994, 282, 663
*
*  Called:
*     iau_S2PV      spherical coordinates to pv-vector
*     iau_PFW06     bias-precession F-W angles, IAU 2006
*     iau_IR        initialize r-matrix to identity
*     iau_RZ        rotate around Z-axis
*     iau_RX        rotate around X-axis
*     iau_RXPV      product of r-matrix and pv-vector
*
*  This revision:   2023 March 20
*
*  SOFA release 2023-10-11
*
*  Copyright (C) 2023 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, PV(3,2)

*  Astronomical unit (m)
      DOUBLE PRECISION DAU
      PARAMETER ( DAU = 149597870.7D3 )

*  JD for J2000.0
      DOUBLE PRECISION DJ00
      PARAMETER (DJ00 = 2451545D0 )

*  Degrees to radians
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*
*  Coefficients for fundamental arguments:
*
*  . Powers of time in Julian centuries
*  . Units are degrees.
*

*  Moon's mean longitude (wrt mean equinox and ecliptic of date)
      DOUBLE PRECISION ELP0, ELP1, ELP2, ELP3, ELP4, ELP, DELP
      PARAMETER ( ELP0 = 218.31665436D0,
     :            ELP1 = 481267.88123421D0,
     :            ELP2 = -0.0015786D0,
     :            ELP3=  1D0 / 538841D0,
     :            ELP4 = -1D0 / 65194000D0 )

*  Moon's mean anomaly
      DOUBLE PRECISION EMP0, EMP1, EMP2, EMP3, EMP4, EMP, DEMP
      PARAMETER ( EMP0 = 134.9633964D0,
     :            EMP1 = 477198.8675055D0,
     :            EMP2 = 0.0087414D0,
     :            EMP3 = 1D0 / 69699D0,
     :            EMP4 = -1D0 / 14712000D0 )

*  Moon's mean elongation
      DOUBLE PRECISION D0, D1, D2, D3, D4, D, DD
      PARAMETER (D0 = 297.8501921D0,
     :           D1 = 445267.1114034D0,
     :           D2 = -0.0018819D0,
     :           D3 = 1D0 / 545868D0,
     :           D4 = 1D0 / 113065000D0 )

*  Sun's mean anomaly
      DOUBLE PRECISION EM0, EM1, EM2, EM3, EM4, EM, DEM
      PARAMETER (EM0 = 357.5291092D0,
     :           EM1 = 35999.0502909D0,
     :           EM2 = -0.0001536D0,
     :           EM3 = 1D0 / 2449000D0,
     :           EM4 = 0D0 )

*  Mean distance of the Moon from its ascending node
      DOUBLE PRECISION F0, F1, F2, F3, F4, F, DF
      PARAMETER (F0 = 93.2720950D0,
     :           F1 = 483202.0175233D0,
     :           F2 = -0.0036539D0,
     :           F3 = 1D0 / 3526000D0,
     :           F4 = 1D0 / 863310000D0 )

*
*  Other arguments
*

*  Meeus A_1, due to Venus (deg)
      DOUBLE PRECISION A10, A11, A1, DA1
      PARAMETER (  A10 = 119.75D0,
     :             A11 = 131.849D0 )

*  Meeus A_2, due to Jupiter (deg)
      DOUBLE PRECISION A20, A21, A2, DA2
      PARAMETER ( A20 = 53.09D0,
     :            A21 = 479264.290D0 )

*  Meeus A_3, due to sidereal motion of the Moon in longitude (deg)
      DOUBLE PRECISION A30, A31, A3, DA3
      PARAMETER ( A30 = 313.45D0,
     :            A31 = 481266.484D0 )

*  Coefficients for Meeus "additive terms" (deg)
      DOUBLE PRECISION AL1, AL2, AL3
      PARAMETER ( AL1 =  0.003958D0,
     :            AL2 =  0.001962D0,
     :            AL3 =  0.000318D0 )
      DOUBLE PRECISION AB1, AB2, AB3, AB4, AB5, AB6
      PARAMETER ( AB1 = -0.002235D0,
     :            AB2 =  0.000382D0,
     :            AB3 =  0.000175D0,
     :            AB4 =  0.000175D0,
     :            AB5 =  0.000127D0,
     :            AB6 = -0.000115D0 )

*  Fixed term in distance (m)
      DOUBLE PRECISION R0
      PARAMETER ( R0 = 385000560D0 )

*  Coefficients for (dimensionless) E factor
      DOUBLE PRECISION E1, E2, E, DE, ESQ, DESQ
      PARAMETER ( E1 = -0.002516D0,
     :            E2 = -0.0000074D0 )

*  Miscellaneous
      INTEGER N, I
      DOUBLE PRECISION T, ELPMF, DELPMF, VEL, VDEL, VR, VDR, A1MF,
     :                 DA1MF, A1PF, DA1PF, DLPMP, SLPMP, VB, VDB, V, DV,
     :                 EMN, EMPN, DN, FN, EN, DEN, ARG, DARG, FARG,
     :                 COEFF, EL, DEL, R, DR, B, DB, GAMB, PHIB, PSIB,
     :                 EPSA, RM(3,3)

*
*  Coefficients for Moon position series (L,B,R)
*
*   TLR(1,N)      =  coefficient of L sine term (deg)
*   TLR(2,N)      =  coefficient of R cosine term (m)
*   TB(N)         =  coefficient B sine term (deg)
*   ITx(1-4,N)    =  coefficients of D, M, M', F in argument
*

      INTEGER NLR, NB
      PARAMETER ( NLR = 60, NB = 60 )
      DOUBLE PRECISION TLR(2,NLR), TB(NB)
      INTEGER ITLR(4,NLR), ITB(4,NB)

*
*  Longitude and distance series
*                                                        D   M   M'  F
*
      DATA (TLR(I, 1),I=1,2) /  6.288774D0, -20905355D0 /,
     :     (ITLR(I, 1),I=1,4) /                          0,  0,  1,  0 /
      DATA (TLR(I, 2),I=1,2) /  1.274027D0,  -3699111D0 /,
     :     (ITLR(I, 2),I=1,4) /                          2,  0, -1,  0 /
      DATA (TLR(I, 3),I=1,2) /  0.658314D0,  -2955968D0 /,
     :     (ITLR(I, 3),I=1,4) /                          2,  0,  0,  0 /
      DATA (TLR(I, 4),I=1,2) /  0.213618D0,   -569925D0 /,
     :     (ITLR(I, 4),I=1,4) /                          0,  0,  2,  0 /
      DATA (TLR(I, 5),I=1,2) / -0.185116D0,     48888D0 /,
     :     (ITLR(I, 5),I=1,4) /                          0,  1,  0,  0 /
      DATA (TLR(I, 6),I=1,2) / -0.114332D0,     -3149D0 /,
     :     (ITLR(I, 6),I=1,4) /                          0,  0,  0,  2 /
      DATA (TLR(I, 7),I=1,2) /  0.058793D0,    246158D0 /,
     :     (ITLR(I, 7),I=1,4) /                          2,  0, -2,  0 /
      DATA (TLR(I, 8),I=1,2) /  0.057066D0,   -152138D0 /,
     :     (ITLR(I, 8),I=1,4) /                          2, -1, -1,  0 /
      DATA (TLR(I, 9),I=1,2) /  0.053322D0,   -170733D0 /,
     :     (ITLR(I, 9),I=1,4) /                          2,  0,  1,  0 /
      DATA (TLR(I,10),I=1,2) /  0.045758D0,   -204586D0 /,
     :     (ITLR(I,10),I=1,4) /                          2, -1,  0,  0 /
      DATA (TLR(I,11),I=1,2) / -0.040923D0,   -129620D0 /,
     :     (ITLR(I,11),I=1,4) /                          0,  1, -1,  0 /
      DATA (TLR(I,12),I=1,2) / -0.034720D0,    108743D0 /,
     :     (ITLR(I,12),I=1,4) /                          1,  0,  0,  0 /
      DATA (TLR(I,13),I=1,2) / -0.030383D0,    104755D0 /,
     :     (ITLR(I,13),I=1,4) /                          0,  1,  1,  0 /
      DATA (TLR(I,14),I=1,2) /  0.015327D0,     10321D0 /,
     :     (ITLR(I,14),I=1,4) /                          2,  0,  0, -2 /
      DATA (TLR(I,15),I=1,2) / -0.012528D0,         0D0 /,
     :     (ITLR(I,15),I=1,4) /                          0,  0,  1,  2 /
      DATA (TLR(I,16),I=1,2) /  0.010980D0,     79661D0 /,
     :     (ITLR(I,16),I=1,4) /                          0,  0,  1, -2 /
      DATA (TLR(I,17),I=1,2) /  0.010675D0,    -34782D0 /,
     :     (ITLR(I,17),I=1,4) /                          4,  0, -1,  0 /
      DATA (TLR(I,18),I=1,2) /  0.010034D0,    -23210D0 /,
     :     (ITLR(I,18),I=1,4) /                          0,  0,  3,  0 /
      DATA (TLR(I,19),I=1,2) /  0.008548D0,    -21636D0 /,
     :     (ITLR(I,19),I=1,4) /                          4,  0, -2,  0 /
      DATA (TLR(I,20),I=1,2) / -0.007888D0,     24208D0 /,
     :     (ITLR(I,20),I=1,4) /                          2,  1, -1,  0 /
      DATA (TLR(I,21),I=1,2) / -0.006766D0,     30824D0 /,
     :     (ITLR(I,21),I=1,4) /                          2,  1,  0,  0 /
      DATA (TLR(I,22),I=1,2) / -0.005163D0,     -8379D0 /,
     :     (ITLR(I,22),I=1,4) /                          1,  0, -1,  0 /
      DATA (TLR(I,23),I=1,2) /  0.004987D0,    -16675D0 /,
     :     (ITLR(I,23),I=1,4) /                          1,  1,  0,  0 /
      DATA (TLR(I,24),I=1,2) /  0.004036D0,    -12831D0 /,
     :     (ITLR(I,24),I=1,4) /                          2, -1,  1,  0 /
      DATA (TLR(I,25),I=1,2) /  0.003994D0,    -10445D0 /,
     :     (ITLR(I,25),I=1,4) /                          2,  0,  2,  0 /
      DATA (TLR(I,26),I=1,2) /  0.003861D0,    -11650D0 /,
     :     (ITLR(I,26),I=1,4) /                          4,  0,  0,  0 /
      DATA (TLR(I,27),I=1,2) /  0.003665D0,     14403D0 /,
     :     (ITLR(I,27),I=1,4) /                          2,  0, -3,  0 /
      DATA (TLR(I,28),I=1,2) / -0.002689D0,     -7003D0 /,
     :     (ITLR(I,28),I=1,4) /                          0,  1, -2,  0 /
      DATA (TLR(I,29),I=1,2) / -0.002602D0,         0D0 /,
     :     (ITLR(I,29),I=1,4) /                          2,  0, -1,  2 /
      DATA (TLR(I,30),I=1,2) /  0.002390D0,     10056D0 /,
     :     (ITLR(I,30),I=1,4) /                          2, -1, -2,  0 /
      DATA (TLR(I,31),I=1,2) / -0.002348D0,      6322D0 /,
     :     (ITLR(I,31),I=1,4) /                          1,  0,  1,  0 /
      DATA (TLR(I,32),I=1,2) /  0.002236D0,     -9884D0 /,
     :     (ITLR(I,32),I=1,4) /                          2, -2,  0,  0 /
      DATA (TLR(I,33),I=1,2) / -0.002120D0,      5751D0 /,
     :     (ITLR(I,33),I=1,4) /                          0,  1,  2,  0 /
      DATA (TLR(I,34),I=1,2) / -0.002069D0,         0D0 /,
     :     (ITLR(I,34),I=1,4) /                          0,  2,  0,  0 /
      DATA (TLR(I,35),I=1,2) /  0.002048D0,     -4950D0 /,
     :     (ITLR(I,35),I=1,4) /                          2, -2, -1,  0 /
      DATA (TLR(I,36),I=1,2) / -0.001773D0,      4130D0 /,
     :     (ITLR(I,36),I=1,4) /                          2,  0,  1, -2 /
      DATA (TLR(I,37),I=1,2) / -0.001595D0,         0D0 /,
     :     (ITLR(I,37),I=1,4) /                          2,  0,  0,  2 /
      DATA (TLR(I,38),I=1,2) /  0.001215D0,     -3958D0 /,
     :     (ITLR(I,38),I=1,4) /                          4, -1, -1,  0 /
      DATA (TLR(I,39),I=1,2) / -0.001110D0,         0D0 /,
     :     (ITLR(I,39),I=1,4) /                          0,  0,  2,  2 /
      DATA (TLR(I,40),I=1,2) / -0.000892D0,      3258D0 /,
     :     (ITLR(I,40),I=1,4) /                          3,  0, -1,  0 /
      DATA (TLR(I,41),I=1,2) / -0.000810D0,      2616D0 /,
     :     (ITLR(I,41),I=1,4) /                          2,  1,  1,  0 /
      DATA (TLR(I,42),I=1,2) /  0.000759D0,     -1897D0 /,
     :     (ITLR(I,42),I=1,4) /                          4, -1, -2,  0 /
      DATA (TLR(I,43),I=1,2) / -0.000713D0,     -2117D0 /,
     :     (ITLR(I,43),I=1,4) /                          0,  2, -1,  0 /
      DATA (TLR(I,44),I=1,2) / -0.000700D0,      2354D0 /,
     :     (ITLR(I,44),I=1,4) /                          2,  2, -1,  0 /
      DATA (TLR(I,45),I=1,2) /  0.000691D0,         0D0 /,
     :     (ITLR(I,45),I=1,4) /                          2,  1, -2,  0 /
      DATA (TLR(I,46),I=1,2) /  0.000596D0,         0D0 /,
     :     (ITLR(I,46),I=1,4) /                          2, -1,  0, -2 /
      DATA (TLR(I,47),I=1,2) /  0.000549D0,     -1423D0 /,
     :     (ITLR(I,47),I=1,4) /                          4,  0,  1,  0 /
      DATA (TLR(I,48),I=1,2) /  0.000537D0,     -1117D0 /,
     :     (ITLR(I,48),I=1,4) /                          0,  0,  4,  0 /
      DATA (TLR(I,49),I=1,2) /  0.000520D0,     -1571D0 /,
     :     (ITLR(I,49),I=1,4) /                          4, -1,  0,  0 /
      DATA (TLR(I,50),I=1,2) / -0.000487D0,     -1739D0 /,
     :     (ITLR(I,50),I=1,4) /                          1,  0, -2,  0 /
      DATA (TLR(I,51),I=1,2) / -0.000399D0,         0D0 /,
     :     (ITLR(I,51),I=1,4) /                          2,  1,  0, -2 /
      DATA (TLR(I,52),I=1,2) / -0.000381D0,     -4421D0 /,
     :     (ITLR(I,52),I=1,4) /                          0,  0,  2, -2 /
      DATA (TLR(I,53),I=1,2) /  0.000351D0,         0D0 /,
     :     (ITLR(I,53),I=1,4) /                          1,  1,  1,  0 /
      DATA (TLR(I,54),I=1,2) / -0.000340D0,         0D0 /,
     :     (ITLR(I,54),I=1,4) /                          3,  0, -2,  0 /
      DATA (TLR(I,55),I=1,2) /  0.000330D0,         0D0 /,
     :     (ITLR(I,55),I=1,4) /                          4,  0, -3,  0 /
      DATA (TLR(I,56),I=1,2) /  0.000327D0,         0D0 /,
     :     (ITLR(I,56),I=1,4) /                          2, -1,  2,  0 /
      DATA (TLR(I,57),I=1,2) / -0.000323D0,      1165D0 /,
     :     (ITLR(I,57),I=1,4) /                          0,  2,  1,  0 /
      DATA (TLR(I,58),I=1,2) /  0.000299D0,         0D0 /,
     :     (ITLR(I,58),I=1,4) /                          1,  1, -1,  0 /
      DATA (TLR(I,59),I=1,2) /  0.000294D0,         0D0 /,
     :     (ITLR(I,59),I=1,4) /                          2,  0,  3,  0 /
      DATA (TLR(I,60),I=1,2) /  0.000000D0,      8752D0 /,
     :     (ITLR(I,60),I=1,4) /                          2,  0, -1, -2 /

*
*  Latitude series
*                                                        D   M   M'  F

      DATA TB( 1) /  5.128122D0 /,
     :     (ITB(I, 1),I=1,4) /                           0,  0,  0,  1 /
      DATA TB( 2) /  0.280602D0 /,
     :     (ITB(I, 2),I=1,4) /                           0,  0,  1,  1 /
      DATA TB( 3) /  0.277693D0 /,
     :     (ITB(I, 3),I=1,4) /                           0,  0,  1, -1 /
      DATA TB( 4) /  0.173237D0 /,
     :     (ITB(I, 4),I=1,4) /                           2,  0,  0, -1 /
      DATA TB( 5) /  0.055413D0 /,
     :     (ITB(I, 5),I=1,4) /                           2,  0, -1,  1 /
      DATA TB( 6) /  0.046271D0 /,
     :     (ITB(I, 6),I=1,4) /                           2,  0, -1, -1 /
      DATA TB( 7) /  0.032573D0 /,
     :     (ITB(I, 7),I=1,4) /                           2,  0,  0,  1 /
      DATA TB( 8) /  0.017198D0 /,
     :     (ITB(I, 8),I=1,4) /                           0,  0,  2,  1 /
      DATA TB( 9) /  0.009266D0 /,
     :     (ITB(I, 9),I=1,4) /                           2,  0,  1, -1 /
      DATA TB(10) /  0.008822D0 /,
     :     (ITB(I,10),I=1,4) /                           0,  0,  2, -1 /
      DATA TB(11) /  0.008216D0 /,
     :     (ITB(I,11),I=1,4) /                           2, -1,  0, -1 /
      DATA TB(12) /  0.004324D0 /,
     :     (ITB(I,12),I=1,4) /                           2,  0, -2, -1 /
      DATA TB(13) /  0.004200D0 /,
     :     (ITB(I,13),I=1,4) /                           2,  0,  1,  1 /
      DATA TB(14) / -0.003359D0 /,
     :     (ITB(I,14),I=1,4) /                           2,  1,  0, -1 /
      DATA TB(15) /  0.002463D0 /,
     :     (ITB(I,15),I=1,4) /                           2, -1, -1,  1 /
      DATA TB(16) /  0.002211D0 /,
     :     (ITB(I,16),I=1,4) /                           2, -1,  0,  1 /
      DATA TB(17) /  0.002065D0 /,
     :     (ITB(I,17),I=1,4) /                           2, -1, -1, -1 /
      DATA TB(18) / -0.001870D0 /,
     :     (ITB(I,18),I=1,4) /                           0,  1, -1, -1 /
      DATA TB(19) /  0.001828D0 /,
     :     (ITB(I,19),I=1,4) /                           4,  0, -1, -1 /
      DATA TB(20) / -0.001794D0 /,
     :     (ITB(I,20),I=1,4) /                           0,  1,  0,  1 /
      DATA TB(21) / -0.001749D0 /,
     :     (ITB(I,21),I=1,4) /                           0,  0,  0,  3 /
      DATA TB(22) / -0.001565D0 /,
     :     (ITB(I,22),I=1,4) /                           0,  1, -1,  1 /
      DATA TB(23) / -0.001491D0 /,
     :     (ITB(I,23),I=1,4) /                           1,  0,  0,  1 /
      DATA TB(24) / -0.001475D0 /,
     :     (ITB(I,24),I=1,4) /                           0,  1,  1,  1 /
      DATA TB(25) / -0.001410D0 /,
     :     (ITB(I,25),I=1,4) /                           0,  1,  1, -1 /
      DATA TB(26) / -0.001344D0 /,
     :     (ITB(I,26),I=1,4) /                           0,  1,  0, -1 /
      DATA TB(27) / -0.001335D0 /,
     :     (ITB(I,27),I=1,4) /                           1,  0,  0, -1 /
      DATA TB(28) /  0.001107D0 /,
     :     (ITB(I,28),I=1,4) /                           0,  0,  3,  1 /
      DATA TB(29) /  0.001021D0 /,
     :     (ITB(I,29),I=1,4) /                           4,  0,  0, -1 /
      DATA TB(30) /  0.000833D0 /,
     :     (ITB(I,30),I=1,4) /                           4,  0, -1,  1 /
      DATA TB(31) /  0.000777D0 /,
     :     (ITB(I,31),I=1,4) /                           0,  0,  1, -3 /
      DATA TB(32) /  0.000671D0 /,
     :     (ITB(I,32),I=1,4) /                           4,  0, -2,  1 /
      DATA TB(33) /  0.000607D0 /,
     :     (ITB(I,33),I=1,4) /                           2,  0,  0, -3 /
      DATA TB(34) /  0.000596D0 /,
     :     (ITB(I,34),I=1,4) /                           2,  0,  2, -1 /
      DATA TB(35) /  0.000491D0 /,
     :     (ITB(I,35),I=1,4) /                           2, -1,  1, -1 /
      DATA TB(36) / -0.000451D0 /,
     :     (ITB(I,36),I=1,4) /                           2,  0, -2,  1 /
      DATA TB(37) /  0.000439D0 /,
     :     (ITB(I,37),I=1,4) /                           0,  0,  3, -1 /
      DATA TB(38) /  0.000422D0 /,
     :     (ITB(I,38),I=1,4) /                           2,  0,  2,  1 /
      DATA TB(39) /  0.000421D0 /,
     :     (ITB(I,39),I=1,4) /                           2,  0, -3, -1 /
      DATA TB(40) / -0.000366D0 /,
     :     (ITB(I,40),I=1,4) /                           2,  1, -1,  1 /
      DATA TB(41) / -0.000351D0 /,
     :     (ITB(I,41),I=1,4) /                           2,  1,  0,  1 /
      DATA TB(42) /  0.000331D0 /,
     :     (ITB(I,42),I=1,4) /                           4,  0,  0,  1 /
      DATA TB(43) /  0.000315D0 /,
     :     (ITB(I,43),I=1,4) /                           2, -1,  1,  1 /
      DATA TB(44) /  0.000302D0 /,
     :     (ITB(I,44),I=1,4) /                           2, -2,  0, -1 /
      DATA TB(45) / -0.000283D0 /,
     :     (ITB(I,45),I=1,4) /                           0,  0,  1,  3 /
      DATA TB(46) / -0.000229D0 /,
     :     (ITB(I,46),I=1,4) /                           2,  1,  1, -1 /
      DATA TB(47) /  0.000223D0 /,
     :     (ITB(I,47),I=1,4) /                           1,  1,  0, -1 /
      DATA TB(48) /  0.000223D0 /,
     :     (ITB(I,48),I=1,4) /                           1,  1,  0,  1 /
      DATA TB(49) / -0.000220D0 /,
     :     (ITB(I,49),I=1,4) /                           0,  1, -2, -1 /
      DATA TB(50) / -0.000220D0 /,
     :     (ITB(I,50),I=1,4) /                           2,  1, -1, -1 /
      DATA TB(51) / -0.000185D0 /,
     :     (ITB(I,51),I=1,4) /                           1,  0,  1,  1 /
      DATA TB(52) /  0.000181D0 /,
     :     (ITB(I,52),I=1,4) /                           2, -1, -2, -1 /
      DATA TB(53) / -0.000177D0 /,
     :     (ITB(I,53),I=1,4) /                           0,  1,  2,  1 /
      DATA TB(54) /  0.000176D0 /,
     :     (ITB(I,54),I=1,4) /                           4,  0, -2, -1 /
      DATA TB(55) /  0.000166D0 /,
     :     (ITB(I,55),I=1,4) /                           4, -1, -1, -1 /
      DATA TB(56) / -0.000164D0 /,
     :     (ITB(I,56),I=1,4) /                           1,  0,  1, -1 /
      DATA TB(57) /  0.000132D0 /,
     :     (ITB(I,57),I=1,4) /                           4,  0,  1, -1 /
      DATA TB(58) / -0.000119D0 /,
     :     (ITB(I,58),I=1,4) /                           1,  0, -1, -1 /
      DATA TB(59) /  0.000115D0 /,
     :     (ITB(I,59),I=1,4) /                           4, -1,  0, -1 /
      DATA TB(60) /  0.000107D0 /,
     :     (ITB(I,60),I=1,4) /                           2, -2,  0,  1 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Centuries since J2000.
      T = ((DATE1 - DJ00) + DATE2) / DJC

*  ---------------------
*  Fundamental arguments
*  ---------------------

*  Arguments (radians) and derivatives (radians per Julian century)
*  for the current date.

*  Moon's mean longitude.
      ELP = DD2R * MOD ( ELP0
     :               + ( ELP1
     :               + ( ELP2
     :               + ( ELP3
     :               +   ELP4 * T ) * T ) * T ) * T, 360D0 )
      DELP = DD2R * (    ELP1
     :               + ( ELP2 * 2D0
     :               + ( ELP3 * 3D0
     :               +   ELP4 * 4D0 * T ) * T ) * T )

*  Moon's mean elongation.
      D = DD2R * MOD ( D0
     :             + ( D1
     :             + ( D2
     :             + ( D3
     :             +   D4 * T ) * T ) * T ) * T, 360D0 )
      DD = DD2R * (    D1
     :             + ( D2 * 2D0
     :             + ( D3 * 3D0
     :             +   D4 * 4D0 * T ) * T ) * T )

*  Sun's mean anomaly.
      EM = DD2R * MOD ( EM0
     :              + ( EM1
     :              + ( EM2
     :              + ( EM3
     :              +   EM4 * T ) * T ) * T ) * T, 360D0 )
      DEM = DD2R * (    EM1
     :              + ( EM2 * 2D0
     :              + ( EM3 * 3D0
     :              +   EM4 * 4D0 * T ) * T ) * T )

*  Moon's mean anomaly.
      EMP = DD2R * MOD ( EMP0
     :               + ( EMP1
     :               + ( EMP2
     :               + ( EMP3
     :               +   EMP4 * T ) * T ) * T ) * T, 360D0 )
      DEMP = DD2R * (    EMP1
     :               + ( EMP2 * 2D0
     :               + ( EMP3 * 3D0
     :               +   EMP4 * 4D0 * T ) * T ) * T )

*  Mean distance of the Moon from its ascending node.
      F = DD2R * MOD ( F0
     :             + ( F1
     :             + ( F2
     :             + ( F3
     :             +   F4 * T ) * T ) * T ) * T, 360D0 )
      DF = DD2R * (    F1
     :             + ( F2 * 2D0
     :             + ( F3 * 3D0
     :             +   F4 * 4D0 * T ) * T ) * T )

*  Meeus further arguments.
      A1 = DD2R * ( A10 + A11*T )
      DA1 = DD2R * AL1
      A2 = DD2R * ( A20 + A21*T )
      DA2 = DD2R * A21
      A3 = DD2R * ( A30 + A31*T )
      DA3 = DD2R * A31

*  E-factor, and square.
      E = 1D0 + ( E1 + E2*T ) * T
      DE = E1 + 2D0*E2*T
      ESQ = E*E
      DESQ = 2D0*E*DE

*  Use the Meeus additive terms (deg) to start off the summations.
      ELPMF = ELP - F;
      DELPMF = DELP - DF;
      VEL = AL1 * SIN(A1)
     :    + AL2 * SIN(ELPMF)
     :    + AL3 * SIN(A2)
      VDEL = AL1 * COS(A1) * DA1
     :     + AL2 * COS(ELPMF) * DELPMF
     :     + AL3 * COS(A2) * DA2

      VR = 0D0
      VDR = 0D0

      A1MF = A1 - F
      DA1MF = DA1 - DF
      A1PF = A1 + F
      DA1PF = DA1 + DF
      DLPMP = ELP - EMP
      SLPMP = ELP + EMP
      VB = AB1 * SIN(ELP)
     :   + AB2 * SIN(A3)
     :   + AB3 * SIN(A1MF)
     :   + AB4 * SIN(A1PF)
     :   + AB5 * SIN(DLPMP)
     :   + AB6 * SIN(SLPMP)
      VDB = AB1 * COS(ELP) * DELP
     :    + AB2 * COS(A3) * DA3
     :    + AB3 * COS(A1MF) * DA1MF
     :    + AB4 * COS(A1PF) * DA1PF
     :    + AB5 * COS(DLPMP) * (DELP-DEMP)
     :    + AB6 * COS(SLPMP) * (DELP+DEMP)

*  -----------------
*  Series expansions
*  -----------------

*  Longitude and distance plus derivatives.
      DO 1 N=NLR,1,-1
         DN = DBLE( ITLR(1,N) )
         I = ITLR(2,N)
         EMN = DBLE( I )
         EMPN = DBLE( ITLR(3,N) )
         FN = DBLE( ITLR(4,N) )
         I = ABS(I)
         IF ( I .EQ. 1 ) THEN
            EN = E
            DEN = DE
         ELSE IF ( I .EQ. 2 ) THEN
            EN = ESQ
            DEN = DESQ
         ELSE
            EN = 1D0
            DEN = 0D0
         END IF
         ARG = DN*D + EMN*EM + EMPN*EMP + FN*F
         DARG = DN*DD + EMN*DEM + EMPN*DEMP + FN*DF
         FARG = SIN(ARG)
         V = FARG * EN
         DV = COS(ARG)*DARG*EN + FARG*DEN
         COEFF = TLR(1,N)
         VEL = VEL + COEFF*V
         VDEL = VDEL + COEFF*DV
         FARG = COS(ARG)
         V = FARG * EN;
         DV = - SIN(ARG)*DARG*EN + FARG*DEN
         COEFF = TLR(2,N)
         VR = VR + COEFF*V
         VDR = VDR + COEFF*DV
 1    CONTINUE
      EL = ELP + DD2R*VEL
      DEL = ( DELP + DD2R*VDEL ) / DJC
      R = ( VR + R0 ) / DAU
      DR = VDR / DAU / DJC

*  Latitude plus derivative.
      DO 2 N=NB,1,-1
         DN = DBLE( ITB(1,N) )
         I = ITB(2,N)
         EMN = DBLE ( I )
         EMPN = DBLE( ITB(3,N) )
         FN = DBLE( ITB(4,N) )
         I = ABS(I)
         IF ( I .EQ. 1 ) THEN
            EN = E
            DEN = DE
         ELSE IF ( I .EQ. 2 ) THEN
            EN = ESQ
            DEN = DESQ
         ELSE
            EN = 1D0
            DEN = 0D0
         END IF
         ARG = DN*D + EMN*EM + EMPN*EMP + FN*F
         DARG = DN*DD + EMN*DEM + EMPN*DEMP + FN*DF
         FARG = SIN(ARG)
         V = FARG * EN
         DV = COS(ARG)*DARG*EN + FARG*DEN
         COEFF = TB(N)
         VB = VB + COEFF*V
         VDB = VDB + COEFF*DV
 2    CONTINUE
      B = VB * DD2R
      DB = VDB * DD2R / DJC

*  ------------------------------
*  Transformation into final form
*  ------------------------------

*  Longitude, latitude to x, y, z (au).
      CALL iau_S2PV ( EL, B, R, DEL, DB, DR, PV )

*  IAU 2006 Fukushima-Williams bias+precession angles.
      CALL iau_PFW06 ( DATE1, DATE2, GAMB, PHIB, PSIB, EPSA )

*  Mean ecliptic coordinates to GCRS rotation matrix.
      CALL iau_IR ( RM )
      CALL iau_RZ ( PSIB, RM )
      CALL iau_RX ( -PHIB, RM )
      CALL iau_RZ ( -GAMB, RM )

*  Rotate the Moon position and velocity into GCRS (Note 7).
      CALL iau_RXPV ( RM, PV, PV )

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2023
*  Standards of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
*  CONDITIONS WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The names of all routines in your derived work shall not
*        include the prefix "iau" or "sofa" or trivial modifications
*        thereof such as changes of case.
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  5. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  6. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  In any published work or commercial product which uses the SOFA
*  software directly, acknowledgement (see www.iausofa.org) is
*  appreciated.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*-----------------------------------------------------------------------

      END
