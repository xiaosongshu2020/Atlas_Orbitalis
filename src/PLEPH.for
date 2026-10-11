c       program testeph
c C
c C     This program reads a binary format file of ploynomial coefficients
c C     fit to the JPL/Calech Planetary and Lunar Ephemerides,
c C     interpolates the coefficients to a set of input times to
c C     evaluate the positions and velocities of the planets at those
c C     times anc compres with a test printout file.
c C
c C$ Disclaimer
c C
c C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
c C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
c C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
c C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
c C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
c C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
c C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
c C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
c C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
c C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
c C
c C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
c C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
c C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
c C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
c C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
c C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
c C
c C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
c C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
c C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
c C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
c C
c C                Version : March 25, 2013
c C
c C     Program TESTEPH
c C     ---------------
c C
c C     TESTEPH tests the JPL ephemeris reading and interpolating routine using
c C     examples computed from the original ephemeris.
c C
c C     TESTEPH contains the reading and interpolating subroutines that are of
c C     eventual interest to the user.  Once TESTEPH is working correctly, the
c C     user can extract those subroutines and the installation process is complete.
c C
c C     You must supply "testpo.XXX" to TESTEPH, via standard input.  "testpo.XXX
c C     is the specially formatted text file that contains the test cases for the
c C     ephmeris, DEXXX.
c C
c C     After the initial identifying text which is concluded by an "EOT" in
c C     columns 1-3, the test file contains the following quantities:
c C
c C         JPL Ephemeris Number
c C         calendar date
c C         Julian Ephemeris Date
c C         target number (1-Mercury, ...,3-Earth, ,,,9-Pluto, 10-Moon, 11-Sun,
c C                        12-Solar System Barycenter, 13-Earth-Moon Barycenter
c C                        14-Nutations, 15-Librations)
c C         center number (same codes as target number)
c C         coordinate number (1-x, 2-y, ... 6-zdot)
c C         coordinate  [au, au/day]
c C
c C     For each test case input, TESTEPH
c C
c C         - computes the corresponding state from data contained
c C           in DExxx,
c C
c C         - compares the two sets,
c C
c C         - writes an error message if the difference between
c C           any of the state components is greater than 10**(-13).
c C
c C         - writes state and difference information for every 100th
c C           test case processed.
c C
c C      This program is written in standard Fortran-77.
c C
c C      HOWEVER, there are two parts which are compiler dependent; both have
c C      to do with opening and reading a direct-access file.  They are dealt
c C      with in the subroutine FSIZERi, i=1,3.  (There are three versions of
c C      this subroutine.
c C
c C      1) The parameter RECL in the OPEN statement is the number of units per
c C         record.  For some compilers, it is given in bytes; in some, it is given
c C         in single precision words.  In the subroutine FSIZER of TESTEPH, the
c C         parameter NRECL must be set to 4 if RECL is given in bytes; NRECL must
c C         be set to 1 if RECL is given in words.  (If in doubt, use 4 for UNIX;
c C         1 for VAX and PC)
c C
c C      2) Also for the OPEN statement, the program needs to know the exact value
c C         of RECL (number of single precision words times NRECL).  Since this
c C         varies from one JPL ephemeris to another, RECL must be determined somehow
c C         and given to the OPEN statement.  There are three methods, depending
c C         upon the compiler.  We have included three versions of the subroutine
c C         FSIZER, one for each method.
c C
c C         a)  Use the INQUIRE statement to find the length of the records
c C             automatically before opening the file.  This works for VAX's;
c C             not in UNIX.
c C
c C         b)  Open the file with an arbitrary value of RECL, read the first record,
c C             and use the information on that record to determine the exact value
c C             of RECL.  Then, close the file and re-open it with the exact value.
c C             This seems to work for UNIX compilers as long as the initial value of
c C             RECL is less than the exact value but large enough to get the required
c C             information from the first file.  (For other compilers, this doesn't
c C             work since you can open a file only with the exact value of RECL.)
c C
c C         c)  Hardwire the value of RECL.  This number is NRECL*1652 for DE200,
c C             NRECL*2036 for DE405, and NRECL*1456 for DE406.
c C
c        INTEGER NMAX
c        PARAMETER (NMAX = 1000)

c        CHARACTER*6  NAMS(NMAX)
c        CHARACTER*3  ALF3

c        DOUBLE PRECISION  DEL
c        DOUBLE PRECISION  ET
c        DOUBLE PRECISION  R(6)
c        DOUBLE PRECISION  SS(3)
c        DOUBLE PRECISION  VALS(NMAX)
c        DOUBLE PRECISION  XI
c        DOUBLE PRECISION  JDEPOC

c        INTEGER  I
c        INTEGER  LINE
c        DATA LINE/0/
c        INTEGER  NVS,NTARG,NCTR,NCOORD
c        INTEGER  NPT
c        DATA NPT/100/

c C      Write a fingerprint to the screen.

c        WRITE(*,*) ' JPL TEST-EPHEMERIS program. ' //
c      .            ' Last modified March 2013.'

c C      Print the ephemeris constants.

c        CALL  CONST (NAMS, VALS, SS, NVS)

c        JDEPOC = 2440400.5d0

c        WRITE (*,'(/3F14.2)') SS

c        DO I=1,NVS
c          IF(NAMS(I) .EQ. 'JDEPOC')JDEPOC = VALS(I)
c          WRITE(6,'(A8,D24.16)')NAMS(I),VALS(I)
c        ENDDO

c C     Skip the file header comments.

c    1  READ(*,'(a3)')ALF3
c       IF(ALF3 .NE. 'EOT') GO TO 1

c       WRITE(*,*)
c      . '  line -- jed --   t#   c#   x#   --- jpl value ---'//
c      . '   --- user value --    -- difference --'
c       WRITE(*,*)

c C     Read a value from the test case; skip if not within the time-range
c C     of the present version of the ephemeris

c    2  READ(*,'(15X,D10.1,3I3,F30.20)',END=9)
c      . ET,NTARG,NCTR,NCOORD,XI

c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c C
c C  NOTE : Over the years, different versions of PLEPH have had a fifth argument:
c C  sometimes, an error return statement number; sometimes, a logical denoting
c C  whether or not the requested date is covered by the ephemeris.  We apologize
c C  for this inconsistency; in this version, we use only the four necessary
c C  arguments and do the testing outside of the subroutine.

c       IF(ET .LT. SS(1)) GO TO 2
c       IF(ET .GT. SS(2)) GO TO 2

c       CALL  PLEPH ( ET, NTARG, NCTR, R )
c       print*, R

c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c c     For testing purposes DEL is the computation of the difference in
c c       read values of the ephemeris converted from the export ascii files
c c       versus those printed form the original integration.
c c     The difference is checked for random parameters, one every 32 day interval.
c c     The agreement is considered okay if DEl is less that 1e-13.
c c     This corresponds to a few cm for body positions, and very small
c c     values for velocities, and angles and their rates.
c c       (A fractional test isn't suitable since sometimes the values will
c c        be near zero for particular components.)
c c     For lunar libration psi, NCOORD=15, the angle accumulates in time
c c     so precision is lost after several centuries. This is handled below
c c     by scaling DEL down for psi, by 100 radians per year from the reference epoch,
c c     as opposed to having different tolerances for different types of quantities.

c       DEL=DABS(R(NCOORD)-XI)

c       IF(NTARG .EQ. 15 .AND. NCOORD .EQ. 3)THEN
c         DEL = DEL/(1.D0+100.D0*ABS(ET-JDEPOC)/365.25d0)
c       ENDIF

c       LINE = LINE+1

c       IF ( MOD(LINE,npt) .EQ. 0 )  WRITE(*,200)
c      & LINE,ET,NTARG,NCTR,NCOORD,XI,R(NCOORD),DEL
c  200  FORMAT(I6,F10.1,3I5,2F25.13,1E13.5)

c C  Print out WARNING if difference greater than tolerance.

c       IF (DEL .GE. 1.D-13) WRITE(*,201)
c      . LINE,ET,NTARG,NCTR,NCOORD,XI,R(NCOORD),DEL
c  201  FORMAT(/ '*****  WARNING : next difference >= 1.D-13  *****'//
c      . I6,F10.1,3I3,2F25.13,1E13.5/' ')

c       GO TO 2

c    9  STOP
c       END


C++++++++++++++++++++++++
C
      SUBROUTINE FSIZER1(NRECL,KSIZE,NRFILE,NAMFIL)
C
C++++++++++++++++++++++++
C
C   Version 1.0 uses the INQUIRE statement to find out the the record length
C   of the direct access file before opening it.  This procedure is non-standard,
C   but seems to work for VAX machines.
C
C  THE SUBROUTINE ALSO SETS THE VALUES OF  NRECL, NRFILE, AND NAMFIL.

C  *****************************************************************
C  *****************************************************************
C
C  THE PARAMETERS NAMFIL, NRECL, AND NRFILE ARE TO BE SET BY THE USER
C
C  *****************************************************************

C  NAMFIL IS THE EXTERNAL NAME OF THE BINARY EPHEMERIS FILE

      CHARACTER*80  NAMFIL

c      NAMFIL='JPLEPH'

C  *****************************************************************

C  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
C  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES
C  (for a VAX, it is probably 1)
C
      NRECL= 4

C  *****************************************************************

C  NRFILE IS THE INTERNAL UNIT NUMBER USED FOR THE EPHEMERIS FILE

c      NRFILE=12

C  *****************************************************************

C   FIND THE RECORD SIZE USING THE INQUIRE STATEMENT


c      IRECSZ=0

      INQUIRE(FILE=NAMFIL,RECL=IRECSZ)

C IF 'INQUIRE' DOES NOT WORK, USUALLY IRECSZ WILL BE LEFT AT 0

      IF(IRECSZ .LE. 0) write(*,*)
     . ' INQUIRE STATEMENT PROBABLY DID NOT WORK'

      KSIZE=IRECSZ/NRECL

      RETURN

      END
C++++++++++++++++++++++++
C
      SUBROUTINE FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL)
C
C++++++++++++++++++++++++
C  THIS SUBROUTINE OPENS THE FILE, 'NAMFIL', WITH A PHONY RECORD LENGTH, READS
C  THE FIRST RECORD, AND USES THE INFO TO COMPUTE KSIZE, THE NUMBER OF SINGLE
C  PRECISION WORDS IN A RECORD.
C
C  THE SUBROUTINE ALSO SETS THE VALUES OF  NRECL, NRFILE, AND NAMFIL.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      SAVE

      INTEGER OLDMAX
      PARAMETER (OLDMAX = 400)
      INTEGER NMAX
      PARAMETER (NMAX = 1000)

      CHARACTER*6 TTL(14,3),CNAM(NMAX)
      CHARACTER*80 NAMFIL

      DIMENSION SS(3)

      INTEGER IPT(3,13)

C  *****************************************************************
C  *****************************************************************
C
C  THE PARAMETERS NRECL, NRFILE, AND NAMFIL ARE TO BE SET BY THE USER
C
C  *****************************************************************

C  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
C  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES
C  (for UNIX, it is probably 4)
C
      NRECL=4

C  NRFILE IS THE INTERNAL UNIT NUMBER USED FOR THE EPHEMERIS FILE

      NRFILE=12

C  NAMFIL IS THE EXTERNAL NAME OF THE BINARY EPHEMERIS FILE

      NAMFIL='data/AtlasData/JPLEPH_DE430_linux'
C      NAMFIL='de430.bsp'

C  *****************************************************************
C  *****************************************************************

C  **  OPEN THE DIRECT-ACCESS FILE AND GET THE POINTERS IN ORDER TO
C  **  DETERMINE THE SIZE OF THE EPHEMERIS RECORD

      MRECL=NRECL*1000

        OPEN(NRFILE,
     *       FILE=NAMFIL,
     *       ACCESS='DIRECT',
     *       FORM='UNFORMATTED',
     *       RECL=MRECL,
     *       STATUS='OLD')

      READ(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT,
     & ((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3)

      CLOSE(NRFILE)

C  FIND THE NUMBER OF EPHEMERIS COEFFICIENTS FROM THE POINTERS

      KMX = 0
      KHI = 0

      DO I = 1,13
         IF (IPT(1,I) .GE. KMX) THEN
            KMX = IPT(1,I)
            KHI = I
         ENDIF
      ENDDO

      ND = 3
      IF (KHI .EQ. 12) ND=2

      KSIZE = 2*(IPT(1,KHI)+ND*IPT(2,KHI)*IPT(3,KHI)-1)

      RETURN

      END
C++++++++++++++++++++++++

      SUBROUTINE FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)
C
C++++++++++++++++++++++++
C
C  THE SUBROUTINE SETS THE VALUES OF  NRECL, KSIZE, NRFILE, AND NAMFIL.

      SAVE

      CHARACTER*80 NAMFIL

C  *****************************************************************
C  *****************************************************************
C
C  THE PARAMETERS NRECL, NRFILE, AND NAMFIL ARE TO BE SET BY THE USER

C  *****************************************************************

C  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
C  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES

       NRECL=4

C  *****************************************************************

C  NRFILE IS THE INTERNAL UNIT NUMBER USED FOR THE EPHEMERIS FILE (DEFAULT: 12)

      NRFILE=12

C  *****************************************************************

C  Set file name.  The record layout (bytes per record) is derived from the
C  file itself in STATE, so a single file serves every platform -- there is
C  no per-OS variant to choose between.

      NAMFIL = 'data/AtlasData/JPLEPH_DE440'

C  KSIZE must be set by the user according to the ephemeris to be read

C  For  de200, set KSIZE to 1652
C  For  de405, set KSIZE to 2036
C  For  de406, set KSIZE to 1456
C  For  de414, set KSIZE to 2036
C  For  de418, set KSIZE to 2036
C  For  de421, set KSIZE to 2036
C  For  de422, set KSIZE to 2036
C  For  de423, set KSIZE to 2036
C  For  de424, set KSIZE to 2036
C  For  de430, set KSIZE to 2036

      KSIZE = 2036

C  *******************************************************************

      RETURN

      END
C++++++++++++++++++++++++++
C
      SUBROUTINE PLEPH ( ET, NTARG, NCENT, RRD )
C
C++++++++++++++++++++++++++
C  NOTE : Over the years, different versions of PLEPH have had a fifth argument:
C  sometimes, an error return statement number; sometimes, a logical denoting
C  whether or not the requested date is covered by the ephemeris.  We apologize
C  for this inconsistency; in this present version, we use only the four necessary
C  arguments and do the testing outside of the subroutine.
C
C     THIS SUBROUTINE READS THE JPL PLANETARY EPHEMERIS
C     AND GIVES THE POSITION AND VELOCITY OF THE POINT 'NTARG'
C     WITH RESPECT TO 'NCENT'.
C
C     CALLING SEQUENCE PARAMETERS:
C
C       ET = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION
C            IS WANTED.
C
C       ** NOTE THE ENTRY DPLEPH FOR A DOUBLY-DIMENSIONED TIME **
C          THE REASON FOR THIS OPTION IS DISCUSSED IN THE
C          SUBROUTINE STATE
C
C     NTARG = INTEGER NUMBER OF 'TARGET' POINT.
C
C     NCENT = INTEGER NUMBER OF CENTER POINT.
C
C            THE NUMBERING CONVENTION FOR 'NTARG' AND 'NCENT' IS:
C
C                1 = MERCURY           8 = NEPTUNE
C                2 = VENUS             9 = PLUTO
C                3 = EARTH            10 = MOON
C                4 = MARS             11 = SUN
C                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
C                6 = SATURN           13 = EARTH-MOON BARYCENTER
C                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
C                            15 = LIBRATIONS, IF ON EPH FILE
C
C             (IF NUTATIONS ARE WANTED, SET NTARG = 14. FOR LIBRATIONS,
C              SET NTARG = 15. SET NCENT=0.)
C
C      RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY
C            OF POINT 'NTARG' RELATIVE TO 'NCENT'. THE UNITS ARE AU AND
C            AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS
C            PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF
C            RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF
C            RADIANS AND RADIANS/DAY.
C
C            The option is available to have the units in km and km/sec.
C            For this, set km=.true. in the STCOMX common block.
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER NMAX
      PARAMETER (NMAX = 1000)

      DIMENSION RRD(6),ET2Z(2),ET2(2),PV(6,13)
      DIMENSION PVST(6,11),PNUT(4)
      DIMENSION SS(3),CVAL(NMAX),PVSUN(6),ZIPS(2)
      DATA ZIPS/2*0.d0/

      LOGICAL BSAVE,KM,BARY
      LOGICAL FIRST
      DATA FIRST/.TRUE./

      INTEGER LIST(12),IPT(39),DENUM

C     JPL ephemeris status (queried through the entry points PLEPH_STATUS /
C     PLEPH_FILENAME at the end of this file -- F90 callers cannot declare
C     a COMMON, so those routines are the only supported interface):
C       IPHINI  0 = not initialized yet, 1 = initialized OK (sticky),
C               -1 = initialization failed (sticky; restarting is required)
C       IPHSTS  outcome of the most recent call:
C               0 = OK, 1 = file/init failure, 2 = epoch outside the span,
C               3 = record read error
      INTEGER IPHINI,IPHSTS
      COMMON/PLEPHSTS/IPHINI,IPHSTS

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT

      COMMON/STCOMX/KM,BARY,PVSUN

C     INITIALIZE ET2 FOR 'STATE' AND SET UP COMPONENT COUNT
C
      ET2(1)=ET
      ET2(2)=0.D0
      GO TO 11

C     ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT
C          (SEE THE DISCUSSION IN THE SUBROUTINE STATE)

      ENTRY DPLEPH(ET2Z,NTARG,NCENT,RRD)

      ET2(1)=ET2Z(1)
      ET2(2)=ET2Z(2)

  11  IF(FIRST) CALL STATE(ZIPS,LIST,PVST,PNUT)
      FIRST=.FALSE.

C     Nothing can be read from an ephemeris that failed to initialize.
C     Return zeros (never garbage) and let the caller consult PLEPH_STATUS.
      IF(IPHINI .NE. 1) THEN
        DO I=1,6
          RRD(I)=0.D0
        ENDDO
        RETURN
      ENDIF

  96  IF(NTARG .EQ. NCENT) RETURN

      DO I=1,12
        LIST(I)=0
      ENDDO

C     CHECK FOR NUTATION CALL

      IF(NTARG.NE.14) GO TO 97
        IF(IPT(35).GT.0) THEN
          LIST(11)=2
          CALL STATE(ET2,LIST,PVST,PNUT)
          DO I=1,4
            RRD(I)=PNUT(I)
          ENDDO
          RRD(5) = 0.d0
          RRD(6) = 0.d0
          RETURN
        ELSE
          DO I=1,4
            RRD(I)=0.d0
          ENDDO
          WRITE(6,297)
  297     FORMAT(' *****  NO NUTATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF

C     CHECK FOR LIBRATIONS

  97  CONTINUE
      DO I=1,6
        RRD(I)=0.d0
      ENDDO

      IF(NTARG.NE.15) GO TO 98
        IF(IPT(38).GT.0) THEN
          LIST(12)=2
          CALL STATE(ET2,LIST,PVST,PNUT)
          DO I=1,6
            RRD(I)=PVST(I,11)
          ENDDO
          RETURN
        ELSE
          WRITE(6,298)
  298     FORMAT(' *****  NO LIBRATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF

C       FORCE BARYCENTRIC OUTPUT BY 'STATE'

  98  BSAVE=BARY
      BARY=.TRUE.

C       SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL

      DO I=1,2
        K=NTARG
        IF(I .EQ. 2) K=NCENT
        IF(K .LE. 10) LIST(K)=2
        IF(K .EQ. 10) LIST(3)=2
        IF(K .EQ. 3) LIST(10)=2
        IF(K .EQ. 13) LIST(3)=2
      ENDDO

C       MAKE CALL TO STATE

      CALL STATE(ET2,LIST,PVST,PNUT)

C     A per-call failure (epoch outside the ephemeris span, record read
C     error) leaves PVST/PVSUN meaningless: return zeros, not garbage.
C     BARY was forced to .TRUE. for this call, so restore it first.
      IF(IPHSTS .NE. 0) THEN
        BARY=BSAVE
        DO I=1,6
          RRD(I)=0.D0
        ENDDO
        RETURN
      ENDIF

      DO I=1,10
        DO J = 1,6
          PV(J,I) = PVST(J,I)
        ENDDO
      ENDDO

      IF(NTARG .EQ. 11 .OR. NCENT .EQ. 11) THEN
      DO I=1,6
        PV(I,11)=PVSUN(I)
      ENDDO
      ENDIF

      IF(NTARG .EQ. 12 .OR. NCENT .EQ. 12) THEN
        DO I=1,6
          PV(I,12)=0.D0
        ENDDO
      ENDIF

      IF(NTARG .EQ. 13 .OR. NCENT .EQ. 13) THEN
        DO I=1,6
          PV(I,13) = PVST(I,3)
        ENDDO
      ENDIF

      IF(NTARG*NCENT .EQ. 30 .AND. NTARG+NCENT .EQ. 13) THEN
        DO I=1,6
          PV(I,3)=0.D0
        ENDDO
        GO TO 99
      ENDIF

      IF(LIST(3) .EQ. 2) THEN
        DO I=1,6
          PV(I,3)=PVST(I,3)-PVST(I,10)/(1.D0+EMRAT)
        ENDDO
      ENDIF

      IF(LIST(10) .EQ. 2) THEN
        DO I=1,6
          PV(I,10) = PV(I,3)+PVST(I,10)
        ENDDO
      ENDIF

  99  DO I=1,6
        RRD(I)=PV(I,NTARG)-PV(I,NCENT)
      ENDDO

      BARY=BSAVE

      RETURN
      END
C+++++++++++++++++++++++++++++++++
C
      SUBROUTINE INTERP(BUF,T,NCF,NCM,NA,IFL,PV)
C
C+++++++++++++++++++++++++++++++++
C
C     THIS SUBROUTINE DIFFERENTIATES AND INTERPOLATES A
C     SET OF CHEBYSHEV COEFFICIENTS TO GIVE POSITION AND VELOCITY
C
C     CALLING SEQUENCE PARAMETERS:
C
C       INPUT:
C
C         BUF   1ST LOCATION OF ARRAY OF D.P. CHEBYSHEV COEFFICIENTS OF POSITION
C
C           T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY
C               COEFFICIENTS AT WHICH INTERPOLATION IS WANTED
C               (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE
C               INTERVAL IN INPUT TIME UNITS.
C
C         NCF   # OF COEFFICIENTS PER COMPONENT
C
C         NCM   # OF COMPONENTS PER SET OF COEFFICIENTS
C
C          NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY
C               (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)
C
C          IFL  INTEGER FLAG: =1 FOR POSITIONS ONLY
C                             =2 FOR POS AND VEL
C
C
C       OUTPUT:
C
C         PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION
C               EXPECTED IS PV(NCM,IFL), DP.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      SAVE
C
      DOUBLE PRECISION BUF(NCF,NCM,*),T(2),PV(NCM,*),PC(18),VC(18)

C
      DATA NP/2/
      DATA NV/3/
      DATA TWOT/0.D0/
      DATA PC(1),PC(2)/1.D0,0.D0/
      DATA VC(2)/1.D0/
C
C       ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET
C       OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME
C       WITHIN THAT SUBINTERVAL.
C
      DNA=DBLE(NA)
      DT1=DINT(T(1))
      TEMP=DNA*T(1)
      L=IDINT(TEMP-DT1)+1

C         TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1)

      TC=2.D0*(DMOD(TEMP,1.D0)+DT1)-1.D0

C       CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,
C       AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.
C       (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE
C       CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)

      IF(TC.NE.PC(2)) THEN
        NP=2
        NV=3
        PC(2)=TC
        TWOT=TC+TC
      ENDIF
C
C       BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED
C       AND ARE STORED IN THE ARRAY 'PC'.
C
      IF(NP.LT.NCF) THEN
        DO 1 I=NP+1,NCF
        PC(I)=TWOT*PC(I-1)-PC(I-2)
    1   CONTINUE
        NP=NCF
      ENDIF
C
C       INTERPOLATE TO GET POSITION FOR EACH COMPONENT
C
      DO 2 I=1,NCM
      PV(I,1)=0.D0
      DO 3 J=NCF,1,-1
      PV(I,1)=PV(I,1)+PC(J)*BUF(J,I,L)
    3 CONTINUE
    2 CONTINUE
      IF(IFL.LE.1) RETURN
C
C       IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH
C       DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.
C
      VFAC=(DNA+DNA)/T(2)
      VC(3)=TWOT+TWOT
      IF(NV.LT.NCF) THEN
        DO 4 I=NV+1,NCF
        VC(I)=TWOT*VC(I-1)+PC(I-1)+PC(I-1)-VC(I-2)
    4   CONTINUE
        NV=NCF
      ENDIF
C
C       INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT
C
      DO 5 I=1,NCM
      PV(I,2)=0.D0
      DO 6 J=NCF,2,-1
      PV(I,2)=PV(I,2)+VC(J)*BUF(J,I,L)
    6 CONTINUE
      PV(I,2)=PV(I,2)*VFAC
    5 CONTINUE
C
      RETURN
C
      END

C+++++++++++++++++++++++++
C
      SUBROUTINE SPLIT_DP(TT,FR)
C
C+++++++++++++++++++++++++
C
C     NOTE: named SPLIT_DP rather than SPLIT.  "SPLIT" is a non-standard
C     Intel Fortran intrinsic subroutine -- CALL SPLIT(string,set,pos[,back])
C     -- so this routine's two-argument call matched a different interface
C     and Intel emitted warning #7532 ("assume 'external'").  Renaming
C     removes the ambiguity, and with it the hazard of a future call with
C     three arguments being silently resolved to the intrinsic.
C
C     THIS SUBROUTINE BREAKS A D.P. NUMBER INTO A D.P. INTEGER
C     AND A D.P. FRACTIONAL PART.
C
C     CALLING SEQUENCE PARAMETERS:
C
C       TT = D.P. INPUT NUMBER
C
C       FR = D.P. 2-WORD OUTPUT ARRAY.
C            FR(1) CONTAINS INTEGER PART
C            FR(2) CONTAINS FRACTIONAL PART
C
C            FOR NEGATIVE INPUT NUMBERS, FR(1) CONTAINS THE NEXT
C            MORE NEGATIVE INTEGER; FR(2) CONTAINS A POSITIVE FRACTION.
C
C       CALLING SEQUENCE DECLARATIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION FR(2)

C       MAIN ENTRY -- GET INTEGER AND FRACTIONAL PARTS

      FR(1)=DINT(TT)
      FR(2)=TT-FR(1)

      IF(TT.GE.0.D0 .OR. FR(2).EQ.0.D0) RETURN

C       MAKE ADJUSTMENTS FOR NEGATIVE INPUT NUMBER

      FR(1)=FR(1)-1.D0
      FR(2)=FR(2)+1.D0

      RETURN

      END


C++++++++++++++++++++++++++++++++
C
      SUBROUTINE STATE(ET2,LIST,PV,PNUT)
C
C++++++++++++++++++++++++++++++++
C
C THIS SUBROUTINE READS AND INTERPOLATES THE JPL PLANETARY EPHEMERIS FILE
C
C     CALLING SEQUENCE PARAMETERS:
C
C     INPUT:
C
C         ET2   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION
C               IS WANTED.  ANY COMBINATION OF ET2(1)+ET2(2) WHICH FALLS
C               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.
C
C                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE
C                   ENTIRE EPOCH IN ET2(1) AND SET ET2(2)=0.
C
C                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET ET2(1) =
C                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION
C                   EPOCH AND SET ET2(2) = FRACTIONAL PART OF A DAY
C                   ELAPSED BETWEEN ET2(1) AND EPOCH.
C
C                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET
C                   ET2(1) = SOME FIXED EPOCH, SUCH AS START OF INTEGRATION,
C                   AND ET2(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH.
C
C        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION
C               IS WANTED FOR EACH OF THE BODIES ON THE FILE.
C
C                         LIST(I)=0, NO INTERPOLATION FOR BODY I
C                                =1, POSITION ONLY
C                                =2, POSITION AND VELOCITY
C
C               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:
C
C                         I = 1: MERCURY
C                           = 2: VENUS
C                           = 3: EARTH-MOON BARYCENTER
C                           = 4: MARS
C                           = 5: JUPITER
C                           = 6: SATURN
C                           = 7: URANUS
C                           = 8: NEPTUNE
C                           = 9: PLUTO
C                           =10: GEOCENTRIC MOON
C                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY
C                           =12: LUNAR LIBRATIONS (IF ON FILE)
C
C     OUTPUT:
C
C          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED
C               QUANTITIES (OTHER THAN NUTATION, STOERD IN PNUT).
C               THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS
C               STATE IN THE ARRAY STARTING AT PV(1,I).
C               (ON ANY GIVEN CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE
C                AFFECTED BY THE  FIRST 10 'LIST' ENTRIES, AND BY LIST(12)
C                IF LIBRATIONS ARE ON THE FILE, ARE SET.
C                THE REST OF THE 'PV' ARRAYIS UNTOUCHED.)
C               THE ORDER OF COMPONENTS STARTING IN PV(1,I) IS: X,Y,Z,DX,DY,DZ.
C
C               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN
C               EQUATOR AND EQUINOX OF J2000 IF THE DE NUMBER IS 200 OR
C               GREATER; OF B1950 IF THE DE NUMBER IS LESS THAN 200.
C
C               THE MOON STATE IS ALWAYS GEOCENTRIC; THE OTHER NINE STATES
C               ARE EITHER HELIOCENTRIC OR SOLAR-SYSTEM BARYCENTRIC,
C               DEPENDING ON THE SETTING OF COMMON FLAGS (SEE BELOW).
C
C               LUNAR LIBRATIONS, IF ON FILE, ARE PUT INTO PV(K,11) IF
C               LIST(12) IS 1 OR 2.
C
C         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,
C               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF
C               QUANTITIES IN NUT IS:
C
C                        D PSI  (NUTATION IN LONGITUDE)
C                        D EPSILON (NUTATION IN OBLIQUITY)
C                        D PSI DOT
C                        D EPSILON DOT
C
C           *   STATEMENT # FOR ERROR RETURN, IN CASE OF EPOCH OUT OF
C               RANGE OR I/O ERRORS.
C
C     COMMON AREA STCOMX:
C
C          KM   LOGICAL FLAG DEFINING PHYSICAL UNITS OF THE OUTPUT
C               STATES. KM = .TRUE., KM AND KM/SEC
C                          = .FALSE., AU AND AU/DAY
C               DEFAULT VALUE = .FALSE.  (KM DETERMINES TIME UNIT
C               FOR NUTATIONS AND LIBRATIONS.  ANGLE UNIT IS ALWAYS RADIANS.)
C
C        BARY   LOGICAL FLAG DEFINING OUTPUT CENTER.
C               ONLY THE 9 PLANETS ARE AFFECTED.
C                        BARY = .TRUE. =\ CENTER IS SOLAR-SYSTEM BARYCENTER
C                             = .FALSE. =\ CENTER IS SUN
C               DEFAULT VALUE = .FALSE.
C
C       PVSUN   DP 6-WORD ARRAY CONTAINING THE BARYCENTRIC POSITION AND
C               VELOCITY OF THE SUN.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      SAVE

      INTEGER OLDMAX
      PARAMETER ( OLDMAX = 400)
      INTEGER NMAX
      PARAMETER ( NMAX = 1000)

      DIMENSION ET2(2),PV(6,11),PNUT(4),T(2),PJD(4),BUF(1500),
     . SS(3),CVAL(NMAX),PVSUN(6)

      INTEGER LIST(12),IPT(3,13)
      INTEGER NREC
      INTEGER*8 FSIZE,RECBYTES
      INTEGER IOS

C     JPL ephemeris status -- must be declared exactly as in PLEPH, CONST
C     and PLEPH_STATUS (same name, same order, same type)
      INTEGER IPHINI,IPHSTS

      LOGICAL FIRST
      DATA FIRST/.TRUE./

      CHARACTER*6 TTL(14,3),CNAM(NMAX)
      CHARACTER*80 NAMFIL

      LOGICAL KM,BARY

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,NUMDE,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL
      COMMON/STCOMX/KM,BARY,PVSUN
      COMMON/PLEPHSTS/IPHINI,IPHSTS

C
C       ENTRY POINT - 1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE
C
      IF(FIRST) THEN
        FIRST=.FALSE.

C       Initialization is attempted exactly once: IPHINI starts at 0, becomes
C       1 on success or -1 on failure, and a failure is sticky (the file will
C       not appear while the program runs).
        IPHINI=0

C ************************************************************************
C ************************************************************************

C THE USER MUST SELECT ONE OF THE FOLLOWING BY DELETING THE 'C' IN COLUMN 1

C ************************************************************************

C        CALL FSIZER1(NRECL,KSIZE,NRFILE,NAMFIL)
C        CALL FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL)
        CALL FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)

      IF(NRECL .EQ. 0) WRITE(*,*)'  ***** FSIZER IS NOT WORKING *****'

C ************************************************************************
C ************************************************************************

C     Read with STREAM I/O instead of ACCESS='DIRECT' + RECL.  RECL units
C     differ between compilers (gfortran counts bytes, Intel counts 4-byte
C     words by default), and a wrong record length usually makes READ
C     succeed while returning all zeros -- which then becomes NaN forces.
C     STREAM addresses by byte offset, with the record length derived from
C     the file itself, so both layouts of DE440 work on both compilers:
C       record 1 always starts at byte 0 -> read header for SS(1:3)
C       records = (SS(2)-SS(1))/SS(3) + 2   (header + constant table)
C       bytes/record = file size / records
C     NOTE: NRECL is unused from here on (FSIZER3 still returns it).

      NCOEFFS=KSIZE/2

        OPEN(NRFILE,
     *       FILE=NAMFIL,
     *       ACCESS='STREAM',
     *       FORM='UNFORMATTED',
     *       STATUS='OLD',IOSTAT=IOS)

C     IOSTAT= on every I/O below keeps a missing or corrupt file from
C     aborting the process inside the runtime library: the failure is
C     reported through IPHINI/IPHSTS instead (see the comment in PLEPH).
      IF(IOS .NE. 0) THEN
        WRITE(0,*)' ***** CANNOT OPEN THE EPHEMERIS FILE *****'
        WRITE(0,*)'   file=',NAMFIL
        IPHINI=-1
        IPHSTS=1
        RETURN
      ENDIF

      READ(NRFILE,POS=1,IOSTAT=IOS)TTL,
     & (CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT,
     & ((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3)
     & ,(CNAM(L),L=OLDMAX+1,NCON)
      IF(IOS .NE. 0) THEN
        WRITE(0,*)' ***** CANNOT READ THE EPHEMERIS HEADER *****'
        WRITE(0,*)'   file=',NAMFIL
        IPHINI=-1
        IPHSTS=1
        RETURN
      ENDIF

      INQUIRE(UNIT=NRFILE,SIZE=FSIZE)
      NREC = IDINT((SS(2)-SS(1))/SS(3)) + 2
      IF(NREC .LE. 2 .OR. FSIZE .LE. 0) THEN
        WRITE(0,*)' ***** EPHEMERIS HEADER HAS BAD SPAN *****'
        WRITE(0,*)'   file=',NAMFIL
        WRITE(0,*)'   SS(1:3)=',SS
        IPHINI=-1
        IPHSTS=1
        RETURN
      ENDIF
      RECBYTES = FSIZE/NREC
      IF(RECBYTES .LT. 8*NCOEFFS) THEN
        WRITE(0,*)' ***** EPHEMERIS RECORD TOO SHORT *****'
        WRITE(0,*)'   file=',NAMFIL
        WRITE(0,*)'   bytes/record=',RECBYTES,' need>=',8*NCOEFFS
        IPHINI=-1
        IPHSTS=1
        RETURN
      ENDIF
      IF(RECBYTES*NREC .NE. FSIZE) THEN
        WRITE(0,*)' warning @STATE: file size is not a multiple of'
        WRITE(0,*)'   the record count; using ',RECBYTES,' B/record'
      ENDIF

      IF(NCON .LE. OLDMAX)THEN
        READ(NRFILE,POS=RECBYTES+1,IOSTAT=IOS)(CVAL(I),I=1,OLDMAX)
      ELSE
        READ(NRFILE,POS=RECBYTES+1,IOSTAT=IOS)(CVAL(I),I=1,NCON)
      ENDIF
      IF(IOS .NE. 0) THEN
        WRITE(0,*)' ***** CANNOT READ THE EPHEMERIS CONSTANTS *****'
        WRITE(0,*)'   file=',NAMFIL
        IPHINI=-1
        IPHSTS=1
        RETURN
      ENDIF

      NRL=0

      IPHINI=1

      ENDIF

C     Guard: never reach the main body with an uninitialized SS/IPT/NRFILE.
C     The guard comes BEFORE the IPHSTS reset, so a sticky init failure is
C     not cleared by later calls.
      IF(IPHINI .NE. 1) RETURN
      IPHSTS=0

C       ********** MAIN ENTRY POINT **********

      IF(ET2(1) .EQ. 0.D0) RETURN

      S=ET2(1)-.5D0
      CALL SPLIT_DP(S,PJD(1))
      CALL SPLIT_DP(ET2(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)+.5D0
      PJD(2)=PJD(2)+PJD(4)
      CALL SPLIT_DP(PJD(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)

C       ERROR RETURN FOR EPOCH OUT OF RANGE

      IF(PJD(1)+PJD(4).LT.SS(1) .OR. PJD(1)+PJD(4).GT.SS(2)) GO TO 98

C       CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL

      NR=IDINT((PJD(1)-SS(1))/SS(3))+3
      IF(PJD(1).EQ.SS(2)) NR=NR-1

        tmp1 = DBLE(NR-3)*SS(3) + SS(1)
        tmp2 = PJD(1) - tmp1
        T(1) = (tmp2 + PJD(4))/SS(3)

C       READ CORRECT RECORD IF NOT IN CORE

      IF(NR.NE.NRL) THEN
        NRL=NR
        READ(NRFILE,POS=(NR-1)*RECBYTES+1,ERR=99)(BUF(K),K=1,NCOEFFS)
      ENDIF

      IF(KM) THEN
      T(2)=SS(3)*86400.D0
      AUFAC=1.D0
      ELSE
      T(2)=SS(3)
      AUFAC=1.D0/AU
      ENDIF

C   INTERPOLATE SSBARY SUN

      CALL INTERP(BUF(IPT(1,11)),T,IPT(2,11),3,IPT(3,11),2,PVSUN)

      DO I=1,6
      PVSUN(I)=PVSUN(I)*AUFAC
      ENDDO

C   CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED

      DO 4 I=1,10
      IF(LIST(I).EQ.0) GO TO 4

      CALL INTERP(BUF(IPT(1,I)),T,IPT(2,I),3,IPT(3,I),
     & LIST(I),PV(1,I))

      DO J=1,6
       IF(I.LE.9 .AND. .NOT.BARY) THEN
       PV(J,I)=PV(J,I)*AUFAC-PVSUN(J)
       ELSE
       PV(J,I)=PV(J,I)*AUFAC
       ENDIF
      ENDDO

   4  CONTINUE

C       DO NUTATIONS IF REQUESTED (AND IF ON FILE)

      IF(LIST(11).GT.0 .AND. IPT(2,12).GT.0)
     * CALL INTERP(BUF(IPT(1,12)),T,IPT(2,12),2,IPT(3,12),
     * LIST(11),PNUT)

C       GET LIBRATIONS IF REQUESTED (AND IF ON FILE)

      IF(LIST(12).GT.0 .AND. IPT(2,13).GT.0)
     * CALL INTERP(BUF(IPT(1,13)),T,IPT(2,13),3,IPT(3,13),
     * LIST(12),PV(1,11))

      RETURN

  98  WRITE(0,*)' *** EPOCH OUTSIDE THE EPHEMERIS SPAN ***'
      WRITE(0,*)'   requested JED=',ET2(1)+ET2(2)
      WRITE(0,*)'   ephemeris span=',SS(1),SS(2)
      IPHSTS=2
      RETURN

   99 WRITE(0,*)' ***** ERROR READING AN EPHEMERIS RECORD *****'
      WRITE(0,*)'   JED=',ET2
      IPHSTS=3
      RETURN

      END
C+++++++++++++++++++++++++++++
C
      SUBROUTINE CONST(NAM,VAL,SSS,N)
C
C+++++++++++++++++++++++++++++
C
C     THIS ENTRY OBTAINS THE CONSTANTS FROM THE EPHEMERIS FILE
C
C     CALLING SEQEUNCE PARAMETERS (ALL OUTPUT):
C
C       NAM = CHARACTER*6 ARRAY OF CONSTANT NAMES
C
C       VAL = D.P. ARRAY OF VALUES OF CONSTANTS
C
C       SSS = D.P. JD START, JD STOP, STEP OF EPHEMERIS
C
C         N = INTEGER NUMBER OF ENTRIES IN 'NAM' AND 'VAL' ARRAYS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      SAVE

      INTEGER NMAX
      PARAMETER (NMAX = 1000)

      CHARACTER*6 NAM(*),TTL(14,3),CNAM(NMAX)

      DOUBLE PRECISION VAL(*),SSS(3),SS(3),CVAL(NMAX),ZIPS(2)
      DOUBLE PRECISION PVST(6,11),PNUT(4)
      DATA ZIPS/2*0.d0/

      INTEGER IPT(3,13),DENUM,LIST(12)
      logical first
      data first/.true./

C     JPL ephemeris status -- same layout as in PLEPH / STATE
      INTEGER IPHINI,IPHSTS

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL
      COMMON/PLEPHSTS/IPHINI,IPHSTS

C  CALL STATE TO INITIALIZE THE EPHEMERIS AND READ IN THE CONSTANTS

      IF(FIRST) CALL STATE(ZIPS,LIST,PVST,PNUT)
      first=.false.

C  An ephemeris that failed to initialize has no constants to report
      IF(IPHINI .NE. 1) THEN
        N=0
        RETURN
      ENDIF

      N=NCON

      DO I=1,3
        SSS(I)=SS(I)
      ENDDO

      DO I=1,N
        NAM(I)=CNAM(I)
        VAL(I)=CVAL(I)
      ENDDO

      RETURN

      END
C+++++++++++++++++++++++++++++
C
      SUBROUTINE PLEPH_STATUS(ISTAT)
C
C+++++++++++++++++++++++++++++
C
C     Outcome of the most recent PLEPH call, for callers that must decide
C     what to do when the JPL ephemeris is unavailable.  This (and
C     PLEPH_FILENAME below) is the only supported interface to the status
C     kept in COMMON/PLEPHSTS, since F90 callers cannot declare a COMMON.
C
C       ISTAT = 0  OK
C               1  the ephemeris file could not be opened or initialized
C                  (sticky: further PLEPH calls return zeros at once)
C               2  the requested epoch is outside the ephemeris time span
C               3  an ephemeris record could not be read
C
C     NOTE: a call with NTARG = NCENT returns before STATE is reached, so it
C     leaves the status of the previous call standing.
C
      INTEGER ISTAT
      INTEGER IPHINI,IPHSTS
      COMMON/PLEPHSTS/IPHINI,IPHSTS

      IF(IPHINI .EQ. 1) THEN
        ISTAT=IPHSTS
      ELSE
        ISTAT=1
      ENDIF

      RETURN
      END
C+++++++++++++++++++++++++++++
C
      SUBROUTINE PLEPH_FILENAME(NAM)
C
C+++++++++++++++++++++++++++++
C
C     Name of the JPL ephemeris file used by PLEPH.  FSIZER3 owns it, so
C     ask rather than duplicating it in the caller.
C
      CHARACTER*80 NAM
      INTEGER NRECL,KSIZE,NRFILE

      CALL FSIZER3(NRECL,KSIZE,NRFILE,NAM)

      RETURN
      END
