!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   GMS :    Time integration of a scalar nonlinear parabolic PDE
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,W,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: W(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION U,V,KAP,TAU

        U = W(1)
        V = W(2)
        TAU = PAR(1)
        KAP = PAR(2)

!      *Set the nonlinear term
        F(1)= -U + U**2/(V*(1+KAP*U**2))
        F(2)= (-V + U**2)/TAU

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,W,PAR,Z)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: W(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: Z

      DOUBLE PRECISION pi,KAP,TAU,Du,Dv,D,EP,L,Tempa, Tempb

      pi=4*DATAN(1.d0)

!      *Set the (constant) parameters
       KAP  = 0.25
       TAU  = 1.
	   EP = 0.02
	   D = 1.
	   L = 1.8
       Du = (EP)**2
	   Dv = D/(TAU)

       PAR(1) = TAU
       PAR(2) = KAP

!      *Set the actual width of the space interval [0,PAR(11)]
       PAR(11) = L

!      *Set the initial data in the (scaled) interval [0,1]
	   Tempa = (exp(50*(Z-.25))-1)/(exp(50*(Z-.25))+1)
	   Tempb = (exp(50*(.75-Z))-1)/(exp(50*(.75-Z))+1)
	   W(1) = 0.346*(Tempa + Tempb)
	   W(2) = 0.287 - 0.03*cos(2*pi*Z)


!      *Also set the space-derivative of the initial data
!      *Note the scaling by PAR(11)
		W(3) = 0.346*(Tempb**2 - Tempa**2)*25/L
		W(4) = 0.06*pi*sin(2*pi*Z)/L


!      *Set the diffusion constants
       PAR(15) = Du
       PAR(16) = Dv

      END SUBROUTINE STPNT

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,W0,W1,FB,IJAC,DBC) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), W0(NDIM), W1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

!      *Define the boundary conditions (Neumann, in this example).
       FB(1)=W0(3)
       FB(2)=W0(4) 
       FB(3)=W1(3) 
       FB(4)=W1(4) 

      END SUBROUTINE BCND

      DOUBLE PRECISION FUNCTION GETU2(U,NDX,NTST,NCOL)
!     ------ --------- -------- -----
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)

        GETU2 = U(2,0)

      END FUNCTION GETU2

      SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST
!---------------------------------------------------------------------- 
! NOTE : 
! Parameters set in this subroutine should be considered as ``solution 
! measures'' and be used for output purposes only.
! 
! They should never be used as `true'' continuation parameters. 
!
! They may, however, be added as ``over-specified parameters'' in the 
! parameter list associated with the AUTO-Constant NICP, in order to 
! print their values on the screen and in the ``p.xxx file.
!
! They may also appear in the list associated with AUTO-constant NUZR.
!
!---------------------------------------------------------------------- 
! For algebraic problems the argument U is, as usual, the state vector.
! For differential equations the argument U represents the approximate 
! solution on the entire interval [0,1]. In this case its values can
! be accessed indirectly by calls to GETP, as illustrated below, or
! by obtaining NDIM, NCOL, NTST via GETP and then dimensioning U as
! U(NDIM,0:NCOL*NTST) in a seperate subroutine that is called by PVLS.
!---------------------------------------------------------------------- 

! Set PAR(3) equal to the L2-norm of U(1)
       PAR(3)=PAR(11)

! Set PAR(4) equal to the minimum of U(2)
       PAR(4)=GETP('MIN',2,U)

! Set PAR(5) equal to the value of U(2) at the left boundary.
       PAR(5)=GETP('BV0',2,U)

! Set PAR(6) equal to the L2-norm of U(2)
       PAR(6)=GETP('NRM',2,U)
!---------------------------------------------------------------------- 
! The first argument of GETP may be one of the following:
!        'NRM' (L2-norm),     'MAX' (maximum),
!        'INT' (integral),    'BV0 (left boundary value),
!        'MIN' (minimum),     'BV1' (right boundary value).
!        'MNT' (t value for minimum)
!        'MXT' (t value for maximum)
!        'NDIM', 'NDX' (effective (active) number of dimensions)
!        'NTST' (NTST from constant file)
!        'NCOL' (NCOL from constant file)
!        'NBC'  (active NBC)
!        'NINT' (active NINT)
!        'DTM'  (delta t for all t values, I=1...NTST)
!        'WINT' (integration weights used for interpolation, I=0...NCOL)
!
! Also available are
!   'STP' (Pseudo-arclength step size used).
!   'FLD' (`Fold function', which vanishes at folds).
!   'BIF' (`Bifurcation function', which vanishes at singular points).
!   'HBF' (`Hopf function'; which vanishes at Hopf points).
!   'SPB' ( Function which vanishes at secondary periodic bifurcations).
!   'EIG' ( Eigenvalues/multipliers, I=1...2*NDIM, alternates real/imag parts).
!   'STA' ( Number of stable eigenvalues/multipliers).
!---------------------------------------------------------------------- 

      END SUBROUTINE PVLS

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
