!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   GRAYSCOTT :    Time integration of a scalar nonlinear parabolic PDE
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
       KAP  = 2.5
       TAU  = 1.
	   EP = 0.02
	   D = 1.
	   L = 2.
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

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
