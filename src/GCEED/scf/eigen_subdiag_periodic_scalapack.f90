!
!  Copyright 2018 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
subroutine eigen_subdiag_periodic(Rmat,evec,iter,ier2)
  use salmon_parallel, only: nproc_size_global
  use scf_data
  implicit none
  
  integer :: iter,ier2
  complex(8) :: Rmat(iter,iter)
  complex(8) :: evec(iter,iter)
  
  character(1) :: JOBZ,UPLO
  real(8) :: W(iter)
  integer :: LWORK
  
  real(8),allocatable :: RWORK(:)
  complex(8),allocatable :: WORK(:)
  
  if(nproc_size_global==1.or.nproc_k>=2)then
    ier2=0
    JOBZ='V'
    UPLO='U'
    LWORK=2*iter-1
    allocate(WORK(LWORK))
    allocate(RWORK(3*iter-2))

    call ZHEEV(JOBZ,UPLO,iter,Rmat,iter,W,WORK,LWORK,RWORK,ier2)
    evec(:,:)=Rmat(:,:)
  
    deallocate(WORK,RWORK)
  else
    call SAMPLE_PZHEEVD_CALL(Rmat,evec,iter)
  end if
  
end subroutine eigen_subdiag_periodic
   
!
      subroutine SAMPLE_PZHEEVD_CALL(Rmat,evec,iter)
!
!
!  -- ScaLAPACK routine (version 1.2) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 10, 1996
!
!     This routine contains a sample call to PDSYEV.
!     When compiled and run, it produces output which can be
!     pasted directly into matlab.
!
!     .. Parameters ..
use salmon_parallel, only: nproc_size_global
use scf_data
integer :: iter
complex(8) :: Rmat(iter,iter)
complex(8) :: evec(iter,iter)
      INTEGER            LWORK, MAXN
      INTEGER            LIWORK
      INTEGER,allocatable :: IWORK(:)
      DOUBLE PRECISION   ZERO
!      PARAMETER          ( LWORK = 264, MAXN = 100, ZERO = 0.0D+0 )
      PARAMETER          ( MAXN = 20000, ZERO = 0.0D+0 )
      INTEGER            LDA
      DOUBLE PRECISION   MONE
      INTEGER            MAXPROCS
!      PARAMETER          ( LDA = MAXN, MONE = -1.0D+0, MAXPROCS = 512 )
      PARAMETER          ( MONE = -1.0D+0, MAXPROCS = 2048 )
      integer :: ii,jj
      character(1) :: SCOPE,TOP
!     ..
!     .. Local Scalars ..
      INTEGER            INFO, N, NB 
!      INTEGER            CONTEXT, I, IAM, INFO, MYCOL, MYROW, N, NB, NPCOL, NPROCS, NPROW

!     ..
!     .. Local Arrays ..
!      INTEGER            DESCA( 50 ), DESCZ( 50 )
!      DOUBLE PRECISION   A( LDA, LDA ), W( MAXN ), WORK( LWORK ), Z( LDA, LDA )
      complex(8)         A( iter, iter ), Z( iter, iter )
      DOUBLE PRECISION   W( iter )
      complex(8), allocatable :: WORK( : )
      real(8), allocatable :: RWORK(:)
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,   &
                         BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,&
                         BLACS_SETUP, DESCINIT, PDLAPRNT,&
                         PDSYEV, NUMROC, INDXG2P
      INTEGER :: NP,NQ
      INTEGER :: NP0,MQ0
      INTEGER :: IAROW, IACOL
      INTEGER :: RSRC_, CSRC_
      PARAMETER  ( RSRC_ = 7, CSRC_ = 8 )
      INTEGER :: RSRC_A, CSRC_A
      INTEGER :: LRWORK
      INTEGER :: NUMROC, INDXG2P
!     ..
!     .. Executable Statements ..
!
!
!     Set up the problem
!
!      N = 4
      N = iter
      NB = 1
!      NPROW = 2
!      NPCOL = 2
      if(nproc_Mxin(1)>1)then
        NPROW = nproc_Mxin(1)
      else if(nproc_Mxin(2)>1)then
        NPROW = nproc_Mxin(2)
      else
        NPROW = nproc_Mxin(3)
      end if
      NPCOL = nproc_size_global/NPROW
      LDA = iter

      RSRC_A = DESCA( RSRC_ )
      CSRC_A = DESCA( CSRC_ )
 
      NP0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPROW )
      MQ0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )

      LWORK=N + ( NP0 + MQ0 + NB ) * NB

      IAROW = INDXG2P( 1, NB, MYROW, RSRC_A, NPROW )
      IACOL = INDXG2P( 1, NB, MYCOL, CSRC_A, NPCOL )
      NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
      NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
      LRWORK=1 + 9*N + 3*NP*NQ

      LIWORK = 7*N + 8*NPCOL + 2

      allocate(WORK(LWORK))
      allocate(RWORK(LRWORK))
      allocate(IWORK(LIWORK))
!
!
      if(iblacsinit==0)then
!     Initialize the BLACS
!
        CALL BLACS_PINFO( IAM, NPROCS2 )
        IF( ( NPROCS2.LT.1 ) ) THEN
           CALL BLACS_SETUP( IAM, NPROW*NPCOL )
        END IF
!
!
!     Initialize a single BLACS context
!
        CALL BLACS_GET( -1, 0, CONTEXT )
        CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
        CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
!
      end if
!     Bail out if this process is not a part of this context.
!
      IF( MYROW.EQ.-1 ) GO TO 20
!
!
!     These are basic array descriptors
!
      if(iblacsinit==0)then
        CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
        CALL DESCINIT( DESCZ, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
!        CALL BLACS_GRIDINFO( DESCA( 2 ), NPROW, NPCOL, MYROW, MYCOL )
      end if
!
!     Build a matrix that you can create with
!     a one line matlab command:  hilb(n) + diag([1:-1/n:1/n])
!
!      CALL PDLAMODHILB( N, A, 1, 1, DESCA, INFO, Rmat, iter )
!
!     Ask PDSYEV to compute the entire eigendecomposition
!
      do jj= 1, N
         do ii= 1, N
           CALL PZELSET( A, ii, jj, DESCA, Rmat(ii,jj) ) 
         end do  
      end do

      CALL PZHEEVD( 'V', 'L', N, A, 1, 1, DESCA, W, Z, 1, 1, &
                    DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

      deallocate(IWORK)
!
!     Print out the eigenvectors
!
!      CALL PDLAPRNT( N, N, Z, 1, 1, DESCZ, 0, 0, 'Z', 6, WORK )


!      evec2=0.d0

!      if(mod(N,NPROW)==0)then
!        maxii=N/NPROW
!      else
!        if(MYROW<mod(N,NPROW))then
!          maxii=N/NPROW+1
!        else
!          maxii=N/NPROW
!        end if
!      end if

!      if(mod(N,NPCOL)==0)then
!        maxjj=N/NPCOL
!      else
!        if(MYCOL<mod(N,NPCOL))then
!          maxjj=N/NPCOL+1
!        else
!          maxjj=N/NPCOL
!        end if
!      end if

!      do jj=0,maxjj-1
!      do ii=0,maxii-1
!        evec2(ii*NPROW+(MYROW+1),jj*NPCOL+(MYCOL+1))=Z(ii+1,jj+1)
!      end do
 !     end do

!      call comm_summation(evec2,evec,N*N,nproc_group_global)

      SCOPE='A'
      TOP=' '
      do jj=1,N
        do ii=1,N
          call PZELGET( SCOPE, TOP, evec(ii,jj), Z, ii, jj, DESCZ )
        end do
      end do

!      CALL BLACS_GRIDEXIT( CONTEXT )
!
   20 CONTINUE
      iblacsinit=1
!
!      CALL BLACS_EXIT( 0 )
!      CALL BLACS_EXIT( 1 )
!
!
!     Uncomment this line on SUN systems to avoid the useless print out
!
!      CALL IEEE_FLAGS( 'clear', 'exception', 'underflow', '')
!
!
! 9999 FORMAT( 'W=diag([', 4D16.12, ']);' )
!
!      STOP

      deallocate(RWORK)
      deallocate(WORK)

      END SUBROUTINE SAMPLE_PZHEEVD_CALL
