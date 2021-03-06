!
!  Copyright 2017 SALMON developers
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
!=======================================================================
!============================ Hartree potential (Solve Poisson equation)
SUBROUTINE Hartree_boundary(trho,wk2)
use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_h, &
                           nproc_id_bound, nproc_size_bound, nproc_group_bound
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
use sendrecvh_sub
use allocate_mat_sub
use deallocate_mat_sub

use omp_lib, only: omp_get_num_threads, omp_get_thread_num, omp_get_max_threads
use misc_routines, only: ceiling_pow2

implicit none
real(8) :: trho(mg_sta(1):mg_end(1),    &
               mg_sta(2):mg_end(2),      &
               mg_sta(3):mg_end(3))

integer,parameter :: maxiter=1000
integer :: ii,jj,kk,ix,iy,iz,lm,LL,icen,pl,cl
integer :: ixbox,iybox,izbox
integer :: k
integer :: istart(0:nproc_size_global-1),iend(0:nproc_size_global-1)
integer :: icount
integer,allocatable :: itrho(:)
integer :: num_center
real(8) :: Ylm
real(8) :: Ylm2(25)
integer :: L2(25)
real(8) :: xx,yy,zz,rr,sum1,xxxx,yyyy,zzzz,rrrr,sumbox1,sumbox2,sumbox3
real(8) :: rholm2box
real(8),allocatable :: rholm(:,:),rholm2(:,:),rholm3(:,:)
integer :: tid
real(8) :: center_trho2(3)
real(8),allocatable :: center_trho(:,:)
real(8),allocatable :: center_trho_nume_deno(:,:)
real(8),allocatable :: center_trho_nume_deno2(:,:)
real(8) :: wk2(ng_sta(1)-Ndh:ng_end(1)+Ndh,    &
               ng_sta(2)-Ndh:ng_end(2)+Ndh,      &
               ng_sta(3)-Ndh:ng_end(3)+Ndh)
real(8) :: xp2,yp2,zp2,xy,yz,xz
real(8) :: deno(25)
real(8) :: rinv
real(8) :: rbox
real(8),allocatable :: Rion2(:,:)

iwk_size=12
call make_iwksta_iwkend

!------------------------- Boundary condition (multipole expansion)

select case( MEO )

case(1)

num_center=1
allocate (rholm((lmax_MEO+1)**2,1))
allocate (rholm2((lmax_MEO+1)**2,1))
allocate(itrho(1))
allocate(center_trho(3,1))
do jj=1,3
  center_trho(jj,1)=0.d0
  center_trho2(jj)=0.d0
end do
do lm=1,(lmax_MEO+1)**2
  rholm(lm,1)=0.d0
  rholm2(lm,1)=0.d0
end do
itrho(1)=1

do LL=0,lmax_MEO
do lm=LL**2+1,(LL+1)**2
  rholm2box=0.d0
!$OMP parallel do reduction ( + : rholm2box)&
!$OMP private(ix,iy,iz,xx,yy,zz,rr,xxxx,yyyy,zzzz,Ylm)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    xx=gridcoo(ix,1)-center_trho(1,1)
    yy=gridcoo(iy,2)-center_trho(2,1)
    zz=gridcoo(iz,3)-center_trho(3,1)
    rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
    call Ylm_sub(xxxx,yyyy,zzzz,lm,Ylm)
    rholm2box=rholm2box+rr**LL*Ylm*trho(ix,iy,iz)*Hvol
  end do
  end do
  end do
  rholm2(lm,1)=rholm2box
end do
end do

case(2)

if(iflag_ps==1)then
  num_center=MI
  allocate (rholm((lmax_MEO+1)**2,MI))
  allocate (rholm2((lmax_MEO+1)**2,MI))
  allocate(itrho(MI))
  allocate(center_trho(3,MI))
  allocate(Rion2(3,MI))
  Rion2(:,:)=Rion(:,:)
end if

!$OMP parallel do private(icen, lm, jj)
do icen=1,num_center
  do lm=1,(lmax_MEO+1)**2
    rholm(lm,icen)=0.d0
    rholm2(lm,icen)=0.d0
  end do
  do jj=1,3
    center_trho(jj,icen)=Rion2(jj,icen)
  end do
  itrho(icen)=1
end do

do icen=1,num_pole_myrank
  do LL=0,lmax_MEO
  do lm=LL**2+1,(LL+1)**2
    rholm2box=0.d0
!$OMP parallel do reduction ( + : rholm2box)&
!$OMP private(jj,ix,iy,iz,xx,yy,zz,rr,xxxx,yyyy,zzzz,Ylm)
    do jj=1,icount_pole(icen)
      ix=icorr_xyz_pole(1,jj,icen)
      iy=icorr_xyz_pole(2,jj,icen)
      iz=icorr_xyz_pole(3,jj,icen)
      xx=gridcoo(ix,1)-Rion2(1,icorr_polenum(icen))
      yy=gridcoo(iy,2)-Rion2(2,icorr_polenum(icen))
      zz=gridcoo(iz,3)-Rion2(3,icorr_polenum(icen))
      rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
      call Ylm_sub(xxxx,yyyy,zzzz,lm,Ylm)
      rholm2box=rholm2box+rr**LL*Ylm*trho(ix,iy,iz)*Hvol
    end do
    rholm2(lm,icorr_polenum(icen))=rholm2box
  end do
  end do
end do

case(3)

num_center=num_pole

allocate (rholm((lmax_MEO+1)**2,num_center))
allocate (rholm2((lmax_MEO+1)**2,num_center))
allocate(itrho(num_center))
allocate(center_trho(3,num_center))
allocate(center_trho_nume_deno(4,num_center))
allocate(center_trho_nume_deno2(4,num_center))

!$OMP parallel do private(icen, jj, lm)
do icen=1,num_center
  do jj=1,4
    center_trho_nume_deno2(jj,icen)=0.d0
  end do
  do lm=1,(lmax_MEO+1)**2
    rholm(lm,icen)=0.d0
    rholm2(lm,icen)=0.d0
  end do
end do

do ii=1,num_pole_myrank
  sum1=0.d0
  sumbox1=0.d0
  sumbox2=0.d0
  sumbox3=0.d0
!$OMP parallel do reduction (+ : sumbox1, sumbox2, sumbox3, sum1) &
!$OMP private(jj,ixbox,iybox,izbox,xx,yy,zz)
  do jj=1,icount_pole(ii)
    ixbox=icorr_xyz_pole(1,jj,ii)
    iybox=icorr_xyz_pole(2,jj,ii)
    izbox=icorr_xyz_pole(3,jj,ii)
    xx=gridcoo(ixbox,1)
    yy=gridcoo(iybox,2)
    zz=gridcoo(izbox,3)
    sumbox1=sumbox1+trho(ixbox,iybox,izbox)*xx
    sumbox2=sumbox2+trho(ixbox,iybox,izbox)*yy
    sumbox3=sumbox3+trho(ixbox,iybox,izbox)*zz
    sum1=sum1+trho(ixbox,iybox,izbox)
  end do
  center_trho_nume_deno2(1,icorr_polenum(ii))=sumbox1
  center_trho_nume_deno2(2,icorr_polenum(ii))=sumbox2
  center_trho_nume_deno2(3,icorr_polenum(ii))=sumbox3
  center_trho_nume_deno2(4,icorr_polenum(ii))=sum1
end do

elp3(201)=get_wtime()
call comm_summation(center_trho_nume_deno2,center_trho_nume_deno,4*num_pole,nproc_group_h)
elp3(202)=get_wtime()
elp3(251)=elp3(251)+elp3(202)-elp3(201)

do ii=1,num_pole
  if(center_trho_nume_deno(4,ii)*Hvol>=1.d-12)then
    itrho(ii)=1
    center_trho(1:3,ii)=center_trho_nume_deno(1:3,ii)/center_trho_nume_deno(4,ii)
  else
    itrho(ii)=0
    center_trho(1:3,ii)=0.d0
  end if
end do

if(omp_get_max_threads() > 16) then
!$omp parallel shared(rholm3,lmax_MEO)
!$omp master
    allocate(rholm3((lmax_MEO+1)**2,0:ceiling_pow2(omp_get_num_threads())-1))
!$omp end master
!$omp end parallel

  rholm2=0.d0
  rholm3=0.d0
  do ii=1,num_pole_myrank
    pl=icorr_polenum(ii)
    cl=icount_pole(ii)
    if(itrho(pl)==1)then
!$omp parallel default(none) &
!$omp          shared(icorr_xyz_pole,gridcoo,center_trho,trho,rholm3) &
!$omp          private(tid,kk,jj,LL,lm,ixbox,iybox,izbox,xx,yy,zz,rr,rinv,xxxx,yyyy,zzzz,Ylm) &
!$omp          firstprivate(ii,pl,cl,lmax_MEO,Hvol)
      tid=omp_get_thread_num()
      rholm3(:,tid)=0.d0

!$omp do
      do jj=1,cl
        ixbox=icorr_xyz_pole(1,jj,ii)
        iybox=icorr_xyz_pole(2,jj,ii)
        izbox=icorr_xyz_pole(3,jj,ii)
        xx=gridcoo(ixbox,1)-center_trho(1,pl)
        yy=gridcoo(iybox,2)-center_trho(2,pl)
        zz=gridcoo(izbox,3)-center_trho(3,pl)
        rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d0-50.d0
        rinv=1.0d0/rr
        xxxx=xx*rinv
        yyyy=yy*rinv
        zzzz=zz*rinv
        do LL=0,lmax_MEO
        do lm=LL**2+1,(LL+1)**2
          call Ylm_sub(xxxx,yyyy,zzzz,lm,Ylm)
          rholm3(lm,tid)=rholm3(lm,tid)+rr**LL*Ylm*trho(ixbox,iybox,izbox)*Hvol
        end do
        end do
      end do
!$omp end do

      kk = ceiling_pow2(omp_get_num_threads())/2
      do while(kk > 0)
        if(tid < kk) then
          rholm3(:,tid) = rholm3(:,tid) + rholm3(:,tid+kk)
        end if
        kk = kk/2
!$omp barrier
      end do
!$omp end parallel
    end if
    rholm2(:,pl)=rholm3(:,0)
  end do
  deallocate(rholm3)
else
  rholm2=0.d0
  do ii=1,num_pole_myrank
    if(itrho(icorr_polenum(ii))==1)then
      rholm=0.d0
      do LL=0,lmax_MEO
      do lm=LL**2+1,(LL+1)**2
        rholm2box=0.d0
!$OMP parallel do reduction ( + : rholm2box)&
!$OMP private(jj,ixbox,iybox,izbox,xx,yy,zz,rr,xxxx,yyyy,zzzz,Ylm)
        do jj=1,icount_pole(ii)
          ixbox=icorr_xyz_pole(1,jj,ii)
          iybox=icorr_xyz_pole(2,jj,ii)
          izbox=icorr_xyz_pole(3,jj,ii)
          xx=gridcoo(ixbox,1)-center_trho(1,icorr_polenum(ii))
          yy=gridcoo(iybox,2)-center_trho(2,icorr_polenum(ii))
          zz=gridcoo(izbox,3)-center_trho(3,icorr_polenum(ii))
          rr=sqrt(xx*xx+yy*yy+zz*zz)+1.d-50 ; xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
          call Ylm_sub(xxxx,yyyy,zzzz,lm,Ylm)
          rholm2box=rholm2box+rr**LL*Ylm*trho(ixbox,iybox,izbox)*Hvol
        end do
        rholm2(lm,icorr_polenum(ii))=rholm2box
      end do
      end do
    end if
  end do
endif

deallocate(center_trho_nume_deno)
deallocate(center_trho_nume_deno2)

end select

if(nproc_size_global==1)then
!$OMP parallel do
  do icen=1,num_center
    rholm(:,icen)=rholm2(:,icen)
  end do
else
  elp3(201)=get_wtime()
  call comm_summation(rholm2,rholm,(lmax_MEO+1)**2*num_center,nproc_group_h)
  elp3(202)=get_wtime()
  elp3(252)=elp3(252)+elp3(202)-elp3(201)
end if

!$OMP parallel do private(iz,iy,ix) collapse(2)
do iz=ng_sta(3)-Ndh,ng_end(3)+Ndh
do iy=ng_sta(2)-Ndh,ng_end(2)+Ndh
do ix=ng_sta(1)-Ndh,ng_end(1)+Ndh
  wk2(ix,iy,iz)=0.d0
end do
end do
end do

do k=1,3

!$OMP parallel do
  do jj=1,lg_num(1)*lg_num(2)*lg_num(3)/minval(lg_num(1:3))*6*Ndh
    wk2bound_h(jj)=0.d0
  end do

  icount=inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)  &
         *inum_Mxin_s(3,nproc_id_global)/inum_Mxin_s(k,nproc_id_global)*2*Ndh
!$OMP parallel do
  do ii=0,nproc_size_bound(k)-1
    istart(ii)=ii*icount/nproc_size_bound(k)+1
    iend(ii)=(ii+1)*icount/nproc_size_bound(k)
  end do
        
  do LL=0,lmax_MEO
    do lm=LL**2+1,(LL+1)**2
      L2(lm)=LL
    end do
  end do

!$OMP parallel do &
!$OMP private(xx,yy,zz,rr,xxxx,yyyy,zzzz,lm,LL,sum1,Ylm2,rrrr,xp2,yp2,zp2,xy,yz,xz,rinv,rbox,deno,icen)
  do jj=istart(nproc_id_bound(k)),iend(nproc_id_bound(k))
    do icen=1,num_center
      if(itrho(icen)==1)then
        xx=gridcoo(icoobox_bound(1,jj,k),1)-center_trho(1,icen)
        yy=gridcoo(icoobox_bound(2,jj,k),2)-center_trho(2,icen)
        zz=gridcoo(icoobox_bound(3,jj,k),3)-center_trho(3,icen)
        rr=sqrt(xx**2+yy**2+zz**2)+1.d-50 
        rinv=1.d0/rr
!        xxxx=xx/rr ; yyyy=yy/rr ; zzzz=zz/rr
        xx=xx*rinv ; yy=yy*rinv ; zz=zz*rinv

        xp2=xx**2
        yp2=yy**2
        zp2=zz**2
        xy=xx*yy
        yz=yy*zz
        xz=xx*zz

        rrrr=(xp2+yp2+zp2)**2

        Ylm2(1)=1.d0
        Ylm2(2)=yy
        Ylm2(3)=zz
        Ylm2(4)=xx
        Ylm2(5)=sqrt(3.d0)*xy                      ! lm=5  (2 -2)
        Ylm2(6)=sqrt(3.d0)*yz                      ! lm=6  (2 -1)
        Ylm2(7)=(2*zp2-xp2-yp2)/2.d0               ! lm=7  (2 0)
        Ylm2(8)=sqrt(3.d0)*xz                      ! lm=8  (2 1)
        Ylm2(9)=sqrt(3.d0/4.d0)*(xp2-yp2)          ! lm=9  (2 2)
        Ylm2(10)=sqrt(5.d0/8.d0)*yy*(3*xp2-yp2)    ! lm=10 (3 -3)
        Ylm2(11)=sqrt(15.d0)*xx*yy*zz               ! lm=11 (3 -2)
        Ylm2(12)=sqrt(3.d0/8.d0)*yy*(4*zp2-xp2-yp2)  ! lm=12 (3 -1)
        Ylm2(13)=zz*(2*zp2-3*xp2-3*yp2)/2.d0         ! lm=13 (3 0)
        Ylm2(14)=sqrt(3.d0/8.d0)*xx*(4*zp2-xp2-yp2)  ! lm=14 (3 1)
        Ylm2(15)=sqrt(15.d0/4.d0)*zz*(xp2-yp2)       ! lm=15 (3 2)
        Ylm2(16)=sqrt(5.d0/8.d0)*xx*(xp2-3*yp2)      ! lm=16 (3 3)

        Ylm2(17)=rrrr*sqrt(35.d0)/2.d0*xy*(xp2-yp2)
        Ylm2(18)=rrrr*sqrt(35.d0/8.d0)*yz*(3*xp2-yp2)
        Ylm2(19)=rrrr*sqrt(5.d0)/2.d0*xy*(7*zp2-1.d0)
        Ylm2(20)=rrrr*sqrt(5.d0/8.d0)*yz*(7*zp2-3.d0)
        Ylm2(21)=rrrr*(35*zp2**2-30*zp2+3.d0)/8.d0
        Ylm2(22)=rrrr*sqrt(5.d0/8.d0)*xz*(7*zp2-3.d0)
        Ylm2(23)=rrrr*sqrt(5.d0)/4.d0*(7*zp2-1)*(xp2-yp2)
        Ylm2(24)=rrrr*sqrt(35.d0/8.d0)*xz*(xp2-3*yp2)
        Ylm2(25)=rrrr*sqrt(35.d0)/8.d0*(xp2**2+yp2**2-6*xp2*yp2)

        deno(1)=rinv
        rbox=rinv*rinv
        deno(2:4)=rbox
        rbox=rbox*rinv
        deno(5:9)=rbox
        rbox=rbox*rinv
        deno(10:16)=rbox
        rbox=rbox*rinv
        deno(17:25)=rbox

        sum1=0.d0
        do lm=1,(lmax_MEO+1)**2
          sum1=sum1+Ylm2(lm)*deno(lm)*rholm(lm,icen)
        end do
        wk2bound_h(jj) = wk2bound_h(jj) + sum1
      end if
    end do
  end do

  if(nproc_size_global==1)then
!$OMP parallel do
    do jj=1,icount
      wkbound_h(jj)=wk2bound_h(jj)
    end do
  else
    elp3(201)=get_wtime()
    call comm_summation( &
      wk2bound_h,              wkbound_h,              icount/2, nproc_group_bound(k), 0                    )
    call comm_summation( &
      wk2bound_h(icount/2+1:), wkbound_h(icount/2+1:), icount/2, nproc_group_bound(k), nproc_size_bound(k)-1)
    elp3(202)=get_wtime()
    elp3(253)=elp3(253)+elp3(202)-elp3(201)
  end if

  if(nproc_id_bound(k)==0) then
!$OMP parallel do
    do jj=1,icount/2
      wk2(icoobox_bound(1,jj,k),icoobox_bound(2,jj,k),icoobox_bound(3,jj,k))=wkbound_h(jj)
    end do
  end if
  if(nproc_id_bound(k)==nproc_size_bound(k)-1) then
!$OMP parallel do
    do jj=icount/2+1,icount
      wk2(icoobox_bound(1,jj,k),icoobox_bound(2,jj,k),icoobox_bound(3,jj,k))=wkbound_h(jj)
    end do
  end if

end do

deallocate(rholm,rholm2)
deallocate(center_trho)

return

END SUBROUTINE Hartree_boundary
