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
! ifunc : +1:  calclation of esp / +0: no calculation of esp
!         +2:  2nd or larger loop (zpsi+htpsi) / +0: 1st loop (tpsi+htpsi)
!         +4:  calculation of rhobox / +0: no calculation of rhobox
!         -1:  special version for N_hamil=4 and nn=3
subroutine add_polynomial_ob(tpsi,htpsi,tpsi_out,iobmax,nn,ifunc,iiik,iiob,iz_s,iz_e,iy_s,iy_e)
use scf_data
implicit none
complex(8) :: tpsi(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,k_sta:k_end)
complex(8) :: htpsi(iwk2sta(1):iwk2end(1)+1,  &
                    iwk2sta(2):iwk2end(2),      &
                    iwk2sta(3):iwk2end(3),     &
                   1:iobnum,k_sta:k_end)
complex(8) :: tpsi_out(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,k_sta:k_end)
integer :: iobmax
integer :: ifunc
integer :: nn
integer :: iob,ix,iy,iz,iik
complex(8) :: cbox
complex(8), parameter :: zi=(0.d0,1.d0)
integer :: iob_allob
integer :: iiik,iiob
integer :: iz_s,iz_e
integer :: iy_s,iy_e

if(ifunc==0)then
!$OMP parallel do collapse(2) private(iz,iy,ix) 
  do iz=iz_s,iz_e
  do iy=iy_s,iy_e
  do ix=mg_sta(1),mg_end(1)
    htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
    tpsi_out(ix,iy,iz,iiob,iiik)=tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
  end do
  end do
  end do
else if(ifunc==1)then
  cbox=0.d0
!$OMP parallel do reduction(+:cbox) collapse(2) private(iz,iy,ix) 
  do iz=iz_s,iz_e
  do iy=iy_s,iy_e
  do ix=mg_sta(1),mg_end(1)
    cbox=cbox+conjg(tpsi(ix,iy,iz,iiob,iiik))*htpsi(ix,iy,iz,iiob,iiik)
    htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
    tpsi_out(ix,iy,iz,iiob,iiik)=tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
  end do
  end do
  end do
  esp2(iiob,iiik)=esp2(iiob,iiik)+dble(cbox)*Hvol
else if(ifunc==2)then
!$OMP parallel do collapse(2) private(iz,iy,ix) 
  do iz=iz_s,iz_e
  do iy=iy_s,iy_e
  do ix=mg_sta(1),mg_end(1)
    htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
    tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
  end do
  end do
  end do
else if(ifunc==3)then
  cbox=0.d0
!$OMP parallel do reduction(+:cbox) collapse(2) private(iz,iy,ix) 
  do iz=iz_s,iz_e
  do iy=iy_s,iy_e
  do ix=mg_sta(1),mg_end(1)
    cbox=cbox+conjg(tpsi(ix,iy,iz,iiob,iiik))*htpsi(ix,iy,iz,iiob,iiik)
    htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
    tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
  end do
  end do
  end do
  esp2(iiob,iiik)=esp2(iiob,iiik)+dble(cbox)*Hvol
else if(ifunc==4)then
  if(ilsda==0)then
    call calc_allob(iiob,iob_allob)
!$OMP parallel do collapse(2) private(iz,iy,ix) 
    do iz=iz_s,iz_e
    do iy=iy_s,iy_e
    do ix=mg_sta(1),mg_end(1)
      htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
      tpsi_out(ix,iy,iz,iiob,iiik)=tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
      rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+  &
        tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
    end do
    end do
    end do
  else
    call calc_allob(iiob,iob_allob)
    if(iob_allob<=MST(1))then
!$OMP parallel do collapse(2) private(iz,iy,ix) 
      do iz=iz_s,iz_e
      do iy=iy_s,iy_e
      do ix=mg_sta(1),mg_end(1)
        htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
        tpsi_out(ix,iy,iz,iiob,iiik)=tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
        rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+  &
          tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
      end do
      end do
      end do
    else
!$OMP parallel do collapse(2) private(iz,iy,ix) 
      do iz=iz_s,iz_e
      do iy=iy_s,iy_e
      do ix=mg_sta(1),mg_end(1)
        htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
        tpsi_out(ix,iy,iz,iiob,iiik)=tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
        rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+  &
          tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
      end do
      end do
      end do
    end if
  end if 
else if(ifunc==5)then
  if(ilsda==0)then
    call calc_allob(iiob,iob_allob)
    cbox=0.d0
!$OMP parallel do reduction(+:cbox) collapse(2) private(iz,iy,ix)
    do iz=iz_s,iz_e
    do iy=iy_s,iy_e
    do ix=mg_sta(1),mg_end(1)
      cbox=cbox+conjg(tpsi(ix,iy,iz,iiob,iiik))*htpsi(ix,iy,iz,iiob,iiik)
      htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
      tpsi_out(ix,iy,iz,iiob,iiik)=tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
      rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+  &
        tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
    end do
    end do
    end do
    esp2(iiob,iiik)=esp2(iiob,iiik)+dble(cbox)*Hvol
  else if(ilsda==1)then
    call calc_allob(iiob,iob_allob)
    cbox=0.d0
    if(iob_allob<=MST(1))then
!$OMP parallel do reduction(+:cbox) collapse(2) private(iz,iy,ix) 
      do iz=iz_s,iz_e
      do iy=iy_s,iy_e
      do ix=mg_sta(1),mg_end(1)
        cbox=cbox+conjg(tpsi(ix,iy,iz,iiob,iiik))*htpsi(ix,iy,iz,iiob,iiik)
        htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
        tpsi_out(ix,iy,iz,iiob,iiik)=tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
        rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+  &
          tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
      end do
      end do
      end do
    else
!$OMP parallel do reduction(+:cbox) collapse(2) private(iz,iy,ix) 
      do iz=iz_s,iz_e
      do iy=iy_s,iy_e
      do ix=mg_sta(1),mg_end(1)
        cbox=cbox+conjg(tpsi(ix,iy,iz,iiob,iiik))*htpsi(ix,iy,iz,iiob,iiik)
        htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
        tpsi_out(ix,iy,iz,iiob,iiik)=tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
        rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+  &
          tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
      end do
      end do
      end do
    end if
    esp2(iiob,iiik)=esp2(iiob,iiik)+dble(cbox)*Hvol
  end if
else if(ifunc==6)then
  if(ilsda==0)then
    call calc_allob(iiob,iob_allob)
!$OMP parallel do collapse(2) private(iz,iy,ix) 
    do iz=iz_s,iz_e
    do iy=iy_s,iy_e
    do ix=mg_sta(1),mg_end(1)
      htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
      tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
      rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+  &
        tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
    end do
    end do
    end do
  else if(ilsda==1)then
    call calc_allob(iiob,iob_allob)
    if(iob_allob<=MST(1))then
!$OMP parallel do collapse(2) private(iz,iy,ix) 
      do iz=iz_s,iz_e
      do iy=iy_s,iy_e
      do ix=mg_sta(1),mg_end(1)
        htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
        tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
        rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+  &
          tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
      end do
      end do
      end do
    else
!$OMP parallel do collapse(2) private(iz,iy,ix) 
      do iz=iz_s,iz_e
      do iy=iy_s,iy_e
      do ix=mg_sta(1),mg_end(1)
        htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
        tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
        rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+  &
          tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
      end do
      end do
      end do
    end if
  end if
else if(ifunc==7)then
  if(ilsda==0)then
    cbox=0.d0
!$OMP parallel do reduction(+:cbox) collapse(2) private(iz,iy,ix) 
    do iz=iz_s,iz_e
    do iy=iy_s,iy_e
    do ix=mg_sta(1),mg_end(1)
      cbox=cbox+conjg(tpsi(ix,iy,iz,iiob,iiik))*htpsi(ix,iy,iz,iiob,iiik)
      htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
      tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
      rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+  &
        tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
    end do
    end do
    end do
    esp2(iiob,iiik)=esp2(iiob,iiik)+dble(cbox)*Hvol
  else if(ilsda==1)then
    call calc_allob(iiob,iob_allob)
    cbox=0.d0
    if(iob_allob<=MST(1))then
!$OMP parallel do reduction(+:cbox) collapse(2) private(iz,iy,ix) 
      do iz=iz_s,iz_e
      do iy=iy_s,iy_e
      do ix=mg_sta(1),mg_end(1)
        cbox=cbox+conjg(tpsi(ix,iy,iz,iiob,iiik))*htpsi(ix,iy,iz,iiob,iiik)
        htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
        tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
        rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+  &
          tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
      end do
      end do
      end do
    else
!$OMP parallel do reduction(+:cbox) collapse(2) private(iz,iy,ix) 
      do iz=iz_s,iz_e
      do iy=iy_s,iy_e
      do ix=mg_sta(1),mg_end(1)
        cbox=cbox+conjg(tpsi(ix,iy,iz,iiob,iiik))*htpsi(ix,iy,iz,iiob,iiik)
        htpsi(ix,iy,iz,iiob,iiik)=-zi*dt*htpsi(ix,iy,iz,iiob,iiik)/dble(nn)
        tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
        rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+  &
          tpsi_out(ix,iy,iz,iiob,iiik)*conjg(tpsi_out(ix,iy,iz,iiob,iiik))*rocc(iob_allob,iiik)*wtk(iiik)
      end do
      end do
      end do
    end if
    esp2(iiob,iiik)=esp2(iiob,iiik)+dble(cbox)*Hvol
  end if
else if(ifunc==-1)then
!$OMP parallel do collapse(2) private(iz,iy,ix) 
  do iz=iz_s,iz_e
  do iy=iy_s,iy_e
  do ix=mg_sta(1),mg_end(1)
    tpsi(ix,iy,iz,iiob,iiik)=-zi*dt*tpsi(ix,iy,iz,iiob,iiik)/dble(nn-1)
    htpsi(ix,iy,iz,iiob,iiik)=(-zi*dt)**2*htpsi(ix,iy,iz,iiob,iiik)/dble(nn-1)/dble(nn)
    tpsi_out(ix,iy,iz,iiob,iiik)=tpsi_out(ix,iy,iz,iiob,iiik)+tpsi(ix,iy,iz,iiob,iiik)+htpsi(ix,iy,iz,iiob,iiik)
  end do
  end do
  end do
end if

end subroutine add_polynomial_ob
