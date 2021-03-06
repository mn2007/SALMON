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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine write_GS_data
  use Global_Variables
  use salmon_global, only: out_dos, &
                         & out_dos_start, &
                         & out_dos_end, &
                         & iout_dos_nenergy, &
                         & out_dos_smearing, &
                         & out_dos_method, &
                         & out_dos_fshift
  use inputoutput, only: unit_length, t_unit_energy_inv, t_unit_energy
  use salmon_parallel, only: nproc_group_global, nproc_id_global, nproc_group_tdks
  use salmon_communication, only: comm_is_root,comm_summation, comm_bcast, comm_sync_all
  use salmon_file, only: open_filehandle
  implicit none
  integer ik,ib,ia,iter,j

  if (comm_is_root(nproc_id_global)) then
    open(403,file=file_GS)
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#grid information-----------------------------------------'
    write(403,*) '#aL =',aL
    write(403,*) '#al(1),al(2),al(3) =',al(1),al(2),al(3)
    write(403,*) '#aLx,aLy,aLz =',aLx,aLy,aLz
    write(403,*) '#bLx,bLy,bLz =',bLx,bLy,bLz
    write(403,*) '#Nd =',Nd
    write(403,*) '#NLx,NLy,NLz=',NLx,NLy,NLz
    write(403,*) '#NL =',NL
    write(403,*) '#Hx,Hy,Hz =',Hx,Hy,Hz
    write(403,*) '#(pi/max(Hx,Hy,Hz))**2 =',(pi/max(Hx,Hy,Hz))**2
    write(403,*) '#(pi/Hx)**2+(pi/Hy)**2+(pi/Hz)**2 =',(pi/Hx)**2+(pi/Hy)**2+(pi/Hz)**2
    write(403,*) '#Hxyz =',Hxyz
    write(403,*) '#NKx,NKy,NKz=',NKx,NKy,NKz
    write(403,*) '#NKxyz =',NKxyz
    write(403,*) '#Sym=',Sym
    write(403,*) '#NK =',NK
    write(403,*) '#NEwald, aEwald =',NEwald, aEwald 
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#GS calc. option------------------------------------------'
    write(403,*) '#FSset_option =',FSset_option
    write(403,*) '#Ncg=',Ncg
    write(403,*) '#Nmemory_MB,alpha_MB =',Nmemory_MB,alpha_MB
    write(403,*) '#NFSset_start,NFSset_every =',NFSset_start,NFSset_every
    write(403,*) '#Nscf=',Nscf
    write(403,*) '#Nscf_conv=',Nscf_conv
    write(403,*) '#NI,NE=',NI,NE
    write(403,*) '#Zatom=',(Zatom(j),j=1,NE)
    write(403,*) '#Lref=',(Lref(j),j=1,NE)
    write(403,*) '#i,Kion(ia)','(Rion(j,a),j=1,3)'
    do ia=1,NI
      write(403,*) '#',ia,Kion(ia)
      write(403,*) '#',(Rion(j,ia),j=1,3)
    end do
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#GS information-------------------------------------------'
    write(403,*) '#NB,Nelec=',NB,Nelec
    write(403,*) '#Eall =',Eall
    write(403,*) '#ddns(iter = Nscf_conv)',ddns(Nscf_conv)
    write(403,*) '#ddns_abs_1e(iter = Nscf_conv)',ddns_abs_1e(Nscf_conv)
    write(403,*) '#esp_var_ave(iter = Nscf_conv)',esp_var_ave(Nscf_conv)
    write(403,*) '#esp_var_max(iter = Nscf_conv)',esp_var_max(Nscf_conv)
    write(403,*) '#NBoccmax is ',NBoccmax
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#band information-----------------------------------------'
    write(403,*) '#Bottom of VB',minval(esp_vb_min(:))
    write(403,*) '#Top of VB',maxval(esp_vb_max(:))
    write(403,*) '#Bottom of CB',minval(esp_cb_min(:))
    write(403,*) '#Top of CB',maxval(esp_cb_max(:))
    write(403,*) '#Fundamental gap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
    write(403,*) '#Fundamental gap[eV]',(minval(esp_cb_min(:))-maxval(esp_vb_max(:)))*2.0*Ry
    write(403,*) '#BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
    write(403,*) '#BG between same k-point[eV]',(minval(esp_cb_min(:)-esp_vb_max(:)))*2.0*Ry
    write(403,*) '#Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
    write(403,*) '#Physicaly upper bound of CB for eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
    write(403,*) '#---------------------------------------------------------'
    write(403,'(a)') ' #iter     total-energy          ddns/nelec         esp_var_ave         esp_var_max'
    do iter=1,Nscf_conv
      !write(403,'(1x,i5,4e20.10)') iter, Eall_GS(iter),ddns(iter),esp_var_ave(iter),esp_var_max(iter)
      write(403,'(1x,i5,4e20.10)') iter, Eall_GS(iter),ddns_abs_1e(iter),esp_var_ave(iter),esp_var_max(iter)
    end do
    close(403)

    !! NOTE: We have changed the SYSNAME_band.out -> SYSNAME_eigen.data
    ! open(405,file=file_band)
    ! write(405,*) '#Bandmap at Ground State'
    ! write(405,*) '#(NK,NB)=','(',NK,NB,')'
    ! do ik=1,NK
    !   do ib=1,NB
    !     write(405,'(1x,2i5,5e26.16e3)') ik,ib,kAc(ik,1),kAc(ik,2),kAc(ik,3) &
    !       ,esp(ib,ik),occ(ib,ik)
    !   enddo
    ! enddo
    ! close(405)

  end if

  
  if(out_dos == 'y')call write_dos_data

  call write_k_data
  call write_eigen_data
  return

  contains


    subroutine write_dos_data
      implicit none
      integer :: fh_dos
      real(8) :: vbmax, cbmin, emax, emin, efermi, eshift
      real(8) :: ww, fk, dw
      integer :: iw
      real(8) :: dos(iout_dos_nenergy), dos_l(iout_dos_nenergy)  
    
      if (comm_is_root(nproc_id_global)) then

        emin = minval(esp(:,:))
        emax = maxval(esp(:,:))
        cbmin = minval(esp_cb_min(:))
        vbmax = maxval(esp_vb_max(:))
      
        if (out_dos_fshift == 'y') then
          efermi = (vbmax + cbmin) * 0.5d0
          eshift = efermi
        else
          eshift = 0d0
        endif
      end if
      call comm_bcast(emin,nproc_group_global)
      call comm_bcast(emax,nproc_group_global)
      call comm_bcast(eshift,nproc_group_global)
      
      out_dos_start = max(out_dos_start, emin - 0.25d0 * (emax - emin))
      out_dos_end = min(out_dos_end, emax + 0.25d0 * (emax - emin))
      
      dos_l = 0d0
      dw = (out_dos_end - out_dos_start) / (iout_dos_nenergy - 1)
      
      select case (out_dos_method)
      case('lorentzian')
        !$omp parallel do private(ik,ib,fk,iw,ww) reduction(+:dos_l) collapse(2)
        do ik = NK_s,NK_e
          do ib = 1,NB
            fk = 2.d0/(NKxyz)*wk(ik)*(out_dos_smearing/pi)          
            do iw = 1, iout_dos_nenergy
              ww =  out_dos_start + (iw-1) * dw + eshift - esp(ib,ik) 
              dos_l(iw) = dos_l(iw) + fk/(ww**2 + out_dos_smearing**2)
            end do
          end do
        end do
        
      case('gaussian')
        !$omp parallel do private(ik,ib,fk,iw,ww) reduction(+:dos_l) collapse(2)
        do ik = NK_s,NK_e
          do ib = 1,NB
            fk = (2.d0 / (NKxyz * sqrt(2.d0*pi) * out_dos_smearing)) * wk(ik)
            do iw = 1, iout_dos_nenergy
              ww =  out_dos_start + (iw-1) * dw + eshift - esp(ib,ik) 
              dos_l(iw) = dos_l(iw) + fk * exp(-(0.5/out_dos_smearing**2)*ww**2)
            end do
          end do
        end do
      end select
      
      call comm_summation(dos_l,dos,iout_dos_nenergy,nproc_group_tdks)

      if (comm_is_root(nproc_id_global)) then
        fh_dos = open_filehandle(file_dos)
        write(fh_dos, '("#",1X,A)') "Density of States"
        
        write(fh_dos, '("#",1X,A,":",1X,A)') "Energy", "Electron energy"
        write(fh_dos, '("#",1X,A,":",1X,A)') "DoS", "Density of States"
        
        write(fh_dos, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "Energy", trim(t_unit_energy%name), &
          & 2, "DoS", trim(t_unit_energy_inv%name)
        do iw = 1, iout_dos_nenergy
          ww =  out_dos_start + (iw-1) * dw 
          write(fh_dos,'(F16.8,99(1X,ES22.14E3))') &
            & ww * t_unit_energy%conv, &
            & dos(iw) * t_unit_energy_inv%conv
        end do
        close(fh_dos)
      end if
      call comm_sync_all
      
    end subroutine write_dos_data

!--------------------------------------------------------------------------------
!! export SYSNAME_k.data file
    subroutine write_k_data()
      implicit none
      integer :: fh_k
      integer :: ik

      if (comm_is_root(nproc_id_global)) then
        fh_k = open_filehandle(file_k_data, status="replace")
        write(fh_k, '("#",1X,A)') "k-point distribution"
        
        write(fh_k, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_k, '("#",1X,A,":",1X,A)') "kx,ky,kz", "Reduced coordinate of k-points"
        write(fh_k, '("#",1X,A,":",1X,A)') "wk", "Weight of k-point"
        
        write(fh_k, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "kx", "none", &
          & 3, "ky", "none", &
          & 4, "kz", "none", &
          & 5, "wk", "none"
        do ik = 1, NK
          write(fh_k, '(I6,99(1X,ES22.14E3))') &
            & ik, & 
            & kAc0(ik,1) / bLx, &
            & kAc0(ik,2) / bLy, &
            & kAc0(ik,3) / bLz, &
            & wk(ik)
        end do !ik
        close(fh_k)
      end if
      call comm_sync_all
    end subroutine write_k_data

  !--------------------------------------------------------------------------------
  !! export SYSNAME_eigen.data file
    subroutine write_eigen_data()
      implicit none
      integer :: fh_eigen
      integer :: ik, ib
    
      if (comm_is_root(nproc_id_global)) then
        fh_eigen = open_filehandle(file_eigen_data, status="replace")
        write(fh_eigen, '("#",1X,A)') "Ground state eigenenergies"
        
        write(fh_eigen, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_eigen, '("#",1X,A,":",1X,A)') "ib", "Band index"
        write(fh_eigen, '("#",1X,A,":",1X,A)') "energy", "Eigenenergy"
        write(fh_eigen, '("#",1X,A,":",1X,A)') "occup", "Occupation"
        
        write(fh_eigen, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "ib", "none", &
          & 3, "energy", trim(t_unit_energy%name), &
          & 4, "occup", "none"
        do ik = 1, NK
          do ib = 1, NB
            write(fh_eigen, '(I6,1X,I6,99(1X,ES22.14E3))') &
              & ik, ib, esp(ib,ik)*t_unit_energy%conv, occ(ib,ik)/wk(ik)*NKxyz
          end do !ib
        end do !ik
        close(fh_eigen)
      end if
      call comm_sync_all
    end subroutine write_eigen_data
    
End Subroutine write_GS_data
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
