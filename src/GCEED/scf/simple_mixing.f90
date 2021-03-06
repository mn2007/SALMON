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
SUBROUTINE simple_mixing(c1,c2)
use scf_data
implicit none

integer :: ix,iy,iz
real(8),intent(IN) :: c1,c2

!rho = c1*rho + c2*matmul( psi**2, occ )
if(ilsda == 0)then
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rho(ix,iy,iz) = c1*rho_in(ix,iy,iz,num_rho_stock) + c2*rho_out(ix,iy,iz,num_rho_stock)
  end do
  end do
  end do
else if(ilsda == 1)then
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rho_s(ix,iy,iz,1) = c1*rho_s_in(ix,iy,iz,1,num_rho_stock) + c2*rho_s_out(ix,iy,iz,1,num_rho_stock)
    rho_s(ix,iy,iz,2) = c1*rho_s_in(ix,iy,iz,2,num_rho_stock) + c2*rho_s_out(ix,iy,iz,2,num_rho_stock)
    rho(ix,iy,iz) = rho_s(ix,iy,iz,1)+rho_s(ix,iy,iz,2)
  end do
  end do
  end do
end if


return

END SUBROUTINE simple_mixing

