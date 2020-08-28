subroutine getpartprop(fieldout,field,np,nz,ny,nx,ix,iy,iz)
  implicit none
  integer :: nx,ny,nz, np
  real :: field(nz,ny,nx)
  real :: ix(np), iy(np), iz(np)
  real, intent(out) :: fieldOut(np)
  integer :: i
  real :: interp
  do i=1,np
     if(ix(i)>=1 .and. iy(i)>=1 .and. iz(i)>=1) then 
        fieldout(i)=interp(field,nz,ny,nx,ix(i),iy(i),iz(i))
     else
        fieldout(i)=-999
     endif
  enddo
end subroutine getpartprop

subroutine getpartprop2(fieldout,field,np,ny,nx,ix,iy,iz)
  implicit none
  integer :: nx,ny,nz, np
  real :: field(ny,nx)
  real :: ix(np), iy(np), iz(np)
  real, intent(out) :: fieldOut(np)
  integer :: i
  real :: interp2
  do i=1,np
     if(ix(i)>=1 .and. iy(i)>=1 .and. iz(i)>=1) then 
        fieldout(i)=interp2(field,ny,nx,ix(i),iy(i))
        !interp2(field,ny,nx,ix,iy)
     else
        fieldout(i)=-999
     endif
  enddo
end subroutine getpartprop2
real function interp(field,nz,ny,nx,ix,iy,iz)
  implicit none
  integer :: nx,ny,nz
  real :: field(nz,ny,nx)
  real :: ix, iy, iz
  real :: s0, s1, sx1, sy1, sz1, sx, sy, sz
  integer :: ix0,iy0,iz0, ix1, iy1, iz1
  
  ix0=int(ix)
  sx=ix-ix0
  iy0=int(iy)
  sy=iy-iy0
  iz0=int(iz)
  sz=iz-iz0
  ix1=ix0+1
  iy1=iy0+1
  iz1=iz0+1
  if (ix1>=nx) ix1=nx
  if (iy1>=ny) iy1=ny
  if (iz1>=nz) iz1=nz
  if (ix1<1) ix1=1
  if (iy1<1) iy1=1
  if (iz1<1) iz1=1
  if (ix0>=nx) ix0=nx
  if (iy0>=ny) iy0=ny
  if (iz0>=nz) iz0=nz
  if (ix0<1) ix0=1
  if (iy0<1) iy0=1
  if (iz0<1) iz0=1
  sx1=1-sx
  sy1=1-sy
  sz1=1-sz
  s0=field(iz0,iy0,ix0)*sx1*sy1+field(iz0,iy0,ix0+1)*(1-sx1)*sy1+&
       field(iz0,iy0+1,ix0)*sx1*(1-sy1)+&
       field(iz0,iy0+1,ix0+1)*(1-sx1)*(1-sy1)
  s1=field(iz0+1,iy0,ix0)*sx1*sy1+field(iz0+1,iy0,ix0+1)*(1-sx1)*sy1+&
       field(iz0+1,iy0+1,ix0)*sx1*(1-sy1)+&
       field(iz0+1,iy0+1,ix0+1)*(1-sx1)*(1-sy1)

  interp=s0*sz1+(1-sz1)*s1
  
end function interp

real function interp2(field,ny,nx,ix,iy)
  implicit none
  integer :: nx,ny
  real :: field(ny,nx)
  real :: ix, iy
  real :: s0, sx1, sy1, sx, sy
  integer :: ix0,iy0,ix1, iy1
  
  ix0=int(ix)
  sx=ix-ix0
  iy0=int(iy)
  sy=iy-iy0
  ix1=ix0+1
  iy1=iy0+1
  if (ix1>=nx) ix1=nx
  if (iy1>=ny) iy1=ny
  if (ix1<1) ix1=1
  if (iy1<1) iy1=1
  if (ix0>=nx) ix0=nx
  if (iy0>=ny) iy0=ny
  if (ix0<1) ix0=1
  if (iy0<1) iy0=1
  sx1=1-sx
  sy1=1-sy
  s0=field(iy0,ix0)*sx1*sy1+field(iy0,ix0+1)*(1-sx1)*sy1+&
       field(iy0+1,ix0)*sx1*(1-sy1)+&
       field(iy0+1,ix0+1)*(1-sx1)*(1-sy1)
  interp2=s0
  
end function interp2

subroutine trace3db(u,v,w,dz,mapf,ix,iy,iz,nz,nx,ny,dx,dy,dt)
  implicit none
  integer :: nx, ny, nz
  real :: ix, iy, iz, dt
  real :: u(nz,ny,nx), v(nz,ny,nx), w(nz,ny,nx), mapf(ny,nx), &
       dz(nz,ny,nx), dx, dy
  real :: u1, v1, w1, dz1
  real :: interp, interp2, m1
  if(ix>=1 .and. ix<=nx .and. &
       iy>=1 .and. iy<=ny .and. &
       iz>=1 .and. iz<=nz) then
     u1=interp(u,nz,ny,nx,ix,iy,iz)
     v1=interp(v,nz,ny,nx,ix,iy,iz)
     w1=interp(w,nz,ny,nx,ix,iy,iz)
     dz1=interp(dz,nz,ny,nx,ix,iy,iz)
     m1=interp2(mapf,ny,nx,ix,iy)
     !print*, u1,v1,w1,ix,iy
     ix=ix+u1*dt/dx*m1
     iy=iy+v1*dt/dy*m1
     iz=iz+w1*dt/dz1
  else
     ix=-999
     iy=-999
     iz=-999
  endif
end subroutine trace3db

subroutine trace2db(u,v,mapf,ix,iy,nx,ny,dx,dy,dt)
  implicit none
  integer :: nx, ny, nz
  real :: ix, iy, iz, dt
  real :: u(ny,nx), v(ny,nx), mapf(ny,nx), &
        dx, dy
  real :: u1, v1, dz1
  real :: interp, interp2, m1
  if(ix>=1 .and. ix<=nx .and. &
       iy>=1 .and. iy<=ny ) then
     u1=interp2(u,ny,nx,ix,iy)
     v1=interp2(v,ny,nx,ix,iy)
     m1=interp2(mapf,ny,nx,ix,iy)
     !print*, u1,v1,w1,ix,iy
     ix=ix+u1*dt/dx*m1
     iy=iy+v1*dt/dy*m1
  else
     ix=-999
     iy=-999
  endif
end subroutine trace2db

subroutine woro_sub(u,v,zm,mapf,&
     lon,lat,iz,nz,nx,ny,dx,dy,dt,tinc,woro)
  implicit none
  integer :: nx, ny, nz, npart, nt
  real :: u(ny,nx), v(ny,nx), &
       mapf(ny,nx), lon(ny,nx), lat(ny,nx),  &
       zm(2,nz,ny,nx)
  real  :: dx,dy
  real, intent(out) :: woro(ny,nx)
  integer ::  iz, im, jm, i, j
  real :: dt, tinc
  real :: ix,iy
  integer :: nt1, it, it1, ip
  real :: interp2, interp
  
  nt1=int((tinc)/dt)
  

  do i=10,nx-10
     do j=10,ny-10
        ix=i
        iy=j
        do it=1,nt1-1
           call trace2db(u,v,mapf,ix,iy,nx,ny,dx,dy,dt)
        end do
        if(ix>=1 .and. ix<=nx .and. iy>=1 .and. iy<=ny) then
           im=int(0.5*(i+ix))
           jm=int(0.5*(j+iy))
           woro(jm,im)=(zm(2,iz,int(iy),int(ix))-&
                zm(1,iz,j,i))/tinc
        endif
     end do
  end do
       
end subroutine woro_sub

subroutine trace3dball(u,v,w,dz,zm,mapf,&
     lon,lat,ix,iy,iz,nz,nx,ny,dx,dy,dt,tinc,nt,&
     npart, posPartOut, lonOut, latOut, hOut,zout)
  implicit none
  integer :: nx, ny, nz, npart, nt
  real :: u(nt,nz,ny,nx), v(nt,nz,ny,nx), w(nt,nz,ny,nx), &
       mapf(ny,nx), lon(ny,nx), lat(ny,nx),  &
       dz(nt,nz,ny,nx), dx, dy, tinc, dt, zm(nt,nz,ny,nx)
  real, intent(out) :: posPartOut(nt,npart,3), lonOut(nt,npart), &
       latOut(nt,npart), hOut(nt,npart), zOut(nt,npart)
  real :: ix(npart), iy(npart), iz(npart)
  integer :: nt1, it, it1, ip
  real :: ut(nz,ny,nx), vt(nz,ny,nx), wt(nz,ny,nx), &
       dzt(nz,ny,nx), zt(nz,ny,nx), st, interp2, interp
  
  nt1=int((nt*tinc)/dt)
  zt=zm(1,:,:,:)
  do ip=1,npart
     posPartOut(1,ip,1)=ix(ip)
     posPartOut(1,ip,2)=iy(ip)
     posPartOut(1,ip,3)=iz(ip)
     lonOut(1,ip)=interp2(lon,ny,nx,ix(ip),iy(ip))
     latOut(1,ip)=interp2(lat,ny,nx,ix(ip),iy(ip))
     hOut(1,ip)=interp(zt,nz,ny,nx,ix(ip),iy(ip),iz(ip))-interp(zt,nz,ny,nx,ix(ip),iy(ip),1.0)
     zout(1,ip)=interp(zt,nz,ny,nx,ix(ip),iy(ip),1.0)
  enddo
  
  do it=1,nt1-1
     it1=int(it*dt/tinc)
     st=it*dt/tinc-it1
     !  print*, st, it1, modulo(it-1,4)
     if (it1+2<=nt ) then
        if(modulo(it-1,1)==0) then
           ut=(1-st)*u(it1+1,:,:,:)+(st)*u(it1+2,:,:,:)
           vt=(1-st)*v(it1+1,:,:,:)+(st)*v(it1+2,:,:,:)
           wt=(1-st)*w(it1+1,:,:,:)+(st)*w(it1+2,:,:,:)
           zt=(1-st)*zm(it1+1,:,:,:)+(st)*zm(it1+2,:,:,:)
           dzt=(1-st)*dz(it1+1,:,:,:)+(st)*dz(it1+2,:,:,:)
           print*, it*dt/3600.
        endif
        if(st<1e-3) then
           do ip=1,npart
              posPartOut(it1+1,ip,1)=ix(ip)
              posPartOut(it1+1,ip,2)=iy(ip)
              posPartOut(it1+1,ip,3)=iz(ip)
              if(ix(ip)>-998 .and. iy(ip)>-998) then
                 lonOut(it1+1,ip)=interp2(lon,ny,nx,ix(ip),iy(ip))
                 latOut(it1+1,ip)=interp2(lat,ny,nx,ix(ip),iy(ip))
                 hOut(it1+1,ip)=interp(zt,nz,ny,nx,ix(ip),iy(ip),iz(ip))-&
                      interp(zt,nz,ny,nx,ix(ip),iy(ip),1.0)
                 zout(it1+1,ip)=interp(zt,nz,ny,nx,ix(ip),iy(ip),1.0)
              else
                 lonOut(it1+1,ip)=-999
                 latOut(it1+1,ip)=-999
                 hOut(it1+1,ip)=-999
                 zout(it1+1,ip)=-999
              endif
           enddo
           
        endif
        !$OMP PARALLEL DO shared(ut,vt,wt,dzt,mapf,ix,iy,iz,nz,nx,ny,dx,dy,dt) private(ip)
        do ip=1,npart
           call trace3db(ut,vt,wt,dzt,mapf,ix(ip),iy(ip),iz(ip),&
                nz,nx,ny,dx,dy,dt)
        enddo
        !$OMP END PARALLEL DO
     endif
  enddo
  
  do ip=1,npart
     posPartOut(nt,ip,1)=ix(ip)
     posPartOut(nt,ip,2)=iy(ip)
     posPartOut(nt,ip,3)=iz(ip)
     if(ix(ip)>-998 .and. iy(ip)>-998) then
        lonOut(nt,ip)=interp2(lon,ny,nx,ix(ip),iy(ip))
        latOut(nt,ip)=interp2(lat,ny,nx,ix(ip),iy(ip))
        hOut(it1+1,ip)=interp(zt,nz,ny,nx,ix(ip),iy(ip),iz(ip))
        hOut(it1+1,ip)=interp(zt,nz,ny,nx,ix(ip),iy(ip),1.0)
     else
        lonOut(nt,ip)=-999
        latOut(nt,ip)=-999
        hOut(it1+1,ip)=-999
        zout(it1+1,ip)=-999
     endif
     

  enddo
  
end subroutine trace3dball
