subroutine compute_grid_integrated_heat(result,tensor, lon,lat, numPdY, numPdX, nlen, npart)

    real *8, intent(in):: lat(numPdY+1,numPdX+1)
    real *8, intent(in) :: lon(numPdY+1,numPdX+1)
    real *8, intent(in) :: tensor(nlen,npart,6)
    real *8, intent(out) :: result(numPdY,numPdX)
    integer :: i,j,k,n, nlen, npart, numPdY, numPdX
  
  
 
   do i=1,numPdY
        do j=1,numPdX
           result(i,j)=0.0
        enddo
    enddo
    
    do i=1,nlen
        do j=1,npart
            do k=1, numPdY
                do n=1, numPdX
                    if (tensor(i,j,1)/=int(-999.9) .and. tensor(i,j,1)/=int(-999.9)) then
                        
                        if (tensor(i,j,1) .gt. lon(k,n)  .and. tensor(i,j,1) .lt. & 
                         lon(k+1,n+1) .and.  tensor(i,j,2) .gt. lat(k, n) .and. tensor(i,j,2) .lt. lat(k+1,n+1)) then
                         result(k,n)=result(k,n)+tensor(i,j,4)*tensor(i,j,6)
                         !write(*,*)tensor(i,j,4),tensor(i,j,5)
                        endif 
                    endif
                enddo
            enddo
        enddo
    enddo
    
    
end subroutine


subroutine compute_grid_integrated_moist(result,tensor, lon,lat, numPdY, numPdX, nlen, npart)

    real *8, intent(in):: lat(numPdY+1,numPdX+1)
    real *8, intent(in) :: lon(numPdY+1,numPdX+1)
    real *8, intent(in) :: tensor(nlen,npart,6)
    real *8, intent(out) :: result(numPdY,numPdX)
    integer :: i,j,k,n, nlen, npart, numPdY, numPdX
  
  
 
   do i=1,numPdY
        do j=1,numPdX
           result(i,j)=0.0
        enddo
    enddo
    
    do i=1,nlen
        do j=1,npart
            do k=1, numPdY
                do n=1, numPdX
                    if (tensor(i,j,1)/=int(-999.9) .and. tensor(i,j,1)/=int(-999.9)) then
                        
                        if (tensor(i,j,1) .gt. lon(k,n)  .and. tensor(i,j,1) .lt. & 
                         lon(k+1,n+1) .and.  tensor(i,j,2) .gt. lat(k, n) .and. tensor(i,j,2) .lt. lat(k+1,n+1)) then
                         result(k,n)=result(k,n)+tensor(i,j,3)*tensor(i,j,6)
                         !write(*,*)tensor(i,j,4),tensor(i,j,5)
                        endif 
                    endif
                enddo
            enddo
        enddo
    enddo
    
    
end subroutine




subroutine search_row(output,matrix,lista,len_lista,numP)
   
    real, intent(in) :: matrix(numP,11)
    integer, intent(in) :: lista(len_lista)
    real, intent(out) :: output(len_lista, 11)
    integer :: i,j, numP
    output(:,:)=-999.9
    do i=1, len_lista
        do j=1, numP
            if (int(matrix(j,1))==lista(i)) then
              output(i,:)=matrix(j,:)
           endif 
        enddo
    enddo
end subroutine

subroutine determined_id(vector, value_mascara,value_mask,len_value_mascara)

    integer, intent(in) :: value_mascara(len_value_mascara)
    integer, intent(out) :: vector(len_value_mascara)
    integer :: i,j, len_value_mascara, value_mask
   
    do i=1,len_value_mascara
       vector(i)=-999
    enddo
   
    do j=1,len_value_mascara
        if (value_mascara(j)==value_mask) then
           vector(j)=j-1
        endif
    enddo
end subroutine


subroutine read_binary_file(output_,filename, nparts,x_l,y_l, x_r,y_r)

    integer,parameter :: rk=kind(1.0)
    real(rk)         :: b1
    integer          :: id1,i,k,idx,index_id, ty
    integer(kind=4)  :: n1,nparts, cant
    integer bytes,aux_bytes
    real *8  :: matrix(nparts,13)
    real *8 :: output(nparts,11)
    !real, parameter :: g0 = 9.780327 !m/s**2
    !real, parameter :: cp=1005.7   !J/kgK
    !real, parameter :: pi = 3.1415927
    real *8, intent(out) :: output_(nparts,11)
    character(500) :: filename
    real, intent(in) :: x_l,y_l, x_r,y_r
    
    open(1,file=filename,access="stream",form="unformatted")
    bytes=(nparts+1)*60 +12
    read(1, pos=1)ty

    if (ty==4)then
       i=17
       aux_bytes=69
    else   
       i=25
       aux_bytes=77
    endif
   
    index_id=1
    do while (i<bytes-60)
        read(1,pos=i) id1 
        matrix(index_id,1)=id1
       idx=2
        do k=i+4,aux_bytes-1,4
            read(1,pos=k) n1
            read(1,pos=k) b1 
            matrix(index_id,idx)=b1
            idx=idx+1
        end do
        i=i+60
        aux_bytes=aux_bytes+60
        index_id=index_id+1
    end do

  
    output(:,1)=matrix(:,1)  ! parcel id
    output(:,2)=matrix(:,2)  ! parcel longitude
    output(:,3)=matrix(:,3)  ! parcel latitude
    output(:,4)=matrix(:,8)  ! parcel specific humidity
    output(:,5)=matrix(:,4)  ! parcel vertical heigh
    output(:,6)=matrix(:,6)  ! topography
    output(:,7)=matrix(:,9)  !  parcel density
    output(:,8)=matrix(:,10)  ! pbl high
    output(:,9)=matrix(:,11)  ! tropopause high
    output(:,10)=matrix(:,12) !  pacel temperature 
    output(:,11)=matrix(:,13) !parcel mass
    !output(:,12)=(cp*matrix(:,11) + g0*(1 + 0.0053024 * (DSIN(matrix(:,3)*pi/180 ))**2  &
    !              - 0.0000058*(DSIN( 2*matrix(:,3)*pi/180))**2)*matrix(:,4))
                  
                  
    output_(:,:)=-999

    cant=1
    do j=1, nparts-1
        if ( output(j,2) .ge. x_l .and. output(j,2) .le. x_r .and. &
           output(j,3) .ge. y_l .and. output(j,3) .le. y_r) then
           output_(cant,:)=output(j,:)
           cant=cant+1
        endif
    enddo

   
return
end subroutine read_binary_file





subroutine len_file(bytes, filename)

    character(500) :: filename
    character(len=256) :: message
    integer :: ios,x
    logical :: foundit
    integer, intent(out) :: bytes
    
    inquire(file=filename,exist=foundit,size=x,iostat=ios,iomsg=message)
    bytes=x
end subroutine



