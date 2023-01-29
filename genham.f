     





          do i1=1,n
          do i2=1,n  
            do j1=1,n
            do j2=1,n
          
          read(*,*)  i1,hhs1(i1,:)
          read(*,*)  i2,hhs2(:,i2)
          
          hhs()=(0.0d0,0.0d0,0.0d0,0.0d0)
          mu(i1,i2,j1,j2)=mu1(i1,i2,:,:)*mu2(:,:,j1,j2)
          hhf(i1,i2,j1,j2)=mu(i1,i2,j1,j2)/dlat**3
          hhcoup(i1,i2,j1,j2)= mu(i1,i2,:,:)*efield+efield*

          hhs(i1,i1,j1,j2)= hhs1(i1,i2,:,:)+hhs(:,:,j1,j2)+
     .                      hhf(i1,i2,j1,j2)+hhcoup((i1,i2,j1,j2)

