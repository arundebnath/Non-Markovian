      program test

      integer*16 i,j
      real*8 k
      real*8 mat(4,4,4)
      integer*16  vec(4)
        
      
c       read(*,*) k
c      read(*,*) k
      open(30,file='test')      
      do i=1,4

c      vec(i)=k
 

c       do j=1,4
       read(*,*) k
       mat(:,i,:)=k
    
      write(30,*) mat(:,i,:)

c     enddo
c      write(30,*) vec(i)
      enddo
      write(*,*) mat(3,:,:)
      close(30)
      end program test    
