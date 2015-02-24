      program realtall
c Real representation
2     continue      
      write(*,*) 'real ='
      read (*,*) a
      call rint(a,b)
      write(*,*) b
      goto 2
      end
c
      subroutine rint(ia,ib)
      dimension in(32)
      do i=1,32
      in(i)=iand(1,ishft(ia,1-i))
      enddo
      write(*,1) (in(i),i=32,1,-1)
1     format(1x,32(i1))
      ib=ia
      return
      end
