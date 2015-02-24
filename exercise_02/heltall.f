      program heltall
c Integer representation
      dimension in(32)
      write(*,*) 'input='
      read (*,*) ia
      do i=1,32
      in(i)=iand(1,ishft(ia,1-i))
      enddo
      write(*,1) (in(i),i=32,1,-1)
1     format(1x,32(i1))
      end
