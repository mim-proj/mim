!
! fft_in��wrapper�ؿ�
! �̾�Ϥ������ƽФ�
!
!  N     : �ǡ�����(Ǥ��)
!  sign  : 1 or -1 (exp�θ������)
!  a     : �����ϥǡ���
!
subroutine fft(N, sign, a)
  implicit none
  integer,intent(in)       :: N
  integer,intent(in)       :: sign      
  complex(8),intent(inout) :: a(0:N-1)

  complex(8) :: t(0:N-1)
  real(8) :: theta
  real(8),parameter :: pi = 3.1415926535898

  theta = sign * 2 * pi / N

  call fft_in(N, theta, a, t)

end subroutine fft


!
! FFT�Ƶ��ؿ�
!
!  N     : �ǡ�����(Ǥ��)
!  theta : +-2 * pi / N (�����Ѵ������Ѵ������)
!  a     : �����ϥǡ���
!  t     : ����ΰ�
! 
! A_k = ��_{j=0}^{N-1} a_j exp(i theta j k / N)
! �����2,3,5...���б�
!
recursive subroutine fft_in(N, theta, a, t)
  implicit none
  integer,intent(in)       :: N
  real(8),intent(in)       :: theta
  complex(8),intent(inout) :: a(0:N-1)
  complex(8),intent(inout) :: t(0:N-1)
  integer :: base, Nb
  integer :: j, k, alpha, p

  if( N == 1 ) return

  ! ʬ���������
  if     ( mod(N, 2) == 0 ) then
     base = 2
  else if( mod(N, 3) == 0 ) then
     base = 3
  else if( mod(N, 5) == 0 ) then
     base = 5
  else if( mod(N, 7) == 0 ) then
     base = 7
  else if( mod(N, 11) == 0 ) then
     base = 11
  else 
     base = N
  end if
     
  Nb = N / base

  ! �����ʳ���ʬ��ǡ��������
  do j=0, Nb-1
     do alpha=0, base-1

        t(j+Nb*alpha) = cmplx(0.0, 0.0)
        do p=0, base-1
           t(j+Nb*alpha) = t(j+Nb*alpha) + a(j+p*Nb) &
                &        * exp( cmplx( 0.0, theta * ( j + p * Nb ) * alpha ) )
        end do

     end do     
  end do

  ! �Ƶ�
  do alpha=0, base-1
     call fft_in(Nb, base*theta, t(Nb*alpha:Nb*(alpha+1)-1), a)
  end do

  ! ��̤��Ǽ
  do k=0, Nb-1
     do p=0, base-1
        a(p+base*k) = t(k+Nb*p)
     end do
  end do

end subroutine fft_in



! N : �ǡ�����(Ǥ�դ�����������ǽ)
! thetain = +-2 * pi (�����Ѵ������Ѵ������)
! A_k = ��_{j=0}^{N-1} a_j exp(i theta j k / N)
subroutine dft(N, theta, fr, fi, ar, ai)
  implicit none
  integer,intent(in)  :: N
  real(8),intent(in)  :: theta
  real(8),intent(in)  :: fr(0:N-1), fi(0:N-1)
  real(8),intent(out) :: ar(0:N-1), ai(0:N-1)
  integer j, k
  real(8) :: temp
  
  do k=0, N-1
     ar(k) = 0
     ai(k) = 0
     do j=0, N-1
        temp = theta * real(k) * j / real(N)
        ar(k) = ar(k) + fr(j) * cos(temp) - fi(j) * sin(temp)
        ai(k) = ai(k) + fi(j) * cos(temp) + fr(j) * sin(temp)
     end do
  end do
end subroutine dft




! FFT�Υ��֥롼����
! ��®�����Ψ�����ʬ����ˤ�����
! thetain = +-2 * pi (�����Ѵ������Ѵ������)
! A_k = ��_{j=0}^{N-1} a_j exp(i theta j k / N)
subroutine fft_quick(N, thetain, ar, ai)
  implicit none
  integer,intent(in) :: N
  real(8),intent(in) :: thetain
  real(8),intent(inout) :: ar(0:N-1), ai(0:N-1)

  integer :: m, mh
  integer :: i, j, k, L, NLOG
  real(8) :: wr, wi, xr, xi
  real(8) :: tmp
  real(8) :: theta
  theta = thetain

! srambler
  i=0
  do j=1, N-1
    k = N / 2
    do while( k <= i )
      i = i - k
      k = k / 2
    end do
    i = i + k
    if( j < i ) then
      tmp = ar(j)
      ar(j) = ar(i)
      ar(i) = tmp
      tmp = ai(j)
      ai(j) = ai(i)
      ai(i) = tmp
    end if
  end do

  NLOG = int( log(real(N+1)) / log(2.0) ) - 1 
!  write(*,*) "N=",N,"NLOG=",NLOG

  do L=0, NLOG
    mh = 2**L
    m = 2**(L+1)
    theta = theta * 0.5
!    write(*,*) "m=",m,"mh=",mh,"theta=",theta
    do i=0, mh-1
      wr = cos( theta*i )
      wi = sin( theta*i )
      do j=i, N-1, M
        k = j + mh
        xr = wr * ar(k) - wi * ai(k)
        xi = wr * ai(k) + wi * ar(k)
        ar(k) = ar(j) - xr
        ai(k) = ai(j) - xi
        ar(j) = ar(j) + xr
        ai(j) = ai(j) + xi
      end do
    end do
  end do
  return
end subroutine fft_quick




