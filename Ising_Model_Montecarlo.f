c 		MARTA XIULAN ARIBÓ HERRERA
c 		PROGRAMA DE MONTECARLO

c--------------------------- DEFINIM VARIABLES   ------------------------

        PROGRAM MONTECARLO
        IMPLICIT NONE
C 		DIMENSIÓ DEL SISTEMA 
        INTEGER*4 L
        PARAMETER (L=8) 

C 		MATRIU DE SPINS
        INTEGER*2 S(1:L,1:L) 

C 		NOM DE L'ARXIU
        CHARACTER*18 NOM

C       COMPTADORS 
        INTEGER*4 I,J,K,N,MCTOT,MCINI,MCD,IMC,IPAS,DE, H,NDAT
        PARAMETER (N = L*L)

        INTEGER*4  NSEED, SEED0, SEED
C       NUMERO DE TEMPERATURES I VECTOR DE TEMPERATURES
        INTEGER*4 NTEMP
        parameter (NTEMP=200)
        REAL*8 T(1:NTEMP), V_E(1:NTEMP), V_DE(1:NTEMP)


        REAL*8 W(-8:8) !tabulem els valors DE

C       FUNCIONS
        REAL*8 genrand_real2, magnetitzacio, energia
        REAL*8 magne, energ !funcions

C       ATRES VARIABLES NECESSÀRIES
        REAL*8 TEMP, x !temperatura
        REAL*8 suma, delta !ENEBIS 
        INTEGER*4 xx,yy

C       PER LES CONDICIONS DE CONTORN
        INTEGER*4 PBC(0:L+1) 

C       PER FER PROMITJOS
        REAL*8 sum,sume,sume2,summ,summ2,sumam,vare,varm !promitjos
        REAL*8 errorm, errore, cv, xi
C 		CALCUL DEL TEMPS DE EXECUCIÓ
        real*4 TIME1, TIME2
        CHARACTER*30 DATE
        call CPU_TIME(TIME1)

        call FDATE(DATE)

        WRITE(*,*) DATE


C-----------------------------------------------------------------------

c---------------------- Definirem les condicions de contorn--------------
c       Es considera una xarxa periòdica 
c-----------------------------------------------------------------------
        PBC(0) = L
        PBC(L+1) = 1
        do i = 1, L
          PBC(i) = i
        enddo
c--------------------------- PROGRAMA PRINCIPAL  -----------------------
        NOM ="MC-L-008-TEMP-0200"
        NSEED = 200
        SEED0 = 123469
        MCTOT = 40000 !pasos de montecarlo
        MCINI = 2000 !saltem algunes passes
        MCD = 20 !agafem configuracions cada cert nombre de passes

c------------------------ VECTOR DE TEMPERATURES -----------------------
        T(1) = 1.4d0
        do K = 2, NTEMP
            T(k)= 1.4d0 + 0.01*k
        enddo

c--------------------------- MONTECARLO --------------------------------
        DO K =1,NTEMP
          TEMP = T(K)

          !per fer els promitjos
          sum = 0.d0 
          sume = 0.d0 
          sume2 = 0.d0
          summ = 0.d0
          summ2 = 0.d0
          sumam = 0.d0
          vare=0.0D0
          varm=0.0D0

c      Podem tabular els valors de l'exponencial, donat que DE pren uns 
c      pocs valors.

        do DE = -8,8
          W(DE) = dexp(-dfloat(DE)/TEMP)
        enddo

        H = 0

        DO SEED = SEED0, SEED0+NSEED-1, 1 !bucle per les diferents llavors
          H = H+1 !comptador extra
          write(*,*) "Llavor = ", SEED, "#Temperatura", K, H, L
          call init_genrand(SEED) 

          !generació de la matriu spins
          do i = 1, L 
              do j = 1, L
                  x = genrand_real2()
                  if (x.LT.0.5d0) then
                      S(i,j) = 1
                  else
                      S(i,j) = -1
                  endif
              enddo
          enddo
                  
          energia = ENERG(S,L,PBC)


        DO IMC = 1, MCTOT !montecarlo

          DO  IPAS = 1, N

            delta = genrand_real2() !per comparar
            xx = int(genrand_real2()*L)+1 !posicions de la matriu random
            yy = int(genrand_real2()*L)+1

            suma=S(xx,PBC(yy+1))+S(xx,PBC(yy-1))+S(PBC(xx+1),yy)
     *+S(PBC(xx-1),yy)

            DE = int(2*S(xx,yy)*suma) !diferencia d'energia

            !condicions
          if (DE.LT.0) then
              S(xx,yy) = - S(xx,yy)
              Energia = Energia + DE
          endif
          if (DE.GT.0) then
              if (delta.LT.w(DE)) then
                S(xx,yy) = - S(xx,yy)
                Energia = Energia + DE
              endif

          endif

          ENDDO !bucle ipas

        if ((IMC.GT.MCINI).AND.(MCD*(IMC/MCD).EQ.IMC)) then
          
          magnetitzacio = MAGNE(S,L)
          sum = sum + 1.d0
          sume = sume + energia
          sume2 = sume2 + energia*energia

          summ = summ + magnetitzacio
          sumam = sumam + abs(magnetitzacio)
          summ2 = summ2 + magnetitzacio*magnetitzacio

        endif

        ENDDO !bucle imc

        ENDDO !bucle seed

        !normalització de  promitjos

        sume = (sume/sum)/real(N) !<e>
        sume2 = (sume2/sum)/(real(N)*real(N)) !<e^2>
        summ = (summ/sum)/real(N) !<m>
        sumam = (sumam/sum)/real(N) !<|m|>
        summ2 = (summ2/sum)/((real(N)*real(N))) !<m^2>
        vare = sume2-sume*sume !var(e)=<e^2>-<e>^2
        varm = summ2-summ*summ !var(m)=<m^2>-<m>^2
        errore = dsqrt(vare/sum)/real(N)
        errorm = dsqrt(varm/sum)/real(N)
        cv = real(N)*((sume2-(sume*sume))/(TEMP*TEMP))
        xi = real(N)*((summ2-(sumam*sumam))/(TEMP))
        V_E(k) = sume
        !surtida de dades

        open (1, file = nom//".dat")

                   !1  2   3   4     5     6    7 
        write(1,*) L,TEMP,sum,sume,sume2,vare,errore,
     *summ,sumam,summ2,varm,errorm, cv, xi 
!       8    9     10   11    12    13 14
        print*, "Temperatura = ", TEMP
        print*, "<E>/N =",sume/N
        print*, "eps_e =", sqrt(sume2-sume*sume)/(N*sqrt(sum))
        print*, "<m> =",summ/N,"<|m|> =",sumam/N,"(<m^2>)^(1/2)",
     *sqrt(summ2/N**2)
        print*, "eps_m =", sqrt(summ2-summ*summ)/(N*sqrt(sum))
        print*, "cv = ", cv, "xi =", xi
C        close(1)

        ENDDO !temperatura bucle
        open (2, file = "derivades_08.dat")
        NDAT = NTEMP
        call derfun(NTEMP,T,V_E,V_DE)

        DO k = 1, NDAT
          write(2,*) T(k), V_DE(k)
        ENDDO

        call CPU_TIME(TIME2)



        WRITE(*,*) DATE
        WRITE(*,*) "CPUTIME = ", TIME2-TIME1


        END PROGRAM 

c--------------------------- Subrutines i funcions ---------------------
c-----------------------------------------------------------------------
c       Funció que calcula la magnetització del sistema
c-----------------------------------------------------------------------
        REAL*8 function magne(S,L)
        implicit NONE
        integer*2 S(1:L,1:L)
        integer*4 I,J,L
        double precision m
        m = 0.d0
        do i = 1,L
            do j = 1, L
                m = m + S(i,j)

            enddo
        enddo
        magne = m
        end 
c-----------------------------------------------------------------------
c       Funció que calcula l'energia del sistema
c-----------------------------------------------------------------------

        REAL*8 function ENERG(S,L,PBC)
        implicit NONE
        integer*2 S(1:L,1:L)
        integer*4 I,J,L
        integer*4 PBC(0:L+1)
        double precision E


        E = 0.d0
        do i = 1, L
          do j = 1,L
            E=E-S(i,j)*S(PBC(i+1),j)-S(i,j)*S(i,PBC(j+1))
            !calcila l'energia de l'espi de la dreta i d'abaix

          enddo
        enddo
        ENERG = E


        end 
c-----------------------------------------------------------------------
c       Subrutina que calcula la derivada
c-----------------------------------------------------------------------
c        NDAT = NTEMP
c        call derfun(NTEMP,T,V_E,V_DE)

        subroutine derfun(ndat,x,fu,dfu)
        implicit none 
        integer ndat, i
        double precision x(ndat), fu(ndat), dfu(ndat), h
        h = x(ndat)-(x(ndat-1)) !intervals constants
        print*,h 
        do i = 1, ndat
          if ((i.LT. ndat).and.(i.GT.1)) then
            dfu(i)= (fu(i+1)-fu(i-1))/(2*h)
          endif
          if (i.EQ.1) then
            dfu(i) = (fu(i+1)- fu(i))/h
          endif
          if (i.EQ.ndat) then
            dfu(i) = (fu(i)-fu(i-1))/h
          endif
        enddo
        return
        end

 
c--------------------------- GENERADOR MT 19937 ------------------------
c  A C-program for MT19937, with initialization improved 2002/1/26.
c  Coded by Takuji Nishimura and Makoto Matsumoto.
c
c  Before using, initialize the state by using init_genrand(seed)  
c  or init_by_array(init_key, key_length).
c
c  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
c  All rights reserved.                          
c  Copyright (C) 2005, Mutsuo Saito,
c  All rights reserved.                          
c
c  Redistribution and use in source and binary forms, with or without
c  modification, are permitted provided that the following conditions
c  are met:
c
c    1. Redistributions of source code must retain the above copyright
c       notice, this list of conditions and the following disclaimer.
c
c    2. Redistributions in binary form must reproduce the above copyright
c       notice, this list of conditions and the following disclaimer in the
c       documentation and/or other materials provided with the distribution.
c
c    3. The names of its contributors may not be used to endorse or promote 
c       products derived from this software without specific prior written 
c       permission.
c
c  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
c  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c
c  Any feedback is very welcome.
c  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
c  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
c
c-----------------------------------------------------------------------
c  FORTRAN77 translation by Tsuyoshi TADA. (2005/12/19)
c
c     ---------- initialize routines ----------
c  subroutine init_genrand(seed): initialize with a seed
c  subroutine init_by_array(init_key,key_length): initialize by an array
c
c     ---------- generate functions ----------
c  integer function genrand_int32(): signed 32-bit integer
c  integer function genrand_int31(): unsigned 31-bit integer
c  double precision function genrand_real1(): [0,1] with 32-bit resolution
c  double precision function genrand_real2(): [0,1) with 32-bit resolution
c  double precision function genrand_real3(): (0,1) with 32-bit resolution
c  double precision function genrand_res53(): (0,1) with 53-bit resolution
c
c  This program uses the following non-standard intrinsics.
c    ishft(i,n): If n>0, shifts bits in i by n positions to left.
c                If n<0, shifts bits in i by n positions to right.
c    iand (i,j): Performs logical AND on corresponding bits of i and j.
c    ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
c    ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
c
c-----------------------------------------------------------------------
c     initialize mt(0:N-1) with a seed
c-----------------------------------------------------------------------
      subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
c
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*
     &          ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
c
      return
      end
c-----------------------------------------------------------------------
c     initialize by an array with array-length
c     init_key is the array for initializing keys
c     key_length is its length
c-----------------------------------------------------------------------
      subroutine init_by_array(init_key,key_length)
      integer init_key(0:*)
      integer key_length
      integer N
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      parameter (N=624)
      integer i,j,k
      integer mt(0:N-1)
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
c
      call init_genrand(19650218)
      i=1
      j=0
      do 100 k=max(N,key_length),1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)
     &           +init_key(j)+j
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        j=j+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
        if(j.ge.key_length)then
          j=0
        endif
  100 continue
      do 200 k=N-1,1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
  200 continue
      mt(0)=TOPBIT_MASK
c
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0xffffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
c
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
c
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
c
      y=mt(mti)
      mti=mti+1
c
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
c
      genrand_int32=y
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0x7fffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int31()
      integer genrand_int31
      integer genrand_int32
      genrand_int31=int(ishft(genrand_int32(),-1))
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1]-real-interval
c-----------------------------------------------------------------------
      function genrand_real1()
      double precision genrand_real1,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real1=r/4294967295.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on (0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real3()
      double precision genrand_real3,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real3=(r+0.5d0)/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1) with 53-bit resolution
c-----------------------------------------------------------------------
      function genrand_res53()
      double precision genrand_res53
      integer genrand_int32
      double precision a,b
      a=dble(ishft(genrand_int32(),-5))
      b=dble(ishft(genrand_int32(),-6))
      if(a.lt.0.d0)a=a+2.d0**32
      if(b.lt.0.d0)b=b+2.d0**32
      genrand_res53=(a*67108864.d0+b)/9007199254740992.d0
      return
      end
c-----------------------------------------------------------------------
c     initialize large number (over 32-bit constant number)
c-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
CC    TOPBIT_MASK = Z'80000000'
CC    ALLBIT_MASK = Z'ffffffff'
CC    UPPER_MASK  = Z'80000000'
CC    LOWER_MASK  = Z'7fffffff'
CC    MATRIX_A    = Z'9908b0df'
CC    T1_MASK     = Z'9d2c5680'
CC    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end
