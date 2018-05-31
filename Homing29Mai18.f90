!Publicação 29/05/18
	program Homing
	IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N) 
	Integer, ALLOCATABLE :: Mp(:,:) , Mpr(:,:) , Mpd(:,:), Msoma(:,:)
	Integer, ALLOCATABLE :: Mmae(:), Mpai(:), Nfilho(:), MX(:)
	character*30 nomerep
	common /c1/ nTpeixe,ngenoma,rmu,nG
	common /c2/ ncomp
	call random_seed() 
	OPEN(UNIT=60,FILE='Input_dam.in',STATUS='old')
	Read(60,*) ntime
	Read(60,*) ncomp, k
	Read(60,*) ngenoma, nG, rmu  
	Read(60,*) rmor,rmig 
	Read(60,*) nsalva
	Read(60,*) nrep
	Read(60,*) nb
	Read(60,*) ntimebar   
	Read(60,*) permjm		
	Read(60,*) permmj		
	Close(60)

Do nr=1, nrep !Replicates
print*, 'Run', nr
	write(nomerep,*) nr

	
	OPEN(UNIT=70,FILE='Fst'//TRIM(adjustl(nomerep))//'.txt',STATUS='UNKNOWN')
	write(70,*)'Gen ', 'N_Trib ','Extant_Pops ', 'Tot_Pop_size ', 'Global_Fst ', 'Fst_Up ', 'Fst_Down'
	OPEN(UNIT=80,FILE='Pairwise_Fst'//TRIM(adjustl(nomerep))//'.txt',STATUS='UNKNOWN')
	write(80,*)'Gen ', 'N_comparisons(N_trib^2) ', 'PopAxPopA, PopAxPopB, AC, AD, ..., BA, BB, BC, BD, ..., CA, CB, CC, ..., NN'
	
!Creation of the Fish

	nTpeixe=k
	
	ALLOCATE (Mp(nTpeixe,ngenoma+1),Mpr(nTpeixe,ngenoma+1),Mpd(nTpeixe,ngenoma+1), Msoma(nTpeixe,ngenoma+1))
	ALLOCATE (Mmae(nTpeixe),Mpai(nTpeixe),Nfilho(ngenoma+1))
	ncont=0
	Mp=0	
	do ind=1,nTpeixe
		call random_number(rand)
		nposicao=1+Int(rand*ncomp)
		Mp(ind,1)= nposicao
		Do j=2,ngenoma+1
			call random_number(rand)
                        ngg=Nint(rand)
			Mp (ind,j)=ngg
		Enddo
	enddo

!Calculating Fst
n=0
call Fst(Mp,n,nvaz, fstt)
call FstUD(Mp,n,nb,Fstup, Fstdown)
		write(70,11) n, ncomp,ncomp-nvaz, nTpeixe, fstt, Fstup, Fstdown
11	FORMAT((I4),(2X,I5),(4X,I5),(7x,I5),(6X, f6.4),(3X, f6.4),(2X, f6.4))
call Pairwise(Mp,n)
	
	!Dynamics
	Do n=1, ntime
			!Damming
			If(n.eq.ntimebar) then
				pabaixo=1.0-Real(nb)/Real(ncomp)
				kbaixo=(k*real(pabaixo))
				kcima=k-kbaixo
				Do i=1, nTpeixe
					call random_number(rand)
					if(rand.le.pabaixo) then
						if(Mp(i,1).le.nb)then
							if(nb.eq.ncomp-1) then
								Mp(i,1)=ncomp
							else
							call random_number(rand)
							if(rand.le.0.5) then
								Mp(i,1)=nb+1
							else
								Mp(i,1)=nb+2
							endif
							endif
						endif
					else
						if(Mp(i,1).gt.nb)then
							if(nb.eq.1) then
								Mp(i,1)=1
							else
								call random_number(rand)
								if(rand.le.0.5) then
									Mp(i,1)=nb
								else
									Mp(i,1)=nb-1
								endif
							endif
						endif	
					endif
				Enddo
			Endif
			
		!Dispersal through the dam
		If(n.ge.ntimebar) then	
			Do i=1, nTpeixe
				if(Mp(i,1).le.nb) then	
						call random_number(rand)
					if(rand.le.permmj) then
						if(nb.eq.ncomp-1) then
							Mp(i,1)=ncomp
						else
							call random_number(rand)
							if(rand.le.0.5) then
								Mp(i,1)=nb+1
							else
								Mp(i,1)=nb+2
							endif
						endif
					endif
				else
				call random_number(rand)
								ntr=Mp(i,1)
					if(rand.le.permjm) then
						if(nb.eq.1) then
							Mp(i,1)=1
						else
							call random_number(rand)
							if(rand.le.0.5) then
								Mp(i,1)=nb
							else
								Mp(i,1)=nb-1
							endif
						endif
					endif
				endif
				ntr=Mp(i,1)
			Enddo
		endif
		
		!Spawning migration
		if(n.lt.ntimebar) then
			Do i=1, nTpeixe
				call random_number(rand)
				if(rand.le.rmig) then
					new=Mp(i,1)
					Do while (new.eq.Mp(i,1))
						call random_number(rand)
						new=1+Int(ncomp*rand)
					enddo
					Mp(i,1)=new	
				endif
			Enddo
		else	
			Do i=1, nTpeixe
				call random_number(rand)
				ntr=Mp(i,1)
				if(rand.le.rmig) then
					new=Mp(i,1)
					if(Mp(i,1).le.nb)then
						if(nb.gt.1) then
							Do while (new.eq.Mp(i,1))
								call random_number(rand)
									new=1+Int(nb*rand)
							Enddo
						endif
					else
						if(nb.lt.ncomp-1) then
							Do while (new.eq.Mp(i,1))
								call random_number(rand)
								new=nb+1+Int((ncomp-nb)*rand)
							Enddo
						endif

					endif
				Mp(i,1)=new
				endif
			Enddo
		Endif
		
		!Mortality
		DeaLLOCATE (Mpd)
		ALLOCATE (Mpd(nTpeixe,ngenoma+1))
		nsobrevive=0
		Do i=1, nTpeixe
			call random_number(rand)
	  		if(rand.gt.rmor) then
				nsobrevive=nsobrevive+1
				Mpd(nsobrevive,:)=Mp(i,:)	
			endif
		end Do
		If(nsobrevive.eq.0) go to 22 !!!!!stop

		DeaLLOCATE (Mp)
		ALLOCATE (Mp(nsobrevive,ngenoma+1))
		Mp=Mpd(1:nsobrevive,:)
		nTpeixe=nsobrevive
		
		!Reproduction
		if(n.lt.ntimebar) then
			if(nTpeixe.lt.k) then
				DeaLLOCATE (Mpr)
				ALLOCATE (Mpr(nTpeixe,ngenoma+1))
				Deallocate(Mmae)
				Allocate(Mmae(nTpeixe))
				call shuffle(Mmae,nTpeixe)
				ncontfilho=0		
				Do i=1, nTpeixe 
					imae=Mmae(i)
					Deallocate(Mpai)
					Allocate(Mpai(nTpeixe))
					call shuffle(Mpai,nTpeixe)
					call offspring(Mp,Mpai,Nfilho,imae,nachou)

					if(nachou.eq.1) then				
						ncontfilho=ncontfilho+1
						Mpr(ncontfilho,:)=Nfilho
						if(nTpeixe+ncontfilho.ge.k) go to 12
					endif 
				Enddo
	12 			continue
				DeaLLocate (Msoma)
				ALLOCATE (Msoma(nTpeixe+ncontfilho,ngenoma+1))
				Msoma(1:nTpeixe,:)=Mp
				Msoma(nTpeixe+1:nTpeixe+ncontfilho,:)=Mpr(1:ncontfilho,:)
				DeaLLocate (Mp)
				ALLOCATE (Mp(nTpeixe+ncontfilho,ngenoma+1))
				nTpeixe=nTpeixe+ncontfilho
				Mp=Msoma
			endif
		else
			nTpeixeB=0
			nTpeixeC=0
			Do i=1, nTpeixe
				if(Mp(i,1).le.nb) nTpeixeC=NtpeixeC+1
			Enddo
			nTpeixeB=nTpeixe-nTpeixeC
			DeaLLOCATE (Mpr)
			ALLOCATE (Mpr(nTpeixe,ngenoma+1))
			Deallocate(Mmae)
			Allocate(Mmae(nTpeixe))
			call shuffle(Mmae,nTpeixe)
			ncontfilho=0		
			Do i=1, nTpeixe
				imae=Mmae(i)
				lmae=Mp(imae,1)
				nrepro=0
				if(lmae.le.nb) then
				   if(nTpeixeC.lt.kcima) nrepro=1
				else
				   if(nTpeixeB.lt.kbaixo) nrepro=1	
				endif
				if(nrepro.eq.1)then
					Deallocate(Mpai)
					Allocate(Mpai(nTpeixe))
					call shuffle(Mpai,nTpeixe)
					call offspring(Mp,Mpai,Nfilho,imae,nachou)
					if(nachou.eq.1) then				
						ncontfilho=ncontfilho+1
						Mpr(ncontfilho,:)=Nfilho
						if(Nfilho(1).le.nb) then
							nTpeixeC=nTpeixeC+1
							else
							nTpeixeB=nTpeixeB+1
						endif
						if(nTpeixeC.eq.kcima.and.nTpeixeB.eq.Kbaixo) go to 13
					endif 
				endif
			Enddo
13 			continue
			DeaLLocate (Msoma)
			ALLOCATE (Msoma(nTpeixe+ncontfilho,ngenoma+1))
			Msoma(1:nTpeixe,:)=Mp
			Msoma(nTpeixe+1:nTpeixe+ncontfilho,:)=Mpr(1:ncontfilho,:)
			DeaLLocate (Mp)
			ALLOCATE (Mp(nTpeixe+ncontfilho,ngenoma+1))
			nTpeixe=nTpeixe+ncontfilho
			Mp=Msoma
		endif
		if(mod(n,nsalva).eq.0) then 
		call Fst(Mp,n,nvaz, fstt)
		call FstUD(Mp,n,nb,Fstup, Fstdown)
		write(70,11) n, ncomp,ncomp-nvaz, nTpeixe, fstt, Fstup, Fstdown
		call Pairwise(Mp,n)		
		endif	
		
	Enddo
	Close(70)
22	continue
	Close(80)
	DeALLOCATE (Mp,Mpr,Mpd, Msoma)
	DeALLOCATE (Mmae,Mpai,Nfilho)

Enddo

	END PROGRAM Homing
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	SUBROUTINE init_random_seed()
        INTEGER :: l, k, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        CALL RANDOM_SEED(size = k)
        ALLOCATE(seed(k))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * (/ (l - 1, l = 1, k) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)
	END SUBROUTINE init_random_seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	subroutine shuffle (NV,Ni)
	IMPLICIT REAL*8 (A-H,O-Z)
        Dimension::NV(Ni)
	  Do n=1,Ni
		NV(n)=n
	  Enddo
	  Do n=1, Ni
		call random_number(h1)
		call random_number(h2)
		i1=1+int(h1*Ni)
		i2=1+int(h2*Ni)
		vh1=NV(i1)
		NV(i1)=NV(i2)
		NV(i2)=vh1
	 Enddo	  	
      endsubroutine shuffle
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	subroutine offspring(Mp,Mpai,Nfilho,imae,nachou)
	IMPLICIT REAL*8 (A-H,O-Z)
	Integer, ALLOCATABLE :: Mrep(:)
	common /c1/ nTpeixe,ngenoma,rmu,nG
	Dimension::Mp(nTpeixe,ngenoma+1)
	Dimension::Mpai(nTpeixe)
	Dimension::Nfilho(ngenoma+1)
	ALLOCATE (Mrep(ngenoma))
	nachou=0
	Mrep=0
	
	lmae=Mp(imae,1)
	
	Do ip=1, nTpeixe
		ipai=Mpai(ip)
		if(ipai.ne.imae) then
			lpai=Mp(ipai,1)
			if(lpai.eq.lmae) then
				kdif=0
			   	Do j=2,ngenoma+1
					kdif=kdif+abs(Mp(imae,j)-Mp(ipai,j))
			  	enddo	
			   	if(kdif.le.nG) then
					nachou=1
					Nfilho(1)=lpai
					call shuffle (Mrep,ngenoma)	
					Do j=1,ngenoma
						if(j.le.(ngenoma/2)) then 
							jj=1+Mrep(j)
							Nfilho(jj)=Mp(ipai,jj)
							call random_number(rand)
							if(rand.le.rmu) Nfilho(jj)=abs(Nfilho(jj)-1)
							
						else
							jj=1+Mrep(j)
							Nfilho(jj)=Mp(imae,jj)
							call random_number(rand)
							if(rand.le.rmu) Nfilho(jj)=abs(Nfilho(jj)-1)
						endif
					Enddo	
					go to 21			     	
				endif
			endif
		endif
	Enddo
 21    continue 	
       endsubroutine offspring
	   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	subroutine Fst(Mp,n,nvaz, fstt)
	IMPLICIT REAL*8 (A-H,O-Z)
	common /c1/ nTpeixe,ngenoma,rmu,nG
	common /c2/ ncomp
	Dimension::Mp(nTpeixe,ngenoma+1)
	Integer, ALLOCATABLE :: NK(:)
	ALLOCATE (NK(ncomp))
	NK=0 
	Do np=1, NTpeixe 
		nposi=Mp(np,1)
		NK(nposi)=NK(nposi)+1
	Enddo
       
	ni2s=0 
	nvaz=0
	Do nr=1, ncomp
		if(NK(nr).gt.0)then 
		        ni=0
			Do np=1, NTpeixe 
				if(Mp(np,1).eq.nr) ni=ni+1	
			Enddo
			ni2s=ni2s+ni**2
		else
		nvaz=nvaz+1
		endif
	Enddo
	frac=real(ni2s)/real(nTpeixe)
	cn=1.0d0/real(ncomp-nvaz-1)*(real(nTpeixe)-frac)
    
	vMSPl=0.0d0				
	vMSGl=0.0d0				
	vSigp=0.0d0				
	Fstl=0.0d0				
	Fstt=0.0d0				
	tMSG=0.0d0				
	tSig=0.0d0			
	Do l=2,ngenoma+1		
		pi2s=0.0d0
		
		npt=0
		Do np=1, nTpeixe
			npt=npt+Mp(np,l)	
		enddo
		pt=real(npt)/real(nTpeixe)
		
		vMSPl=0.0d0			
		p2s=0.0d0 			
		p2=0.0d0 			
		Do nr=1, ncomp
			if(NK(nr).gt.0) then
				npi=0
				ncont=0
				Do np=1, nTpeixe
					if(Mp(np,1).eq.nr)then
						npi=npi+Mp(np,l)
						ncont=ncont+1	
					endif
				Enddo
				pi=real(npi)/real(ncont)
				pi2s=pi2s+real(ncont)*(pi-pt)**2	
				
				p2=real(ncont)*pi*(1.0d0-pi)			
				p2s=p2s+p2							
			endif		

		Enddo
		vMSPl=1.0d0/real(ncomp-nvaz-1)*pi2s
		vMSGl=1.0d0/(nTpeixe-ncomp+nvaz)*p2s				
		vSigp=1.0d0/cn*(vMSPl-vMSGl)				
		Fstl=vSigp/(vSigp+vMSGl)					
		tSig=tSig+vSigp
		tMSG=tMSG+vMSGl
	Enddo

Fstt=tSig/(tSig+tMSG)
        endsubroutine Fst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FstUD(Mp,n,nb, Fstup, Fstdown)
		IMPLICIT REAL*8 (A-H,O-Z)
		common /c1/ nTpeixe,ngenoma,rmu,nG
		common /c2/ ncomp
		Dimension::Mp(nTpeixe,ngenoma+1), Mpup(nTpeixe,ngenoma+1), Mpdown(nTpeixe,ngenoma+1), NK(ncomp)
		Nup=0
		Ndown=0
	
	Do nt=1, nTpeixe
		if (Mp(nt,1).le.nb) then
			Nup=Nup+1
			Mpup(Nup,:)=Mp(nt,:)
		else
			Ndown=Ndown+1
			Mpdown(Ndown,:)=Mp(nt,:)
		endif
	enddo
	
	NK=0 
	Do np=1,ntpeixe
		nposi=Mp(np,1)
		NK(nposi)=NK(nposi)+1
	Enddo

	ni2s=0 
	nvaz=0
	Do nr=1, nb
		if(NK(nr).gt.0)then 
		        ni=0
			Do np=1, Nup 
				if(Mpup(np,1).eq.nr) ni=ni+1	
			Enddo
			ni2s=ni2s+ni**2
		else
		nvaz=nvaz+1
		endif
	Enddo
	frac=real(ni2s)/real(Nup)
	cn=1.0d0/real(nb-nvaz-1)*(real(Nup)-frac)
		
	vMSPl=0.0d0				
	vMSGl=0.0d0				
	vSigp=0.0d0				
	Fstl=0.0d0				
	Fstt=0.0d0				
	tMSG=0.0d0				
	tSig=0.0d0
	
	Do l=2,ngenoma+1	
		pi2s=0.0d0
		npt=0
		Do np=1, Nup
			npt=npt+Mpup(np,l)	
		enddo
		pt=real(npt)/real(Nup)

		vMSPl=0.0d0			
		p2s=0.0d0 			
		p2=0.0d0 			
		Do nr=1, nb
			if(NK(nr).gt.0) then
				npi=0
				ncont=0
				Do np=1, Nup
					if(Mpup(np,1).eq.nr)then
						npi=npi+Mpup(np,l)
						ncont=ncont+1	
					endif
				Enddo
				pi=real(npi)/real(ncont)
				pi2s=pi2s+real(ncont)*(pi-pt)**2	

				p2=real(ncont)*pi*(1.0d0-pi)			
				p2s=p2s+p2
			endif		

		Enddo
		vMSPl=1.0d0/real(nb-nvaz-1)*pi2s
		vMSGl=1.0d0/(Nup-nb+nvaz)*p2s				
		vSigp=1.0d0/cn*(vMSPl-vMSGl)				
		Fstl=vSigp/(vSigp+vMSGl)					!
		tSig=tSig+vSigp
		tMSG=tMSG+vMSGl
	Enddo

Fstup=tSig/(tSig+tMSG)

	ni2s=0 
	nvaz=0
	Do nr=nb+1, ncomp
		if(NK(nr).gt.0)then 
		        ni=0
			Do np=1, Ndown 
				if(Mpdown(np,1).eq.nr) ni=ni+1	
			Enddo
			ni2s=ni2s+ni**2
		else
		nvaz=nvaz+1
		endif
	Enddo
	frac=real(ni2s)/real(Ndown)
	cn=1.0d0/real(ncomp-nb-nvaz-1)*(real(Ndown)-frac)
	
	vMSPl=0.0d0				
	vMSGl=0.0d0				
	vSigp=0.0d0				
	Fstl=0.0d0				
	Fstt=0.0d0				
	tMSG=0.0d0				
	tSig=0.0d0			
	Do l=2,ngenoma+1		
		pi2s=0.0d0
		npt=0
		Do np=1, Ndown
			npt=npt+Mpdown(np,l)	
		enddo
		pt=real(npt)/real(Ndown)

		vMSPl=0.0d0			
		p2s=0.0d0 			
		p2=0.0d0 			
		Do nr=nb+1, ncomp
			if(NK(nr).gt.0) then
				npi=0
				ncont=0
				Do np=1, Ndown
					if(Mpdown(np,1).eq.nr)then
						npi=npi+Mpdown(np,l)
						ncont=ncont+1	
					endif
				Enddo
				pi=real(npi)/real(ncont)
				pi2s=pi2s+real(ncont)*(pi-pt)**2	
				p2=real(ncont)*pi*(1.0d0-pi)			
				p2s=p2s+p2
			endif		

		Enddo

		vMSPl=1.0d0/real(ncomp-nb-nvaz-1)*pi2s
		vMSGl=1.0d0/(Ndown-ncomp+nb+nvaz)*p2s				
		vSigp=1.0d0/cn*(vMSPl-vMSGl)				
		Fstl=vSigp/(vSigp+vMSGl)					!
		tSig=tSig+vSigp
		tMSG=tMSG+vMSGl
		

	Enddo

Fstdown=tSig/(tSig+tMSG)

	endsubroutine FstUD
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	subroutine Pairwise(Mp,n)
	IMPLICIT REAL*8 (A-H,O-Z)
	common /c1/ nTpeixe,ngenoma,rmu,nG
	common /c2/ ncomp
	Dimension::Mp(nTpeixe,ngenoma+1)
	Integer, ALLOCATABLE :: NK(:)
	Real, ALLOCATABLE :: PFst(:)
	ALLOCATE (NK(ncomp))
	npair=int(ncomp**2)
	ALLOCATE (PFst(npair))
	npairpop=0
	NK=0 
	Do np=1, NTpeixe 
		nposi=Mp(np,1)
		NK(nposi)=NK(nposi)+1
	Enddo

	Do npop=1,ncomp
		Do npopp=1, ncomp
		ni2s=0 
			If (npop.eq.npopp) then
			fstt=0.0d0
			go to 14
			endif			
			If (NK(npop).gt.0.and.NK(npopp).gt.0) then
					nipop=0
					nipopp=0
					nipop=NK(npop)
					nipopp=NK(npopp)
			else
				Fstt=2.0d0
				go to 14
			endif
			ni2s=nipop**2+nipopp**2
			nTpops=NK(npop)+NK(npopp)
			frac=real(ni2s)/real(nTpops)
			cn=real(nTpops)-frac
			
			vMSPl=0.0d0				
			vMSGl=0.0d0				
			vSigp=0.0d0				
			Fstl=0.0d0				
			Fstt=0.0d0				
			tMSG=0.0d0				
			tSig=0.0d0			
			Do l=2,ngenoma+1		
				pi2s=0.0d0
				npt=0
				Do np=1, nTpeixe
					if(MP(np,1).eq.npop.or.MP(np,1).eq.npopp) then 
						npt=npt+Mp(np,l)
					endif
				enddo
				if(npop.eq.npopp) npt=2*npt
				pt=real(npt)/real(nTpops)

				vMSPl=0.0d0			
				p2s=0.0d0 			
				p2=0.0d0 			
				Do nr=1, ncomp
					if (nr.eq.npop.or.nr.eq.npopp) then
						npi=0
						ncont=0
						Do np=1, nTpeixe
							if(Mp(np,1).eq.nr)then
								npi=npi+Mp(np,l)
								ncont=ncont+1	
							endif
						Enddo
						pi=real(npi)/real(ncont)
						pi2s=pi2s+real(ncont)*(pi-pt)**2	
						p2=real(ncont)*pi*(1.0d0-pi)			
						p2s=p2s+p2
					endif
				Enddo
				if (npop.eq.npopp) p2s=2*p2s
				vMSPl=real(pi2s)
				vMSGl=real(p2s)/int(nTpops-2)				
				vSigp=1.0d0/cn*(vMSPl-vMSGl)				
				Fstl=vSigp/(vSigp+vMSGl)					
				tSig=tSig+vSigp
				tMSG=tMSG+vMSGl
				
			Enddo
			fstt=tSig/(tSig+tMSG)
14	continue	
			npairpop=int(npairpop+1)	
			PFst(npairpop)=fstt		
		Enddo
	Enddo

write(80,15) n, npair, PFst
  15	FORMAT( I5, 2X, I3, 2X, 400(2x, f6.4))
        endsubroutine Pairwise
