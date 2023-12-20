!This code simulates an SIR model with a quarantine strategy on ER networks with cliques
!At time t=0, there is only one infected individual and the rest of the population is susceptible.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Input
!N_nodeI: Total number of individuals
!kminI,kmaxI:   Minimum and maximum number of cliques per individual.
!kminC,kmaxC:   Minimum and maximum clique size
!lambdaC,lambdaI: parameter of the Poisson distribution

!nrea  : number of realizations
!tb : Isolation time
!tr : Recovery time
!beta : probability of infection
!ff : probability of detection
!nn : maximum number of times a susceptible individual will obey a quarantine order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output
!File: Generates a table where the first column represents the detection probability f, and the second column represents the corresponding fraction of recovered people R at the final stage. Each row corresponds to an individual realization


module globals
  implicit none

  integer kminC,kmaxC,kminI,kmaxI
  integer N_node,N_nodeC,N_nodeI
  integer Rtot
  integer tr,tb,nn

  real(8) lambdaC,lambdaI
  real(8) kmedC,kmedI  
  real(8) r  !random number
  real(8) ff,beta,I0

  integer,allocatable,dimension(:)::edge,edgep
  integer,allocatable,dimension(:)::node,nodep
  integer,allocatable,dimension(:)::kk,kkp
  integer,allocatable,dimension(:)::State
  
  real(8),allocatable,dimension(:)::PkC,PkI
    
end module globals

module random
   integer::q1=1,q2=104
   integer ir(670)
   integer::sem(670)
   real(8)::nmax=2*(2**30-1)+1
   real(8) sig
end module random


Program MainProgram
  use globals
  use random

  implicit none

  integer i,j,rea,nrea
  integer nptosf
  real(8) product
  
  character(8) Var01
  character(5) Var02
  character(3) Var03
  character(3) Var04
  character(3) Var05
  character(5) Var06

  call initialize_random
  
  !Parameters

  N_nodeI=1000000

  nrea   =20

  tr     =1
  tb     =1
  nn     =1
  beta   =0.2d0
 
  

  kminC  =0
  kmaxC  =20
  kminI  =0
  kmaxI  =20  
  lambdaC=3d0
  lambdaI=3d0
  
  
  
  nptosf  =100
  
  write(Var01,'(I8)') N_nodeI
  call zeros(Var01)
  write(Var02,'(F5.3)') ff
  call zeros(Var02)
  write(Var03,'(I3)') tr
  call zeros(Var03)
  write(Var04,'(I3)') tb
  call zeros(Var04)
  write(Var05,'(I3)') nn
  call zeros(Var05)  
  write(Var06,'(F5.3)') beta
  call zeros(Var06)

  allocate(PkC(kminC:kmaxC),PkI(kminI:kmaxI))
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Degree distribution Cliques
  PkC           =0d0
  PkC(0)        =exp(-lambdaC)

  product      =1d0
  do i=1,kminC-1
     product=product*i
  enddo  
  do i=kminC,kmaxC
     if(i==0)cycle
     product=product*i         
     PkC(i)=1.d0*exp(-lambdaC)*lambdaC**(i)/dble(product)
  enddo  

  PkC  =PkC/sum(PkC)
  
  kmedC=0d0
  do i=kminC,kmaxC
     kmedC=kmedC+i*PkC(i)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !Degree distribution Individuals
  PkI           =0d0
  PkI(0)        =exp(-lambdaI)

  product      =1d0
  do i=1,kminI-1
     product=product*i
  enddo  
  do i=kminI,kmaxI
     if(i==0)cycle
     product=product*i         
     PkI(i)=1.d0*exp(-lambdaI)*lambdaI**(i)/dble(product)
  enddo  

  PkI  =PkI/sum(PkI)
     
  kmedI=0d0
  do i=kminI,kmaxI
     kmedI=kmedI+i*PkI(i)
  enddo
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  N_nodeC=kmedI*(N_nodeI/kmedC)
  N_node=N_nodeC+N_nodeI
    
  print*,"Number of cliques", N_nodeC,"Number of individuals",N_nodeI, "N_nodeC+N_nodeI",N_node


  allocate(kk(N_node),node(N_node),State(N_nodeI),kkp(N_nodeI),nodep(N_nodeI))  
  

  open(2,file='Scatter_ff_'//Var02//'_bet_'//Var06//'_NI_'//Var01//'_tr_'//Var03//'_tb_'//Var04//'_nn_'//Var05//'.dat')
  write(2,*)'#',rea
  write(2,*)'#','kminC,kmaxC,kminI,kmaxI,lambdaC,lambdaI'
  write(2,*)'#',kminC,kmaxC,kminI,kmaxI,lambdaC,lambdaI     
  do rea=1,nrea
     if(mod(rea,10)==0)print*,rea
     call ConfModel     
     do j=1,nptosf
        ff  =0.01d0*j
        State =1
        call SIR
        write(2,*) ff,Rtot/dble(N_nodeI)
     enddo
    
     deallocate(edge,edgep)          
  enddo
  close(2)
  

end Program MainProgram

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   SIR
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SIR
  use globals
  use random
  implicit none
  integer i,j
  integer sel,nod,vec
  integer Ninf,NisoI,NisoS
  integer Ninfnew,NisoInew,NisoIindirecnew,NisoSnew
  integer tt,Ntot
  integer, dimension(N_nodeI):: ListIndInf,ListIndIsolInf,ListIndIsolSus,ListIndIsolInfIndirectlyNew
  integer, dimension(N_nodeI):: ListIndInfNew,ListIndIsolInfNew,ListIndIsolSusNew
  integer, dimension(N_nodeI):: Clock,ClockIso,NumberIso

  !NumberIso(i) is the number of times that an individual "i" has already been placed under quarantine.

  
  !State=0  :susceptible
  !State=1  :infected
  !State=2  :recovered
  !State=3  :isolatedI
  !State=4  :isolatedS
  
  
  State         =0
  Ninf          =0   !number of infected individuals
  NisoI         =0   !number of isolated individuals who are infected
  NisoS         =0   !number of isolated individuals who are susceptible

  ListIndInf    =0   !List of infected individuals
  ListIndIsolInf=0   !List of isolated individuals
  ListIndIsolSus=0   !List of isolated individuals who are susceptible


  Rtot          =0   !number of recovered individuals

  Clock         =0
  ClockIso      =0
  NumberIso     =0

  tt            =0

  !INDEX-CASES
  call rand
  sel=N_nodeI*r+1
  State(sel)=1

  Clock(sel)=0    
  Ninf =Ninf+1
  ListIndInf(Ninf) =sel
  
  
  do while(Ninf>0.or.NisoI>0)

     do i=1,NisoS
        nod  =ListIndIsolSus(i)
        ClockIso(nod)=ClockIso(nod)+1
     enddo
     do i=1,Ninf
        nod  =ListIndInf(i)
        Clock(nod)=Clock(nod)+1
     enddo    

     !STEP 1: INFECTION
     !Those infected individuals who have not been isolated, will now transmit the disease to their susceptible neighbors.
     
     Ninfnew       =0
     ListIndInfNew =0
     
     do i=1,Ninf
        nod =ListIndInf(i)
        if(State(nod)==3)cycle  !infected people who have been placed under quarantine in the previous step cannot transmit the disease further
        do j=1,kkp(nod)
           vec =edgep(nodep(nod)+j-1)
           if(State(vec)==0)then
              call rand
              if(r<beta)then
                 State(vec)  =1
                 Ninfnew     =Ninfnew+1
                 ListIndInfNew(Ninfnew)=vec
              endif
           endif
        enddo     
     enddo

    
     !STEP 2: DETECTION AND QUARANTINE
     !The newly infected individuals are detected and isolated with probability ff
     NisoInew         =0
     ListIndIsolInfNew=0
     do i=1,Ninf
        nod  =ListIndInf(i)
        if(NumberIso(nod)>=nn)cycle !people with fatigue refuse further testing
        call rand
        if(r<ff)then
           NisoInew =NisoInew+1
           ListIndIsolInfNew(NisoInew)=nod
           State(nod)=3
        endif
     enddo
     
     !STEP 3: CONTACT TRACING AND QUARANTINE
     NisoSnew           =0
     NisoIindirecnew    =0
     ListIndIsolSusNew  =0
     ListIndIsolInfIndirectlyNew=0
     do i=1,NisoInew
        nod =ListIndIsolInfNew(i)
        do j=1,kkp(nod)
           vec =edgep(nodep(nod)+j-1)
           if(NumberIso(vec)>=nn)cycle !neighbors with fatigue do not follow quarantine orders
           if(State(vec)==0)then
              !a susceptible individual is placed under quarantine
              NisoSnew      =NisoSnew+1
              State(vec)    =4
              ClockIso(vec) =0
              NumberIso(vec)=NumberIso(vec)+1
              ListIndIsolSusNew(NisoSnew) =vec
           else if(State(vec)==1)then
              !an infected individual is placed under quarantine.
              State(vec) =3
              NisoIindirecnew  =NisoIindirecnew+1
              ListIndIsolInfIndirectlyNew(NisoIindirecnew)=vec
           endif
        enddo
     enddo

     ListIndIsolInfNew(NisoInew+1:NisoInew+NisoIindirecnew)=ListIndIsolInfIndirectlyNew(1:NisoIindirecnew)
     NisoInew=NisoInew+NisoIindirecnew
     
     ListIndIsolSus(1+NisoS:NisoS+NisoSnew)=ListIndIsolSusNew(1:NisoSnew)
     NisoS=NisoS+NisoSnew    
     
     !STEP 4: RECOVERY
     do i=1,Ninf
        nod =ListIndInf(i)
        if(State(nod)==3)cycle
        if(Clock(nod)<tr)then                     
           Ninfnew   =Ninfnew+1
           ListIndInfNew(Ninfnew)=nod
        else
           State(nod)=2
           Rtot      =Rtot+1
        endif
     enddo
     Ninf        =Ninfnew
     ListIndInf  =ListIndInfNew     
     ! Update the list of infected nodes, removing those that have been isolated.
     Ninfnew=0
     ListIndInfNew=0
     do i=1,Ninf
        nod =ListIndInf(i)
        if(State(nod)==1)then
           Ninfnew=Ninfnew+1
           ListIndInfNew(Ninfnew)=nod  
        endif
     enddo
     Ninf        =Ninfnew
     ListIndInf  =ListIndInfNew
     
     ! Now, we update the state of isolated infected individuals.
     do i=1,NisoI
        nod =ListIndIsolInf(i)
        if(Clock(nod)<tr)then
           Clock(nod)=Clock(nod)+1 
           NisoInew  =NisoInew+1
           ListIndIsolInfNew(NisoInew) =nod
        else
           State(nod)=2
           Rtot      =Rtot+1           
        endif
     enddo
     NisoI          =NisoInew     
     ListIndIsolInf =ListIndIsolInfNew


     !STEP 5: REINTEGRATION
    
     NisoSnew=0
     ListIndIsolSusNew=0
     
     do i=1,NisoS
        nod  =ListIndIsolSus(i)
        if(ClockIso(nod)==tb)then
           ClockIso(nod)  =0
           State(nod)     =0
        else
           NisoSnew       =NisoSnew+1
           ListIndIsolSusNew(NisoSnew) =nod
        endif
     enddo    
     
     NisoS           =NisoSnew
     ListIndIsolSus  =ListIndIsolSusNew     

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !New time step
     tt    =tt+1
  enddo

end subroutine SIR




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Configuration Model
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine creates three arrays: kk, edge, and node, which contain all the network topology information.
!kk(i) is the degree of node "i".
!The subset edge(node(i):node(i)+kk(i)-1) contains the list of neighbors of node "i".

subroutine ConfModel
  use globals
  use random
  implicit none
  integer i,j,k,ll
  integer nod1,lug,nod2,vec
  integer Aux,flag
  integer m,n,stubmaxpos,stubminpos
  integer counting, nstubs,nstubsC,nstubsI
  integer nstubAux,stubposC,stubposI

  real(8) ws

  integer, allocatable, dimension(:)::listStubC,listStubI,kkAux

  real(8),dimension(kminC-1:kmaxC)::CumulativePkC
  real(8),dimension(kminI-1:kmaxI)::CumulativePkI
  
  
  allocate(kkAux(n_node))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Cumulative distribution group 1 (cliques)
  
  CumulativePkC   =0d0
  ws             =0d0
  do i=kminC,kmaxC
     ws                = ws+PkC(i)
     CumulativePkC(i)   = ws
  enddo
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Cumulative distribution group 2 (Individuals)
  
  CumulativePkI   =0d0
  ws             =0d0
  do i=kminI,kmaxI
     ws                = ws+PkI(i)
     CumulativePkI(i)   = ws
  enddo  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!Assigning Connectivity to all nodes
2 kk=0 
  do i=1,N_nodeC
     call rand
     do j=kminC,kmaxC
        if(CumulativePkC(j-1)<r.and.r<=CumulativePkC(j)) then
           kk(i)    = j
           exit
        endif
     enddo
  enddo
  do i=N_nodeC+1,N_node
     call rand
     do j=kminI,kmaxI
        if(CumulativePkI(j-1)<r.and.r<=CumulativePkI(j)) then
           kk(i)    = j
           exit
        endif
     enddo
  enddo  

  nstubsC =sum(kk(1:N_nodeC))
  nstubsI =sum(kk(N_nodeC+1:N_node))
  nstubs  =nstubsC+nstubsI
  
  
  if(abs(nstubsC-nstubsI)>0.01d0*(kmedI*N_nodeI)) goto 2
  
  allocate(listStubC(nstubsC),listStubI(nstubsI))  
  allocate(edge(nstubs))

  listStubC     =0
  listStubI     =0
  nstubAux      =0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!List of stubs

  do i=1,N_nodeC
     if(kk(i)==0) cycle
     do k=1,kk(i)
        nstubAux           =nstubAux+1
        listStubC(nstubAux)=i
     enddo
  enddo
  
  nstubAux      =0
  do i=N_nodeC+1,N_node
     if(kk(i)==0) cycle
     do k=1,kk(i)
        nstubAux           =nstubAux+1
        listStubI(nstubAux)=i
     enddo
  enddo  

  node(1)       =1
  do i=1,n_node-1
     node(i+1)  = node(i)+kk(i)
  enddo

  edge        =0
  kkAux       =0
  counting    =0
  do while(nstubsC>0.and.nstubsI>0)
3    if(counting==100) then
	deallocate(listStubC,listStubI,edge)
	go to 2
     end if
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Randomly choose two stubs
     
     call rand
     stubposC    =r*nstubsC+1
     call rand
     stubposI    =r*nstubsI+1

     m        =listStubC(stubposC)   !the first stub corresponds to node "m"
     n        =listStubI(stubposI)   !the second stub corresponds to node "n"

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Checking whether these two stubs can be joined or not

     if(m.ne.n) then
	do k=1,kkAux(m)
           if(edge(node(m)+k-1)==n)then
              counting     =counting+1
              go to 3   !if nodes "m" and "n" are already connected, we have to choose another pair of stubs
           endif
	enddo
	edge(node(m)+kkAux(m))   =n
	edge(node(n)+kkAux(n))   =m

	kkAux(m)                 =kkAux(m)+1
	kkAux(n)                 =kkAux(n)+1
	listStubC(stubposC)      =listStubC(nstubsC)
	listStubI(stubposI)      =listStubI(nstubsI)
	
	nstubsC             =nstubsC-1
	nstubsI             =nstubsI-1
	counting            =0
     else
        counting           =counting+1
        go to 3 !self-loops are forbidden
     endif
  enddo
  
  kk   =kkAux
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!One-mode Projection
  
  do i=N_nodeC+1,N_node
     nod1 =i-N_nodeC
     Aux  =0
     do j=1,kk(i)
        lug      =edge(node(i)+j-1)
        Aux      =Aux+kk(lug)
     enddo
     kkp(nod1)  =Aux-kk(i)
  enddo    
  
  nodep    =0
  nodep(1) =1

  do i=2,N_nodeI
     nodep(i)=nodep(i-1)+kkp(i-1)
  enddo
  
    
  allocate(edgep(sum(kkp)))
  
  kkp    =0
  edgep  =0

  do i=N_nodeC+1,N_node
     nod1 =i-N_nodeC
     
     do j=1,kk(i)
        lug      =edge(node(i)+j-1)
        do k=1,kk(lug)
           nod2  =edge(node(lug)+k-1)-N_nodeC
           if(nod2==nod1)cycle
           flag  =0
           do ll=1,kkp(nod1)
              vec   =edgep(nodep(nod1)+ll-1)
              if(vec==nod2)then
                 flag =1
                 exit
              endif
           enddo
           if(flag==0)then
              edgep(nodep(nod1)+kkp(nod1))   =nod2
              edgep(nodep(nod2)+kkp(nod2))   =nod1
              kkp(nod1)  =kkp(nod1)+1
              kkp(nod2)  =kkp(nod2)+1
           endif
        enddo
     enddo  
  enddo    
  
  
  deallocate(kkAux,listStubC,listStubI)
end subroutine ConfModel



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Seed for the Pseudorandom number generator
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_random
  use globals
  use random
  implicit none

  integer i,see(33)
  integer hour

  CALL SYSTEM_CLOCK(hour)
  !hour=26896048
  see=hour
  CALL RANDOM_SEED(put=see)		
  CALL RANDOM_SEED(get=see)

  do i=1,670
     call random_number(r)
     sem(i)   =r*nmax
     ir(i)    =i+1
  enddo
  ir(670)=1
  return
end subroutine initialize_random
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Pseudo-random number generator
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rand
  use globals
  use random
  implicit none

 1 q1=ir(q1)
  q2=ir(q2)
  sem(q1)= IEOR(sem(q1),sem(q2))
  r=dfloat(sem(q1))/nmax
  if(r==1d0) go to 1
  return
end subroutine rand


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   Zero Character
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zeros(b)
 implicit none
 
 character(*):: b
 integer i
 
 do i=1,len(b)
   if (b(i:i)==' ') then  
      b(i:i)='0'
   endif
 end do
end subroutine zeros


