!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  USER FUNCTIONS FOR THE BERING GLACIER PROBLEM
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!
! User Function MinZsBottom
! Return the bedrock altitude for a given nodes
!!------------------------------------------------------------------------------!!
FUNCTION MinZsBottom ( Model, nodenumber, znode) RESULT(Zbed)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber  
   REAL(KIND=dp) :: znode, Zbed 
   INTEGER :: imin, Npt, t
   INTEGER :: NMAX, i, j,Nb, Nbx, Nby, ib, ix, iy
   REAL(KIND=dp) :: x, y, z, xb0, yb0, x1, x2, y1, y2, zi(2,2) 
   REAL(KIND=dp) :: R, Rmin, dbx, dby, lbx, lby
   REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE xb, yb, zb

   Nbx = 861
   Nby = 321
   xb0 = 336000.0d0
   yb0 = 6660000.0d0      
   lbx = 172000.0
   lby = 64000.0

   Nb = Nbx * Nby
   dbx = lbx / (Nbx-1.0)
   dby = lby / (Nby-1.0)

   IF (FirstTime) THEN
        FirstTime = .False.
        
! Open the DEM
! Grid specificity can be read in mesh_input.dat
! 
        OPEN(10,file="bed_JPL_v7_fourptsall.dat")
        
        ALLOCATE(xb(Nb), yb(Nb), zb(Nb))
        READ(10,*)(xb(i), yb(i), zb(i), i=1,Nb)
        CLOSE(10)
   END IF
        
! Compute zbed for that point (x,y)

        x = Model % Nodes % x (nodenumber)
        y = Model % Nodes % y (nodenumber)
         
! Find zbed for that point from the Bedrock MNT 

        ix = INT((x-xb0)/dbx)+1
        iy = INT((y-yb0)/dbx)+1
        ib = Nbx * (iy - 1) + ix
        
        x1 = xb(ib)
        x2 = xb(ib+1)
        y1 = yb(ib)
        y2 = yb(ib + Nbx)
        
        zi(1,1) = zb(ib)
        zi(2,1) = zb(ib+1)
        zi(2,2) = zb(ib + Nbx + 1)
        zi(1,2) = zb(ib + Nbx)
        
        
        IF ((zi(1,1)<-9990.0).OR.(zi(1,2)<-9990.0).OR.(zi(2,1)<-9990.0).OR.(zi(2,2)<-9990.0)) THEN
           IF ((zi(1,1)<-9990.0).AND.(zi(1,2)<-9990.0).AND.(zi(2,1)<-9990.0).AND.(zi(2,2)<-9990.0)) THEN
           ! Find the nearest point avalable
             Rmin = 9999.0
             DO i=1, Nb
               IF (zb(i)>0.0) THEN
                 R = SQRT((x-xb(i))**2.0+(y-yb(i))**2.0)
                 IF (R<Rmin) THEN
                   Rmin = R
                   imin = i
                 END IF
               END IF
             END DO
            zbed = zb(imin)
                        
           ELSE
            ! Mean value over the avalable data
             zbed = 0.0
             Npt = 0
             DO i=1, 2
               DO J=1, 2
                  IF (zi(i,j) > 0.0) THEN 
                     zbed = zbed + zi(i,j)
                     Npt = Npt + 1
                  END IF   
               END DO
             END DO
             zbed = zbed / Npt
             
           END IF
        ELSE
          zbed = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dbx*dby)      
        END IF

END FUNCTION MinZsBottom

! Function InterpolateDEM
!!------------------------------------------------------------------------------!!
FUNCTION InterpolateDEM (x, y, xb, yb, zb, Nbx, Nby, xb0, yb0, lbx, lby) RESULT(zbed)
   USE DefUtils
   IMPLICIT NONE
   INTEGER :: imin, Npt, t
   INTEGER :: NMAX, i, j, Nb, Nbx, Nby, ib, ix, iy
   REAL(KIND=dp) :: x, y, zbed, xb0, yb0, x1, x2, y1, y2, zi(2,2) 
   REAL(KIND=dp) :: R, Rmin, lbx, lby, dbx, dby
   REAL(KIND=dp) :: xb(Nbx*Nby), yb(Nbx*Nby), zb(Nbx*Nby)       

! Find zbed for that point from the Bedrock MNT 
        dbx = lbx / (Nbx-1.0)
        dby = lby / (Nby-1.0)
        Nb = Nbx*Nby

        ix = INT((x-xb0)/dbx)+1
        iy = INT((y-yb0)/dbx)+1
        ib = Nbx * (iy - 1) + ix
        
        x1 = xb(ib)
        x2 = xb(ib+1)
        y1 = yb(ib)
        y2 = yb(ib + Nbx)
        
        zi(1,1) = zb(ib)
        zi(2,1) = zb(ib+1)
        zi(2,2) = zb(ib + Nbx + 1)
        zi(1,2) = zb(ib + Nbx)
        
        
        IF ((zi(1,1)<-9990.0).OR.(zi(1,2)<-9990.0).OR.(zi(2,1)<-9990.0).OR.(zi(2,2)<-9990.0)) THEN
           IF ((zi(1,1)<-9990.0).AND.(zi(1,2)<-9990.0).AND.(zi(2,1)<-9990.0).AND.(zi(2,2)<-9990.0)) THEN
           ! Find the nearest point avalable
             Rmin = 9999.0
             DO i=1, Nb
               IF (zb(i)>0.0) THEN
                 R = SQRT((x-xb(i))**2.0+(y-yb(i))**2.0)
                 IF (R<Rmin) THEN
                   Rmin = R
                   imin = i
                 END IF
               END IF
             END DO
            zbed = zb(imin)
                        
           ELSE
            ! Mean value over the avalable data
             zbed = 0.0
             Npt = 0
             DO i=1, 2
               DO J=1, 2
                  IF (zi(i,j) > 0.0) THEN 
                     zbed = zbed + zi(i,j)
                     Npt = Npt + 1
                  END IF   
               END DO
             END DO
             zbed = zbed / Npt
             
           END IF
        ELSE
          zbed = (zi(1,1)*(x2-x)*(y2-y)+zi(2,1)*(x-x1)*(y2-y)+zi(1,2)*(x2-x)*(y-y1)+zi(2,2)*(x-x1)*(y-y1))/(dbx*dby)      
        END IF
END FUNCTION InterpolateDEM

!!!!!!!!!!!!!!!!!!!
! User Function TopSurface
! Return the TopfaceElevation from the DEM 
! Make sure the top surface is above the bedrock...
!!------------------------------------------------------------------------------!!
FUNCTION TopSurface ( Model, nodenumber, znode) RESULT(Zsurf)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber  
   REAL(KIND=dp) :: znode, Zbed, zsurf 
   INTEGER :: imin, Npt, t
   INTEGER :: NMAX, i, j,Nb, Nbx, Nby, ib, ix, iy
   INTEGER :: Ns, Nsx, Nsy
   REAL(KIND=dp) :: x, y, z, xb0, yb0, x1, x2, y1, y2, zi(2,2) 
   REAL(KIND=dp) :: xs0, ys0
   REAL(KIND=dp) :: R, Rmin, lbx, lby
   REAL(KIND=dp) :: lsx, lsy, InterpolateDEM
   REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:)       
   REAL(KIND=dp), ALLOCATABLE :: xs(:), ys(:), zs(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE xb, yb, zb, xs, ys, zs

   Nsx = 861
   Nsy = 321
   xs0 = 336000.0d0
   ys0 = 6660000.0d0      
   lsx = 172000.0
   lsy = 64000.0

   Ns = Nsx * Nsy

   Nbx = 861
   Nby = 321
   xb0 = 336000.0d0
   yb0 = 6660000.0d0      
   lbx = 172000.0
   lby = 64000.0

   Nb = Nbx * Nby

   IF (FirstTime) THEN
        FirstTime = .False.
        OPEN(10,file="s0_BBGS_v2_59959.dat")
        ALLOCATE(xs(Ns), ys(Ns), zs(Ns))
        READ(10,*)(xs(i), ys(i), zs(i), i=1,Ns)
        CLOSE(10)

        OPEN(10,file="bed_JPL_v7_fourptsall.dat")
        ALLOCATE(xb(Nb), yb(Nb), zb(Nb))
        READ(10,*)(xb(i), yb(i), zb(i), i=1,Nb)
        CLOSE(10)
   END IF

        
! Compute zbed for that point (x,y)

        x = Model % Nodes % x (nodenumber)
        y = Model % Nodes % y (nodenumber)
         
        zbed = InterpolateDEM(x,y,xb,yb,zb,Nbx,Nby,xb0,yb0,lbx,lby)
        zsurf = InterpolateDEM(x,y,xs,ys,zs,Nsx,Nsy,xs0,ys0,lsx,lsy)
        IF (zsurf-zbed < EPSILON(zsurf)) zsurf = zbed + 1.0
END FUNCTION TopSurface 

!!!!!!!!!!!!!!!!!!!
! User Function BottomSurface
! Return the BottomSurfaceElevation from the DEM 
!!------------------------------------------------------------------------------!!
FUNCTION BottomSurface ( Model, nodenumber, znode) RESULT(Zbed)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber  
   REAL(KIND=dp) :: znode, Zbed 
   INTEGER :: imin, Npt, t
   INTEGER :: NMAX, i, j, Nb, Nbx, Nby
   REAL(KIND=dp) :: x, y, z, xb0, yb0 
   REAL(KIND=dp) :: R, Rmin, lbx, lby, InterpolateDEM
   REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE xb, yb, zb

   Nbx = 861
   Nby = 321
   xb0 = 336000.0d0
   yb0 = 6660000.0d0      
   lbx = 172000.0
   lby = 64000.0

   Nb = Nbx * Nby

   IF (FirstTime) THEN
        FirstTime = .False.
        OPEN(10,file="bed_JPL_v7_fourptsall.dat")
        ALLOCATE(xb(Nb), yb(Nb), zb(Nb))
        READ(10,*)(xb(i), yb(i), zb(i), i=1,Nb)
        CLOSE(10)
   END IF
        
! Compute zbed for that point (x,y)

        x = Model % Nodes % x (nodenumber)
        y = Model % Nodes % y (nodenumber)
         
        zbed = InterpolateDEM(x,y,xb,yb,zb,Nbx,Nby,xb0,yb0,lbx,lby)
END FUNCTION BottomSurface 

!!!!!!!!!!!!!!!!!!!
! User Function SurgeWave
! Track propagation of surge wave and assign basal sliding coeff
!!------------------------------------------------------------------------------!!
FUNCTION SurgeWave( model, n, time ) RESULT(result)
	USE DefUtils
	IMPLICIT None
	TYPE(Model_t) :: model
	TYPE(Solver_t), Target :: Solver
	TYPE(Element_t) ::  Element
	TYPE(Element_t), POINTER :: BoundaryElement
	TYPE(Nodes_t) :: ElementNodes
	TYPE(ValueList_t), POINTER :: BC
	INTEGER :: n, i, insize, starti, vel, tdI, fronti, w, length, nboundary 
	INTEGER :: effectarea, disti, temp(1), f, BoundaryElementNode
	REAL(KIND=dp) :: time, result
	REAL(KIND=dp) :: x, y, ref
	! Define variable for start distance along track
	REAL(KIND=dp) :: dstart
	! Define variables for location of surge front at time td 
	REAL(KIND=dp) :: xt, yt, td
	! Variables for flowline data and basal fric param along flowline and distance measure
	REAL(KIND=dp), ALLOCATABLE :: xf(:), yf(:), df(:), dbeta(:), dist(:)
	LOGICAL :: FirstTime=.True.

	SAVE FirstTime
        SAVE xf, yf, df

	! Time in days (add 0.1 so it is converted to the correct integer)
	td = (365 * time) + 0.1
	tdI = INT(td)

	! -----------------------
	! get element information
	! -----------------------
	BoundaryElement => model % CurrentElement
	IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
	      CALL FATAL('SurgeWave','No boundary element found')
	END IF
	BC => GetBC()
	IF ( .NOT. ASSOCIATED(BC) ) THEN
	      CALL FATAL('SurgeWave','No boundary condition found')
	END IF

	! Get locations of nodes
	nboundary = GetElementNOFNodes(BoundaryElement)
	DO BoundaryElementNode=1,nboundary
		IF (  n == BoundaryElement % NodeIndexes(BoundaryElementNode) ) EXIT
	END DO

	! get element nodes
	! -----------------
	CALL GetElementNodes(ElementNodes,BoundaryElement)
	x = ElementNodes % x(BoundaryElementNode)
	y = ElementNodes % y(BoundaryElementNode)

	! Define size of flowline file
	insize = 15730

	! Define start location as index in flowline data (rift i=9000, 2011 i = 12300)
	starti = 12300

	!Define surge wave velocity (divide by the seperation distance, so 50m/day => 5)
	vel = 3.8

	! Deifine effected area width (divided by 10)
	w = 2000

	! Define surge front width f (divided by 10)
	f = 60

	! Load flowline data (along track flow seperated by 10m)
	IF (FirstTime) THEN
		FirstTime = .False.
		OPEN(10,file="flowline_BBGS.dat")
		ALLOCATE(xf(insize), yf(insize), df(insize))
		READ(10,*)(xf(i), yf(i), df(i), i=1,insize)
		CLOSE(10)
	END IF

	! Find start location along track
	dstart = df(starti)

	! Calculate legth of effected area and location of surge front (and index)
	length = vel * tdI
	fronti = starti + length ! Add to go downglacier

	! Actual length of effected area is min of width an dlengthi
	effectarea = MIN(length,w)
	
	! put default pointer back to where it belongs
        ! --------------------------------------------
        model % CurrentElement => BoundaryElement

	! Calculate basal sliding coeff for each along flowline point
	ALLOCATE(dbeta(insize))
	ALLOCATE(dist(insize))
	DO i=1, insize
        	IF (i > (fronti+f)) THEN
        		dbeta(i) = 0.0001
		ELSE IF (i <=(fronti +f) .AND. i >= fronti) THEN
			dbeta(i) = 0.000001 * (1 + 9*(i - fronti)/f)
        	ELSE IF (i < (fronti - effectarea)) THEN !Subtract to go upglacier
        		dbeta(i) = 0.0001
        	ELSE ! Effected area
        		dbeta(i) = 0.000001 * (1 + 9 * (fronti - i)/w)
        	END IF
		!dbeta(i) = 0.0001
		! find reference to flowline distance point (squared distance)
		dist(i) = (x-xf(i))**2 + (y-yf(i))**2
        END DO

	temp = MINLOC(dist)
	disti = temp(1)
	!IF (dbeta(disti) <0.00005) THEN
	!	write(*,*)'Nodenumber',n,'x',x,'y',y,'Beta',dbeta(disti)
	!END IF

	! give friction coeff
	result = dbeta(disti)


END FUNCTION SurgeWave
