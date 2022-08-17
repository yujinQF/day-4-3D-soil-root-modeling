!> \file RootGrowthNeighb
!> \brief Contains functions orginating from RootGrowth
Module RootGrowthNeighb

Contains
!***************************************************************************
!> ###finds 8 corner points of cube (corner) and the voxel ID (imin)
!> that contains root tip or segment originating or ending point###
  Subroutine Neighb(x,y,z,corner,imin,lgrids)
    Use typedef
    Use GridData
    USE GridData2
    Use DomData
    Use MPIutils, Only: stop_program
    Implicit None

    Real(dp), Intent(in) :: x,y,z
    Real(dp):: delta1
    Integer(ap):: ixmin,ixymin
    Integer(ap), Intent(out) :: imin,corner(8)
    Logical :: lgrids
    Integer(ap):: nel_, nex_, ney_, nez_, nelm_
    Real(dp):: dxgrid_, dygrid_, dzgrid_, xmin_, ymin_, zmax_
    INTEGER(ap), ALLOCATABLE,DIMENSION (:,:) :: elmnod_
    !> \param x root x coordinate
    !> \param y root y coordinate
    !> \param z root z coordinate
    !> \param corner soil corner coordinates of voxel surrounding root coords
    !> \param nElm = total number of elements
    !>  ne* = number of elements (half-cuboids) in '*'-direction
    !> \param nel = nex*ney (number of elements per horizontal element-layer)
    !> \param imin = element-ID# of odd-numbered half of the cuboid containing the root
    !>        tip or segment

	
    IF (lgrids) THEN
        dxgrid_ = dxgrid2
        dygrid_ = dygrid2
        dzgrid_ = dzgrid2
        nel_ = nel2
        nex_ = nex2
        ney_ = ney2
        nez_ = nez2
        xmin_ = xmin2
        ymin_ = ymin2
        zmax_ = zmax2
        nelm_ = nelm2
        ALLOCATE(elmnod_(1:8,1:nElm_))
        elmnod_ = elmnod2

   Else
        dxgrid_ = dxgrid
        dygrid_ = dygrid
        dzgrid_ = dzgrid
        nel_ = nel
        nex_ = nex
        ney_ = ney
        nez_ = nez
        xmin_ = xmin
        ymin_ = ymin
        zmax_ = zmax
        nelm_ = nelm
        ALLOCATE(elmnod_(1:8,1:nElm_))
        elmnod_ = elmnod
    END IF

    ! find imin by checking distance between tip coord´s and cuboid center-points
    delta1 = x-(xmin_)
    If (delta1.Eq.nex_*dxgrid_) Then
       ixmin = nex_
    Elseif (delta1.Gt.nex_*dxgrid_) Then
       Print*, delta1, nex_*dxgrid_, nex_, xmin_,dxgrid_, x
       Call stop_program('Root node out of the soil domain (Xdirection)')
    Else
       ixmin = 1+Floor(delta1/dxgrid_)
    End If

    ! search along y-axis (start at ixmin):
    delta1 = y-(ymin_)
    If (delta1.Eq.ney_*dygrid_) Then
       ixymin = ixmin + (Floor(delta1/dygrid_)-1)*ney_
    Elseif (delta1.Gt.ney_*dygrid_) Then
       Print*, delta1, ney_*dygrid_, ney_, ymin_,dygrid_, y
       Call stop_program('Root node out of the soil domain (Ydirection)')
    Else
       ixymin = ixmin+Floor(delta1/dygrid_)*nex_
    End If

    ! search along z-axis (start at ixymin):
    delta1 = Abs(z-zmax_)
    If (delta1.Eq.nez_*dzgrid_) Then
       imin = ixymin + (Floor(delta1/dzgrid_)-1)*nel_
    Elseif (delta1.Gt.nez_*dzgrid_) Then
       Print*, delta1, nez_*dzgrid_, nez_, dzgrid_, z
       Call stop_program('Root node out of the soil domain (Zdirection)')
    Else
       imin = ixymin+Floor(delta1/dzgrid_)*nel_
    End If
    If(imin.Gt.nelm_) Print*, x,y,z

    ! assign cuboid corner nodes:
    corner(1) = elmnod_(1,imin)
    corner(2) = elmnod_(2,imin)
    corner(3) = elmnod_(3,imin)
    corner(4) = elmnod_(4,imin)
    corner(5) = elmnod_(5,imin)
    corner(6) = elmnod_(6,imin)
    corner(7) = elmnod_(7,imin)
    corner(8) = elmnod_(8,imin)

    Return
  End Subroutine Neighb
  !--------------------------------------------------------------
  !>  neighbouring nodes of each node of an element
  Subroutine FlowNeighb(corner,Flow_corner)
    Use typedef
    Use GridData
    Implicit None

    Integer(ap) :: i,c
    Integer(ap), Intent(in) :: corner(8)
    Integer(ap), Intent(out) :: Flow_corner(1:8,1:6)
    !> \param corner nodes of the centre element
    !> \param Flow_corner neighbouring nodes of the centre element

    Flow_corner=0
    Do c=1,8
       Flow_corner(c,1) = corner(c) + 1
       Flow_corner(c,2) = corner(c) - 1
       Flow_corner(c,3) = corner(c) + 1 + ney
       Flow_corner(c,4) = corner(c) - 1 - ney
       Flow_corner(c,5) = corner(c) + (nPt/(nez+1))
       Flow_corner(c,6) = corner(c) - (nPt/(nez+1))
       Do i=1,6
          If (Flow_corner(c,i) .Lt. 0) Flow_corner(c,i) = 0
       Enddo
    End Do

    Return
  End Subroutine FlowNeighb
  !--------------------------------------------------------------
  !>  neighbouring elements of an element iElmx
  Subroutine ElemNeighb(iElm,NeighElm)
    Use typedef
    Use GridData
    Implicit None

    Integer(ap) :: ix,iy,iz,i
    Integer(ap), Intent(in) :: iElm
    Integer(ap), Intent(out) :: NeighElm(1:6)
    !> \param iElm centre element
    !> \param NeighElm neighbouring elements (6)

    NeighElm = 0
    !> z-direction
    iz = Int(Ceiling(Real(iElm)/Real(nex)/Real(ney)))
    If (iz.Le.1) Then
       NeighElm(5) = 0
    Else
       NeighElm(5) = iElm + (nELm/nez)
    End If
    If (iz.Ge.nez) Then
       NeighElm(6) = 0
    Else
       NeighElm(6) = iElm - (nElm/nez)
    End If

    !> y-direction
    iy = Int(Ceiling(Real(iElm-(iz-1)*nex*ney) /Real(nex)))
    If (iy.Le.1) Then
       NeighElm(3) = 0
    Else
       NeighElm(3) = iElm - nex
    End If
    If (iy.Ge.ney) Then
       NeighElm(4) = 0
    Else
       NeighElm(4) = iElm + nex
    End If

    !> x-direction
    ix = iElm - (iz-1)*nex*ney - (iy-1)*nex
    If (ix.Le.1) Then
       NeighElm(1) = 0
    Else
       NeighElm(1) = iElm - 1
    End If
    If (ix.Ge.nex) Then
       NeighElm(2) = 0
    Else
       NeighElm(2) = iElm + 1
    End If
    Do i = 1,6
       If (NeighElm(i).Lt.0 .Or. NeighElm(i).Gt.nElm) NeighElm(i)=0
    Enddo

    Return
  End Subroutine ElemNeighb
  
End Module RootGrowthNeighb
