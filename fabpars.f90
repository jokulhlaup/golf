!f2py subroutine for fabric parameters. Supplemental to opil. 
Module fp

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Finds euler angles for 'n' sets of eigenvectors)          !
!!   Input: evs are eigenvectors for 'n' fabric samples, 3 3x1 !
!!   vectors arranged as a 9x1 array.                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine fangles(evs,angles,n)


  Double Precision, Intent(in), Dimension(n,9) :: evs
  Double Precision, Intent(out), Dimension(3)  :: angles

  Do i=1,n
    angles[i,1]=acos(-evs[i,6]/sqrt(1-eigvec[i,9]^2))
    angles[i,2]=acos(eigvec[i,9])
    angles[i,3]=acos(eigvec[i,6]/sqrt(1-eigvec[3,3]^2))
  End Do

End Subroutine fangles


