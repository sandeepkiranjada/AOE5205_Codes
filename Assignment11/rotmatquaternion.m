function R = rotmatquaternion(q)
%
%  Copyright (c) 2019 Mark L. Psiaki.  All rights reserved.  
%
%  This function computes the 3-by-3 direction
%  cosines rotation matrix for the transformation from
%  one right-handed coordinate system to another
%  right-handed coordinate system as parameterized by
%  the rotation's unit-normalized quaternion:
%
%    R = R(q)
%  
%
%  Inputs:
%
%    q                     The 4-by-1 unit-normalized quaternion
%                          parameterization of the rotation.
%                          This input is non-dimensional.
%
%                          Note:  Be careful to enter a q
%                          such that q'*q = sum(q.^2) = 1.
%                          Otherwise, the columns and rows
%                          of R will not be unit-normalized.
%
%  Outputs:
%
%    R                     The 3-by-3 direction cosines matrix
%                          for the rotation that is parameterized
%                          by q.  It is non-dimensional.
% 

%
%  Retrieve the vector and scalar portions of the quaternion.
%  
   qvec = q(1:3,1);
   q4 = q(4,1);
%
%  Compute the cross-product-equivalent matrix for the vector
%  part of the quaternion.
%
   CPEQmat_qvec = [          0,   -qvec(3),   qvec(2);...
                       qvec(3),          0,  -qvec(1);...
                      -qvec(2),    qvec(1),         0];
%
%  Compute the rotation matrix using the formula that
%  employs the vector part and the scalar part of the 
%  quaternion and that includes a term involving the outer 
%  product of the vector part and another term involving
%  the cross-product equivalent matrix for the vector part.
%
%  If one feels unsure of how to implement this formula,
%  it is acceptable to implement the equivalent
%  formula that gives element-by-element formulas for the
%  9 elements of R as functions of the 4 elements of q.
%  
   R = (q4^2 - qvec'*qvec).*eye(3) + 2.*(qvec*qvec') - (2*q4).*CPEQmat_qvec;