function [QDPR]=qrtmul(QPR,Q)
%
%  Copyright (c) 2000 Mark L. Psiaki.  All rights reserved.  
%
%  "And whatever you do, whether in word or deed, do it all in the name of 
%   the Lord Jesus, giving thanks to God the Father through Him." 
%   (Colossians 3:17).
%
%  This program was originally developed by M.L.Psiaki in FORTRAN('qrtmul.for') 
%  (final version: 01/29/88) 
%
%  Modified to MATLAB version on 1999 by Hee Jung
%
%  Further optimized by Mark Psiaki in Feb. 2000. 
% 
%  This function does quarternion multiplication to give QDPR as a function of
%  QPR and Q where the coordinate transformations for these quarternions are 
%  related by:
%                     A(QDPR) = A(QPR)*A(Q)
%
%  The quarternions are assumed to be normalized upon entering this subroutine.
%  See. WERTZ pp. 415,416.
% 
%  Inputs
%   QPR :
%   Q :
% 
% Outputs 
%   QDPR: 
%

	QDPR = QPR*Q(4);
	QDPR(1:3) = QDPR(1:3) + QPR(4)*Q(1:3);
	QDPR(4) = QDPR(4) - QPR(1:3)'*Q(1:3);
  	QDPR(1) = QDPR(1) + QPR(3)*Q(2) - QPR(2)*Q(3);
	QDPR(2) = QDPR(2) + QPR(1)*Q(3) - QPR(3)*Q(1);
	QDPR(3) = QDPR(3) + QPR(2)*Q(1) - QPR(1)*Q(2);
