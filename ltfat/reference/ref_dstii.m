function c=ref_dstii(f)
%REF_DSTII  Reference Discrete Sine Transform type II
%   Usage:  c=ref_dstii(f);
%
%   The transform is computed by an FFT of 4 times the length of f.
%
%
%   Url: http://ltfat.github.io/doc/reference/ref_dstii.html

% Copyright (C) 2005-2023 Peter L. Soendergaard <peter@sonderport.dk> and others.
% This file is part of LTFAT version 2.6.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

L=size(f,1);
W=size(f,2);

% Create weights.
w=ones(L,1);
w(L)=1/sqrt(2);
w=w*sqrt(2/L);

% Create transform matrix.
F=zeros(L);

for m=0:L-1
  for n=0:L-1
    F(m+1,n+1)=w(n+1)*sin(pi*(n+1)*(m+.5)/L);
  end;
end;

% Compute coefficients.
c=F'*f;



