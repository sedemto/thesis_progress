%TIME_8
%
%   The purpose of this test is to determine whether it pays of to do
%   special optimizations for p=1 q=2 and p=2 q=3%
%
%   Url: http://ltfat.github.io/doc/timing/time_8.html

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


Lr=[480000*sf^2,262144*sf^2,900*sf^2];
ar=[     600*sf,        512,       2];
Mr=[     800*sf,       1024,  600*sf];

for ii=1:length(Lr)

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);
  
  g=rand(L,1);
  f=rand(L,1);
  
  [L, a, M]
  
  gf=comp_wfac(g,a,M);
  
  c1=mex_dgt_fac_1(f,gf,a,M,1);

  % This is an attempt of cacheflushing
  f=rand(L,1); gf=rand(size(gf));

  c2=mex_dgt_fac_7(f,gf,a,M,1);

end;



