function r=corr(vect1, vect2, lag);

if nargin < 3;
  lag = 0;
end;

[m,n]=size(vect1);
if n == 1;
  vect1 = vect1';
end
[l,k]=size(vect2);
if l == 1
  vect2 = vect2';
end
if lag < 0;
  tem = vect1';
  vect1 = vect2';
  vect2 = tem;
  clear tem;
  lag = abs(lag);
end
%
if nargin == 2;
  v1 = (vect1 - mean(vect1));
  v2 = (vect2 - mean(vect2));
elseif nargin == 3;
  v1 = (vect1((lag+1):length(vect1)) - mean(vect1((lag+1):length(vect1))));
  v2 = (vect2(1:(length(vect2)-lag)) - mean(vect2(1:(length(vect2)-lag))));
end
%
r = (v1 * v2) / (sqrt(v1 * v1') * sqrt(v2' * v2));

