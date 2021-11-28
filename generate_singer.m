% This function returns a singer matrix of order pxp^2 that satisfies the RTP
% with specified values for k, t, delta.

function [ A , p ] = generate_singer ( n , k , t , delta )

if nargin < 4
delta = 0.5 ;
end

if nargin < 3
t = 1.5 ;
end

%% Checking for consistency and generating error messages.

if delta > sqrt((t-1)/t)
disp('Error: delta must be less than the square root of (t-1)/t')
end

tk = ceil(t*k) ;

p_bound = max( ( (tk-1)/delta )^2 , sqrt(n) ) ;

% Find all primes less than or equal to p_bound + 100.
% The presumption is that, at the kinds of small numbers encountered in compressed sensing,
% there will be at least one prime number between p_bound and p_bound + 100. 

p_list=primes( p_bound + 100 ) ;

p_index=find((p_bound<p_list));

p=p_list(p_index(1));

% Then we compute the chirp matrix.

s = ceil(n/p) ;
I = sqrt(-1) ;

for y=1:s
for x=1:p
r = (y-1)*p + x ;
xb = x-1 ; yb = y-1 ;
for t=1:p
tb = t-1 ;
C(t,r) = (-1)^(tb*yb) * exp ( ( I*pi*(2*xb + yb*tb)*tb ) / p);
end
end
end

A = (1/sqrt(p))*C(:,1:n) ;

end