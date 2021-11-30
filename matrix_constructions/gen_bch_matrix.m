function Phi = gen_bch_matrix(n , k , p)

% Phi = gen_bch_matrix(n , k , p)
%
% This function, based on the p-ary BCH codes, constructs an m*n matrix 
% with a coherence value less than 1/(k-1) which is a sufficient condition 
% for establishing the RIP of order "k" and also perfect reconstrcution 
% of (k/2)-sparse vectors. The details of the design can be found in: 
%
% [*] A. Amini, V. Montazerhodjat and F. Marvasti, "Matrices with small 
%     Coherence using p-ary Block Codes," 
%     http://bigwww.epfl.ch/amini/Papers/GenRIP%20(two-column).pdf
%
% The elements of the matrix are on the unit circle (prior to column 
% normalization). The matrix is used to sample (k/2)-sparse vectors of the 
% size n*1. In fact, "k" and "n" are the standard notations of the CS
% parameters. The output is a m*n complex (except for p=2) matrix where "m"
% stands for the number of samples.
%
% "p" is the charachteristic of the field and should be a prime integer.
%
% The function is written by Arash Amini

Circ_En         = 0;    % Circ_En = 1 implies that the columns of the matrix 
                        % should be sorted in a way that the circular shifted 
                        % versions are adjacent. This procedure considerably 
                        % increases the computational time.

if nargin == 2
    p       = 2;
end

% determining \tilde{m} for the requested n
i_alg           = ceil(log2(k * p / (p-1) / 2) / log2(p));
tildem          = max(i_alg - 1 , 1);
hx              = [];

while p ^ (length(hx) - 2) < n
    tildem      = tildem + 1;
    [hx , gx , hxtext]      = PolyBCHp(p , tildem , i_alg);
end
m               = p ^ tildem  -  1;

% if p ^ (length(hx) - 2) > n
%     disp(['!!!! Note: The same number of samples can be used (the same "k") for n <= ' , num2str(p ^ (length(hx) - 2)) , ' !!!!'])
% end

[hx , gx , hxtext]      = PolyBCHp(p , tildem , i_alg);



% modifing the code for the even parity
gx              = mod(-[gx  0]  + [0  gx]  ,  p);


% generating the matrix
Phi             = zeros(m  ,  n);
if Circ_En == 1     % the columns of the matrix are put together such that the circular shifted versions are neighbors
    ind_flag    = [1  zeros(1  ,  p ^ (length(hx) - 2) - 1)];
    Col_Num     = 1;
    uncoded     = de2bi(0  ,  length(hx) - 2 , p , 'right-msb');
    
    aa          = bi2de(Phi([m , 1 : length(hx) - 3]  ,  Col_Num).' , p , 'right-msb');
    while aa < p ^ (length(hx) - 2)
        if ind_flag(aa + 1) == 0
            uncoded     = de2bi(aa  ,  length(hx) - 2 , p , 'right-msb');
            Col_Num     = Col_Num + 1;
            [quet , remnd]  = deconv([uncoded  zeros(1 , length(gx)-1)]  ,  gx(end : -1 : 1));
            Phi(:  ,  Col_Num)  = mod([uncoded  zeros(1 , length(gx)-1)] - remnd  ,  p).';
            ind_flag(aa + 1)    = 1;
            
            aa          = bi2de(Phi([m , 1 : length(hx) - 3]  ,  Col_Num).' , p , 'right-msb');
        else
            aa          = aa + 1;
        end
    end
    
    Phi         = Phi(: , 1 : n);
        
    
else
    for aa = 0 : n - 1
        %     aa
        uncoded     = de2bi(aa , length(hx) - 2 , p , 'right-msb');
        Phi(:  ,  aa + 1) = mod(conv(uncoded , gx) ,  p).';
    end
end


Phi         = exp(2 * j * pi / p * Phi) / sqrt(m);

end








function [hx , gx , hxtext] = PolyBCHp(p , tildem , i)

% [hx , gx , hxtext] = PolyBCHp(p , tildem , i)
%
% This function finds the binary parity check ("hx") and code generating 
% ("gx") polynomials of the BCH code used for CS purposes.
%
% "p" is the charachteristic of the field and should be a prime integer.
% The field order is p^tildem where "tildem" is the second input variable.
% ( GF(p^tildem) )
% 
% The roots of these binary polynomials belong to GF(p^tildem); moreover, 
% the roots of h(x) are restricted to lie in the set 
% {1, alpha, alpha^2, ..., alpha^( p^(tildem-1) + p^(tildem - i) )}
% where "alpha" is one of the primitive roots of the field.
%
% Remark:   The outputs "hx" and "gx" are p-ary row vectors which show the 
%           polynomial coefficients with the ascending order with respect 
%           to the powers of 'x'; i.e., the first elements are constant 
%           terms of the polynomial. In addition, "hxtext" is a string
%           which shows the parity check polynomial.


% setting the field and its primitive root
q               = p ^ tildem;
field           = gftuple((-1 : q - 2).' , tildem , p);
alpha           = 2;        % the primitive root



% finding the root powers of the parity check polynomial ( h(x) )
Hbin            = binHseq(tildem , i);
Hseq            = bi2de(Hbin , p , 'left-msb');


% evaluating the parity check polynomial by its roots
Hpoly           = 0;
for bb = 1 : min(length(Hseq) , q-1)
    root_inv    = gfsub(-1 , Hseq(bb) , field);
    Hpoly       = gfconv(Hpoly  , [root_inv   0] , field);
end
hx              = gftuple(Hpoly.' , tildem , p);
hx              = hx(: , 1).';
deg_h           = length(Hpoly) - 1;


% producing the string that represents the parity check polynomial
hxtext          = ['     h(x) = ' , num2str(hx(1))];
for bb = 2 : deg_h + 1
    if hx(bb) ~= 0
        hxtext  = [hxtext , ' + ' , num2str(hx(bb)) , 'x^' , num2str(bb - 1)];
    end
end

% finding g(x) by dividing x^(p^tildem - 1) - 1 by h(x)
BigPoly         = [gfsub(-1 , 0 , field)  -ones(1 , q-2)  0];
[Gpoly , chk]   = gfdeconv(BigPoly , Hpoly , field);
gx              = gftuple(Gpoly.' , tildem , p);
gx              = gx(: , 1).';

end








function Hbin = binHseq(tildem , i)

% Hbin = binHseq(tildem , i)
%
% This function finds all binary sequences of length "tildem" such that 1s
% are circularly spaced with at least $i$ zeros. These binary sequences are
% put together as rows of a matrix ("H")

Allseq          = de2bi([0 : 2^tildem - 1]  ,  'left-msb');
isOK            = zeros(1 , 2^tildem);
AllInd          = 1 : tildem;

for aa = 1 : 2^tildem
    seq         = Allseq(aa , :);
    loc1        = AllInd(seq == 1);
    
    
    if length(loc1) > 1         
        diffs   = loc1(2 : end) - loc1(1 : end-1);          % spacing between the 1s in the middle
        diffs   = [loc1(1) + tildem - loc1(end)    diffs];  % including the circular spacing between 
                                                            % the first and the last 1
        if min(diffs) > i       % whether there exists at least "i" zeros in between
            isOK(aa)    = 1;
        end
        
    else
        isOK(aa)= 1;
    end
end

Hbin            = Allseq(isOK == 1  ,  :);
end
