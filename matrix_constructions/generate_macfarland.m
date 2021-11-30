% Function to generate the Macfarland matrix
function seq = generate_macfarland(root,N) 
    if mod(N,2)==1
    for n=0:(N-1)
        seq(n+1)=exp(-j*(pi*root*n*(n+1))/N);
    end
    elseif mod(N,2)==0
        for n=0:(N-1)
            seq(n+1)=exp(-j*(pi*root*n*(n))/N);
        end
    end
end