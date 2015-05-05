module Ising

export ising

function configuracion(L)
    sigma=zeros(L,L)
    for i=1:L
        for j=1:L
            sigma[i,j]=int(rand(-1:2:1))
        end
    end
    return sigma
end

function energia(sigma)
    L=size(sigma,1)
    E=0
    for i=1:L
        for j=1:L
            E+=-sigma[i,j]*(sigma[mod1(i+1,L),j]+sigma[i,mod1(j+1,L)])
        end
    end
    return E
end

function dE(sigma,i,j)
    L=size(sigma,1)
    E1=2*sigma[i,j]*(sigma[mod1(i+1,L),j]+sigma[mod1(i-1,L),j]+sigma[i,mod1(j+1,L)]+sigma[i,mod1(j-1,L)])
    return E1
end
function ising(L,T,t1,t2)
    sigma = int(configuracion(L))
    E=zeros(t2)
    M=zeros(t2)
    E[1]=energia(sigma)
    M[1]=sum(sigma)
    for i=2:t2
        E[i]=E[i-1]
        M[i]=M[i-1]
        a1=rand(1:L)
        a2=rand(1:L)
        E1=dE(sigma,a1,a2)
        b=min(e^(-E1/T),1)
        c=rand()
        if c<b
            sigma[a1,a2]=-sigma[a1,a2]
            E[i]=E[i]+E1
            M[i]=M[i]+2*sigma[a1,a2]
        end
    end
    Ep=sum(E[t1+1:t2])/(t2-t1)
    Mp=sum(M[t1+1:t2])/(t2-t1)
    return Ep,E,M,Mp
end

end
