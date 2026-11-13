from itertools import combinations
from functools import cache
from bridge import table, folded
from scsmall import sc

def evaluate(T,m,odd):
    n=len(T); V=(1<<n)-1
    if n==1:
        return 1,1
    rho=lambda i:(i+m)%(2*m) if i<2*m else i
    def reflect(S):
        return sum(1<<rho(i) for i in range(n) if S>>i&1)
    def arc(i,j):
        return T[i]>>j&1
    @cache
    def E(S,i):
        R=S^(1<<i)
        return sum(E(R,j) for j in range(n) if R>>j&1 and arc(j,i)) if R else 1
    H=F=0
    for cut in combinations(range(n),m):
        S=sum(1<<i for i in cut)
        for x in (range(n) if odd else [-1]):
            if odd and S>>x&1:
                continue
            R=V^S^(1<<x if odd else 0)
            mirror=(reflect(R),rho(x) if odd else -1)
            if (S,x)>mirror:
                continue
            e=[E(S,i) for i in cut]
            if (S,x)==mirror:
                if odd:
                    f=sum(v*arc(i,x) for i,v in zip(cut,e))
                    H+=f*f
                    F+=f
                else:
                    H+=sum(u*v*arc(i,rho(j)) for i,u in zip(cut,e) for j,v in zip(cut,e))
                    F+=sum(u*arc(i,rho(i)) for i,u in zip(cut,e))
            else:
                J=[j for j in range(n) if R>>j&1]
                if odd:
                    value=sum(v*arc(i,x) for i,v in zip(cut,e))*sum(arc(x,j)*E(mirror[0],rho(j)) for j in J)
                else:
                    value=sum(v*arc(i,j)*E(mirror[0],rho(j)) for i,v in zip(cut,e) for j in J)
                H+=2*value
    return H,F

checked=0
for n in range(2,8):
    for _,T in sc(n//2,n%2):
        m=n//2; odd=n%2
        if evaluate(T,m,odd)!=(sum(table(T)[-1]),folded(T,m,odd)):
            raise ValueError((n,T))
        checked+=1
print('exact paired full-count / fixed-path cuts:',checked)
