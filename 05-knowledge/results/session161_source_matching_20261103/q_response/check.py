from array import array
from math import comb
from fractions import Fraction

def L(q,r): return ((r+1)**q+(r-1)**q)//2

def charge(P,r):
    b=bin(P)[2:]; k=len(b)
    return sum(int(bit)*r**i*L(k-1-i,r) for i,bit in enumerate(b))

def delta(j,q,r):
    P=2**j+2**q-1; R=P+2; k=j+1
    return 2*r**(k+1)*(charge(R,r)-charge(P,r))

def T(q,r):
    return L(q,r)-sum(r**a*L(q-a,r) for a in range(1,q))

def literal(j,q,plus):
    P=2**j+2**q-1; R=P+2
    p=bin(P)[2:]; s=''.join('1' if v=='0' else '0' for v in bin(R)[2:])[::-1]
    w=p+'01'+s; u=p+'10'+s; m=len(w); n=m+2
    edges=[0]*n
    for i in range(m):
        for z in range(i+1,m):edges[i] |= 1<<z
    a,b=m,m+1
    for i in range(m):
        if w[i]=='1':edges[a] |=1<<i
        else:edges[i]|=1<<a
        if u[i]=='1':edges[b] |=1<<i
        else:edges[i]|=1<<b
    if plus:edges[a]|=1<<b
    else:edges[b]|=1<<a
    return edges

def wrap(adj):
    n=len(adj); a,b=n,n+1; out=adj.copy()
    for i in range(n):out[i]|=1<<a
    out += [1<<b, (1<<n)-1]
    return out

def H(adj):
    n=len(adj); total=1<<n; inc=[0]*n
    for u in range(n):
        for v in range(n):
            if (adj[u]>>v)&1:inc[v]|=1<<u
    dp=array('Q',[0])*(total*n)
    for v in range(n): dp[(1<<v)*n+v]=1
    for mask in range(1,total):
        if not mask&(mask-1):continue
        bs=mask
        while bs:
            bit=bs & -bs; v=bit.bit_length()-1; bs-=bit
            prev=mask^bit; ps=inc[v]&prev; totalv=0
            while ps:
                q=ps & -ps; u=q.bit_length()-1; ps-=q
                totalv+=dp[prev*n+u]
            dp[mask*n+v]=totalv
    return sum(dp[(total-1)*n+v] for v in range(n))

def main():
    checks=0
    for j in range(2,15):
        for q in range(1,j):
            for r in range(1,8):
                d=delta(j,q,r)
                expected=2*r**(2*j+2-q)*T(q,r)
                closed=2*r**(2*j+2-q)*(r**q+((1-r)*(r+1)**q+(r+1)*(r-1)**q)//2)
                if d!=expected or d!=closed:raise ValueError((j,q,r,d,expected,closed))
                checks+=1
            wrapped=delta(j,q,1)+delta(j,q,2)
            formula=2+2**(2*j+2-q)*(2**(q+1)-3**q+3)
            if wrapped!=formula:raise ValueError((j,q,wrapped,formula))
            if (wrapped<0)!=(q>=3):raise ValueError((j,q,wrapped))
    print('charge controls:', checks, 'including q=1..j-1, j=2..14, r=1..7')
    for q in (1,2,3):
        j=4; am=literal(j,q,False); ap=literal(j,q,True)
        hm,hp=H(am),H(ap); wm,wp=H(wrap(am)),H(wrap(ap))
        d=delta(j,q,1); wd=delta(j,q,1)+delta(j,q,2)
        if hp-hm!=d or wp-wm!=wd:raise ValueError((j,q,hm,hp,wm,wp,d,wd))
        print('literal j,q, old minus/plus, wrapped minus/plus:',j,q,hm,hp,wm,wp,flush=True)
if __name__=='__main__':main()
