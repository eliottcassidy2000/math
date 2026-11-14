"""Literal certificates for 99*2^m + [13,15,...,35] and the independent overlaps."""
from itertools import product

def need(q):
    if not q:
        raise RuntimeError('incorrect corridor certificate')

def tt(n): return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def ordinal(A,B):return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)
def add(A,mask): return tuple(x|((1<<len(A)) if not mask>>i&1 else 0) for i,x in enumerate(A))+(mask,)
def star(I,U,W,ic,wc):return add(ordinal(ordinal(I,U),W),ic|(wc<<(len(I)+len(U))))
def reach(T,i):
    seen={i}; todo=[i]
    while todo:
        a=todo.pop()
        for j in range(len(T)):
            if T[a]>>j&1 and j not in seen:
                seen.add(j); todo.append(j)
    return seen
def components(T):
    R=[reach(T,i) for i in range(len(T))]
    return {tuple(j for j in sorted(ri) if i in R[j]) for i,ri in enumerate(R)}
def induce(T,vs):return tuple(sum(1<<k for k,j in enumerate(vs) if T[i]>>j&1) for i in vs)
def validate(T):
    need(all(not (r>>i&1) and r>>len(T)==0 for i,r in enumerate(T)))
    need(all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(len(T)) for j in range(i+1,len(T))))
def endtable(T):
    n=len(T); M=(1<<n)-1; E=[[0]*n for _ in range(M+1)]
    for i in range(n): E[1<<i][i]=1
    for S in range(1,M+1):
        for i,a in enumerate(E[S]):
            if a:
                r=T[i]&~S
                while r:
                    b=r&-r; r-=b; E[S|b][b.bit_length()-1]+=a
    return E
def profile(T,mask):
    M=(1<<len(T))-1
    e=endtable(T); p=endtable(tuple(M^r^(1<<i) for i,r in enumerate(T)))
    h=[1]+[sum(r) for r in e[1:]]
    W=[sum(k for i,k in enumerate(r) if mask>>i&1) for r in p]
    F=[sum(k for i,k in enumerate(r) if mask>>i&1) for r in e]
    E=[sum(k for i,k in enumerate(r) if not(mask>>i&1)) for r in e]
    return (h[M],
            sum(W[s]*h[M^s] for s in range(1,M+1)),
            sum(F[s]*W[M^s] for s in range(1,M)),
            h[M]+sum(E[s]*h[M^s] for s in range(1,M+1)),
            sum(endtable(add(T,mask))[-1]))
def count(T): return sum(endtable(T)[-1])
def seam(A,v,j):
    n=len(A); B=tt(j+1)
    return (tuple(r|(((1<<(j+1))-1)<<n if i!=v else ((1<<j)-1)<<n)
                  for i,r in enumerate(A))+
            tuple(r<<n|(1<<v if i==j else 0) for i,r in enumerate(B)))
def busch(m):
    P=tt(m+4); bits=(1,0,1,0)+(1,)*m+(0,)
    U=add(ordinal(P,(0,)),sum(x<<i for i,x in enumerate(bits)))
    A=add(U,1<<(len(U)-1))
    return seam(A,len(A)-1,1)

def main():
    from json import loads
    from pathlib import Path
    bank=loads(Path(__file__).with_name('corridor_witnesses.json').read_text())
    last=bank['last_offset']
    rails={int(c):row for c,row in bank['rails'].items()}
    initial={int(n):tuple(T) for n,T in bank['initial'].items()}
    need(set(rails)==set(range(1,last+1,2)))
    used=set()
    def pair(I,ic,W,wc):
        h,a,l,_,_=profile(I,ic); k,_,_,b,d=profile(W,wc)
        return (h,a,l),(k,b,d)
    for c,r in rails.items():
        I,ic,W,wc=tuple(r['I']),r['ic'],tuple(r['W']),r['wc']
        f,j=pair(I,ic,W,wc); h,a,l=f; k,b,d=j
        need(list(f)==r['Iprofile'] and list(j)==r['Wprofile'])
        shift=r['shift']
        need((a*b,h*d-k*l)==(99*2**shift,c))
        T=star(I,tt(1),W,ic,wc); validate(T)
        C=[S for S in components(T) if len(S)>1]
        need(len(C)==1 and len(I) in C[0])
        core=induce(T,C[0])
        need(count(core)==99*2**(shift+1)+c)
        for m in range(1,shift+1):
            n=99*2**m+c; used.add(n);T0=initial[n];validate(T0)
            need(count(T0)==n and len(components(T0))==1)
    need(used==set(initial))
    C3=(2,4,1)
    orig=(tt(2),3,(0,9,3,5,15),22)
    other=(tt(2),3,(8,1,3,6,15,31),62)
    for I,ic,W,wc in (orig,other):
        (h,a,l),(k,b,d)=pair(I,ic,W,wc)
        need((a*b,h*d-k*l)==(99,23))
    for U in (C3,tt(2)):
        x=count(U); h=[1]+[sum(e) for e in endtable(U)[1:]]
        Q=sum(a*b for a,b in zip(h,h[::-1]))
        for I,ic,W,wc in (orig,other):
            T=star(I,U,W,ic,wc); validate(T)
            need(count(T)==99*Q+23*x and len(components(T))==1)
    for m in (1,2,3):
        A=busch(m); validate(A)
        need(len(components(A))==1 and count(A)==99*2**m+23)
    for gamma,u,j,nq in ((2,1,2,6),(6,1,3,7),(2,3,5,6)):
        need(gamma*3**u+6*2**u-3==3*(2**j+1))
    print(f'99 corridor: {len(rails)} exact offsets, {len(initial)} strong phase bridges, independent overlap verified')
if __name__=='__main__':main()
