"""Discriminant-glue certificates for pure-A rank 24 and an exact Golay tape."""
from itertools import permutations, product
from collections import Counter
from math import prod

M={
 3:((0,1,1,1,1,1),(1,0,1,2,2,1),(1,1,0,1,2,2),
    (1,2,1,0,1,2),(1,2,2,1,0,1),(1,1,2,2,1,0)),
 5:((1,2,2),(2,1,3),(2,3,1)),
 7:((2,3),(4,2)),
 13:((5,),)
}

def span(p,G):
    n=len(G[0])
    return {tuple(sum(a*g[j] for a,g in zip(x,G))%p for j in range(n))
            for x in product(range(p),repeat=len(G))}

def check_glue(p,G):
    C=span(p,G);r=len(G[0])
    assert len(C)**2==p**r and r*(p-1)==24
    weights=[sum(x*(p-x) for x in v) for v in C if any(v)]
    assert min(weights)==4*p and all(w%(2*p)==0 for w in weights)
    return C

def signed_swap(v,p,pairing,signs,twist):
    z=[0]*len(v)
    for (i,j),a in zip(pairing,signs):
        z[i]=a*v[j]%p
        z[j]=twist*a*v[i]%p
    return tuple(z)

def glues():
    signs={3:(1,1,-1,-1,-1,1),5:(1,1,1),7:(1,1),13:(1,)}
    for p,A in M.items():
        m=len(A)
        G=[tuple(int(i==j) for j in range(m))+A[i] for i in range(m)]
        C=check_glue(p,G)
        E=tuple((i,i+1) for i in range(0,2*m,2))
        assert all(signed_swap(v,p,E,signs[p],-1) in C for v in G)
        # Obstruction criterion: m must be even and its half-dimensional
        # hyperbolic discriminant must equal that of a square determinant.
        assert m%2 or pow((-1)**(m//2),(p-1)//2,p)!=1
        if p==7:
            image={}
            for P in permutations(range(4)):
                mult=0
                for s in product((1,-1),repeat=4):
                    if all(tuple(s[i]*v[P[i]]%7 for i in range(4)) in C
                           for v in G):mult+=1
                if mult:image[P]=mult
            even={P for P in permutations(range(4))
                  if sum(P[i]>P[j] for i in range(4) for j in range(i+1,4))%2==0}
            assert set(image)==even and set(image.values())=={2}
    G=[]
    for i in range(4):
        g=[0]*8;g[i:i+4]=[3,1,2,1];g[7]=-sum(g)%4
        G.append(tuple(g))
    C=check_glue(4,G)
    E=((0,1),(2,5),(3,7),(4,6))
    assert all(signed_swap(v,4,E,[1]*4,1) in C for v in G)

def golay():
    # Binary QR generator on 23 finite positions, parity at infinity.
    G=[0xae3<<i for i in range(12)]
    G=[g|((g.bit_count()%2)<<23) for g in G]
    C=[0]
    for g in G:C += [x^g for x in C]
    C=set(C)
    assert len(C)==2**12 and min(x.bit_count() for x in C if x)==8
    assert all(x.bit_count()%4==0 for x in C)
    # An explicit fixed-point-free projective involution.
    P=[23 if i==0 else 0 if i==23 else -pow(i,-1,23)%23 for i in range(24)]
    assert all(P[i]!=i and P[P[i]]==i for i in range(24))
    def f(x):return sum(((x>>i)&1)<<P[i] for i in range(24))
    assert all(f(g) in C for g in G)
    fixed={x for x in C if f(x)==x}
    assert {x^f(x) for x in C}==fixed
    assert Counter(x.bit_count() for x in fixed)=={0:1,8:15,12:32,16:15,24:1}
    twist=[x for x in C if f(x)==x^((1<<24)-1)]
    assert len(twist)==64 and {x.bit_count() for x in twist}=={12}
    d=min(twist)
    profile=Counter(((c&d).bit_count(),(c&~d).bit_count()) for c in twist)
    assert profile=={(0,12):1,(4,8):15,(6,6):32,(8,4):15,(12,0):1}
    k=36
    def marked(a,b):
        return 3**(a+b)*(4*2**(2*k-a-b)+3*2**(k-a)+2**(k-b)-1)
    assert len({marked(a,b) for a,b in profile})==5
    assert sum(v for (a,b),v in profile.items() if marked(a,b)==marked(b,a))==32
    print('twisted dodecad',hex(d),'physical split profile:',sorted(profile.items()))
    out=[]
    for k in (22,23,24):
        W=Counter((x&((1<<k)-1)).bit_count() for x in C)
        h=lambda a,b:2**(a+b+2)+3*2**a+2**b-1
        heights={h(1+w,1+k-w):m for w,m in W.items()}
        assert len(heights)==len(W)
        assert all(2**(k+4)<=x+1<2**(k+5) for x in heights)
        out.append(sorted(W.items()))
    assert [len(x) for x in out]==[11,8,5]
    assert h(13,13)==268468223
    print('Golay tape histograms:',out)

def incidence():
    n=13; D={0,1,3,9}
    A=[[int((j-i)%n in D) for j in range(n)] for i in range(n)]
    assert all(sum(v)==4 for v in A)
    assert all(sum(A[i][j] for i in range(n))==4 for j in range(n))
    assert all(sum(a*b for a,b in zip(A[i],A[j]))==1+3*(i==j)
               for i in range(n) for j in range(n))
    assert all(A[i][j]==A[-j%n][-i%n] for i in range(n) for j in range(n))
    for u,v in product(range(n),repeat=2):
        if v==5*u%n:
            assert -4*u%n==5*4*v%n
    for v,p in ((7,7),(21,3),(57,3)):
        assert v%p==0 and v%(p*p)!=0
        assert all((x*x+y*y)%p for x,y in product(range(p),repeat=2) if x or y)
    print('Singer polarity and glue similarities checked')

if __name__=='__main__':
    assert {h for h in range(2,26) if 24%(2*(h-1))==0}=={4,2,5,3,7,13}
    assert {h for h in range(2,26) if 24%(h-1)==0}=={2,3,4,5,7,9,13,25}
    glues()
    golay()
    incidence()
    print('All signed-swap, codimension and pitch checks passed.')
