"""Literal checks for signed modular path-cover transfer and its corollaries."""
from collections import defaultdict
from itertools import permutations, product
from math import comb, factorial

def tt(n):
    return tuple(tuple(i<j for j in range(n)) for i in range(n))

def circulant(n, out):
    return tuple(tuple((j-i)%n in out for j in range(n)) for i in range(n))

def conv(T):
    return tuple(zip(*T))

def edge_ok(T, p, anti=False):
    return all(T[i][j] == (T[p[j]][p[i]] if anti else T[p[i]][p[j]])
               for i in range(len(T)) for j in range(i))

def profiles(T, tau):
    n = len(T)
    O, Z = [0]*(n+1), [0]*(n+1)
    assert tau is None or (edge_ok(T,tau,True) and all(tau[tau[i]]==i for i in range(n)))
    for perm in permutations(range(n)):
        bad = sum(1<<i for i in range(n-1) if not T[perm[i]][perm[i+1]])
        need = bad.bit_count()
        for free in range(n-need):
            O[need+1+free] += comb(n-1-need,free)
        if tau is None or any(perm[i] != tau[perm[-1-i]] for i in range(n)):
            continue
        # For a fixed concatenated permutation, a palindromic ordered cover
        # is exactly a cut set invariant under index reflection.
        coeff = [1]
        for i in range(n//2):
            j = n-2-i
            assert bool(bad>>i&1)==bool(bad>>j&1)
            if bad>>i&1:
                continue
            weight = 1 if i==j else 2
            coeff = coeff + [0]*weight
            for r in range(len(coeff)-1,weight-1,-1):
                coeff[r]+=coeff[r-weight]
        for extra,c in enumerate(coeff):
            Z[need+1+extra]+=c
    return O,Z

def aut_data(T,tau):
    parts=defaultdict(list)
    for i,row in enumerate(T):
        parts[sum(row)].append(i)
    p=list(range(len(T))); g=c=0
    for blocks in product(*(permutations(a) for a in parts.values())):
        for src,dst in zip(parts.values(),blocks):
            for i,j in zip(src,dst):p[i]=j
        if edge_ok(T,p):
            g+=1
            c+=(tau is not None and all(p[tau[i]]==tau[p[i]] for i in range(len(T))))
    return g,c

def inflate(Q,Ts,sigma,local):
    addresses=[(i,j) for i,T in enumerate(Ts) for j in range(len(T))]
    book={v:i for i,v in enumerate(addresses)}
    T=tuple(tuple(Q[i][k] if i!=k else Ts[i][j][l] for k,l in addresses)
            for i,j in addresses)
    tau=tuple(book[sigma[i],local[i][j]] for i,j in addresses)
    return T,tau

def kernel(Q, sigma, factors, local):
    """Enumerate literal ordered quotient walk tuples and filter symmetry."""
    cap=tuple(map(len,factors))
    n=sum(cap); d=len(cap)
    prof=[profiles(T,local[i] if i==sigma[i] else None)
          for i,T in enumerate(factors)]
    O=[0]*(n+1); Z=[0]*(n+1)
    def visit(words,counts):
        k=len(words)
        if all(counts):
            value=1
            for i,c in enumerate(counts):value*=prof[i][0][c]
            O[k]+=value
            if tuple(tuple(sigma[z] for z in reversed(s)) for s in reversed(words))==words:
                val=1
                for i,c in enumerate(counts):
                    j=sigma[i]
                    if i==j:val*=prof[i][1][c]
                    elif i<j:val*=prof[i][0][c]
                Z[k]+=val
        for i in range(d):
            if counts[i]==cap[i]:continue
            counts[i]+=1
            if words and Q[words[-1][-1]][i]:
                visit(words[:-1]+(words[-1]+(i,),),counts)
            visit(words+((i,),),counts)
            counts[i]-=1
    visit((),[0]*d)
    return O,Z

def rev(n):return tuple(reversed(range(n)))

def ident(n):return tuple(range(n))

def tricol():
    Q=circulant(3,{1});s=(1,0,2)
    for m in (1,2,3):
        yield Q,[tt(m)]*3,s,[rev(m)]*3,"equal cyclic"
    yield Q,[Q]*3,s,[s]*3,"equal cyclic"
    for M in (tt(2),tt(3),Q):
        m=len(M);yield Q,[M,conv(M),tt(1)],s,[ident(m),ident(m),(0,)],"mirror norm"

def controls():
    yield from tricol()
    yield tt(2),[tt(3)]*2,(1,0),[rev(3)]*2,"no fixed quotient module"
    yield tt(3),[tt(2),circulant(3,{1}),tt(2)],(2,1,0),[rev(2),rev(3),rev(2)],"fixed middle"
    q=circulant(5,{1,2});s=tuple((-i)%5 for i in range(5))
    yield q,[tt(m) for m in (1,2,1,1,2)],s,[rev(m) for m in (1,2,1,1,2)],"prime quotient"


def color_moments(T):
    n=len(T)
    H={}
    for bits in range(1<<n):
        V=[i for i in range(n) if bits>>i&1]
        U=tuple(tuple(T[i][j] for j in V) for i in V)
        H[bits]=profiles(U,None)[0][1] if V else 1
    sums=[0]*3
    for j in (1,2,3):
        for colors in product(range(j),repeat=n):
            bins=[0]*j
            for u,i in enumerate(colors): bins[i]|=1<<u
            weight=1
            for b in bins: weight*=H[b]
            sums[j-1]+=weight
    return tuple(sums)

def stack(Ts):
    addresses=[(i,j) for i,T in enumerate(Ts) for j in range(len(T))]
    return tuple(tuple(i<k if i!=k else Ts[i][j][l] for k,l in addresses) for i,j in addresses)

def check_residue_palette():
    C3=circulant(3,{1})
    S4=[list(row) for row in tt(4)]
    S4[0][3]=False;S4[3][0]=True;S4=tuple(map(tuple,S4))
    T5=circulant(5,{1,2})
    C=color_moments(C3); A=color_moments(S4); B=color_moments(T5)
    assert C==(3,12,33) and A[0]==5 and B[0]==15
    moment_seeds=[tuple(x**t for x in C) for t in range(8)]
    for seed in (A,B):
        moment_seeds += [tuple(x*c**(2*t) for x,c in zip(seed,C)) for t in range(4)]
    assert { (4*c-11*h)%32 for h,b,c in moment_seeds } == set(range(1,32,2))
    for r in range(1,13):
        modulus=1<<r;period=1<<max(0,r-5)
        got={ (4*c*pow(9,ell,modulus)+5*b*pow(4,ell,modulus)-11*h)%modulus
              for h,b,c in moment_seeds for ell in range(r,r+period) }
        assert got==set(range(1,modulus,2)),(r,sorted(got))
    for D in ((),tt(1),C3):
        h,b,c=color_moments(D)
        # Independent direct computation of the literal inflated tournament.
        center=stack([tt(1),D,tt(1)])
        factors=[tt(2),tt(2),center]
        addresses=[(i,j) for i,T in enumerate(factors) for j in range(len(T))]
        T=tuple(tuple(C3[i][k] if i!=k else factors[i][j][l] for k,l in addresses) for i,j in addresses)
        assert profiles(T,None)[0][1]==4*c*9+5*b*4-11*h
    print("16-seed dyadic palette complete modulo every 2^r through r=12; 3 literal inflated path-count checks.")

def main():
    checks=0
    for Q,M,s,l,kind in controls():
        T,tau=inflate(Q,M,s,l)
        O,Z=profiles(T,tau)
        assert (O,Z)==kernel(Q,s,M,l),(kind,len(T))
        g,c=aut_data(T,tau)
        assert all(v%g==0 for v in O)
        assert all(v%c==0 for v in Z)
        assert all(O[k]%factorial(k)==0 and Z[k]%(2**(k//2)*factorial(k//2))==0 for k in range(1,len(T)+1))
        if kind=="equal cyclic":
            S=M[0]
            o,z=profiles(S,l[0]);h,k=aut_data(S,l[0])
            a=[x//h for x in o]+[0]
            b=[x//k for x in z]+[0]
            m=len(S)
            assert (g,c)==(3*h**3,h*k)
            assert O[1]//g==sum(a[i]**3+a[i+1]*a[i]**2+a[i+1]**2*a[i] for i in range(1,m+1))
            assert Z[1]//c==sum(a[i]*(b[i-1]+b[i]+b[i+1]) for i in range(1,m+1))
        if kind=="mirror norm":
            o,_=profiles(M[0],None); ng,_=aut_data(M[0],None)
            h,u=o[1],o[2]
            assert (g,c)==(ng**2,ng)
            assert O[1]==(h+u)**2+2*h*h and Z[1]==h+u
            assert O[1]%8==3
        print(kind,len(T),f"H={O[1]}",f"signed_fixed={Z[1]}",f"Aut={g}",f"centralizer={c}",flush=True)
        checks+=len(T)
    print(f"{checks} complete profile coordinates agree in each of the ordinary and signed comparisons.")
    check_residue_palette()

if __name__=="__main__":main()
