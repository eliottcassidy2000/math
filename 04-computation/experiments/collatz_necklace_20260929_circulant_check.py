"""Fair consecutive splits of cycle necklaces make the cycle system circulant:
for T_{q,d}(y)=y/2 (even), (qy+d)/2 (odd), a cycle of shape (K,X) with a j-fold fair split
(j arcs of length m=K/j, each with x=X/j odd letters) at cut elements n_0..n_{j-1} and
share carries c_i satisfies 2^m n_{i+1} = q^x n_i + d c_i, hence for zeta=exp(2 pi i/j):
   N_k (2^m zeta^k - q^x) = d C_k,   N_k = sum_i zeta^{-ik} n_i,  C_k = sum_i zeta^{-ik} c_i,
and 2^K - q^X = prod_k (2^m - zeta^k q^x) = prod_{e|j} Phi_e(2^m, q^x).
This script verifies these identities exactly on real integer cycles."""
import sys
from math import gcd
from sympy import symbols, cyclotomic_poly, Poly, ZZ, rem, expand, Symbol, factorint
z=Symbol('z')

def T(y,q,d): return y//2 if y%2==0 else (q*y+d)//2
def carry(w,q):
    c=0
    for j,b in enumerate(w):
        if b: c=q*c+(1<<j)
    return c
def cycle_from(y0,q,d,maxlen=10**6):
    ys=[y0]; y=T(y0,q,d)
    while y!=y0:
        ys.append(y); y=T(y,q,d)
        if len(ys)>maxlen: return None
    return ys
def word(ys): return [y%2 for y in ys]
def fair_splits(w,j):
    K=len(w); m=K//j; X=sum(w); x=X//j; ww=w+w; out=[]
    for r in range(m):
        if all(sum(ww[r+i*m:r+(i+1)*m])==x for i in range(j)): out.append(r)
    return out
def zmod(expr,j):
    """reduce a polynomial in z modulo Phi_j(z) -> canonical sympy expr"""
    return rem(Poly(expand(expr),z,domain=ZZ),Poly(cyclotomic_poly(j,z),z,domain=ZZ)).as_expr()

def verify_cycle(ys,q,d,label):
    w=word(ys); K=len(w); X=sum(w); g=gcd(K,X)
    c=carry(w,q); D=(1<<K)-q**X
    assert ys[0]*D==d*c, "periodic point formula failed"
    print(f"\n{label}: least element {min(ys)}, shape (K,X)=({K},{X}), gcd={g}, clock 2^K-q^X={D if abs(D)<10**30 else 'big'}, min/max ratio {min(ys)/max(ys):.4f}")
    # primitivity of the word
    assert len(set(tuple(w[i:]+w[:i]) for i in range(K)))==K, "non-primitive word?"
    for j in range(2,g+1):
        if g%j: continue
        m=K//j; x=X//j
        # exact cyclotomic factorization of the clock
        prod=1
        for e in range(1,j+1):
            if j%e==0:
                # homogenized Phi_e(2^m, q^x) = (q^x)^phi(e) * Phi_e(2^m/q^x)
                Pe=Poly(cyclotomic_poly(e,z),z,domain=ZZ); deg=Pe.degree()
                val=sum(int(co)*(2**m)**k*(q**x)**(deg-k) for (k,),co in Pe.terms())
                prod*=val
        assert prod==D, "cyclotomic factorization failed"
        rs=fair_splits(w,j)
        print(f"  j={j}: (m,x)=({m},{x}); clock = prod_(e|j) Phi_e(2^m,q^x) verified; fair {j}-splits at rotations {rs if len(rs)<8 else str(rs[:8])+'...'} ({len(rs)} of {m})")
        for r in rs[:2]:
            ns=[ys[(r+i*m)%K] for i in range(j)]
            ww=w+w; cs=[carry(ww[r+i*m:r+(i+1)*m],q) for i in range(j)]
            for i in range(j):
                assert (1<<m)*ns[(i+1)%j]==q**x*ns[i]+d*cs[i], "affine share step failed"
            # DFT identities in Z[zeta_j]: N_k (2^m z^k - q^x) - d C_k == 0 mod Phi_j
            for k in range(j):
                Nk=sum(ns[i]*z**((-i*k)%j) for i in range(j))
                Ck=sum(cs[i]*z**((-i*k)%j) for i in range(j))
                assert zmod(Nk*((1<<m)*z**k-q**x)-d*Ck,j)==0, "DFT identity failed"
            s=sum(ns); sc=sum(cs)
            assert s*((1<<m)-q**x)==d*sc
            line=f"    r={r}: cut elements {ns if j<=4 else str(ns[:4])+'...'}; sum n_i = d*sum c_i/(2^m-q^x) ok"
            if j%2==0:
                a=sum((-1)**i*ns[i] for i in range(j)); ac=sum((-1)**i*cs[i] for i in range(j))
                assert a*((1<<m)+q**x)==-d*ac
                line+=f"; alt sum n_i = -d*alt sum c_i/(2^m+q^x) = {a} ok"
            print(line+"; all j DFT identities exact in Z[zeta_j]")
    return (K,X,g)

# ---- 1. small 3x+d census, both sheets (d odd, prime to 3, |d|<=Dmax), least |element| <= B|d|
Dmax=int(sys.argv[1]) if len(sys.argv)>1 else 100
B=int(sys.argv[2]) if len(sys.argv)>2 else 100
q=3
print("=== 1. primitive cycles of 3x+d, |d|<=%d, least element <= %d|d|, with gcd(K,X)>1"%(Dmax,B))
hits=[]
for ad in range(1,Dmax+1,2):
    if ad%3==0: continue
    for d in (ad,-ad):
        seen=set()
        for y0 in range(1,B*ad+1):
            if y0 in seen: continue
            y=y0; steps=0; path=[]
            while True:
                path.append(y); y=T(y,q,d); steps+=1
                if y<y0 or steps>4000: break
                if y==y0:
                    ys=cycle_from(y0,q,d); seen.update(ys)
                    if True:
                        w=word(ys); K=len(w); X=sum(w); g=gcd(K,X)
                        # primitive (word) and not a multiple of a cycle of a smaller d
                        prim= (gcd(gcd(*ys),ad)==1)
                        if g>1 and prim:
                            hits.append((d,min(ys),K,X,g))
                    break
print("d, least, K, X, gcd:", hits[:40], "... total", len(hits))
for (d,y0,K,X,g) in hits[:6]:
    verify_cycle(cycle_from(y0,3,d),3,d,f"3x{'+' if d>0 else ''}{d} cycle")

# ---- 2. Fermat-Catalan free cycles
print("\n=== 2. Fermat-Catalan / perfect-power clocks and their free cycles")
verify_cycle(cycle_from(65,3,-49),3,-49,"3x-49, free shape (5,4) from 2^5+7^2=3^4")
# 7x+169 shape (9,3): all necklaces free; list cycles and 3-fold fair splits
from itertools import combinations
necks=set()
for ones in combinations(range(9),3):
    w=[0]*9
    for p in ones: w[p]=1
    necks.add(min(tuple(w[i:]+w[:i]) for i in range(9)))
print(f"7x+169: {len(necks)} necklaces of shape (9,3); clock 2^9-7^3 = {2**9-7**3} = 13^2; 2^3-7 = 1")
n3=0
for wn in sorted(necks):
    x0=169*carry(list(wn),7)//(2**9-7**3)
    ys=cycle_from(x0,7,169)
    assert ys is not None and word(ys)==list(wn)[:len(ys)]
    if len(ys)<9:
        print(f"  word {''.join(map(str,wn))}: NON-primitive word (period {len(ys)}): cycle {ys} = 169 x cycle of 7x+1")
        continue
    prim = gcd(gcd(*ys),169)==1
    rs=fair_splits(list(wn),3)
    n3+= bool(rs)
    print(f"  word {''.join(map(str,wn))}: cycle min {min(ys)} max {max(ys)} {'primitive' if prim else 'NON-primitive (13 x cycle of 7x+13)'}; 3-fair splits at {rs}")
    if rs: verify_cycle(ys,7,169,"   7x+169 3-split instance")
print(f"  necklaces with a 3-fold fair split: {n3} of {len(necks)}")
# Eisenstein square behind 7^3+13^2=2^9
print("  2^3 - 7*z == (3 - z)^2 mod Phi_3(z):", zmod((8-7*z)-(3-z)**2,3)==0)
verify_cycle(cycle_from(343*carry([1,0,0,0,0,1,0,0,0],13)//(2**9-13**2),13,343),13,343,"13x+343, free shape (9,2) from 7^3+13^2=2^9 (word 100001000)")
verify_cycle(cycle_from(-4913*carry([1,0,0,0,1,0,0],71)//(2**7-71**2),71,-4913),71,-4913,"71x-4913, free shape (7,2) from 2^7+17^3=71^2 (word 1000100)")
verify_cycle(cycle_from(9*carry([1,1,0,0],5)//(2**4-5**2)*-1,5,-9),5,-9,"5x-9, free shape (4,2) from 2^4+3^2=5^2 (word 1100)")

# ---- 3. Belaga-Mignotte long cycle d=17021 least element 5: shape (2140,1088), gcd 4
print("\n=== 3. the Belaga-Mignotte long cycle of 3x+17021 (least element 5)")
ys=cycle_from(5,3,17021)
verify_cycle(ys,3,17021,"3x+17021 long cycle")
ys=cycle_from(101,3,14303)
verify_cycle(ys,3,14303,"3x+14303 long cycle")

# ---- 4. all free cycles of the perfect-power clocks with X>=2
print("\n=== 4. every necklace of a free perfect-power shape is an integer cycle (THM-4484 shift criterion)")
def all_necklaces(K,X):
    s=set()
    for ones in combinations(range(K),X):
        w=[0]*K
        for p in ones: w[p]=1
        s.add(min(tuple(w[i:]+w[:i]) for i in range(K)))
    return sorted(s)
for (q,d,K,X,idn) in [(3,-49,5,4,"2^5+7^2=3^4"),(7,169,9,3,"7^3+13^2=2^9"),(13,343,9,2,"7^3+13^2=2^9"),(71,-4913,7,2,"2^7+17^3=71^2"),(5,-9,4,2,"2^4+3^2=5^2"),(9,-49,5,2,"2^5+7^2=9^2"),(17,-225,6,2,"2^6+15^2=17^2")]:
    D=(1<<K)-q**X; assert abs(D)==abs(d)
    out=[]
    for wn in all_necklaces(K,X):
        x0=d*carry(list(wn),q)//D
        ys=cycle_from(x0,q,d); assert ys is not None
        prim=(len(ys)==K) and gcd(gcd(*ys),abs(d))==1
        out.append((''.join(map(str,wn)),min(ys),len(ys),'prim' if prim else f'period {len(ys)}, gcd {gcd(gcd(*ys),abs(d))}'))
    print(f" {q}x{'+' if d>0 else ''}{d} [{idn}], shape ({K},{X}), clock {D}: {len(out)} necklaces ->", out)
