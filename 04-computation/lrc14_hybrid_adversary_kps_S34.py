"""kps-S34: can a SMART HYBRID adversary (divisibility below + COVER the free moduli via
CRT-set residues) push q_min above the first free modulus? Tests the double-alignment cost:
divisibility costs log q (1 runner @ residue 0); covering costs ~13 log q (all 13 residues).
Measures how far q_min can be pushed, and whether it stays O(log M)."""
from math import gcd, prod
import random
random.seed(2)
def lcm(a,b): return a*b//gcd(a,b)
def danger(v,a,q): r=(v%q)*a%q; return 14*min(r,q-r)<q
def is_lonely(sp,a,q): return not any(danger(v,a,q) for v in sp)
def blocked(sp,q): return not any(is_lonely(sp,a,q) for a in range(1,q))
def is_cov(sp): return all(any(v%q==0 for v in sp) for q in range(2,15))
def q_min(sp,Qm=1200):
    for q in range(2,Qm+1):
        for a in range(1,q):
            if is_lonely(sp,a,q): return q
    return None
def first_free(sp,Qm=1200):
    for q in range(2,Qm+1):
        if all(v%q for v in sp): return q
    return None
def danger_set(r,q):
    return set(a for a in range(q) if 14*min((r*a)%q,q-(r*a)%q)<q)
def cover_residues(q,restarts=200):
    """13 nonzero residues covering Z/q, or None."""
    ds={r:danger_set(r,q) for r in range(1,q)}
    for _ in range(restarts):
        unc=set(range(1,q)); ch=[]
        for step in range(13):
            r=(random.choice(range(1,q)) if step==0 and random.random()<0.5
               else max(range(1,q),key=lambda r:len(ds[r]&unc)))
            ch.append(r); unc-=ds[r]
            if not unc: return ch
    return None

def crt(a1,m1,a2,m2):
    # x ≡ a1 mod m1, x ≡ a2 mod m2  (gcd(m1,m2)=1)
    g=gcd(m1,m2); assert g==1
    inv=pow(m1,-1,m2)
    x=(a1 + m1*((a2-a1)*inv % m2))%(m1*m2)
    return x, m1*m2

def divisibility_blocker(N,Q):
    bins=[]
    for q in range(2,Q+1):
        pl=False
        for b in bins:
            l=lcm(b[0],q)
            if l<=2*N: b[0]=l; pl=True; break
        if not pl:
            if q<=2*N: bins.append([q])
            else: return None
    if len(bins)>13: return None
    Ds=[b[0] for b in bins]
    while len(Ds)<13: Ds.append(1)
    return Ds[:13]

def build_hybrid(N, cover_upto_free=3, spread=200):
    """divisibility-block {2..Q}, then COVER the first `cover_upto_free` free moduli by
    assigning CRT-set residues. Returns family + info."""
    Q=14
    while divisibility_blocker(N,Q+1) is not None and Q<400: Q+=1
    Ds=divisibility_blocker(N,Q)
    if Ds is None: return None
    base=[max(1,D) for D in Ds]
    # find the first few free moduli of the plain divisibility family (approx): q>Q
    fam0=[]
    for D in base:
        v=D*max(1,(N+D-1)//D)
        fam0.append(v)
    frees=[]
    q=Q+1
    while len(frees)<cover_upto_free and q< Q+120:
        if all(v%q for v in fam0): frees.append(q)
        q+=1
    # assign covering residues at each free modulus
    covsets=[]
    for qf in frees:
        cs=cover_residues(qf, restarts=300)
        if cs is None: covsets=None; break
        covsets.append((qf,cs))
    if not covsets: return ('no_coverset', fam0, Q)
    # build each runner: v_i ≡ 0 mod D_i, v_i ≡ cs_j[i] mod qf_j  (CRT), land in [N, spread*N]
    fam=[]
    ok=True
    for i,D in enumerate(base):
        a,m=0,max(1,D)
        good=True
        for (qf,cs) in covsets:
            if gcd(m,qf)!=1: good=False; break
            a,m=crt(a%m,m, cs[i]%qf, qf)
        if not good: ok=False; break
        # smallest v>=N with v≡a mod m
        v=a + m*((N-a + m-1)//m if a<N else 0)
        if v<N: v+= m*((N-v)//m +1)
        fam.append(v)
    if not ok or len(set(fam))<13: return ('crt_fail', fam0, Q)
    return ('hybrid', sorted(set(fam)), Q, frees)

print(f"{'N':>10} {'Q_div':>6} {'first_free':>10} {'q_min(div)':>10} | {'q_min(hybrid)':>13} {'frees_covered':>14} {'maxratio':>9}")
import math
for e in range(3,9):
    N=10**e
    Ds=divisibility_blocker(N,14)
    Q=14
    while divisibility_blocker(N,Q+1) is not None and Q<400: Q+=1
    Ds=divisibility_blocker(N,Q)
    fam0=sorted(set(D*max(1,(N+D-1)//D) for D in Ds))
    if len(fam0)<13: continue
    ff=first_free(fam0); qm0=q_min(fam0)
    res=build_hybrid(N, cover_upto_free=4, spread=10**6)
    if res and res[0]=='hybrid':
        _,famh,_,frees=res
        qmh=q_min(famh, 1500)
        # which of the targeted frees are actually blocked
        covered=[qf for qf in frees if blocked(famh,qf)]
        mr=max(famh)/min(famh)
        print(f"{N:>10,} {Q:>6} {str(ff):>10} {str(qm0):>10} | {str(qmh):>13} {str(covered):>14} {mr:>9.1f}")
    else:
        tag=res[0] if res else 'None'
        print(f"{N:>10,} {Q:>6} {str(ff):>10} {str(qm0):>10} | hybrid={tag}")
print()
print("READING: if hybrid q_min stays ~ C*log(M) (a few free moduli beyond first_free),")
print("the census is Theta(log M): the adversary CAN cover a FEW free moduli (paying 13x)")
print("but runs out of CRT capacity -> q_min = O(log M). No fixed-q census, but O(log M) one.")
