"""kps-S34 (corrected): the hybrid adversary can only COVER free PRIMES (composite free q
have gcd(D_i,q)>1 -- runners already carry the small factors). Measure how far covering free
primes pushes q_min, and whether it stays O(log M)."""
from math import gcd
import random
random.seed(3)
def lcm(a,b): return a*b//gcd(a,b)
def danger(v,a,q): r=(v%q)*a%q; return 14*min(r,q-r)<q
def is_lonely(sp,a,q): return not any(danger(v,a,q) for v in sp)
def blocked(sp,q): return not any(is_lonely(sp,a,q) for a in range(1,q))
def q_min(sp,Qm=2000):
    for q in range(2,Qm+1):
        for a in range(1,q):
            if is_lonely(sp,a,q): return q
    return None
def first_free(sp,Qm=2000):
    for q in range(2,Qm+1):
        if all(v%q for v in sp): return q
    return None
def is_prime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
def danger_set(r,q): return set(a for a in range(q) if 14*min((r*a)%q,q-(r*a)%q)<q)
def cover_residues(q,restarts=400):
    ds={r:danger_set(r,q) for r in range(1,q)}
    for _ in range(restarts):
        unc=set(range(1,q)); ch=[]
        for step in range(13):
            r=(random.choice(range(1,q)) if step==0 and random.random()<0.6
               else max(range(1,q),key=lambda r:len(ds[r]&unc)))
            ch.append(r); unc-=ds[r]
            if not unc: return ch
    return None
def crt(a1,m1,a2,m2):
    inv=pow(m1,-1,m2); x=(a1+m1*((a2-a1)*inv%m2))%(m1*m2); return x,m1*m2
def div_blocker(N,Q):
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
    while len(Ds)<13: Ds.append(2)   # padding: multiples of 2 (already covered)
    return Ds[:13]

def hybrid(N, want_cover, cap_mult):
    Q=14
    while div_blocker(N,Q+1) is not None and Q<600: Q+=1
    Ds=div_blocker(N,Q)
    if Ds is None: return None
    fam0=sorted(set(D*max(1,(N+D-1)//D) for D in Ds))
    # free PRIMES above Q
    frees=[]; q=Q+1
    while len(frees)<want_cover and q<Q+300:
        if is_prime(q) and all(v%q for v in fam0): frees.append(q)
        q+=1
    covsets=[]
    for qf in frees:
        cs=cover_residues(qf)
        if cs is None: break
        covsets.append((qf, cs + [1]*(13-len(cs))))
    fam=[]
    for i,D in enumerate(Ds):
        a,m=0,D
        bad=False
        for (qf,cs) in covsets:
            if gcd(m,qf)!=1: bad=True; break
            a,m=crt(a%m,m,cs[i]%qf,qf)
        if bad: return ('gcd_fail',fam0,Q,frees,[])
        # land v>=N, v≡a mod m, within [N, cap_mult*N]
        if a<N: v=a+m*((N-a+m-1)//m)
        else: v=a
        fam.append(v)
    fam=sorted(set(fam))
    if len(fam)<13: return ('collide',fam0,Q,frees,[])
    covered=[qf for (qf,_) in covsets if blocked(fam,qf)]
    return ('ok',fam,Q,frees,covered)

import math
print(f"{'N':>10} {'Q_div':>6} {'first_free':>10} {'qmin_div':>9} | {'want':>4} {'qmin_hyb':>9} {'#covered':>9} {'ratio':>7} {'clogM':>7}")
for e in range(3,9):
    N=10**e
    Q=14
    while div_blocker(N,Q+1) is not None and Q<600: Q+=1
    Ds=div_blocker(N,Q)
    fam0=sorted(set(D*max(1,(N+D-1)//D) for D in Ds))
    ff=first_free(fam0); qm0=q_min(fam0)
    for want in [6]:
        res=hybrid(N, want, cap_mult=10**9)
        if res and res[0]=='ok':
            _,fam,_,frees,covered=res
            qmh=q_min(fam,2000); mr=max(fam)/min(fam)
            print(f"{N:>10,} {Q:>6} {str(ff):>10} {str(qm0):>9} | {want:>4} {str(qmh):>9} {len(covered):>9} {mr:>7.1f} {qmh/math.log(N):>7.1f}" if qmh else f"{N:>10,} qmh=None")
        else:
            print(f"{N:>10,} {Q:>6} {str(ff):>10} {str(qm0):>9} | hybrid={res[0] if res else None}")
print()
print("READING: #covered = free primes the adversary blocked by residue-alignment.")
print("q_min(hybrid) vs c*log(M): if the multiplier stays bounded, census = Theta(log M).")
print("Covering costs 13x (all 13 residues fixed); divisibility 1x -> adversary reach O(log M).")
