"""kps-S34: THE CLOSURE TEST. q_min = first_free_modulus was EXACT for the divisibility
blocker. Is `q_min <= first_free` UNIVERSAL? i.e. does the FIRST FREE MODULUS always carry a
lonely witness -- even against an adversary that SACRIFICES divisibility to ALIGN danger sets
at the free modulus? If yes, LRC(14)-compressed reduces to: modular LRC(14) holds at the
first free modulus (a specific q* ~ log M). Try HARD to break it."""
from math import gcd
import random
random.seed(1)
def lcm(a,b): return a*b//gcd(a,b)

def danger(v, a, q):
    r = (v % q) * a % q; m = min(r, q-r); return 14*m < q      # dist(va/q) < 1/14

def is_lonely(speeds, a, q): return not any(danger(v,a,q) for v in speeds)
def blocked(speeds, q):    return not any(is_lonely(speeds,a,q) for a in range(1,q))
def is_covering(speeds):   return all(any(v%q==0 for v in speeds) for q in range(2,15))
def q_min(speeds, Qmax=600):
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if is_lonely(speeds,a,q): return q
    return None
def first_free(speeds, Qmax=600):
    for q in range(2,Qmax+1):
        if all(v%q for v in speeds): return q
    return None

# ---- Adversary 1: divisibility blocker (baseline, q_min=first_free expected) ----
def blocker(N, Q):
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
    fam=[]
    for b in bins:
        l=b[0]; r=l*max(1,(N+l-1)//l)
        if r>2*N: r=l if N<=l<=2*N else None
        if r is None: return None
        fam.append(r)
    k=0
    while len(fam)<13: fam.append(N+ k +1); k+=1
    return fam[:13]

# ---- Adversary 2: TRY to block the first free modulus q* by choosing runners
#      whose residues mod q* make the 13 danger-sets cover Z/q*. ----
def try_block_free(N, Q, tries=4000):
    """block {2..Q} by divisibility AND try to cover the first-free q* by residue alignment.
    Uses cofactor freedom: runner = base * c, c chosen so runner mod q* hits a target."""
    base = blocker(N,Q)
    if base is None: return None
    fam = sorted(set(base))
    if len(fam)<13: return None
    qs = first_free(fam)
    if qs is None: return None
    # can we perturb runners (keep divisibility, keep in [N,2N], change residue mod qs)
    # to make qs blocked? search cofactor multiples.
    best = fam
    for _ in range(tries):
        cand=[]
        ok=True
        for v in fam:
            # v = g * (v//g) ; keep divisibility by all d|v with d<qs; nudge by +/- multiples of that g'
            # simplest: try v + t*L where L = product of small primes dividing v (keeps covering-ish)
            # here: try nearby multiples of the largest divisor <= 2N structure -- approximate:
            # pick a divisor d of v that is the "blocking part", replace cofactor
            d = 1
            for p in range(2,Q+1):
                if v%p==0 and lcm(d,p)<=v and lcm(d,p) < N: d=lcm(d,p)
            # runner must stay divisible by d, land in [N,2N]
            lo=(N+d-1)//d; hi=(2*N)//d
            if hi<lo: cand.append(v); continue
            c=random.randint(lo,hi)
            cand.append(d*c)
        cand=sorted(set(cand))
        if len(cand)<13: continue
        if not is_covering(cand): continue
        if first_free(cand)!=qs: continue          # keep same first free modulus
        if blocked(cand, qs):                       # SUCCESS: free modulus blocked!
            return ('BLOCKED_FREE', cand, qs)
    return ('free_has_witness', best, qs)

print("=== Test 1: is q_min <= first_free for the divisibility blocker + random covering? ===")
viol=0; tested=0
for e in range(3,10):
    N=10**e
    fam=blocker(N, 14)
    for Q in range(14, 300):
        f=blocker(N,Q)
        if f is None: break
        fam=f
    fam=sorted(set(fam))
    if len(fam)<13: continue
    qm=q_min(fam); ff=first_free(fam); tested+=1
    flag = "" if (qm is not None and ff is not None and qm<=ff) else "  <-- VIOLATION"
    if flag: viol+=1
    print(f"N=10^{e}: q_min={qm}  first_free={ff}{flag}")
# random compressed covering
print("--- random compressed covering families (ratio<=2) ---")
rv=0; rt=0
for _ in range(3000):
    N=random.choice([500,2000,10000,50000])
    fam=sorted(set(random.randint(N,2*N) for _ in range(13)))
    if len(fam)<13 or not is_covering(fam): continue
    qm=q_min(fam,200); ff=first_free(fam,200); rt+=1
    if qm is not None and ff is not None and qm>ff:
        rv+=1
        if rv<=5: print(f"  VIOLATION fam={fam} q_min={qm} first_free={ff}")
print(f"random compressed covering: {rt} tested, {rv} with q_min>first_free")

print()
print("=== Test 2: can an adversary BLOCK the first free modulus (align danger sets)? ===")
blocked_free=0; attempts=0
for e in range(3,9):
    N=10**e
    Q=14
    while blocker(N,Q+1) is not None and Q<300: Q+=1
    res=try_block_free(N, Q, tries=3000)
    if res is None: continue
    attempts+=1
    tag,fam,qs=res
    print(f"N=10^{e}: first_free q*={qs}  ->  {tag}")
    if tag=='BLOCKED_FREE': blocked_free+=1
print(f"\nadversary blocked the first free modulus in {blocked_free}/{attempts} magnitudes")
print()
print("READING: q_min<=first_free universal + no adversary blocks q* => THE CLOSURE:")
print("LRC(14)-compressed reduces to `modular LRC(14) holds at the first free modulus q*~log M`.")
