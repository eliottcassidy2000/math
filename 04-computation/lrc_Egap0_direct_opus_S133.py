import numpy as np, random
random.seed(133)
GRID=300000; xs=(np.arange(GRID)+0.5)/GRID
def Egap0(E):
    # gap containing 0 for phases {frac(e_i x)} (0 NOT a phase) = min(frac) + 1 - max(frac) = wrap gap
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0)
    return float((ph.min(axis=1) + 1 - ph.max(axis=1)).mean())
def Emaxgap(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
    return float(g.max(axis=1).mean())
thr=1/7
print("=== E[gap containing 0] (room at origin, 0 not a phase) -- direct-bound candidate (opus-S133) ===")
print(f"    E[maxgap] >= E[gap0]; if E[gap0] > 1/7={thr:.5f} for ALL E => density floor (clean fixed-phi)\n")
tests={"AP{1..13}":list(range(1,14)),"{2..14}":list(range(2,15)),
       "minEmaxgap shape":[0,2,3,4,5,6,7,8,9,10,12,17,28],  # note has 0 -- shift
       "{1,2,3,4,10..18}sat":[1,2,3,4,10,11,12,13,14,15,16,17,18],
       "primes":[2,3,5,7,11,13,17,19,23,29,31,37,41]}
minv=1.0; worst=None
for name,E in tests.items():
    E=[e for e in E if e!=0] or E
    if len(set(E))<13: E=list(range(1,14))
    g0=Egap0(E); mg=Emaxgap(E)
    if g0<minv: minv=g0; worst=(name,E)
    print(f"  {name:<20} E[gap0]={g0:.5f} {'>1/7' if g0>thr else '<1/7 !!'}   E[maxgap]={mg:.5f}")
# adversarial: MINIMIZE E[gap0] over 13-clusters
print("\n  adversarial minimization of E[gap0]:")
best=Egap0(list(range(1,14)))
for _ in range(40):
    E=random.sample(range(1,40),13); cur=Egap0(E)
    for _ in range(50):
        i=random.randrange(13); c=E.copy(); c[i]=random.randint(1,60)
        if len(set(c))<13: continue
        cv=Egap0(c)
        if cv<cur: E,cur=c,cv
    if cur<best: best,worstE=cur,sorted(E)
print(f"  MIN E[gap0] found = {best:.5f}  ({'>1/7 -- fixed-phi bound HOLDS' if best>thr else '<1/7 -- fixed-phi FAILS, need the max after all'})")
if best<=thr: print(f"    (min at {worstE})")
print(f"\n  1/7={thr:.5f}. => {'E[gap0]>1/7 uniformly => CLEAN direct bound E[maxgap]>=E[gap0]>1/7!' if best>thr else 'E[gap0] dips below 1/7 => the origin gap alone is not enough; the max-over-phi is needed.'}")
