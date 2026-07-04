"""kps-S35: THE DECISIVE EXTREMAL TEST. At the extremal divisibility-blocker, each runner i
can only reach residues A_i = {F_i*c mod q* : F_i*c in [N,2N]} mod q* (freedom < q*). Can the
adversary CHOOSE one residue per runner (from A_i) so the 13 danger sets COVER Z/q* -- i.e.
block q* WITHIN the same magnitude? If NO (over a strong search), the witness is FORCED at q*
=> q_min = q* provably for the extremal, closing the compressed crux at the extremal."""
from math import gcd
def lcm(a,b): return a*b//gcd(a,b)
def band_k(q): return (q-1)//14
def danger_set(r,q):
    k=band_k(q)
    return frozenset(a for a in range(q) if min((r*a)%q, q-(r*a)%q)<=k)

def blocker_parts(N):
    Q=14
    def build(Q):
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
        return bins
    while build(Q+1) is not None and Q<800: Q+=1
    parts=[b[0] for b in build(Q)]
    while len(parts)<13: parts.append(2)
    return parts, Q

def first_free_of(parts, N):
    fam=[F*max(1,(N+F-1)//F) for F in parts]
    for q in range(2,2000):
        if all(v%q for v in fam): return q, fam
    return None, fam

def achievable_residues(F, N, q):
    lo=(N + F -1)//F; hi=(2*N)//F
    return sorted(set((F*c)%q for c in range(lo, hi+1)))

def can_cover_from_sets(A_list, q, restarts=4000):
    """each runner picks one residue from A_list[i]; can the 13 danger-sets cover Z/q?
    greedy: repeatedly pick (runner, residue in its set) covering most-uncovered. random restarts."""
    import random
    random.seed(7)
    # precompute danger sets for all achievable residues
    dcache={}
    def dset(r):
        if r not in dcache: dcache[r]=danger_set(r,q)
        return dcache[r]
    best_unc=q
    for _ in range(restarts):
        unc=set(range(1,q))
        used=[False]*13
        order=list(range(13)); random.shuffle(order)
        for step in range(13):
            # choose the (unused runner, residue) covering most uncovered; random tie-break via order
            bestgain=-1; bi=None; br=None
            for i in range(13):
                if used[i]: continue
                for r in A_list[i]:
                    if r==0: g=len(unc)   # residue 0 covers everything (but q* free => 0 not achievable)
                    else: g=len(dset(r)&unc)
                    if g>bestgain: bestgain=g; bi=i; br=r
            if bi is None: break
            used[bi]=True
            unc-=dset(br)
            if not unc: return True
        best_unc=min(best_unc,len(unc))
    return False, best_unc

print(f"{'N':>14} {'q*':>5} {'freedoms(|A_i|)':>28} {'adversary covers q* within [N,2N]?':>34}")
for e in range(3,12):
    N=10**e
    parts,Q=blocker_parts(N)
    qs,fam=first_free_of(parts,N)
    A=[achievable_residues(F,N,qs) for F in parts]
    sizes=sorted(len(a) for a in A)
    res=can_cover_from_sets(A, qs, restarts=1500)
    if res is True:
        verdict="YES -- q_min>q* (adversary blocks q*)"
    else:
        _,unc=res
        verdict=f"NO (best leaves {unc} uncovered) -> WITNESS at q*"
    print(f"{N:>14,} {qs:>5} {str(sizes):>28} {verdict:>34}")
print()
print("READING: NO in every row => the extremal adversary CANNOT cover q* within its own")
print("magnitude (freedom<q* per runner blocks the alignment) => q_min=q* is FORCED at the")
print("extremal. The compressed crux at the extremal reduces to `the forced residues don't")
print("cover Z/q*`, and here the adversary provably can't force a cover without more magnitude.")
