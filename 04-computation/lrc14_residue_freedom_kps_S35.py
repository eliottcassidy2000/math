"""kps-S35: THE EXTREMAL ARGUMENT. At q* ~ 13 ln M (max divisibility), the runners are HEAVY
(divisibility part F_i ~ M), so their residue freedom mod q* -- the # of multiples of F_i that
are <= M, = floor(M/F_i) -- collapses toward 1. With ~1 choice, the residue mod q* is FORCED
(not adversarially settable) => generic => witness (f_q -> 0). Quantify the freedom collapse."""
from math import gcd, log
def lcm(a,b): return a*b//gcd(a,b)
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

def blocker_with_parts(N):
    """strongest compressed divisibility-blocker in [N,2N]; return (family, parts F_i, Q)."""
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
    bins=build(Q)
    parts=[b[0] for b in bins]
    while len(parts)<13: parts.append(2)
    fam=[F*max(1,(N+F-1)//F) for F in parts]   # smallest multiple of F_i in [N,2N]
    return fam, parts, Q

print(f"{'N':>16} {'q*':>5} {'min F_i':>10} {'max F_i':>16} {'min freedom':>11} {'median freedom':>14} {'#forced(freedom<q*)':>19}")
for e in range(3,16):
    N=10**e
    fam,parts,Q=blocker_with_parts(N)
    M=max(fam)
    qs=first_free(fam)
    # residue freedom of runner i mod q* = number of multiples of F_i in [1,M] = floor(M/F_i)
    freedoms=sorted(M//F for F in parts)
    forced=sum(1 for f in freedoms if f < qs)     # freedom < q* => can't reach all residues mod q*
    print(f"{N:>16,} {qs:>5} {min(parts):>10,} {max(parts):>16,} {freedoms[0]:>11} {freedoms[len(freedoms)//2]:>14} {forced:>19}/13")
print()
print("READING: freedom_i = floor(M/F_i) = # settable residues mod q* for runner i.")
print("If freedom_i < q* the runner CANNOT hit an arbitrary residue mod q* (its residue is")
print("partly forced). #forced/13 near 13 => the adversary has ~no room to ALIGN the danger")
print("sets at q* => residues generic => witness (f_q->0). This is `alignment costs magnitude`")
print("at the residue level: max divisibility (heavy F_i) leaves no cofactor room to cover q*.")
