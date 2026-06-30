"""
The construction's THREE-DISTANCE structure (mac-mini's tool, concrete). At witness zeta_6, speeds map to
the AP; compute the gaps (Steinhaus <=3 lengths). And the antipodal +1: M = n/(n(n-1)+1) = block 1/(n-1)
perturbed by the lcm-killer adding 1 to the denominator.
"""
from collections import Counter
def constr_gaps(n):
    D=n*n-n+1
    S=list(range(1,n-1))+[(n-1)*n]
    pts=sorted((s*n)%D for s in S)  # speeds * zeta_6 mod D
    gaps=[]
    for i in range(len(pts)):
        nxt=pts[(i+1)%len(pts)]
        g=(nxt-pts[i])%D
        gaps.append(g)
    # the gap containing 0 = the wrap gap (from max pt to min pt through 0)
    zero_gap=(pts[0]-pts[-1])%D
    return pts,gaps,D,zero_gap
print("CONSTRUCTION three-distance gaps (Steinhaus): at zeta_6 the speeds are an AP; gaps <=3 lengths.")
print(f"{'n':>3} {'distinct gaps':>20} {'0-gap':>6} {'M=n/D':>9} {'min dist from 0':>15}")
for n in [7,8,10,14,20]:
    pts,gaps,D,zg=constr_gaps(n)
    dg=sorted(set(gaps))
    cnt=Counter(gaps)
    md=min(min(p,D-p) for p in pts)
    print(f"{n:>3} {str(dg):>20} {zg:>6} {n}/{D}  min-dist={md} (=n)  gaps={dict(cnt)}")
print("\n  => THREE distinct gaps {1, n, 2n}: the AP step n (x n-3), the antipodal-closure gap 1 (x1),")
print("     and the 0-gap 2n (x1, = n+n, contains 0). M = n/D (0 sits centered in the 2n-gap).")
print("     The '1' gap is the killer closing the AP antipodally (n-2)n -> (n-1)^2 = D-n, distance 1.")
print()
print("THE ANTIPODAL +1 (why n/Phi_6, building on klein's lcm-killer + antipodal binding):")
for n in [7,10,14,20]:
    D=n*n-n+1
    print(f"  n={n}: block {{1..n-2}} alone has M=1/(n-1)=1/{n-1}={1/(n-1):.5f}; killer n(n-1)={n*(n-1)} = lcm(n-1,n) = -1 mod D")
    print(f"        => denominator n(n-1) -> n(n-1)+1 = {n*(n-1)+1}=Phi_6, so M: n/(n(n-1)) -> n/(n(n-1)+1) = {n}/{D}={n/D:.5f}")
print("  => M(construction) = the dense block's 1/(n-1) with the lcm-killer adding +1 to the denominator")
print("     (the antipode, klein's (1,-1) binding). The killer kills t=1/(n-1) (gap 0) and shifts to zeta_6.")
print()
print("klein's DENSE-CORE partial bound: block {1..n-2} has lonely gap 1/(n-1) > 1/n; killer perturbs to")
print("  n/Phi_6 still > 1/n (margin (n-1)/(n.Phi_6) ~ 1/n^2). => conjecture PROVABLE on dense-core coverings.")
