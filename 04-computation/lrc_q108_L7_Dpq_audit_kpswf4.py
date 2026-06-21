from fractions import Fraction as Fr

# Audit the L7 proof: Var(h_j) and the exact D_{p,q} vs 14/p.
# h_j(a) = Leb{ s in [0,L) : sec(a+s)=j } where sec(y)=floor(7*frac(y)).
# = overlap of [a,a+L) with the union of strips [j/7+n, (j+1)/7+n).
# Periodic in a with period 1. Trapezoid. Var over one period.

def hj(a, L, j):
    # a in [0,1) say; integrate length of [a,a+L) landing in sector j (mod 1, mod 7-sectors)
    # exact: sample the strip structure. Strip I_j = [j/7,(j+1)/7) mod 1.
    # overlap of [a,a+L) with the 1-periodic strip family {[j/7+n,(j+1)/7+n)}.
    lo = a; hi = a + L
    total = Fr(0)
    n0 = int(lo) - 1
    n1 = int(hi) + 2
    for n in range(n0, n1):
        s0 = Fr(j,7) + n
        s1 = Fr(j+1,7) + n
        left = max(lo, s0); right = min(hi, s1)
        if right > left:
            total += right - left
    return total

def total_variation(L, j, denom=2):
    # compute TV of a -> hj(a,L,j) over a in [0,1). Breakpoints where slope changes:
    # at a = (j/7 + n) - {0 or L} mod 1. Use a fine exact grid on the breakpoints.
    # breakpoints: a such that lo or hi crosses a strip edge: a = j/7+n or a=(j+1)/7+n
    #   or a+L = j/7+n or a+L=(j+1)/7+n.
    bps = set()
    for n in range(-2,4):
        for edge in (Fr(j,7)+n, Fr(j+1,7)+n):
            bps.add(edge % 1)
            bps.add((edge - L) % 1)
    bps = sorted(b for b in bps if 0 <= b < 1)
    bps.append(Fr(1))
    # sample midpoints, compute hj at the breakpoints themselves (piecewise linear -> TV = sum |jumps in value| at vertices)
    vals = [hj(b % 1 if b<1 else Fr(0), L, j) for b in bps]
    # Actually hj is continuous piecewise-linear; TV = sum over consecutive bps of |hj(b_{i+1})-hj(b_i)|
    tv = Fr(0)
    pts = bps
    hvals = [hj(p, L, j) for p in pts[:-1]] + [hj(Fr(0),L,j)]
    for i in range(len(pts)-1):
        tv += abs(hj(pts[i+1] if pts[i+1]<1 else Fr(0), L, j) - hj(pts[i], L, j))
    return tv

# test for several p/q in (1, 43/20]
print("Audit Var(h_j) = 2/7 for L=p/(7q) with p/q in (1,43/20]:")
for (p,q) in [(2,1),(3,2),(5,3),(7,4),(7,5),(43,20),(5,4),(4,3)]:
    L = Fr(p,7*q)
    tvs = [total_variation(L,j) for j in range(7)]
    print(f"  p/q={p}/{q}: L={L}={float(L):.4f}  TV(h_j) for j=0..6: {[str(t) for t in tvs]}  all==2/7: {all(t==Fr(2,7) for t in tvs)}")

# Now directly compute D_{p,q} and compare to 14/p
def sec(y):
    f = y % 1
    return int(7*f) if int(7*f)<7 else 6
def mu(i,j,p,q):
    # Leb{v in [0,1): sec(qv)=i and sec(pv)=j}
    bps = {Fr(0),Fr(1)}
    for c in (p,q):
        for t in range(0,7*c+1):
            bps.add(Fr(t,7*c))
    bps = sorted(b for b in bps if 0<=b<=1)
    tot = Fr(0)
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid=(a+b)/2
        if sec(q*mid)==i and sec(p*mid)==j:
            tot += b-a
    return tot
def D(p,q):
    tot=Fr(0)
    for i in range(7):
        for j in range(7):
            tot += abs(mu(i,j,p,q)-Fr(1,49))
    return tot
print("\nDirect D_{p,q} vs 14/p:")
maxDp=Fr(0)
for (p,q) in [(2,1),(3,2),(5,3),(7,4),(7,5),(43,20),(5,4),(4,3),(9,5),(11,6),(13,7)]:
    if __import__('math').gcd(p,q)!=1: continue
    d=D(p,q); bound=Fr(14,p)
    print(f"  p/q={p}/{q}: D={float(d):.5f}  14/p={float(bound):.5f}  D<=14/p: {d<=bound}  D*p={float(d*p):.4f}")
    if d*p>maxDp: maxDp=d*p
print(f"  sup D*p over these = {float(maxDp):.4f}  (proof claims <=14; verified sup ~20/7={float(Fr(20,7)):.4f})")
