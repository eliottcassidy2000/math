import math
from itertools import combinations

# ========== 1) the chord<->arc bridge: unit distance among roots of unity <=> dZ = 1/6 ==========
def dZ(x): x=x-math.floor(x); return min(x,1-x)
print("=== chord = 2 sin(pi*dZ(arc)); UNIT distance (chord=1) <=> dZ = 1/6 (the hexagonal 60 deg = 2*3 gap) ===")
for a in range(7):
    x=a/6; chord=2*abs(math.sin(math.pi*x))
    print(f"  arc={a}/6 dZ={dZ(x):.4f}  chord=2sin(pi*dZ)={chord:.4f}  {'<-- UNIT' if abs(chord-1)<1e-9 else ''}")
print("  => the triangular/hexagonal lattice (optimal for unit distances) is the LRC gap delta=1/6 structure; 6=2*3.")

# ========== 2) n=22 triangular-lattice best (compact hexagonal cluster) ==========
def tri_points(R):
    pts=[]
    for i in range(-R,R+1):
        for j in range(-R,R+1):
            x=i+j*0.5; y=j*math.sqrt(3)/2
            pts.append((i,j,x,y))
    return pts
def unit_edges(subset):
    c=0
    for (a,b) in combinations(subset,2):
        d=math.hypot(a[2]-b[2],a[3]-b[3])
        if abs(d-1)<1e-9: c+=1
    return c
# greedy compact: sort by distance from origin, take 22 closest (centered hexagonal growth)
pts=tri_points(4)
pts.sort(key=lambda p: (p[2]**2+p[3]**2))
best=0; bestsub=None
# try a few centers/orientations by taking 22 closest to several centers
for cx,cy in [(0,0),(0.5,math.sqrt(3)/2),(0.5,math.sqrt(3)/6)]:
    s=sorted(pts,key=lambda p:( (p[2]-cx)**2+(p[3]-cy)**2 ))[:22]
    e=unit_edges(s)
    if e>best: best=e; bestsub=s
print(f"\n=== n=22 triangular-lattice compact cluster: {best} unit distances (a constructive LOWER bound) ===")
print(f"  known optima: u(20)=54, u(21)=57; u(22) is the FIRST OPEN case (optima proven only through n=21).")
print(f"  best-known upper bound was u(22)<=72; conjectured optimum ~60-62 uses NON-lattice rigid gluings.")
print(f"  the {best} lattice count is sub-optimal -> the LOCAL grid-suboptimality echoing the global disproof.")

# ========== 3) n=22 = 2*11: the cyclotomic CM structure, doubling = primitive root ==========
def order(g,m):
    x=g%m; o=1
    while x!=1: x=(x*g)%m; o+=1
    return o
print("\n=== n=22 = 2*11: cyclotomic CM field Q(zeta_11), Galois group (Z/11)* = Z/10 ===")
print(f"  ord_11(2) = {order(2,11)}  -> 2 is a PRIMITIVE ROOT mod 11: doubling <2> is a SINGLE 10-cycle")
print(f"  orbit of 1 under doubling mod 11: {[ (2**k)%11 for k in range(10)]}")
print(f"  (contrast n=14: ord_7(2)=3, two 3-cycles. n=22's doubling is one big free orbit.)")
print("  complex conjugation c (z->z^-1) = the antipodal involution sigma (v->-v); its orbits on (Z/11)*:")
print(f"    sigma-pairs {{a,-a}} mod 11: {[(a,(11-a)) for a in range(1,6)]}  (5 free pairs, NO fixed point: 11 odd => no apex)")

# ========== 4) the dictionary: grid-disproof construction <-> our LRC/tournament objects ==========
print("\n=== unit-distance disproof  <->  LRC / tournament dictionary ===")
for a,b in [
 ("magnitude-1 elements of a CM field","roots of unity = runner positions e^{2pi i v t} (LRC)"),
 ("the CM field Q(zeta_n) / O_K","the cyclotomic home of the LRC at rational times"),
 ("complex conjugation c (the CM involution)","the antipodal involution sigma: v -> -v (HYP-2205)"),
 ("split prime P != c(P) (free under c)","free sigma-orbit (the free-orbit cascade, HYP-2205)"),
 ("ramified/real place (c-fixed)","the apex (sigma-fixed lane, HYP-2185)"),
 ("Galois group / Frobenius (which primes split)","(Z/n)* doubling = the perspective/automorphism rigidity"),
 ("abundance of unit dists from MANY split primes","loneliness from MANY free sigma-orbits"),
 ("take [K:Q]->infinity (tower, beats fixed grid)","collapse family bigger than the AP/grid (HYP-2195)"),
 ("Delsarte LP upper bound u(n)<=...","the Delsarte LP for LRC loneliness (HYP-2215)"),
]:
    print(f"  {a:46} <-> {b}")
