import math, itertools
print("=== HADWIGER-NELSON ↔ LRC ↔ UNIT DISTANCE: one object, three questions ===")
print("  the UNIT-DISTANCE GRAPH G of a homogeneous space (plane / circle-torus):")
for q,who in [
 ("max EDGES of G on n points", "UNIT DISTANCE problem"),
 ("chromatic number chi(G) of the whole space", "HADWIGER-NELSON"),
 ("max INDEPENDENT (1-avoiding / lonely) set alpha(G)", "LONELY RUNNER (the lonely set / p_0)"),
]:
    print(f"   {q:48} = {who}")
print("  BOUND: chi(G) * alpha(G) >= n  (a k-coloring's largest class is an independent set of size >= n/k).")
print("  => dense unit distances (small alpha, LRC hard) FORCE large chi (HN); a big lonely set (LRC) CAPS chi.")

print("\n=== the measurable bridge: chi_m(R^2) >= 1 / m_1, m_1 = max density of a 1-AVOIDING set ===")
print("  HN measurable chromatic number chi_m(plane) >= 1/m_1; m_1 in [0.229, 0.259] (Croft/Keleti) => chi_m >= 5 (Falconer).")
print("  LRC analogue: the lonely set is a 'forbidden-avoiding' set on the circle; its measure p_0 = the m_1-analogue.")
print("  => '1-avoiding density' (HN) and 'lonely measure p_0' (LRC) are THE SAME quantity on plane vs circle; both Delsarte-LP (S621).")

# the KEYS: insights of one apply to the other
print("\n=== the KEYS (insight transfer) ===")
keys=[
("HN -> LRC", "the Delsarte/Witsenhausen-Oler Fourier LP bounding the 1-avoiding density = the LP bounding p_0 (HYP-2215); the autocorrelation of the forbidden region."),
("LRC -> HN", "the resonance / additive-chain collapse (where p_0=0) = where the 1-avoiding set is forced small = the unit-distance-graph clusters that force chi up (the Moser/spindle rigidity)."),
("UD -> both", "the EXTREMAL configs: triangular lattice (dense unit, 3-colorable=hexagonal, delta=1/6) = the TAME corner of all three; CM-tower/non-lattice = the RECORDS (UD disproof, HN de Grey chi>=5) = leaving the lattice."),
("the pi/3 thread", "triangular UDG: chi=3 (hexagonal 3-coloring), densest unit distances, delta=1/6 chord; the cube-root/Phi_3 is the shared 3 = the lattice's tame value on all three."),
("tie/coloring (S633)", "G<=H => chi(G)<=chi(H) (tie induction); adding unit ties raises chi (HN) and shrinks alpha (LRC) and counts edges (UD) -- one monotone family."),
]
for a,b in keys: print(f"  {a:12}: {b}")

# verify chi*alpha>=n on the Moser spindle and triangular patches
def isE(P,i,j): return abs(math.hypot(P[i][0]-P[j][0],P[i][1]-P[j][1])-1)<1e-9
def edges(P): return [(i,j) for i,j in itertools.combinations(range(len(P)),2) if isE(P,i,j)]
def chi(n,E):
    adj=[set() for _ in range(n)]
    for a,b in E: adj[a].add(b); adj[b].add(a)
    for k in range(1,n+1):
        col=[-1]*n
        def bt(v):
            if v==n: return True
            for c in range(k):
                if all(col[u]!=c for u in adj[v]):
                    col[v]=c
                    if bt(v+1): return True
                    col[v]=-1
            return False
        if bt(0): return k
    return n
def alpha(n,E):
    adj=[set() for _ in range(n)]
    for a,b in E: adj[a].add(b); adj[b].add(a)
    best=0
    for r in range(n,0,-1):
        for S in itertools.combinations(range(n),r):
            if all(b not in adj[a] for a,b in itertools.combinations(S,2)):
                return r
    return 0
t=math.acos(5/6)
def rot(p,a): return (p[0]*math.cos(a)-p[1]*math.sin(a), p[0]*math.sin(a)+p[1]*math.cos(a))
O=(0,0);p1=(1,0);p2=rot(p1,math.pi/3);p3=(p1[0]+p2[0],p1[1]+p2[1]);q1=rot(p1,t);q2=rot(p2,t);q3=rot(p3,t)
spindle=[O,p1,p2,p3,q1,q2,q3]
def trip(R): return [(i+j*0.5,j*math.sqrt(3)/2) for i in range(-R,R+1) for j in range(-R,R+1)]
def comp(p,n): return sorted(p,key=lambda x:x[0]**2+x[1]**2)[:n]
print("\n=== verify chi * alpha >= n on unit-distance graphs ===")
for name,P in [("Moser spindle n=7",spindle),("tri patch n=7",comp(trip(3),7)),("tri patch n=12",comp(trip(4),12))]:
    n=len(P); E=edges(P); c=chi(n,E); a=alpha(n,E)
    print(f"  {name:18}: n={n} edges={len(E)} chi={c} alpha={a} -> chi*alpha={c*a} >= n={n}? {c*a>=n}")
