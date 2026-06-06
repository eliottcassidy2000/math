import math, itertools
# ================= REFRAME 1: signed speeds (directions) in LRC =================
def dZ(x): x=x-math.floor(x); return min(x,1-x)
print("=== REFRAME 1: negative speed = sigma (v->-v); loneliness is SIGN-INVARIANT (gauge) ===")
V=[1,2,3,4,5]; delta=1/6
import random; random.seed(0)
ok=True
for _ in range(5000):
    t=random.random()
    for eps in [[1]*5, [(-1)**i for i in range(5)], [random.choice([-1,1]) for _ in range(5)]]:
        d1=min(dZ(v*t) for v in V)              # all positive
        d2=min(dZ(e*v*t) for e,v in zip(eps,V)) # signed
        if abs(d1-d2)>1e-12: ok=False
print(f"  min_i ||v_i t|| invariant under arbitrary signs (every-other / random): {ok}  -> loneliness depends only on |speeds|")
print("  => the SIGN is a GAUGE; the physical object is the sigma-quotient. 'Every other clockwise' = the even/odd 2-coloring gauge (the 2-adic, S629/S630).")
# the RELATIVE-motion trick: from runner k's frame, others move at v_i - v_k (same dir) or v_i + v_k (opposite dir)
print("\n  the relative-motion trick (view from runner k): same-direction relative speed = v_i - v_k (DIFFERENCE);")
print("  opposite-direction = v_i + v_k (SUM). Directions turn the difference-set into a sum-set:")
k=3
print(f"   from runner {k}: same-dir rel speeds {sorted(set(abs(v-k) for v in V if v!=k))}; opp-dir rel speeds {sorted(set(v+k for v in V if v!=k))}")
print("   -> Goldbach(sum)/Lemoine vs difference; choosing directions selects the relative-speed lattice (the resonance gauge).")

# ================= REFRAME 2: the distance tournament (orient by </=/>1) =================
print("\n=== REFRAME 2: distance tournament — orient each pair by distance </=/>1 (ties = unit edges = the equator) ===")
def tri(R):
    return [(i+j*0.5, j*math.sqrt(3)/2) for i in range(-R,R+1) for j in range(-R,R+1)]
def compact(pts,n): return sorted(pts,key=lambda p:p[0]**2+p[1]**2)[:n]
def trichotomy(P):
    lt=eq=gt=0
    for i,j in itertools.combinations(range(len(P)),2):
        d=math.hypot(P[i][0]-P[j][0],P[i][1]-P[j][1])
        if d<1-1e-9: lt+=1
        elif d<1+1e-9: eq+=1
        else: gt+=1
    return lt,eq,gt
big=tri(4)
print(f"  {'n':>3} {'d<1 (arrows A)':>14} {'d=1 (TIES=unit)':>16} {'d>1 (arrows B)':>15}")
for n in [4,7,12,19]:
    P=compact(big,n); lt,eq,gt=trichotomy(P)
    print(f"  {n:>3} {lt:>14} {eq:>16} {gt:>15}")
print("  => the triangular lattice is UNIT-MINIMAL: NO pairs at distance <1 (arrow-class A is EMPTY); all non-unit pairs are >1 (class B).")
print("     so the lattice 'distance tournament' is degenerate: ties (unit) + one arrow class. Maximizing unit distances = maximizing the EQUATOR (ties).")

# the CM-coherence point: |beta|=1 in ONE embedding <=> ALL embeddings (CM property) => the tie-set is the unit group O_K^1 (algebraic), not a metric accident.
print("\n=== does the tournament framing simplify the grid disproof? the CM-coherence answer ===")
print("  unit distance = |z_i - z_j| = 1 = a TIE in the distance order = |beta|=1 for the lattice difference beta.")
print("  CM property (S623): |beta|=1 in one embedding <=> |beta|=1 in ALL embeddings. So the tie-set = O_K^1 (the unit group), an ALGEBRAIC subvariety,")
print("  not a metric accident. => the distance tournament's EQUATOR (ties) is algebraically COHERENT for the CM lattice but INcoherent for a generic lattice.")
print("  That coherence IS the disproof's mechanism: the CM tower inflates the coherent tie-set (split primes => many norm-1 units, S623).")
print("  SIMPLIFICATION (honest): the tournament reframe CLARIFIES (not trivializes) — it recasts 'max unit distances' as 'max the algebraically-coherent equator',")
print("  pinpointing WHY lattice (incoherent equator, ~6 ties/point) loses to CM-tower (coherent equator, n^0.014 ties/point). The orientation (</>) is the inside/outside the unit circle = the CM involution's sign.")
