import math, itertools
# ============ EVERYTHING AS GRAPH COLORINGS — the reframe table ============
print("=== reframe: every problem is a COLORING / a deletion-contraction (Tutte) invariant ===")
print("  unifying object: the TUTTE polynomial T(G;x,y); deletion-contraction T(G)=T(G-e)+T(G/e) = TIE INDUCTION")
print("  specializations: chromatic polynomial (proper colorings), independence polynomial (H, my partition fn), flow/reliability.")
for prob, color in [
 ("tournament H", "color = vertex ORDER (Ham path); 3-cycles = ODD cycles = the chromatic obstruction (chi>=3)"),
 ("unit distance", "color the PLANE so unit-apart differ = HADWIGER-NELSON, chi(R^2) in {5,6,7}; UDG=tie-graph, maximize edges => force chi up"),
 ("LRC", "color the CLOCK by nearest runner / by depth; lonely = an uncovered color; covering = a coloring"),
 ("Collatz", "color by PARITY (2-coloring) = the parity vector; even/odd = bipartite 2-coloring"),
 ("Goldbach/Lemoine", "color primes; even=node-pair (unordered), odd=arc (ordered) = edge 2-coloring (S630)"),
 ("forbidden H = 7,21", "3-cycle = triangle = K3 = chi 3 (odd cycle); resonance = non-bipartite = chromatic obstruction"),
]:
    print(f"  {prob:20}: {color}")

# "some nodes, some edges, some both" = the simplicial / Tutte structure (0-cells, 1-cells, 2-cells)
print("\n=== 'some nodes, some edges, some both' = 0/1/2-cells (vertex/edge/triangle colorings) ===")
print("  nodes (0-cells) = points/runners/numbers; edges (1-cells) = distances/arcs/relations; both (2-cells) = 3-cycles/unit-triangles (the resonance).")
print("  the resonance lives at the 2-cells (filled triangles) = the odd-cycle/chromatic obstruction = the pi/3 cube-root.")

# ============ TIE INDUCTION: chromatic number monotone; the 3-cycle/odd-cycle obstruction ============
print("\n=== the 3-cycle (resonance atom) = the minimal ODD cycle = the obstruction to 2-coloring ===")
print("  K3 (triangle) chi=3; even cycle chi=2 (bipartite); odd cycle chi=3. resonance = odd cycle = non-bipartite.")
print("  tie induction: G<=H (more ties = fewer edges) => chi(G)<=chi(H). adding ties can only LOWER chromatic number.")

# ============ HADWIGER-NELSON: the unit-distance chromatic number; the Moser spindle ============
def isprime(n):
    if n<2: return False
    for i in range(2,int(n**.5)+1):
        if n%i==0: return False
    return True
def chromatic_number(n, edges):
    adj=[set() for _ in range(n)]
    for a,b in edges: adj[a].add(b); adj[b].add(a)
    # greedy lower-bound via brute k-coloring search for small graphs
    for k in range(1,n+1):
        # backtracking k-coloring
        color=[-1]*n
        def bt(v):
            if v==n: return True
            for c in range(k):
                if all(color[u]!=c for u in adj[v]):
                    color[v]=c
                    if bt(v+1): return True
                    color[v]=-1
            return False
        if bt(0): return k
    return n
# Moser spindle (n=7): the famous 4-chromatic unit-distance graph
t=math.acos(5/6)
def rot(p,a): return (p[0]*math.cos(a)-p[1]*math.sin(a), p[0]*math.sin(a)+p[1]*math.cos(a))
O=(0,0); p1=(1,0); p2=rot(p1,math.pi/3); p3=(p1[0]+p2[0],p1[1]+p2[1])
q1=rot(p1,t); q2=rot(p2,t); q3=rot(p3,t)
spindle=[O,p1,p2,p3,q1,q2,q3]
def uedges(P):
    E=[]
    for i,j in itertools.combinations(range(len(P)),2):
        if abs(math.hypot(P[i][0]-P[j][0],P[i][1]-P[j][1])-1)<1e-9: E.append((i,j))
    return E
E=uedges(spindle)
print(f"\n=== HADWIGER-NELSON: chromatic number of unit-distance graphs ===")
print(f"  Moser spindle (n=7): {len(E)} unit edges, chromatic number = {chromatic_number(7,E)}  (forces chi(plane) >= 4; de Grey 2018: >=5)")
def tri(R): return [(i+j*0.5,j*math.sqrt(3)/2) for i in range(-R,R+1) for j in range(-R,R+1)]
def compact(pts,n): return sorted(pts,key=lambda p:p[0]**2+p[1]**2)[:n]
for n in [7,12,19]:
    P=compact(tri(4),n); E=uedges(P)
    print(f"  triangular patch n={n}: {len(E)} unit edges, chromatic number = {chromatic_number(n,E)} (triangular UDG is 3-chromatic = the hexagonal 3-coloring)")
