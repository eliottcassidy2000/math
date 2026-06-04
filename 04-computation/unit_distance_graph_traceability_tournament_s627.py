import itertools, math
def tri_lattice(R):
    pts=[]
    for i in range(-R,R+1):
        for j in range(-R,R+1):
            pts.append((i+j*0.5, j*math.sqrt(3)/2))
    return pts
def unit_edges(P):
    n=len(P); E=set()
    for i in range(n):
        for j in range(i+1,n):
            if abs(math.hypot(P[i][0]-P[j][0],P[i][1]-P[j][1])-1)<1e-9: E.add((i,j))
    return E
def has_ham_path(n, adj):
    # adj: set of frozensets/edges; DP over subsets for a directed-agnostic Hamiltonian path
    nb=[0]*n
    for (i,j) in adj:
        nb[i]|=(1<<j); nb[j]|=(1<<i)
    full=(1<<n)-1
    # dp[mask][v] reachable path covering mask ending at v
    from functools import lru_cache
    dp=[[False]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=True
    for mask in range(1<<n):
        for v in range(n):
            if not dp[mask][v]: continue
            avail=nb[v] & ~mask
            w=avail
            while w:
                u=(w & -w).bit_length()-1; w&=w-1
                dp[mask|(1<<u)][u]=True
    return any(dp[full][v] for v in range(n))

def compact_subset(pts, n):
    return sorted(pts, key=lambda p: p[0]**2+p[1]**2)[:n]

print("=== triangular-lattice compact patches: UDG unit-edge count & Hamiltonian-path (traceable?) ===")
print(f"{'n':>3} {'unit_edges':>11} {'UDG_traceable':>14} {'~3n-sqrt(12n)':>14}")
big=tri_lattice(5)
for n in range(3,21):
    P=compact_subset(big,n); E=unit_edges(P)
    trace=has_ham_path(n,E)
    approx=3*n-math.sqrt(12*n)
    print(f"{n:>3} {len(E):>11} {str(trace):>14} {approx:>14.1f}")

# Moser spindle (n=7), famous unit-distance graph (chromatic 4)
print("\n=== special graphs ===")
# Moser spindle coordinates
a=math.acos(5/6)  # rhombus angle
def rot(p,t): return (p[0]*math.cos(t)-p[1]*math.sin(t), p[0]*math.sin(t)+p[1]*math.cos(t))
O=(0,0); 
# two rhombi sharing O, rotated by arccos(5/6)
p1=(1,0); p2=rot(p1, math.pi/3); p3=( (p1[0]+p2[0]), (p1[1]+p2[1]) )  # rhombus O,p1,p3,p2
t=math.acos(5/6)
q1=rot(p1,t); q2=rot(p2,t); q3=rot(p3,t)
spindle=[O,p1,p2,p3,q1,q2,q3]
E=unit_edges(spindle)
print(f"Moser spindle n=7: unit_edges={len(E)} UDG_traceable={has_ham_path(7,E)}")

# the tournament mapping: order points (canonical: by angle around centroid), base transitive, FLIP unit edges
def to_tournament(P, E):
    n=len(P); cx=sum(p[0] for p in P)/n; cy=sum(p[1] for p in P)/n
    order=sorted(range(n), key=lambda i: math.atan2(P[i][1]-cy, P[i][0]-cx))
    pos={v:k for k,v in enumerate(order)}
    beats=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j: continue
            a,b=min(pos[i],pos[j]),max(pos[i],pos[j])  # base: earlier beats later
            ei,ej=(i,j) if i<j else (j,i)
            unit=(ei,ej) in E
            # base transitive: order[a]->order[b]; flip if unit
            u,v=order[a],order[b]  # u earlier; base u->v; if unit flip to v->u
            if unit: beats[v][u]=True
            else: beats[u][v]=True
    return beats
def ham_paths_count(n,beats):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c: continue
            for nx in range(n):
                if not (mask>>nx)&1 and beats[last][nx]: dp[mask|(1<<nx)][nx]+=c
    return sum(dp[full][v] for v in range(n))
def cyclic_triples_unit_parity(P,E):
    # among all triples, count by #unit edges (0,1,2,3); 3 = equilateral unit triangle
    n=len(P); cnt={0:0,1:0,2:0,3:0}; eqtri=0
    for tri in itertools.combinations(range(n),3):
        u=sum(1 for x,y in itertools.combinations(tri,2) if (min(x,y),max(x,y)) in E)
        cnt[u]+=1
        if u==3: eqtri+=1
    return cnt, eqtri
print("\n=== tournament H (angular order, flip-unit) + equilateral-unit-triangle count ===")
print(f"{'n':>3} {'H(#hampaths)':>12} {'#eq_unit_triangles':>18} {'triple #unit-edge dist':>26}")
for n in range(3,12):
    P=compact_subset(big,n); E=unit_edges(P)
    beats=to_tournament(P,E); H=ham_paths_count(n,beats)
    cnt,eq=cyclic_triples_unit_parity(P,E)
    print(f"{n:>3} {H:>12} {eq:>18}   {dict(cnt)}")
