"""
tournament_spectrum_discovery_kps.py  (kind-pasteur-2026-06-27-S31ah)

The programmatic TECHNIQUE GENERATOR: enumerate small tournaments, collect the
ACHIEVED spectrum of each invariant, and report the GAPS -- each persistent gap
is a candidate IMPOSSIBILITY CERTIFICATE (mechanizing "find a new {7,21}").

Discovers:
  (1) H-spectrum per n + gaps (should rediscover {7,21} and transient gaps).
  (2) Which SMALL graphs are realizable as Omega(T) (complement = forbidden
      Omega substructures, generalizing THM-201 K3 / THM-202 P4).
  (3) alpha-vector (cycle-conflict) spectrum + gaps.
"""
import itertools, sys, random
from tournament_certificate_engine_kps import (
    all_tournaments, random_tournament, H_value, alpha_vector,
    conflict_graph, odd_cycle_counts, independence_poly, scores)

# ---- canonical form for small graphs (<=7 vertices) ----
def graph_canon(m, E):
    """canonical adjacency tuple over all vertex permutations (m<=7)."""
    best=None
    for perm in itertools.permutations(range(m)):
        adj=tuple(1 if E[perm[i]][perm[j]] else 0
                  for i in range(m) for j in range(i+1,m))
        if best is None or adj<best: best=adj
    return (m, best)

def graph_name(m, canon_bits):
    """small human label for a canonical small graph."""
    edges=sum(canon_bits)
    if m==0: return "empty"
    if m==1: return "K1"
    full=m*(m-1)//2
    if edges==0: return f"{m}*K1 (no edges)"
    if edges==full: return f"K{m}"
    if m==2: return "K2"
    if m==3:
        return {1:"K2+K1",2:"P3"}.get(edges,"?3")
    if m==4:
        return {1:"K2+2K1",2:"(P3 or 2K2)",3:"(P4 or K3+K1 or star)",
                4:"C4/paw",5:"K4-e",6:"K4"}.get(edges,f"4v-{edges}e")
    return f"{m}v-{edges}e"

def all_small_graphs(m):
    """yield canonical forms of all graphs on m vertices."""
    seen=set(); pairs=list(itertools.combinations(range(m),2))
    for bits in range(1<<len(pairs)):
        E=[[False]*m for _ in range(m)]
        for k,(i,j) in enumerate(pairs):
            if bits>>k&1: E[i][j]=E[j][i]=True
        c=graph_canon(m,E)
        if c not in seen: seen.add(c); yield c
    return

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(31)

    print("="*70); print(" (1) H-SPECTRUM per n + gaps"); print("="*70)
    Hsets={}
    for n in range(3,7):
        S=set()
        for adj in all_tournaments(n):
            S.add(H_value(adj))
        Hsets[n]=S
        mx=max(S); gaps=[v for v in range(1,mx+1,2) if v not in S]
        print(f" n={n}: H achieved (sorted)={sorted(S)}")
        print(f"        odd gaps below max({mx}) = {gaps}")
    # n=7 sample
    S7=set()
    for _ in range(60000):
        S7.add(H_value(random_tournament(7, random)))
    mx=max(S7); gaps7=[v for v in range(1,mx+1,2) if v not in S7]
    print(f" n=7 (60k sample): gaps below {mx} = {sorted(gaps7)[:30]}{' ...' if len(gaps7)>30 else ''}")
    persistent=[v for v in (7,21) if all(v not in Hsets[n] for n in Hsets) and v not in S7]
    print(f" PERSISTENT gaps (forbidden all n): {persistent}  <- the impossibility certificates")

    print("\n"+"="*70); print(" (2) which SMALL Omega graphs are REALIZABLE as Omega(T)?"); print("="*70)
    # collect realized small-Omega iso classes (tournaments with few odd cycles => small Omega)
    realized={}  # m -> set of canon graphs
    for n in range(3,7):
        for adj in all_tournaments(n):
            m,E=conflict_graph(adj)
            if m<=5:  # only small Omega
                c=graph_canon(m,E)
                realized.setdefault(m,set()).add(c[1])
    for m in range(1,6):
        allg=set(c[1] for c in all_small_graphs(m)) if m<=5 else set()
        real=realized.get(m,set())
        forb=allg-real
        print(f" Omega on {m} vertices: {len(real)}/{len(allg)} realizable; "
              f"FORBIDDEN={len(forb)}")
        for fb in sorted(forb):
            print(f"     FORBIDDEN Omega: {graph_name(m,fb)}  bits={fb}")

    print("\n"+"="*70); print(" (3) alpha-vector (conflict-census) spectrum n<=6"); print("="*70)
    avset={}
    for n in range(3,7):
        s=set()
        for adj in all_tournaments(n):
            s.add(tuple(alpha_vector(adj)))
        avset[n]=s
        print(f" n={n}: #distinct alpha-vectors={len(s)}; "
              f"alpha1 values={sorted(set(a[1] if len(a)>1 else 0 for a in s))}")
    a1_all=set()
    for n in avset:
        for a in avset[n]: a1_all.add(a[1] if len(a)>1 else 0)
    a1gaps=[v for v in range(max(a1_all)+1) if v not in a1_all]
    print(f" alpha_1 (total odd cycles) gaps over n<=6: {a1gaps}  (note: alpha1=3 the H=7 driver)")
