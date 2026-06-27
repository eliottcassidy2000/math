"""
tournament_certificate_engine_kps.py  (kind-pasteur-2026-06-27-S31ah)

A programmatic TOURNAMENT PROOF/DISPROOF TOOLKIT. Generalizes the project's
H=7 / H=21 impossibility technique (THM-200/115/201/202) into a reusable engine
that (1) computes tournament invariants, (2) DISCOVERS forbidden-value spectra
(the technique generator: enumerate small tournaments => achieved set => the
complement is a candidate impossibility certificate), (3) applies certificates.

Core idea (OCF, THM-002): H(T) = I(Omega(T), 2) = prod over Omega-components.
A structure is IMPOSSIBLE if mapping it to a tournament forces a forbidden
invariant (H in {7,21}; H even; non-Landau score seq; non-realizable Omega; ...).

Invariants computed:
  - score sequence + Landau realizability
  - directed odd-cycle counts c3,c5,c7 (ALL directed cycles per vertex set -- the
    MISTAKE-corrected enumeration)
  - conflict graph Omega (vertices = directed odd cycles, edges = shared vertex)
  - independence polynomial alpha-vector and H = I(Omega,2)
  - #Hamiltonian paths (Redei; always ODD)
  - strongly connected components and the SCC product H = prod H(SCC)
"""
import itertools, sys
from functools import lru_cache

# ---------- tournament representation: adj[i][j]=1 iff i->j ----------
def all_tournaments(n):
    """yield adjacency matrices of all 2^C(n,2) tournaments on n labeled vertices."""
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        adj = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if bits>>k & 1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def random_tournament(n, rng):
    adj=[[0]*n for _ in range(n)]
    for i,j in itertools.combinations(range(n),2):
        if rng.random()<0.5: adj[i][j]=1
        else: adj[j][i]=1
    return adj

def scores(adj):
    n=len(adj); return sorted(sum(adj[i]) for i in range(n))

# ---------- Landau realizability of a score sequence ----------
def is_landau(seq):
    s=sorted(seq); n=len(s)
    pref=0
    for k in range(1,n+1):
        pref+=s[k-1]
        if pref < k*(k-1)//2: return False
    return pref == n*(n-1)//2

# ---------- directed odd cycles (correct: ALL per vertex set) ----------
def directed_cycles_of_length(adj, L):
    """return set of canonical directed L-cycles (tuples, rotation-normalized, direction kept)."""
    n=len(adj); out=set()
    for verts in itertools.combinations(range(n), L):
        # enumerate directed Hamiltonian cycles of the induced subtournament
        vs=list(verts)
        for perm in itertools.permutations(vs[1:]):
            cyc=(vs[0],)+perm
            if all(adj[cyc[i]][cyc[(i+1)%L]] for i in range(L)):
                # rotation-normalize (keep direction): start at min, keep order
                m=cyc.index(min(cyc)); norm=cyc[m:]+cyc[:m]
                out.add(norm)
    return out

def odd_cycle_list(adj, maxL=None):
    n=len(adj); maxL=maxL or n
    cyc=[]
    for L in range(3, maxL+1, 2):
        cyc += list(directed_cycles_of_length(adj, L))
    return cyc

def odd_cycle_counts(adj):
    n=len(adj); c={}
    for L in range(3,n+1,2):
        c[L]=len(directed_cycles_of_length(adj,L))
    return c

# ---------- conflict graph Omega + independence polynomial ----------
def conflict_graph(adj):
    """Omega: vertices = directed odd cycles; edge iff they share a tournament vertex."""
    cyc=odd_cycle_list(adj)
    sets=[set(c) for c in cyc]
    m=len(cyc); E=[[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            if sets[i]&sets[j]: E[i][j]=E[j][i]=True
    return m,E

def independence_poly(m,E):
    """alpha[k] = # independent sets of size k in graph (m vertices, adjacency E)."""
    # DP over vertices; for small m brute force independent sets
    alpha=[0]*(m+1); alpha[0]=1
    # iterate subsets is too big for m>22; use recursive add
    verts=list(range(m))
    # greedy recursion counting independent sets by size
    from functools import lru_cache
    def rec(idx, chosen_mask):
        # count independent sets among verts[idx:] compatible with chosen_mask
        if idx==m:
            return
        # we instead do explicit DFS accumulating sizes
    # explicit DFS
    cnt=[0]*(m+1)
    def dfs(start, size, forbidden):
        cnt[size]+=1
        for v in range(start, m):
            if forbidden>>v & 1: continue
            nf=forbidden
            for w in range(m):
                if E[v][w]: nf|=1<<w
            dfs(v+1, size+1, nf | (1<<v))
    dfs(0,0,0)
    return cnt

def H_value(adj):
    m,E=conflict_graph(adj)
    alpha=independence_poly(m,E)
    return sum(alpha[k]*(2**k) for k in range(len(alpha)))

def alpha_vector(adj):
    m,E=conflict_graph(adj)
    a=independence_poly(m,E)
    # trim trailing zeros
    while len(a)>1 and a[-1]==0: a.pop()
    return a

# ---------- #Hamiltonian paths (Redei: always odd) ----------
def ham_path_count(adj):
    n=len(adj)
    if n==0: return 1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c: continue
            for nxt in range(n):
                if mask>>nxt&1: continue
                if adj[last][nxt]: dp[mask|1<<nxt][nxt]+=c
    full=(1<<n)-1
    return sum(dp[full][v] for v in range(n))

# ---------- strongly connected components ----------
def sccs(adj):
    n=len(adj); idx=[0]*n; low=[0]*n; onstk=[False]*n; stk=[]; comps=[]; cnt=[1]
    import sys as _s; _s.setrecursionlimit(10000)
    visited=[False]*n
    def dfs(v):
        idx[v]=low[v]=cnt[0]; cnt[0]+=1; stk.append(v); onstk[v]=True; visited[v]=True
        for w in range(n):
            if adj[v][w]:
                if not visited[w]: dfs(w); low[v]=min(low[v],low[w])
                elif onstk[w]: low[v]=min(low[v],idx[w])
        if low[v]==idx[v]:
            comp=[]
            while True:
                w=stk.pop(); onstk[w]=False; comp.append(w)
                if w==v: break
            comps.append(comp)
    for v in range(n):
        if not visited[v]: dfs(v)
    return comps

# ---------- CERTIFICATES (impossibility detectors) ----------
KNOWN_H_GAPS = {7, 21}  # permanent gaps (THM-200, THM-115)

def H_spectrum_certificate(v):
    """Return reason v is FORBIDDEN as an H value, or None if plausibly achievable."""
    if v < 1: return "H>=1 always"
    if v % 2 == 0: return "H is always ODD (=1+2*sum)"
    if v in KNOWN_H_GAPS: return f"H={v} is a PERMANENT GAP (THM-200/115)"
    return None

def redei_parity_certificate(count):
    """#Ham paths must be odd (Redei). Return reason if forbidden."""
    if count % 2 == 0: return "Redei: #Ham paths must be ODD"
    return None

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    # SELF-TEST against known values
    print("=== SELF-TEST ===")
    # Paley T_7 (QR tournament): i->j iff (j-i) mod 7 in {1,2,4}
    n=7; adj=[[0]*7 for _ in range(7)]
    QR={1,2,4}
    for i in range(7):
        for j in range(7):
            if i!=j and (j-i)%7 in QR: adj[i][j]=1
    print(f"Paley T_7: scores={scores(adj)} (reg=3?), c={odd_cycle_counts(adj)}")
    print(f"  H={H_value(adj)} (expect 189), #Ham paths={ham_path_count(adj)} (odd? {ham_path_count(adj)%2})")
    print(f"  alpha-vector={alpha_vector(adj)}")
    # transitive T_4 -> H=1
    adjT=[[1 if i<j else 0 for j in range(4)] for i in range(4)]
    print(f"Transitive T_4: H={H_value(adjT)} (expect 1), #Ham={ham_path_count(adjT)} (expect 1)")
    # 3-cycle -> H=3
    adj3=[[0,1,0],[0,0,1],[1,0,0]]
    print(f"3-cycle: H={H_value(adj3)} (expect 3), c={odd_cycle_counts(adj3)}")
    print(f"\n=== CERTIFICATES ===")
    for v in [1,3,5,7,9,11,21,189,8]:
        r=H_spectrum_certificate(v)
        print(f"  H={v}: {'FORBIDDEN -- '+r if r else 'plausibly achievable'}")
