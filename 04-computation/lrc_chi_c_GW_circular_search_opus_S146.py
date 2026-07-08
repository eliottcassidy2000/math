"""
lrc_chi_c_GW_circular_search_opus_S146.py  (opus-2026-07-07-S146, subagent of HYP-5197)

CIRCULAR CHROMATIC NUMBER of  G_GW = Cay(Z, +-GW),  GW = {1..11, 13, 24}.

KNOWN (THM-652 / S144-S145):
    chi_f(G_GW) = 1/mu = 13   (mu = 1/13, rigid optimum {0,12} mod 26)
    chi(G_GW)   = 14 = 1/M    (integer rung faithful; no 13-coloring: odd-cycle matching obstruction)
    ==> chi_c(G_GW) in (13, 14].    THIS SCRIPT probes where in (13,14] it lies.

A (p,q)-coloring: c: V -> Z_p with circular distance |c(x)-c(y)|_p >= q on every edge
   (|a|_p := min(a mod p, -a mod p)).  chi_c = inf{p/q : a (p,q)-coloring exists}.
   Circular-distance >= q  <=>  q <= |c(x)-c(y)| <= p-q  (ordinary abs diff), since p >= 2q.

TWO RIGOROUS SIDES (both machine-checkable):
  * UPPER bound chi_c <= p/q : exhibit + verify a PERIODIC (p,q)-coloring c: Z_N -> Z_p
        (valid on the infinite graph iff no s in GW is 0 mod N -- no self loop; then the
         circulant Cay(Z_N, +-(GW mod N)) constraints are EXACTLY the infinite-graph
         constraints under the periodic lift.  No "N>2max" needed; only no-self-loop.)
  * LOWER bound chi_c >  p/q : exhibit a FINITE segment [0..L] (induced subgraph of G_GW)
        that is NOT (p,q)-colorable (SAT proves UNSAT rigorously) -> G has no (p,q)-coloring.

KEY STRUCTURAL FACT (why this is subtle): a LINEAR coloring c(i)=a*i mod p needs
   min_{s in GW} ||a s / p|| >= q/p.   sup_t min_s ||s t|| = M(GW) = 1/14 (tight LRC(14)
   instance), so LINEAR colorings achieve p/q = 14 but NEVER below.  Any chi_c < 14 must
   come from a NON-LINEAR periodic coloring.  (Section A verifies the linear wall.)

Engines: (L) linear arithmetic scan; (S) exact SAT (pysat/cadical153) for circulant &
segment feasibility incl. UNSAT proofs; (H) min-conflicts local search (fast finder).

Usage:  python lrc_chi_c_GW_circular_search_opus_S146.py [mode]
  modes: linear controls circulant segment ladder all   (default: all)
"""
import sys, time, math, random, itertools
from fractions import Fraction as F
from threading import Timer

GW = list(range(1, 12)) + [13, 24]          # {1,...,11,13,24}
GWSET = set(GW)
WMAX = max(GW)                               # 24

# ----------------------------------------------------------------- basic helpers
def cdist(a, b, p):
    """circular distance |a-b|_p on Z_p."""
    d = (a - b) % p
    return min(d, p - d)

def gcd(a, b):
    return math.gcd(a, b)

# ----------------------------------------------------------------- verifiers (independent)
def verify_periodic(coloring, p, q, N, S=GW):
    """Rigorously verify c:Z_N->Z_p is a proper (p,q)-coloring of the infinite graph.
       Returns (ok, reason)."""
    if len(coloring) != N:
        return False, f"len {len(coloring)} != N {N}"
    if any((c < 0 or c >= p) for c in coloring):
        return False, "color out of range"
    for s in S:
        if s % N == 0:
            return False, f"self-loop: s={s} == 0 mod N={N}"
    for r in range(N):
        for s in S:
            d = cdist(coloring[r], coloring[(r + s) % N], p)
            if d < q:
                return False, f"edge r={r} s={s}: cdist={d} < q={q}"
    return True, "ok"

def verify_segment(coloring, p, q, L, S=GW):
    """Verify c on vertices 0..L (segment) is a proper (p,q)-coloring. coloring len L+1."""
    for i in range(L + 1):
        for s in S:
            j = i + s
            if j <= L:
                if cdist(coloring[i], coloring[j], p) < q:
                    return False, f"edge {i}-{j}: cdist {cdist(coloring[i],coloring[j],p)} < {q}"
    return True, "ok"

# ================================================================= (A) LINEAR scan
def linear_scan(pmax=3000, verbose=True):
    """For c(i)=a*i mod p: feasible (p,q) with q=min_s |a s|_p.  Best ratio p/q over all
       a,p<=pmax.  Confirms the linear wall = 14 (no linear coloring beats 14)."""
    import numpy as np
    Sarr = np.array(GW, dtype=np.int64)
    best_ratio = None; best_tuple = None
    # also: for each q, smallest p admitting a linear (p,q)-coloring
    best_p_for_q = {}
    for p in range(14, pmax + 1):
        a = np.arange(1, p, dtype=np.int64)[:, None]           # (p-1,1)
        prod = (a * Sarr[None, :]) % p                          # (p-1,13)
        cd = np.minimum(prod, p - prod)                         # circular dist
        m = cd.min(axis=1)                                      # (p-1,) best q per a
        qmax = int(m.max())
        if qmax >= 1:
            ai = int(m.argmax()); q = qmax; aval = int(a[ai, 0])
            r = F(p, q)
            if best_ratio is None or r < best_ratio:
                best_ratio = r; best_tuple = (p, q, aval)
            for qq in range(1, qmax + 1):
                if qq not in best_p_for_q or p < best_p_for_q[qq][0]:
                    # smallest p with some a achieving min_s|a s|_p >= qq
                    idxs = np.nonzero(m >= qq)[0]
                    if len(idxs):
                        best_p_for_q[qq] = (p, int(a[idxs[0], 0]))
    if verbose:
        print("=== (A) LINEAR coloring scan  c(i)=a*i mod p,  p<=%d ===" % pmax)
        print(f"  best linear ratio p/q = {best_tuple[0]}/{best_tuple[1]} = "
              f"{best_tuple[0]/best_tuple[1]:.6f}  (a={best_tuple[2]})")
        print("  smallest p admitting a linear (p,q)-coloring, per q:")
        for q in sorted(best_p_for_q):
            p, a = best_p_for_q[q]
            print(f"    q={q:2d}: p={p:4d}  ratio={p/q:.5f}  (a={a})  [={F(p,q)}]")
        print("  => LINEAR WALL: best linear ratio is %s. Linear colorings cannot beat"
              % best_ratio)
        print("     14 because sup_t min_s||s t|| = M(GW) = 1/14. Sub-14 needs NON-LINEAR.")
    return best_ratio, best_tuple, best_p_for_q

# ================================================================= SAT machinery
_SOLVER = "cadical153"

def _build_cnf(nvert, edges, p, q):
    """CNF for: proper (p,q)-coloring of graph (nvert vertices, edges list of (u,v)).
       Returns (clauses, var(v,c) function via dict, nvars)."""
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    pool = IDPool()
    def V(v, c):
        return pool.id(('x', v, c))
    # reserve all color vars first
    for v in range(nvert):
        for c in range(p):
            V(v, c)
    clauses = []
    # exactly-one color per vertex
    for v in range(nvert):
        lits = [V(v, c) for c in range(p)]
        clauses.append(list(lits))                              # at-least-one
        amo = CardEnc.atmost(lits, bound=1, vpool=pool, encoding=EncType.seqcounter)
        clauses.extend(amo.clauses)                             # at-most-one
    # edge constraints: forbid |c(u)-c(v)|_p < q
    win = range(-(q - 1), q)                                    # deltas with cdist<q
    for (u, v) in edges:
        for a in range(p):
            for dl in win:
                b = (a + dl) % p
                clauses.append([-V(u, a), -V(v, b)])
    return clauses, V, pool.top

def sat_feasible(nvert, edges, p, q, budget=3_000_000, solver=_SOLVER):
    """Return (status, coloring|None). status in {'SAT','UNSAT','UNKNOWN'}.
       Deterministic conflict-budget limiting (the Timer/interrupt path does NOT fire
       reliably on this platform -- verified; conf_budget does).  UNKNOWN = budget hit."""
    from pysat.solvers import Solver
    clauses, V, _ = _build_cnf(nvert, edges, p, q)
    s = Solver(name=solver, bootstrap_with=clauses)
    if budget:
        s.conf_budget(int(budget))
        res = s.solve_limited()
    else:
        res = s.solve()
    if res is None:
        s.delete(); return 'UNKNOWN', None
    if not res:
        s.delete(); return 'UNSAT', None
    model = set(l for l in s.get_model() if l > 0)
    coloring = [None] * nvert
    for v in range(nvert):
        for c in range(p):
            if V(v, c) in model:
                coloring[v] = c; break
    s.delete()
    return 'SAT', coloring

def circulant_edges(N, S=GW):
    """edges of Cay(Z_N, +-(S mod N)); None if a self-loop occurs (infeasible N)."""
    for s in S:
        if s % N == 0:
            return None
    E = set()
    for r in range(N):
        for s in S:
            u, v = r, (r + s) % N
            if u != v:
                E.add((min(u, v), max(u, v)))
    return list(E)

def segment_edges(L, S=GW):
    E = []
    for i in range(L + 1):
        for s in S:
            if i + s <= L:
                E.append((i, i + s))
    return E

def sat_circulant(p, q, N, budget=3_000_000):
    E = circulant_edges(N)
    if E is None:
        return 'SELFLOOP', None
    st, col = sat_feasible(N, E, p, q, budget)
    if st == 'SAT':
        ok, why = verify_periodic(col, p, q, N)
        if not ok:
            return 'BUG(' + why + ')', col
    return st, col

def sat_segment(p, q, L, budget=3_000_000):
    E = segment_edges(L)
    st, col = sat_feasible(L + 1, E, p, q, budget)
    if st == 'SAT':
        ok, why = verify_segment(col, p, q, L)
        if not ok:
            return 'BUG(' + why + ')', col
    return st, col

# ================================================================= (H) local search finder
def local_search_periodic(p, q, N, restarts=40, iters=4000, seed=0, S=GW):
    """min-conflicts finder for a periodic (p,q)-coloring on Z_N. Returns coloring|None."""
    E = circulant_edges(N)
    if E is None:
        return None
    adj = [[] for _ in range(N)]
    for (u, v) in E:
        adj[u].append(v); adj[v].append(u)
    rng = random.Random(seed * 1000003 + p * 131 + q * 17 + N)
    def confs(col):
        c = 0
        for (u, v) in E:
            if cdist(col[u], col[v], p) < q:
                c += 1
        return c
    for _ in range(restarts):
        col = [rng.randrange(p) for _ in range(N)]
        for _ in range(iters):
            # find a conflicted vertex
            bad = [u for u in range(N) if any(cdist(col[u], col[w], p) < q for w in adj[u])]
            if not bad:
                ok, _why = verify_periodic(col, p, q, N)
                return col if ok else None
            u = rng.choice(bad)
            # choose color for u minimizing conflicts (min-conflicts, random tie-break)
            best = []; bestc = None
            for c in range(p):
                cc = sum(1 for w in adj[u] if cdist(c, col[w], p) < q)
                if bestc is None or cc < bestc:
                    bestc = cc; best = [c]
                elif cc == bestc:
                    best.append(c)
            col[u] = rng.choice(best)
        if confs(col) == 0:
            ok, _ = verify_periodic(col, p, q, N)
            if ok:
                return col
    return None

# ================================================================= experiments
RATIO_LADDER = []   # list of (p,q) with 13 < p/q < 14, reduced, sorted by value
def build_ladder(qmax=9):
    seen = set(); out = []
    for q in range(1, qmax + 1):
        for p in range(13 * q + 1, 14 * q):        # strictly between 13 and 14
            if gcd(p, q) != 1:
                continue
            r = F(p, q)
            if r in seen:
                continue
            seen.add(r); out.append((p, q))
    out.sort(key=lambda pq: pq[0] / pq[1])
    return out

def exp_controls(budget=8_000_000):
    print("\n=== (B) CONTROLS (sanity) ===")
    # (14,1) MUST succeed:  c(i)=i mod 14
    col = [i % 14 for i in range(14)]
    ok, why = verify_periodic(col, 14, 1, 14)
    print(f"  (14,1) explicit c(i)=i mod 14 on N=14: verify={ok} ({why}) -> chi_c <= 14  [expected OK]")
    # (13,1): chi=14 => no 13-coloring of the INFINITE graph. The density bound never
    # forbids a FINITE segment's 13-coloring, so segments stay SAT until (perhaps huge) L.
    for L in (30, 60, 120):
        st, _ = sat_segment(13, 1, L, budget)
        print(f"  (13,1) segment [0..{L}] : {st}"
              + ("  (13-colorable finite window; chi=14 is asymptotic, THM-652)"
                 if st == 'SAT' else ""))
        if st == 'UNSAT':
            print(f"    -> finite segment [0..{L}] is not 13-colorable => chi(G_GW)>=14 (matches THM-652)")
            break

def lcm(a, b):
    return a * b // math.gcd(a, b)

def _periods_for(p, q):
    """Natural period candidates: multiples of 26 (rigid optimum period), p, and small mults."""
    base = [26 * k for k in range(1, 9)]                 # 26..208
    base += [p, 2 * p, 3 * p, 13 * q, 26 * q]
    Ns = sorted({N for N in base if 25 <= N <= 320 and circulant_edges(N) is not None})
    return Ns

def exp_sat_ladder(qmax=8, budget=4_000_000):
    """DECISIVE per-period SAT over the whole ratio ladder (13<p/q<14), ascending.
       For each ratio, SAT is complete for each period N: SAT => periodic coloring exists
       (upper bound); UNSAT for all tested N => no coloring of those periods.  First SAT
       found = best rigorous upper bound.  Local search pre-probes larger N cheaply."""
    print("\n=== (C) SAT circulant feasibility over the ratio ladder (13 < p/q < 14) ===")
    ladder = build_ladder(qmax)
    print("  ladder:", [f"{p}/{q}={p/q:.3f}" for p, q in ladder])
    best_upper = None
    table = {}
    for (p, q) in ladder:
        r = F(p, q)
        # cheap N for large q to keep SAT tractable; full N for small q
        Ns = _periods_for(p, q)
        if q >= 6:
            Ns = [N for N in Ns if N <= 156]
        found = None; verds = []
        for N in Ns:
            t0 = time.time()
            st, col = sat_circulant(p, q, N, budget)
            dt = time.time() - t0
            verds.append((N, st, round(dt, 1)))
            if st == 'SAT':
                found = (N, col); break
            if st.startswith('BUG'):
                found = ('BUG', st); break
        table[(p, q)] = (found, verds)
        if found and found[0] != 'BUG':
            N, col = found
            ok, why = verify_periodic(col, p, q, N)
            print(f"  {p}/{q}={p/q:.4f}: *** SAT period N={N}, verify={ok} -> chi_c <= {r} ***")
            if ok and (best_upper is None or r < best_upper[0]):
                best_upper = (r, p, q, N, col)
        elif found and found[0] == 'BUG':
            print(f"  {p}/{q}={p/q:.4f}: BUG {found[1]}  verds={verds}")
        else:
            allun = all(s == 'UNSAT' for _, s, _ in verds)
            summ = " ".join(f"N{N}:{s}" for N, s, _ in verds)
            print(f"  {p}/{q}={p/q:.4f}: {'UNSAT(all tested N)' if allun else 'not found'}  [{summ}]")
    if best_upper:
        print(f"  >>> best SAT upper bound: chi_c <= {best_upper[0]} = "
              f"{best_upper[1]/best_upper[2]:.4f}  (period N={best_upper[3]})")
    else:
        print("  >>> NO periodic (p,q)-coloring with p/q<14 found on any tested period.")
    return best_upper, table

def exp_localsearch_deep(qmax=8):
    """Deep min-conflicts over larger periods to catch any coloring SAT's period set missed."""
    print("\n=== (D) DEEP local-search over larger periods (finder only) ===")
    ladder = build_ladder(qmax)
    best = None
    for (p, q) in ladder:
        hit = None
        Ns = sorted({26 * k for k in range(1, 13)} | {p, 2 * p, 3 * p, 4 * p})
        Ns = [N for N in Ns if 26 <= N <= 340 and circulant_edges(N) is not None]
        for N in Ns:
            col = local_search_periodic(p, q, N, restarts=14, iters=1600, seed=7)
            if col is not None:
                hit = (N, col); break
        if hit:
            ok, _ = verify_periodic(hit[1], p, q, hit[0])
            print(f"  {p}/{q}={p/q:.4f}: LS FOUND N={hit[0]} verify={ok}")
            if ok and (best is None or F(p, q) < best[0]):
                best = (F(p, q), p, q, hit[0], hit[1])
        else:
            print(f"  {p}/{q}={p/q:.4f}: LS none")
    if best:
        print(f"  best LS upper bound: chi_c <= {best[0]} (period {best[3]})")
    else:
        print("  LS found NO periodic coloring with p/q<14.")
    return best

def chi_c_segment(L, ladder, budget=4_000_000):
    """Bracket chi_c([0..L]) = inf feasible ratio, by binary search over the ladder.
       Feasibility is MONOTONE in the ratio (K_{p/q} -> K_{p'/q'} iff p/q<=p'/q'),
       so there is a clean SAT/UNSAT threshold.  Returns (lo_unsat, hi_sat, cache)."""
    lo, hi = 0, len(ladder) - 1      # index into ascending ladder
    best_unsat = None; best_sat = None; cache = {}
    while lo <= hi:
        mid = (lo + hi) // 2
        p, q = ladder[mid]
        st, _ = sat_segment(p, q, L, budget)
        cache[(p, q)] = st
        if st == 'SAT':
            best_sat = (p, q); hi = mid - 1
        elif st == 'UNSAT':
            best_unsat = (p, q); lo = mid + 1
        else:  # UNKNOWN: treat as barrier, stop refining downward
            break
    return best_unsat, best_sat, cache

def exp_segment_trajectory(qmax=8, budget=4_000_000):
    """Rigorous lower bounds via chi_c([0..L]) trajectory (increases to chi_c(G))."""
    print("\n=== (E) chi_c([0..L]) trajectory -> rigorous lower bounds ===")
    ladder = build_ladder(qmax)          # ascending, includes 14 boundary? no (strictly<14)
    # append (14,1) as top so a fully-14 segment is representable
    ladder_ext = ladder + [(14, 1)]
    best_lb = None
    for L in (26, 40, 52, 78, 100, 130, 160, 200, 240):
        bu, bs, _ = chi_c_segment(L, ladder_ext, budget)
        lo_s = f"{bu[0]}/{bu[1]}={bu[0]/bu[1]:.4f}" if bu else "none<"
        hi_s = f"{bs[0]}/{bs[1]}={bs[0]/bs[1]:.4f}" if bs else ">=14"
        print(f"  L={L:3d}: chi_c([0..L]) in ( {lo_s} , {hi_s} ]   "
              f"(largest UNSAT ratio -> chi_c(G) > {bu[0]/bu[1]:.4f})" if bu else
              f"  L={L:3d}: chi_c([0..L]) <= {hi_s} (no sub-14 obstruction yet)")
        if bu and (best_lb is None or F(bu[0], bu[1]) > best_lb):
            best_lb = F(bu[0], bu[1])
    if best_lb is not None:
        print(f"  >>> best rigorous lower bound: chi_c(G_GW) > {best_lb} = {float(best_lb):.5f}")
    else:
        print("  >>> no finite-segment obstruction below 14 within tested L.")
    return best_lb

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    print("#" * 78)
    print("# chi_c(G_GW) circular chromatic search  --  GW =", GW)
    print("# chi_f=13, chi=14 (THM-652); locating chi_c in (13,14].")
    print("#" * 78)
    if mode in ("linear", "all"):
        linear_scan(pmax=2000)
    if mode in ("controls", "all"):
        exp_controls()
    if mode in ("circulant", "sat", "all"):
        exp_sat_ladder()
    if mode in ("ladder", "ls", "all"):
        exp_localsearch_deep()
    if mode in ("segment", "all"):
        exp_segment_trajectory()
    print(f"\n[total {time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
