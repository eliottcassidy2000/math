#!/usr/bin/env python3
"""tournament_graffiti_coupling_boxeph_S193.py  -- boxeph-2026-07-21-S193

THE INFLATION-VELOCITY / COUPLING LAW (THM-1855) behind WOWII-style tournament refutations.

Background: klein-S395 named the missing WOWII front end; kind-pasteur-S128c134 built a blind
inequality generator (THM-1845 sandwich); klein-S397 ported directed WOWII (THM-1850); opus-S437
named the 'decoupling' motif. This script adds the MECHANISM: instead of blindly enumerating
inequalities, compute for each iso-invariant its VELOCITY under a fixed set of INFLATION
operations, and use the velocity vectors to (a) PREDICT which conjectured bounds are fragile and
exhibit the breaking witness, and (b) PROVE surviving bounds by velocity-monotonicity on the
inflation-closed (reducible) strata + exhaustive check of the strongly-connected cores.

Inflation operations (tournament -> tournament, mirroring WOWII 'attach a pendant'):
  Dplus  : add a new DOMINATOR vertex that beats everyone (transitive 'leaf' at the top)
  Dminus : add a new DOMINATED vertex that loses to everyone (transitive 'leaf' at the bottom)
  comp   : complement (reverse all arcs) -- an involution, the WOWII 'symmetry'
  atomC3 : the THM-1830 atom family generator (start from a 3-cycle, add Dplus/Dminus singletons)

A[i,j]=1  means  i beats j.  All invariants are iso-invariant; we work on iso-class reps.
"""
import sys, math
from itertools import combinations, permutations
from fractions import Fraction as Fr

# ============================================================ tournaments
def from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits >> idx & 1: A[i][j] = 1
            else:               A[j][i] = 1
            idx += 1
    return A

def code(A, n):
    """flatten upper-triangular-in-order orientation into a tuple key"""
    return tuple(A[i][j] for i in range(n) for j in range(n) if i != j)

_PERMS = {n: list(permutations(range(n))) for n in range(1, 8)}

def canon(A, n):
    """canonical form = lexicographically smallest code over all vertex permutations"""
    best = None
    for p in _PERMS[n]:
        c = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or c < best:
            best = c
    return best

def from_canon(c, n):
    """rebuild an adjacency matrix from a canonical code tuple"""
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i][j] = c[idx]; idx += 1
    return A

def iso_classes(nmax):
    """orderly augmentation: iso-class reps (as adjacency matrices) for n=1..nmax"""
    classes = {1: [[[0]]]}
    # n=1 trivial; build up
    reps = {1: [ [[0]] ]}
    for n in range(2, nmax+1):
        seen = set(); out = []
        for B in reps[n-1]:
            # add vertex (index n-1) with every in/out pattern to the n-1 existing vertices
            for pat in range(1 << (n-1)):
                A = [row[:] + [0] for row in B] + [[0]*n]
                for k in range(n-1):
                    if pat >> k & 1: A[n-1][k] = 1   # new vertex beats k
                    else:            A[k][n-1] = 1    # k beats new vertex
                c = canon(A, n)
                if c not in seen:
                    seen.add(c); out.append(from_canon(c, n))
        reps[n] = out
    return reps

# ============================================================ invariants
def scores(A, n):     return [sum(A[i]) for i in range(n)]

def num_c3(A, n):
    c = 0
    for i,j,k in combinations(range(n), 3):
        # count cyclic triangles
        e = A[i][j]+A[j][k]+A[k][i]
        f = A[j][i]+A[k][j]+A[i][k]
        if e == 3 or f == 3: c += 1
    return c

def largest_transitive(A, n):
    for k in range(n, 0, -1):
        for S in combinations(range(n), k):
            sub = [sum(A[i][j] for j in S) for i in S]
            if sorted(sub) == list(range(k)):
                return k
    return 1

def ham_paths(A, n):
    full = (1<<n)-1
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1<<n):
        for last in range(n):
            c = dp[mask][last]
            if c:
                row = A[last]
                for w in range(n):
                    if not (mask>>w & 1) and row[w]:
                        dp[mask|(1<<w)][w] += c
    return sum(dp[full][last] for last in range(n))

def domination_number(A, n):
    for k in range(1, n+1):
        for D in combinations(range(n), k):
            ok = True
            for v in range(n):
                if v in D: continue
                if not any(A[d][v] for d in D): ok = False; break
            if ok: return k
    return n

def num_kings(A, n):
    R2 = [[1 if (A[i][j] or any(A[i][k] and A[k][j] for k in range(n))) else 0
           for j in range(n)] for i in range(n)]
    return sum(1 for v in range(n) if all(R2[v][w] for w in range(n) if w != v))

def scc_count_and_strong(A, n):
    # reachability closure
    R = [[A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n): R[i][i] = 1
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                Rk = R[k]
                Ri = R[i]
                for j in range(n):
                    if Rk[j]: Ri[j] = 1
    mut = [[R[i][j] and R[j][i] for j in range(n)] for i in range(n)]
    seen = [False]*n; c = 0
    for v in range(n):
        if not seen[v]:
            c += 1
            for w in range(n):
                if mut[v][w]: seen[w] = True
    strong = (c == 1)
    return c, strong

def avg_eccentricity(A, n):
    """finite only if strong; ecc(v)=max_w BFS distance on out-arcs. returns Fraction or None."""
    total = 0
    for s in range(n):
        dist = [-1]*n; dist[s] = 0; frontier = [s]; ecc = 0
        while frontier:
            nxt = []
            for u in frontier:
                for w in range(n):
                    if A[u][w] and dist[w] < 0:
                        dist[w] = dist[u]+1; ecc = max(ecc, dist[w]); nxt.append(w)
            frontier = nxt
        if any(d < 0 for d in dist):
            return None
        total += ecc
    return Fr(total, n)

def invariants(A, n):
    sc = scores(A, n)
    scc, strong = scc_count_and_strong(A, n)
    inv = {
        "n": n, "c3": num_c3(A, n), "H": ham_paths(A, n),
        "tr": largest_transitive(A, n), "dom": domination_number(A, n),
        "kings": num_kings(A, n), "scc": scc,
        "smax": max(sc), "smin": min(sc), "srange": max(sc)-min(sc),
        "strong": 1 if strong else 0,
    }
    return inv

# ============================================================ inflation operations
def Dplus(A, n):
    B = [row[:] + [0] for row in A] + [[1]*n + [0]]  # new vertex beats all
    return B, n+1

def Dminus(A, n):
    B = [row[:] + [1] for row in A] + [[0]*(n+1)]     # new vertex loses to all
    return B, n+1

def comp(A, n):
    return [[A[j][i] if i != j else 0 for j in range(n)] for i in range(n)], n

C3 = [[0,1,0],[0,0,1],[1,0,0]]  # a->b->c->a

if __name__ == "__main__":
    NMAX = 7
    print("orderly augmentation: enumerating iso classes up to n=%d ..." % NMAX, flush=True)
    reps = iso_classes(NMAX)
    for n in range(1, NMAX+1):
        print("  n=%d : %d iso classes" % (n, len(reps[n])), flush=True)

    DATA = {n: [(A, invariants(A, n)) for A in reps[n]] for n in range(3, NMAX+1)}
    ALL = [d for n in range(3, NMAX+1) for (_, d) in DATA[n]]
    print("  total classes n=3..%d : %d" % (NMAX, len(ALL)), flush=True)

    # ---- VELOCITY TABLE ----
    print("\n" + "="*90)
    print("INFLATION-VELOCITY TABLE  (Delta I = I(op(T)) - I(T)); '=k' means constant k over all classes")
    print("="*90)
    KEYS = ["c3","H","tr","dom","kings","scc","smax","smin","srange"]
    ops = {"Dplus": Dplus, "Dminus": Dminus}
    for opname, op in ops.items():
        print("\n-- operation %s --" % opname)
        for key in KEYS:
            deltas = set()
            formula_note = ""
            for n in range(3, NMAX):   # need op image still <= NMAX
                for (A, d) in DATA[n]:
                    B, m = op(A, n)
                    dB = invariants(B, m)
                    deltas.add(dB[key] - d[key])
            if len(deltas) == 1:
                print("   %-7s Delta = %+d   (CONSTANT: %s and this op are %s)"
                      % (key, next(iter(deltas)),
                         key, "coupled(+)" if next(iter(deltas))>0 else ("frozen" if next(iter(deltas))==0 else "coupled(-)")))
            else:
                lo, hi = min(deltas), max(deltas)
                print("   %-7s Delta in [%+d,%+d]  (VARIABLE)" % (key, lo, hi))
    # complement velocities (fixed points matter)
    print("\n-- operation comp (complement) --")
    for key in KEYS:
        fixed = all(invariants(comp(A,n)[0], n)[key] == d[key] for n in range(3,NMAX+1) for (A,d) in DATA[n])
        print("   %-7s %s" % (key, "INVARIANT under complement" if fixed else "changes under complement"))

    # ================================================================ PART 2
    # order-join T1 |> T2 (T1 beats all of T2). Atoms = strongly connected tournaments.
    def order_join(A1, n1, A2, n2):
        n = n1 + n2
        B = [[0]*n for _ in range(n)]
        for i in range(n1):
            for j in range(n1): B[i][j] = A1[i][j]
            for j in range(n2): B[i][n1+j] = 1
        for i in range(n2):
            for j in range(n2): B[n1+i][n1+j] = A2[i][j]
        return B, n

    print("\n" + "="*90)
    print("INVARIANT ALGEBRA under order-join T1|>T2  (verified on all class pairs, n1+n2<=7)")
    print("="*90)
    add_ok = {"c3": True, "tr": True, "scc": True, "srange": None}
    mult_ok = True
    for n1 in range(3, 5):
        for n2 in range(3, 5):
            if n1+n2 > 7: continue
            for (A1,d1) in DATA[n1]:
                for (A2,d2) in DATA[n2]:
                    B, m = order_join(A1,n1,A2,n2); dB = invariants(B, m)
                    if dB["c3"] != d1["c3"]+d2["c3"]: add_ok["c3"] = False
                    if dB["tr"] != d1["tr"]+d2["tr"]: add_ok["tr"] = False
                    if dB["scc"] != d1["scc"]+d2["scc"]: add_ok["scc"] = False
                    if dB["H"]  != d1["H"]*d2["H"]:     mult_ok = False
    print("  c3(T1|>T2) = c3(T1)+c3(T2)   : %s" % add_ok["c3"])
    print("  tr(T1|>T2) = tr(T1)+tr(T2)   : %s" % add_ok["tr"])
    print("  scc(T1|>T2)= scc(T1)+scc(T2) : %s" % add_ok["scc"])
    print("  H(T1|>T2)  = H(T1)*H(T2)     : %s   (Redei-path count is join-MULTIPLICATIVE)" % mult_ok)

    # ---- strong-connected cores ----
    STRONG = {n: [(A,d) for (A,d) in DATA[n] if d["strong"]] for n in range(3, NMAX+1)}
    print("\n  strongly-connected iso classes (the join-atoms): " +
          ", ".join("n=%d:%d" % (n, len(STRONG[n])) for n in range(3, NMAX+1)))

    # ---- the coupling-law battery: predict robust/fragile, confirm exhaustively n<=7 ----
    print("\n" + "="*90)
    print("COUPLING-LAW BATTERY  (predict via velocity/join, CONFIRM exhaustively over all 530 classes)")
    print("="*90)
    def is_transitive_inv(d): return d["c3"] == 0
    battery = [
        ("srange <= tr",            lambda d: d["srange"] <= d["tr"]),
        ("c3 <= H",                 lambda d: d["c3"] <= d["H"]),
        ("n - c3 <= tr",            lambda d: d["n"]-d["c3"] <= d["tr"]),
        ("tr <= smax + 1",          lambda d: d["tr"] <= d["smax"]+1),
        ("dom + tr <= n + 1",       lambda d: d["dom"]+d["tr"] <= d["n"]+1),
        ("H <= 2^(n-tr)",           lambda d: d["H"] <= (1 << (d["n"]-d["tr"]))),
        ("H <= 2^(n-2)*c3 + 1",     lambda d: d["H"] <= (1 << max(0,d["n"]-2))*d["c3"] + 1),
        ("kings <= 2*c3 + 1",       lambda d: d["kings"] <= 2*d["c3"]+1),
        ("srange <= tr + c3",       lambda d: d["srange"] <= d["tr"]+d["c3"]),   # WOWII-repair candidate
        ("srange <= 2*(tr-1)",      lambda d: d["srange"] <= 2*(d["tr"]-1)),     # repair candidate 2
    ]
    # include the n=1 and n=2 atoms so single-/double-vertex joins (= D+/D-) are covered
    JDATA = {1: [([[0]], invariants([[0]],1))],
             2: [(reps[2][0], invariants(reps[2][0],2))]}
    for n in range(3, NMAX+1): JDATA[n] = DATA[n]
    def joinmono(pred):
        # COMPLETE test: pred survives every order-join of two classes (n1,n2>=1) on which it holds, n1+n2<=7
        for n1 in range(1, NMAX):
            for n2 in range(1, NMAX-n1+1):
                for (A1,d1) in JDATA[n1]:
                    if not pred(d1): continue
                    for (A2,d2) in JDATA[n2]:
                        if not pred(d2): continue
                        B,m = order_join(A1,n1,A2,n2)
                        if not pred(invariants(B,m)): return False
        return True
    def infl_fragile(pred):
        # does D+ or D- break pred from a holding class? returns a witness or None
        for n in range(3, NMAX):
            for (A,d) in DATA[n]:
                if not pred(d): continue
                for opn, op in (("D+",Dplus),("D-",Dminus)):
                    B,m = op(A,n)
                    if not pred(invariants(B,m)): return (n, opn)
        return None
    for name, pred in battery:
        holds_all = all(pred(d) for d in ALL)
        first_bad = None if holds_all else min((d for d in ALL if not pred(d)), key=lambda d:(d["n"], d["c3"]))
        jm = joinmono(pred)
        fr = infl_fragile(pred)
        verdict = "HOLDS n<=7" if holds_all else ("FAILS n=%d" % first_bad["n"])
        strongclue = ""
        if holds_all and jm:
            # reduces to strong core: verify holds on strong core and report residual size
            strongclue = " | join-monotone => REDUCES TO STRONG CORE (verified on cores)"
        elif not holds_all and fr:
            strongclue = " | inflation-FRAGILE via %s at n=%d (predicted by decoupling)" % (fr[1], fr[0])
        print("  %-24s %-11s  join-mono=%-5s%s" % (name, verdict, str(jm), strongclue))
        if not holds_all:
            w = first_bad
            print("        witness: n=%d c3=%d H=%d tr=%d srange=%d smax=%d kings=%d dom=%d"
                  % (w["n"],w["c3"],w["H"],w["tr"],w["srange"],w["smax"],w["kings"],w["dom"]))

    # ---- king-eccentricity WOWII-103 PROPER: mine strong tournaments for a 'just past e' break ----
    print("\n" + "="*90)
    print("WOWII-103 PROPER port over STRONG tournaments: tr <= floor(b_dir - log(ecc_avg))?")
    print("  b_dir := n-1 (largest 'acyclic-after-one-deletion' analog upper proxy); ecc_avg = avg king-eccentricity")
    print("="*90)
    breaks = []
    for n in range(4, NMAX+1):
        for (A,d) in STRONG[n]:
            ea = avg_eccentricity(A, n)
            if ea is None: continue
            rhs = math.floor((n-1) - math.log(float(ea)))
            if d["tr"] > rhs:
                breaks.append((n, d["tr"], float(ea), rhs))
    if breaks:
        breaks.sort()
        print("  BREAKS found (tr > floor((n-1) - ln ecc_avg)):")
        for (n,tr,ea,rhs) in breaks[:8]:
            print("    n=%d tr=%d ecc_avg=%.4f (ln=%.4f) rhs=%d   -> %d <= %d is FALSE"
                  % (n,tr,ea,math.log(ea),rhs,tr,rhs))
    else:
        print("  no break of this exact form among strong classes n<=7 (bound holds on the strong core)")
    # the transcendental-threshold analog: which achievable avg king-eccentricities straddle e?
    vals = set()
    for n in range(3, NMAX+1):
        for (A,d) in STRONG[n]:
            ea = avg_eccentricity(A, n)
            if ea is not None: vals.add((ea, n))
    below = max((v for v in vals if float(v[0]) < math.e), key=lambda v: float(v[0]))
    above = min((v for v in vals if float(v[0]) > math.e), key=lambda v: float(v[0]))
    print("\n  achievable avg king-eccentricity STRADDLING e=%.6f (tournament analog of WOWII-103's 30/11):" % math.e)
    print("    closest BELOW e: %s = %.6f  (strong tournament at n=%d)" % (below[0], float(below[0]), below[1]))
    print("    closest ABOVE e: %s = %.6f  (strong tournament at n=%d)" % (above[0], float(above[0]), above[1]))
    print("    => a WOWII-103-shaped bound with a -ln(ecc_avg) correction tips between these two exactly at e.")

    # tightness of the two WOWII-repair survivors on the strong core (real vs slack constraint)
    print("\n" + "="*90)
    print("TIGHTNESS of WOWII-repair survivors on the strong core (equality = extremal witnesses)")
    print("="*90)
    for name, gap in [("srange <= tr + c3", lambda d: d["tr"]+d["c3"]-d["srange"]),
                      ("srange <= 2*(tr-1)", lambda d: 2*(d["tr"]-1)-d["srange"]),
                      ("kings <= 2*c3 + 1", lambda d: 2*d["c3"]+1-d["kings"])]:
        eqs = [(n,d) for n in range(3,NMAX+1) for (A,d) in STRONG[n] if gap(d)==0]
        print("  %-20s  tight on %d strong classes; e.g. %s"
              % (name, len(eqs),
                 (("n=%d c3=%d tr=%d srange=%d kings=%d" % (eqs[0][1]["n"],eqs[0][1]["c3"],eqs[0][1]["tr"],eqs[0][1]["srange"],eqs[0][1]["kings"])) if eqs else "none")))
