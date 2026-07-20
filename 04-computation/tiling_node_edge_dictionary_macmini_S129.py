#!/usr/bin/env python3
"""
The tiling <-> node <-> edge dictionary, and the death of linear invariants
                                                        (mac-mini-2026-07-20-S129)
===============================================================================
Owner asked four things:
  (1) how do tiling SETS map to iso-class NODES in the merged metagraph, EXACTLY;
  (2) compute tilings / nodes / edges from each other efficiently -- look for tricks;
  (3) a descending star-type invariant has to come from a BASE-PATH-INDEPENDENT
      subgroup, the natural candidate being  \\cap Gamma_P  over all spanning paths;
  (4) (implicitly, from THM-1405) why is the star quotient transverse to isomorphism?

The punchline of (3)+(4) is a one-line theorem, verified in Part D:
  the affine S_n-action on the arc space has  c(tau_k) = e_{(k,k+1)}  for an adjacent
  transposition, and S_n is edge-transitive, so the "isomorphism-difference" subgroup
  is ALL of F_2^E.  Hence NO proper subgroup has iso-invariant orbits -- the search
  for a descending star-type invariant is provably empty, and THM-1405's transversality
  is forced rather than accidental.  Over A_n exactly one dimension survives.
"""
import numpy as np
from math import factorial
from itertools import permutations, combinations

# ------------------------------------------------------------------ scaffolding
def scaffold(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    tiles = [(i, j) for (i, j) in pairs if j - i >= 2]     # K_n minus the base path
    return pairs, idx, tiles, len(pairs), len(tiles)

def perm_table(p, pairs, idx):
    """A_new[t] = A[src[t]] ^ fl[t]  for the relabeling p (affine: fl is c(p))."""
    E = len(pairs)
    src = np.empty(E, dtype=np.int64); fl = np.zeros(E, dtype=np.uint8)
    for e, (i, j) in enumerate(pairs):
        a, b = p[i], p[j]
        t = idx[(min(a, b), max(a, b))]
        src[t] = e; fl[t] = 1 if a > b else 0
    return src, fl

def canon_batch(A, pairs, idx, n):
    """canonical code = min over all n! relabelings; A is (N,E) uint8."""
    E = len(pairs); pow2 = (1 << np.arange(E, dtype=np.int64))
    best = None
    for p in permutations(range(n)):
        src, fl = perm_table(p, pairs, idx)
        code = (A[:, src] ^ fl) @ pow2
        best = code if best is None else np.minimum(best, code)
    return best

def code_to_A(codes, E):
    c = np.asarray(codes, dtype=np.int64).reshape(-1, 1)
    return ((c >> np.arange(E, dtype=np.int64)) & 1).astype(np.uint8)

def beats(A_row, pairs, n):
    """adjacency: w[i][j] = 1 iff i -> j.  bit=1 means the pair is REVERSED."""
    w = np.zeros((n, n), dtype=np.int8)
    for e, (i, j) in enumerate(pairs):
        if A_row[e]: w[j][i] = 1
        else:        w[i][j] = 1
    return w

def ham_paths(w, n):
    """H(T) = number of directed Hamiltonian paths, by subset DP."""
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c or not (S >> v & 1): continue
            for u in range(n):
                if S >> u & 1 or not w[v][u]: continue
                dp[S | 1 << u][u] += c
    return sum(dp[(1 << n) - 1])

def aut_order(A_row, pairs, idx, n):
    E = len(pairs); pow2 = (1 << np.arange(E, dtype=np.int64))
    me = int(A_row @ pow2); k = 0
    for p in permutations(range(n)):
        src, fl = perm_table(p, pairs, idx)
        if int((A_row[src] ^ fl) @ pow2) == me: k += 1
    return k

def classes_upto(n):
    """iso-class representative codes, grown by vertex extension."""
    reps = {0}                                    # n = 1
    for k in range(2, n + 1):
        pairs, idx, _, E, _ = scaffold(k)
        prevE = k * (k - 1) // 2 - (k - 1)        # edges among the first k-1 vertices
        oldpairs, oldidx, _, _, _ = scaffold(k - 1)
        cand = []
        for r in reps:
            base = np.zeros(E, dtype=np.uint8)
            for e, (i, j) in enumerate(oldpairs):
                if r >> e & 1: base[idx[(i, j)]] = 1
            newe = [idx[(i, k - 1)] for i in range(k - 1)]
            for mask in range(1 << (k - 1)):
                v = base.copy()
                for b, e in enumerate(newe): v[e] = (mask >> b) & 1
                cand.append(v)
        A = np.array(cand, dtype=np.uint8)
        reps = set(int(c) for c in canon_batch(A, pairs, idx, k).tolist())
    return sorted(reps)

# =================================================================== PART A
print("=" * 78)
print("PART A -- the EXACT tiling-set -> iso-class-node map:  t(C) = H(C) / |Aut(C)|")
print("=" * 78)
print(f"{'n':>3} {'classes':>8} {'sum t':>10} {'2^m':>10} {'match':>6} {'all t odd':>10} {'|Aut| odd':>10}")
CLASSDATA = {}
for n in range(3, 8):
    pairs, idx, tiles, E, m = scaffold(n)
    reps = classes_upto(n)
    rows = []
    for r in reps:
        Ar = code_to_A([r], E)[0]
        w = beats(Ar, pairs, n)
        H = ham_paths(w, n); a = aut_order(Ar, pairs, idx, n)
        assert H % a == 0, "H not divisible by |Aut| -- orbit-stabilizer broken"
        rows.append((r, H, a, H // a))
    CLASSDATA[n] = rows
    tot = sum(t for _, _, _, t in rows)
    print(f"{n:>3} {len(rows):>8} {tot:>10} {1<<m:>10} {str(tot == 1<<m):>6} "
          f"{str(all(t % 2 for _, _, _, t in rows)):>10} "
          f"{str(all(a % 2 for _, _, a, _ in rows)):>10}")

print()
print("  Consequences, all forced by Redei (H odd) + |Aut(T)| odd:")
for n in range(3, 8):
    rows = CLASSDATA[n]
    print(f"    n={n}: every fibre |t| is ODD, #classes={len(rows)} is "
          f"{'EVEN' if len(rows)%2==0 else 'ODD'}  (sum of {len(rows)} odds = 2^m)")

# =================================================================== PART B
print()
print("=" * 78)
print("PART B -- computing nodes / tilings / edges from each other (the tricks)")
print("=" * 78)
for n in range(4, 8):
    pairs, idx, tiles, E, m = scaffold(n)
    rows = CLASSDATA[n]
    # trick 1: metagraph edges from Aut-ORBITS on arcs, never from 2^E enumeration
    flips = 0; nbr = {}
    for r, H, a, t in rows:
        Ar = code_to_A([r], E)[0]
        # orbits of Aut(T) on arcs
        seen = set(); orbs = []
        stab = []
        pow2 = (1 << np.arange(E, dtype=np.int64)); me = int(Ar @ pow2)
        for p in permutations(range(n)):
            src, fl = perm_table(p, pairs, idx)
            if int((Ar[src] ^ fl) @ pow2) == me: stab.append(src)
        for e in range(E):
            if e in seen: continue
            o = set()
            frontier = {e}
            while frontier:
                x = frontier.pop(); o.add(x)
                for src in stab:
                    y = int(np.where(src == x)[0][0])
                    if y not in o: frontier.add(y)
            seen |= o; orbs.append(sorted(o))
        flips += len(orbs)
        # one canonicalization per ORBIT instead of per arc
        cand = []
        for o in orbs:
            v = Ar.copy(); v[o[0]] ^= 1; cand.append(v)
        cc = canon_batch(np.array(cand, dtype=np.uint8), pairs, idx, n)
        for o, c in zip(orbs, cc.tolist()):
            nbr[(r, int(c))] = nbr.get((r, int(c)), 0) + len(o)
    brute = (1 << E) * E
    smart = flips
    print(f"  n={n}: arc-metagraph via Aut-orbits -> {smart:>7} canonicalizations "
          f"(brute force over labelled tournaments: {brute:>12})  speedup {brute//max(smart,1):>10}x")
    # symmetry check: flip-count C->C' must equal C'->C  (flipping is an involution)
    ok = True
    for (c1, c2), k in nbr.items():
        A1 = dict(((r, (H, a, t)) for r, H, a, t in rows))
        n1 = factorial(n) // A1[c1][1]
        n2 = factorial(n) // A1[c2][1] if c2 in A1 else None
        if n2 is None: ok = False; break
        if n1 * k != n2 * nbr.get((c2, c1), 0): ok = False; break
    print(f"        involution consistency  (n!/|Aut(C)|)*N(C,C') == (n!/|Aut(C')|)*N(C',C): {ok}")

print()
print("  THE DICTIONARY (each computed from the others, no 2^m enumeration needed):")
print("    node -> tilings :  t(C) = H(C)/|Aut(C)|          [orbit-stabilizer on the path fibration]")
print("    node -> merged  :  fibre = (2 - [C is SC]) * t(C)  [H and Aut are complement-invariant]")
print("    node -> edges   :  flip ONE arc per Aut(C)-orbit, weight by orbit size")
print("    global check    :  sum_C H(C)/|Aut(C)| = 2^C(n-1,2)")

# =================================================================== PART C
print()
print("=" * 78)
print("PART C -- the base-path-independent subgroup:  \\cap Gamma_P over spanning paths")
print("=" * 78)
for n in range(4, 8):
    pairs, idx, tiles, E, m = scaffold(n)
    # Gamma_P lives on E \ E(P); as a subspace of F_2^E it is supported off E(P).
    # Every edge lies on SOME Hamiltonian path, so the intersection is supported on
    # the empty set.  Verify the covering claim directly.
    covered = set()
    for p in permutations(range(n)):
        for k in range(n - 1):
            i, j = p[k], p[k + 1]
            covered.add(idx[(min(i, j), max(i, j))])
    print(f"  n={n}: edges covered by some spanning path: {len(covered)}/{E}"
          f"   =>  \\cap_P Gamma_P = {{0}}  ({'trivial' if len(covered)==E else 'NONTRIVIAL!'})")

print()
print("  But Gamma_P is the RESTRICTION of a base-path-independent group:")
print("     Cut(K_n minus P) = Cut(K_n) restricted to E \\ E(P)   [incidence rows restrict]")
print("  and Cut(K_n) = the SEIDEL SWITCHING group (switch at S reverses all S:S^c arcs).")
for n in range(3, 9):
    E = n * (n - 1) // 2
    print(f"  n={n}: switching group dim n-1 = {n-1};  switching classes = 2^(C(n,2)-n+1) = "
          f"{1 << (E - n + 1):>10};  (all-ones NOT a cut, so complementation is independent)")

# =================================================================== PART D
print()
print("=" * 78)
print("PART D -- THE THEOREM: no proper subgroup has iso-invariant orbits")
print("=" * 78)
print("  Gamma_min := smallest subgroup with  P.x + x in Gamma  for all x, all P.")
print("  Orbits of Gamma are unions of iso classes  <=>  Gamma contains Gamma_min.")
print()

def span_dim(vecs, E):
    basis, piv = [], []
    for v in vecs:
        for b, p in zip(basis, piv):
            if v >> p & 1: v ^= b
        if v:
            basis.append(v); piv.append(v.bit_length() - 1)
    return basis, piv

def parity(p):
    n = len(p); s = 0
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]: s ^= 1
    return s

print(f"{'n':>3} {'C(n,2)':>7} {'dim Gamma_min (S_n)':>20} {'quotient':>9} "
      f"{'dim (A_n)':>10} {'quotient':>9}")
for n in range(3, 8):
    pairs, idx, _, E, _ = scaffold(n)
    genS, genA = [], []
    for p in permutations(range(n)):
        src, fl = perm_table(p, pairs, idx)
        c = int(sum(int(b) << t for t, b in enumerate(fl)))          # the affine part c(P)
        vs = [c] + [(1 << t) ^ (1 << int(src[t])) for t in range(E) if src[t] != t]
        genS += vs
        if parity(p) == 0: genA += vs
    bS, _ = span_dim(genS, E); bA, _ = span_dim(genA, E)
    print(f"{n:>3} {E:>7} {len(bS):>20} {E-len(bS):>9} {len(bA):>10} {E-len(bA):>9}")

print()
print("  ONE-LINE PROOF of the S_n column.  For the adjacent transposition tau_k the ONLY")
print("  pair whose endpoint order reverses is (k,k+1) itself, so  c(tau_k) = e_{(k,k+1)}.")
print("  Gamma_min is S_n-invariant and S_n is edge-transitive on K_n, hence Gamma_min")
print("  contains every e_f, i.e. Gamma_min = F_2^E.  So the only subgroup whose orbits are")
print("  unions of iso classes is the WHOLE space -- no star group, no switching group, no")
print("  \\cap Gamma_P, no subgroup whatsoever, descends.  THM-1405 transversality is FORCED.")
print()
print("  Over A_n one dimension survives: the total parity  sum_e x_e  = the parity of the")
print("  feedback arc set w.r.t. the reference order.  It is A_n-invariant and flips by")
print("  sgn(P), so it is a CHIRALITY, not an invariant.")

# every |Aut| is odd => Aut <= A_n => every iso class splits into exactly 2 A_n-classes
print()
print("  Corollary (tournament chirality).  |Aut(T)| is odd, so every element of Aut(T) has")
print("  odd order, so Aut(T) <= A_n.  Hence every iso class splits into EXACTLY TWO")
print("  A_n-classes, separated by FAS-parity.  Check:")
print(f"{'n':>3} {'S_n classes':>12} {'A_n classes':>12} {'ratio':>7} {'parity separates':>17}")
for n in range(3, 7):
    pairs, idx, _, E, _ = scaffold(n)
    pow2 = (1 << np.arange(E, dtype=np.int64))
    N = 1 << E
    A = code_to_A(np.arange(N), E)
    bestS = bestA = None
    for p in permutations(range(n)):
        src, fl = perm_table(p, pairs, idx)
        code = (A[:, src] ^ fl) @ pow2
        bestS = code if bestS is None else np.minimum(bestS, code)
        if parity(p) == 0:
            bestA = code if bestA is None else np.minimum(bestA, code)
    nS, nA = len(set(bestS.tolist())), len(set(bestA.tolist()))
    # does FAS-parity separate the two A_n-classes inside each S_n-class?
    par = A.sum(axis=1) % 2
    sep = True
    byS = {}
    for i in range(N): byS.setdefault(int(bestS[i]), set()).add((int(bestA[i]), int(par[i])))
    for s, pieces in byS.items():
        ap = {}
        for a, q in pieces: ap.setdefault(a, set()).add(q)
        if len(ap) != 2 or any(len(v) != 1 for v in ap.values()) or \
           len({list(v)[0] for v in ap.values()}) != 2: sep = False
    print(f"{n:>3} {nS:>12} {nA:>12} {nA/nS:>7.1f} {str(sep):>17}")
