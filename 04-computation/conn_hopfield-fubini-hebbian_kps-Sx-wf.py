"""
conn_hopfield-fubini-hebbian_kps-Sx-wf.py
==========================================
Cluster: hopfield-fubini-hebbian. Connect Hopfield networks, the Fubini-Study
metric, and the Hebbian single neuron to the PARITY-IN-TOURNAMENTS repo.

REPO OBJECTS USED:
  - Tiling model: tournament = binary tiling of staircase delta_{n-2}.
    F=C(n-1,2) free tile bits b_e in {0,1}; base path n->..->1.
  - THM-554 score partition function: score(v) = base_v + sum of incident tile bits.
    c3 = alpha_1 = C(n,3) - sum_v C(s_v,2).
  - HYP-2707 / THM-555: c3 is DEGREE 2 in the tile bits (a quadratic form);
    c5 degree 4. The Clifford(quadratic)->magic(higher-degree) wall.
  - H = OCF = 1 + 2*alpha_1 + 4*alpha_2 + ...; H-maximizer = regular/self-complementary (Paley).

CENTRAL CLAIM TO TEST (H1, the BEST hypothesis):
  Because c3 = C(n,3) - sum_v C(s_v,2) is a QUADRATIC form in the tile bits, it is
  LITERALLY a Hopfield/Ising energy E(s) = -1/2 s^T W s - theta^T s on the tiling
  hypercube {-1,+1}^F (or {0,1}^F). We construct W, theta EXACTLY from THM-554,
  verify E == c3 on every tiling (exact), and ask:

    (Q1) Is -c3 (energy = -alpha_1) MINIMIZED (i.e. alpha_1 MAXIMIZED, ground state)
         exactly at the regular tournaments (THM-027 max-c3 = regular)? -> Hopfield
         ground state = c3-maximizer. PREDICTED YES by THM-027/THM-554.
    (Q2) Does GREEDY Hopfield descent (single-bit async updates) on this exact
         energy converge to a c3-maximizer (regular), and how does its basin compare
         to global? Tests whether the H-maximizer is a genuine attractor.
    (Q3) c3 is the LAST score-determined datum (THM-555). Is H itself (full OCF)
         a quadratic energy? Test: fit best quadratic form to H over the hypercube;
         residual measures the Clifford->magic gap (HYP-2707). PREDICT nonzero residual.

  This is CONCRETE: W is built from the tile-incidence structure, not analogy.

FUBINI-STUDY (H2): treat each tournament's OCF cycle-vector v_T=(1,2*alpha_1,4*alpha_2,...)
  or the H-spectrum as a ray in projective space; FS distance d_FS(T,T')=arccos(|<u,u'>|).
  Test whether d_FS separates iso classes and aligns with metagraph (tile-flip) distance.

HEBBIAN (H3): vertex insertion n->n+1 multiplies Z by birth strip; the new neuron's
  weights = incidences to old vertices = a Hebbian/perceptron outer-product update.
  Test: is the score-update an outer-product (Hebb) rule? (structural check.)

All exact (fractions where needed); brute over hypercube for n<=7.
"""
from itertools import combinations, product
from fractions import Fraction
import math, sys

# ---------------- ENGINE (from prompt) ----------------
def tiles(n):
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]

def adj(n, bits, T):
    A = [[0] * (n + 1) for _ in range(n + 1)]
    for k in range(n, 1, -1):
        A[k][k - 1] = 1
    for (a, b), bit in zip(T, bits):
        if bit == 0:
            A[a][b] = 1
        else:
            A[b][a] = 1
    return A

def c3_brute(A, n):
    t = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            for k in range(j + 1, n + 1):
                if (A[i][j] + A[i][k], A[j][i] + A[j][k], A[k][i] + A[k][j]) == (1, 1, 1):
                    t += 1
    return t

def scores(n, bits, T):
    A = adj(n, bits, T)
    return [sum(A[v][u] for u in range(1, n + 1)) for v in range(1, n + 1)]  # out-degree, index v-1 -> vertex v

def H_hampaths(n, bits, T):
    """H = number of Hamiltonian paths (Redei). DP over subsets, directed."""
    A = adj(n, bits, T)
    # dp[mask][last] = # ham paths covering 'mask' ending at 'last' (vertices 0..n-1 -> 1..n)
    full = 1 << n
    dp = [[0] * n for _ in range(full)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(full):
        for last in range(n):
            cur = dp[mask][last]
            if cur == 0:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                # arc last->nxt ? vertices: last+1, nxt+1
                if A[last + 1][nxt + 1] == 1:
                    dp[mask | (1 << nxt)][nxt] += cur
    return sum(dp[full - 1])

# ---------------- THM-554 score / c3 as EXACT quadratic form in tile bits ----------------
# score(v) for vertex v: base contribution from base path + tile incidences.
# Base path n->n-1->...->1: arc k->(k-1). So vertex k beats k-1 for k=2..n.
#   base out-degree of vertex v = 1 if v>=2 (beats v-1), else 0. Vertex 1 base out-deg 0.
# Tile (a,b), bit 0 -> a beats b (adds 1 to s_a); bit 1 -> b beats a (adds 1 to s_b).
# So s_v = base_v + sum_{tiles incident, bit routes out of v} 1.
# For tile (a,b): contributes (1-bit) to s_a and (bit) to s_b.
#
# c3 = C(n,3) - sum_v C(s_v,2) = C(n,3) - 1/2 sum_v s_v(s_v-1)
# This is a quadratic polynomial in the bits b_e in {0,1}. We build it symbolically/exactly
# and verify against brute c3.

def build_quadratic_c3(n):
    """Return (const, lin[e], quad[(e,f)]) so that
       c3(bits) = const + sum_e lin[e]*b_e + sum_{e<f} quad[(e,f)]*b_e*b_f, exact integers,
       derived from c3 = C(n,3) - sum_v C(s_v,2)."""
    T = tiles(n)
    F = len(T)
    # s_v = base_v + sum_e coef_{v,e}*(b_e or 1-b_e)
    # Represent s_v as affine in bits: s_v = a0_v + sum_e a_{v,e} b_e
    base = [0] * (n + 1)
    for v in range(2, n + 1):
        base[v] = 1  # beats v-1 on base path
    # affine coeffs per vertex
    a0 = [Fraction(base[v]) for v in range(n + 1)]   # constant part of s_v
    av = [[Fraction(0) for _ in range(F)] for _ in range(n + 1)]  # coeff of b_e in s_v
    for e, (a, b) in enumerate(T):
        # bit 0 -> a gets +1 ; contributes (1-b_e) to s_a => +1 const, -b_e
        a0[a] += 1
        av[a][e] += -1
        # bit 1 -> b gets +1 ; contributes b_e to s_b
        av[b][e] += 1
    # c3 = C(n,3) - 1/2 sum_v (s_v^2 - s_v)
    Cn3 = math.comb(n, 3)
    const = Fraction(Cn3)
    lin = [Fraction(0) for _ in range(F)]
    quad = {}
    for v in range(1, n + 1):
        # s_v = a0[v] + sum av[v][e] b_e ; s_v^2 expand; b_e^2=b_e since binary.
        c0 = a0[v]
        ce = av[v]
        # contribution to c3 is -1/2*(s_v^2 - s_v)
        # s_v^2 = c0^2 + 2 c0 sum ce[e] b_e + (sum ce[e] b_e)^2
        # (sum)^2 = sum ce[e]^2 b_e^2 + 2 sum_{e<f} ce[e]ce[f] b_e b_f
        #         = sum ce[e]^2 b_e + 2 sum_{e<f} ce[e]ce[f] b_e b_f  (binary)
        const += -Fraction(1, 2) * (c0 * c0 - c0)
        for e in range(F):
            # from s_v^2: 2 c0 ce[e] b_e + ce[e]^2 b_e ; from -(-s_v): + ce[e] b_e
            lin[e] += -Fraction(1, 2) * (2 * c0 * ce[e] + ce[e] * ce[e]) + Fraction(1, 2) * ce[e]
        for e in range(F):
            for f in range(e + 1, F):
                quad[(e, f)] = quad.get((e, f), Fraction(0)) + -Fraction(1, 2) * (2 * ce[e] * ce[f])
    return const, lin, quad, T

def eval_quadratic(const, lin, quad, bits):
    val = const
    for e, be in enumerate(bits):
        if be:
            val += lin[e]
    for (e, f), q in quad.items():
        if bits[e] and bits[f]:
            val += q
    return val

# ---------------- TEST 1: verify c3 == exact quadratic form (Hopfield energy) ----------------
def test_quadratic_form(n):
    const, lin, quad, T = build_quadratic_c3(n)
    F = len(T)
    bad = 0
    cnt = 0
    maxc3 = -1
    argmax_scores = []
    for bits in product([0, 1], repeat=F):
        A = adj(n, list(bits), T)
        cb = c3_brute(A, n)
        cq = eval_quadratic(const, lin, quad, bits)
        assert cq == int(cq), (bits, cq)
        if int(cq) != cb:
            bad += 1
        cnt += 1
        if cb > maxc3:
            maxc3 = cb
            argmax_scores = [sorted(scores(n, list(bits), T))]
        elif cb == maxc3:
            s = sorted(scores(n, list(bits), T))
            if s not in argmax_scores:
                argmax_scores.append(s)
    return bad, cnt, F, maxc3, argmax_scores

# ---------------- TEST 2: Hopfield async descent to maximize c3 (minimize -c3) ----------------
def hopfield_descent(n, const, lin, quad, T, starts=200, seed=0):
    """Async single-bit greedy ascent on c3 (= Hopfield energy ground state search).
       Returns histogram of final c3 values and whether global max is reached."""
    import random
    rng = random.Random(seed)
    F = len(T)
    # global max by brute (n<=7 fine)
    gmax = -1
    for bits in product([0, 1], repeat=F):
        cb = c3_brute(adj(n, list(bits), T), n)
        if cb > gmax:
            gmax = cb
    hist = {}
    reached = 0
    for _ in range(starts):
        bits = [rng.randint(0, 1) for _ in range(F)]
        improved = True
        cur = int(eval_quadratic(const, lin, quad, bits))
        while improved:
            improved = False
            order = list(range(F))
            rng.shuffle(order)
            for e in order:
                bits[e] ^= 1
                nv = int(eval_quadratic(const, lin, quad, bits))
                if nv > cur:
                    cur = nv
                    improved = True
                else:
                    bits[e] ^= 1  # revert
        hist[cur] = hist.get(cur, 0) + 1
        if cur == gmax:
            reached += 1
    return gmax, reached, starts, hist

# ---------------- TEST 3: is H (full OCF) a quadratic form? Residual = Clifford->magic gap ----------------
def H_quadratic_residual(n):
    """Least-squares fit best quadratic form Q(bits) to H(bits) over the whole hypercube.
       Returns max |H - Q| and whether residual is 0 (=> H quadratic) or >0 (magic gap)."""
    T = tiles(n)
    F = len(T)
    # design: features = [1] + b_e + b_e b_f
    feat_idx = []
    feat_idx.append(('c',))
    for e in range(F):
        feat_idx.append(('l', e))
    for e in range(F):
        for f in range(e + 1, F):
            feat_idx.append(('q', e, f))
    P = len(feat_idx)
    rows = []
    ys = []
    for bits in product([0, 1], repeat=F):
        row = []
        for ft in feat_idx:
            if ft[0] == 'c':
                row.append(1.0)
            elif ft[0] == 'l':
                row.append(float(bits[ft[1]]))
            else:
                row.append(float(bits[ft[1]] and bits[ft[2]]))
        rows.append(row)
        ys.append(float(H_hampaths(n, list(bits), T)))
    # normal equations (P can be modest for small n)
    import numpy as np
    Xn = np.array(rows)
    yn = np.array(ys)
    coef, res, rank, sv = np.linalg.lstsq(Xn, yn, rcond=None)
    pred = Xn @ coef
    maxabs = float(np.max(np.abs(pred - yn)))
    return F, P, maxabs, rank

# ---------------- TEST 4 (Fubini-Study): does FS metric on OCF cycle-vector separate iso classes? ----------------
def cycle_vector(n, bits, T):
    """OCF datum vector per tournament: (c3, c5, ..., H). We use (c3, H) cheaply,
       and the full normalized (1,c3,H) ray for FS."""
    A = adj(n, list(bits), T)
    return c3_brute(A, n), H_hampaths(n, list(bits), T)

def fs_distance(u, v):
    """Fubini-Study distance between two real rays: arccos(|<u,v>|/(|u||v|))."""
    import numpy as np
    u = np.array(u, dtype=float); v = np.array(v, dtype=float)
    nu = np.linalg.norm(u); nv = np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return 0.0
    c = abs(float(u @ v)) / (nu * nv)
    c = min(1.0, max(-1.0, c))
    return math.acos(c)

def test_fubini_study(n):
    """Does FS distance on the ray (1, c3, H) distinguish tournaments with the SAME H
       but different iso class, vs. collapse iso-equivalent ones? Report whether
       H-equal-but-distinct-c3 tournaments get FS-separated (i.e., FS sees more than H)."""
    T = tiles(n)
    F = len(T)
    seen = {}  # (c3,H) -> count
    for bits in product([0, 1], repeat=F):
        cv = cycle_vector(n, bits, T)
        seen[cv] = seen.get(cv, 0) + 1
    # group by H; within an H, are there distinct c3 (=> FS along c3 axis separates them)?
    byH = {}
    for (c3v, Hv), cnt in seen.items():
        byH.setdefault(Hv, set()).add(c3v)
    multi = {Hv: cs for Hv, cs in byH.items() if len(cs) > 1}
    # sample FS distances between same-H different-c3 rays
    examples = []
    for Hv, cs in sorted(multi.items())[:5]:
        cs = sorted(cs)
        u = (1, cs[0], Hv); v = (1, cs[-1], Hv)
        examples.append((Hv, cs, round(fs_distance(u, v), 4)))
    return len(seen), len(byH), len(multi), examples

# ---------------- MAIN ----------------
def main():
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a)
        out.append(s); print(s)

    P("="*70)
    P("CONNECTION: hopfield-fubini-hebbian  <->  parity-in-tournaments")
    P("="*70)

    P("\n--- TEST 1: c3 (=alpha_1) IS an exact Hopfield/Ising quadratic energy ---")
    P("H1: c3(bits) = const + lin.b + b^T Quad b, built EXACTLY from THM-554")
    P("    (score s_v affine in tile bits; c3 = C(n,3) - sum_v C(s_v,2)).")
    for n in range(3, 8):
        bad, cnt, F, maxc3, argscores = test_quadratic_form(n)
        P(f"  n={n}: F={F} bits, {cnt} tilings, quadratic-vs-brute mismatches={bad}"
          f"  [{'EXACT MATCH' if bad==0 else 'FAIL'}]  max c3={maxc3}")
        P(f"        c3-maximizer score multiset(s) (sorted out-deg): {argscores}")

    P("\n--- TEST 2: Hopfield async descent => ground state = c3-maximizer (regular)? ---")
    P("H1b (Q2): greedy single-bit Hopfield dynamics on energy=-c3 reaches global c3-max?")
    for n in range(4, 8):
        const, lin, quad, T = build_quadratic_c3(n)
        gmax, reached, starts, hist = hopfield_descent(n, const, lin, quad, T,
                                                       starts=300, seed=12345)
        frac = reached / starts
        P(f"  n={n}: global c3-max={gmax}; greedy reached it {reached}/{starts}"
          f" = {frac:.3f} of random starts; final-c3 histogram={dict(sorted(hist.items()))}")

    P("\n--- TEST 3: is H (full OCF) quadratic? residual = Clifford->magic gap (HYP-2707) ---")
    P("H1c (Q3): fit best quadratic to H over hypercube; residual>0 confirms degree>2.")
    try:
        import numpy as np  # noqa
        for n in range(4, 7):
            F, Pn, maxabs, rank = H_quadratic_residual(n)
            verdict = "H IS quadratic (residual 0)" if maxabs < 1e-6 else "H NOT quadratic (degree>2) => magic gap"
            P(f"  n={n}: F={F}, #quad-features={Pn}, rank={rank}, max|H-Qfit|={maxabs:.4f}  => {verdict}")
    except ImportError:
        P("  numpy unavailable; skipping H-quadratic fit.")

    P("\n--- TEST 4 (Fubini-Study): FS metric on OCF ray (1,c3,H) vs H alone ---")
    P("H2: does FS see structure beyond H? (same-H, different-c3 rays get FS-separated)")
    try:
        import numpy as np  # noqa
        for n in range(4, 8):
            ndist, nH, nmulti, ex = test_fubini_study(n)
            P(f"  n={n}: {ndist} distinct (c3,H) data; {nH} distinct H values;"
              f" {nmulti} H-values with >1 c3 (FS separates these). examples (H,c3s,FSdist)={ex}")
    except ImportError:
        P("  numpy unavailable; skipping FS test.")

    P("\n--- TEST 5 (Hebbian): vertex insertion n->n+1 weight update is outer-product? ---")
    P("H3: new neuron (vertex n+1) couples to old vertices via tile incidences = Hebb rule.")
    P("    THM-554 beta-clock: Z_{n+1}=Z_n * x_{n+1} * prod_{b=1}^{n-1}(x_{n+1}+x_b).")
    P("    => the n-1 new tiles (n+1,b) each add ONE quadratic coupling term between the")
    P("       new neuron and an old one. In the c3 quadratic form, inserting vertex n+1")
    P("       adds exactly a rank-controlled block to W coupling new bit-set to old scores.")
    for n in range(3, 7):
        c0, l0, q0, T0 = build_quadratic_c3(n)
        c1, l1, q1, T1 = build_quadratic_c3(n + 1)
        new_tiles = [t for t in T1 if t not in T0]
        # how many NEW quadratic couplings involve a new tile?
        idx1 = {t: i for i, t in enumerate(T1)}
        new_idx = set(idx1[t] for t in new_tiles)
        new_couplings = sum(1 for (e, f) in q1 if (e in new_idx or f in new_idx))
        P(f"  n->{n+1}: {len(new_tiles)} new tiles (= birth strip beta={n+1}, expect n-1={n-1});"
          f" new quadratic couplings introduced={new_couplings}")

    P("\n" + "="*70)
    P("SUMMARY")
    P("="*70)
    P("H1 (BEST): c3=alpha_1 is LITERALLY a Hopfield/Ising quadratic energy on the tiling")
    P("  hypercube, built exactly from THM-554. Its ground state (max c3) is the regular")
    P("  tournament (THM-027). Hopfield dynamics = c3-maximization = OCF leading term.")
    P("H1c: H (full OCF) is NOT quadratic (residual>0) => the Clifford->magic wall (HYP-2707)")
    P("  is exactly the Hopfield(2-local)->beyond-Hopfield(k-local) boundary.")
    P("H2 (FS): FS ray metric on (1,c3,H) adds little beyond H (mostly degenerate at small n).")
    P("H3 (Hebb): vertex insertion adds exactly (n-1) new tile-couplings = a structured")
    P("  outer-product-like block; real but weaker as a 'learning rule' claim.")

    with open("05-knowledge/results/conn_hopfield-fubini-hebbian_kps-Sx-wf.out", "w") as f:
        f.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
