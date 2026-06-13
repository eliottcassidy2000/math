#!/usr/bin/env python3
"""
skew_doubling_core_kps1.py — kind-pasteur-2026-06-09-S1

Core library + base verification for the SKEW-SYLVESTER DOUBLING of tournaments
(THM-447, T767, HYP-2332..2337).

D(T) on 2n vertices, dominance blocks  M' = [[M, M+I], [M-I, -M]]
  = three copies of T's arc set (within copy 1, both cross directions)
    + one NEGATED copy (copy 2 = T^op) + twin arcs i -> i'.
The skew analog of Sylvester/Walsh doubling [[H,H],[H,-H]].

Compare:
  T[K_2]  (OPEN-Q-045): [[M, M+I], [M-I,  M]]   (all-positive doubling)
  SC blowup (OPEN-Q-044): [[M, -M+I], [-M-I, M]] (negated CROSS blocks)

Checks here (exhaustive over labeled tournaments n=3..5):
  C1  D(T) is a tournament; scores = {2 s_i + 1} u {n-1 (x n)}
  C2  M'^2 == I_2 (x) (2 M^2 - I)  exactly (integer matrices)
  C3  S' = M'+I is skew-type; skew-Hadamard preserved (bordered C_3 seed)
  C4  canonical Hamiltonian path n,...,1,1',2',...,n' is valid in D(T)
  C5  spectral map lambda -> sqrt(2 lambda^2 + 1) (numerical)
  C6  family laws: squares of T[K_2] and SC-blowup matrices (NOT block-diagonal)
  C7  H on iso classes n=3..5 and their doubles (data for HYP-2334)
  C8  D(T)^op ?= D(T^op) up to iso (SC preservation probe)

Output: 05-knowledge/results/skew_doubling_core_kps1.out
"""
import itertools, sys
import numpy as np

# ---------- tournament basics ----------

def all_tournaments(n):
    """Yield 0/1 adjacency matrices of all labeled tournaments on n vertices."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        yield A

def M_of(A):
    return A - A.T

def A_of(M):
    return (M > 0).astype(np.int64)

def scores(A):
    return tuple(int(x) for x in A.sum(axis=1))

def canon(A):
    """Canonical form: lexicographically minimal adjacency bitstring over all relabelings."""
    n = A.shape[0]
    best = None
    for p in itertools.permutations(range(n)):
        P = A[np.ix_(p, p)]
        key = P[np.triu_indices(n, 1)].tobytes() + P[np.tril_indices(n, -1)].tobytes()
        if best is None or key < best[0]:
            best = (key, P.copy())
    return best[1]

def iso_classes(n):
    seen = {}
    for A in all_tournaments(n):
        C = canon(A)
        seen.setdefault(C.tobytes(), C)
    return list(seen.values())

def is_iso(A, B):
    n = A.shape[0]
    if scores(A) != scores(B) and sorted(scores(A)) != sorted(scores(B)):
        return False
    for p in itertools.permutations(range(n)):
        if np.array_equal(A[np.ix_(p, p)], B):
            return True
    return False

# ---------- Hamiltonian path count (held-karp DP) ----------

def H_count(A):
    """Number of directed Hamiltonian paths."""
    n = A.shape[0]
    adj = [int(''.join(str(A[i, j]) for j in range(n - 1, -1, -1)), 2) for i in range(n)]
    # adj[i] bitmask of out-neighbors
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if not c:
                continue
            avail = adj[last] & ~mask
            while avail:
                b = avail & -avail
                nxt = b.bit_length() - 1
                dp[mask | b][nxt] += c
                avail ^= b
    full = (1 << n) - 1
    return sum(dp[full])

# ---------- the doubling family ----------

def D_skew(A):
    """Skew-Sylvester double: [[M, M+I],[M-I, -M]]."""
    n = A.shape[0]
    M = M_of(A)
    I = np.eye(n, dtype=np.int64)
    Mp = np.block([[M, M + I], [M - I, -M]])
    return A_of(Mp), Mp

def D_blowup(A):
    """T[K_2] (OPEN-Q-045): [[M, M+I],[M-I, M]]."""
    n = A.shape[0]
    M = M_of(A)
    I = np.eye(n, dtype=np.int64)
    Mp = np.block([[M, M + I], [M - I, M]])
    return A_of(Mp), Mp

def D_scblowup(A):
    """SC blowup (OPEN-Q-044): [[M, -M+I],[-M-I, M]]."""
    n = A.shape[0]
    M = M_of(A)
    I = np.eye(n, dtype=np.int64)
    Mp = np.block([[M, -M + I], [-M - I, M]])
    return A_of(Mp), Mp

# ---------- skew-Hadamard machinery ----------

def border(A):
    """Border tournament: new vertex 0 dominating all. Returns S matrix of order n+1."""
    n = A.shape[0]
    M = M_of(A)
    Mb = np.zeros((n + 1, n + 1), dtype=np.int64)
    Mb[0, 1:] = 1
    Mb[1:, 0] = -1
    Mb[1:, 1:] = M
    return Mb + np.eye(n + 1, dtype=np.int64)

def is_skew_hadamard(S):
    m = S.shape[0]
    return (np.abs(S) == 1).all() and np.array_equal(S + S.T, 2 * np.eye(m, dtype=np.int64)) \
        and np.array_equal(S @ S.T, m * np.eye(m, dtype=np.int64))

def normalize_first_row(S):
    """Negate column j AND row j wherever S[0,j] = -1 (preserves skew-type and Hadamard)."""
    S = S.copy()
    for j in range(1, S.shape[0]):
        if S[0, j] == -1:
            S[:, j] *= -1
            S[j, :] *= -1
    return S

def core_tournament(S):
    """Drop border row/col of a normalized skew-Hadamard; return core adjacency."""
    Sc = S[1:, 1:]
    M = Sc - np.eye(Sc.shape[0], dtype=np.int64)
    return A_of(M)

def is_DRT(A):
    """Doubly regular: regular and |N+(u) cap N+(v)| constant over pairs."""
    n = A.shape[0]
    s = set(scores(A))
    if len(s) != 1:
        return False
    vals = set()
    for u in range(n):
        for v in range(u + 1, n):
            vals.add(int((A[u] * A[v]).sum()))
    return len(vals) == 1

# ---------- verification harness ----------

def check_all(n, out):
    I2n = np.eye(2 * n, dtype=np.int64)
    c1 = c2 = c3 = c4 = total = 0
    for A in all_tournaments(n):
        total += 1
        n_ = A.shape[0]
        Ad, Mp = D_skew(A)
        # C1: tournament + scores
        ok_t = np.array_equal(Ad + Ad.T, np.ones((2 * n_, 2 * n_), dtype=np.int64) - I2n)
        s = scores(A)
        sd = scores(Ad)
        ok_s = sorted(sd) == sorted([2 * x + 1 for x in s] + [n_ - 1] * n_)
        c1 += ok_t and ok_s
        # C2: spectral law
        M = M_of(A)
        law = np.block([[2 * M @ M - np.eye(n_, dtype=np.int64), np.zeros((n_, n_), dtype=np.int64)],
                        [np.zeros((n_, n_), dtype=np.int64), 2 * M @ M - np.eye(n_, dtype=np.int64)]])
        c2 += np.array_equal(Mp @ Mp, law)
        # C3: skew-type
        Sp = Mp + I2n
        c3 += np.array_equal(Sp + Sp.T, 2 * I2n)
        # C4: canonical Ham path (vertices n-1..0 then 0'..(n-1)' in our 0-indexing:
        # base path of T assumed n-1 -> n-2 -> ... -> 0 IF T contains it; here instead
        # check the constructive claim: for ANY Ham path p of T, p, twin, reversed p' is a path of D)
        # use any Ham path found by greedy insertion (every tournament has one)
        p = ham_path(A)
        ok_p = all(Ad[p[i], p[i + 1]] for i in range(n_ - 1))
        ok_p &= bool(Ad[p[-1], n_ + p[-1]])  # twin arc at the path's END vertex
        rp = [n_ + v for v in reversed(p)]
        ok_p &= all(Ad[rp[i], rp[i + 1]] for i in range(n_ - 1))
        c4 += ok_p
    out.write(f"n={n}: total={total}  C1 tournament+scores={c1}  C2 M'^2 law={c2}"
              f"  C3 skew-type={c3}  C4 canonical-path={c4}\n")
    return c1 == total and c2 == total and c3 == total and c4 == total

def ham_path(A):
    """Insertion algorithm: every tournament has a Hamiltonian path."""
    n = A.shape[0]
    path = [0]
    for v in range(1, n):
        placed = False
        for k in range(len(path) + 1):
            before_ok = (k == 0) or A[path[k - 1], v]
            after_ok = (k == len(path)) or A[v, path[k]]
            if before_ok and after_ok:
                path.insert(k, v)
                placed = True
                break
        assert placed
    return path

def main():
    out = open("05-knowledge/results/skew_doubling_core_kps1.out", "w", encoding="utf-8")
    def w(s=""):
        out.write(s + "\n"); out.flush(); print(s)

    w("=== skew_doubling_core_kps1 — THM-447 base verification ===")
    w("")
    # C1-C4 exhaustive
    allok = True
    for n in (3, 4, 5):
        allok &= check_all(n, out)
    w(f"C1-C4 exhaustive n=3..5: {'ALL PASS' if allok else 'FAILURES — see above'}")
    w("")

    # C3b: skew-Hadamard preservation along the bordered tower (seed: 1-vertex)
    w("--- skew-Hadamard tower from order-1 seed [1] ---")
    S = np.array([[1]], dtype=np.int64)
    for k in range(1, 6):
        n = S.shape[0]
        I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
        w(f"order {S.shape[0]}: skew-Hadamard = {is_skew_hadamard(S)}")
    w("")

    # bordered C_3 -> order 4 -> doubled -> order 8 -> core on 7 vertices: DRT? Paley?
    w("--- C_3 doubling tower (bordered) ---")
    C3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    S4 = border(C3)
    w(f"bordered C_3: order 4 skew-Hadamard = {is_skew_hadamard(S4)}")
    I4 = np.eye(4, dtype=np.int64)
    S8 = np.block([[S4, S4], [S4 - 2 * I4, 2 * I4 - S4]])
    w(f"doubled: order 8 skew-Hadamard = {is_skew_hadamard(S8)}")
    S8n = normalize_first_row(S8)
    w(f"normalized order 8 still skew-Hadamard = {is_skew_hadamard(S8n)}")
    T7 = core_tournament(S8n)
    w(f"core on 7 vertices: scores={scores(T7)}  DRT={is_DRT(T7)}  H={H_count(T7)}")
    # Paley T_7: i -> j iff (j - i) is a QR mod 7 (QRs: 1,2,4)
    P7 = np.zeros((7, 7), dtype=np.int64)
    for i in range(7):
        for j in range(7):
            if (j - i) % 7 in (1, 2, 4):
                P7[i, j] = 1
    w(f"Paley T_7: H={H_count(P7)}  iso(core, Paley) = {is_iso(T7, P7)}")
    w("")

    # C5 numerical spectral map on a sample
    w("--- C5: eigenvalue map lambda -> sqrt(2 lambda^2 + 1) (sample n=5) ---")
    rng = np.random.default_rng(7)
    A = next(itertools.islice(all_tournaments(5), 333, None))
    M = M_of(A)
    _, Mp = D_skew(A)
    lam = np.sort(np.abs(np.linalg.eigvals(M).imag))
    lamp = np.sort(np.abs(np.linalg.eigvals(Mp.astype(float)).imag))
    pred = np.sort(np.concatenate([np.sqrt(2 * lam ** 2 + 1)] * 2))
    w(f"lambda(T)      = {np.round(lam, 6)}")
    w(f"lambda(D(T))   = {np.round(lamp, 6)}")
    w(f"predicted      = {np.round(pred, 6)}")
    w(f"match = {np.allclose(np.sort(lamp), pred, atol=1e-8)}")
    w("")

    # C6: family squared laws (symbolic check via random tournaments)
    w("--- C6: family M'^2 structure (n=5 sample) ---")
    for name, fn in (("D_skew", D_skew), ("T[K2]", D_blowup), ("SC-blowup", D_scblowup)):
        _, Mp = fn(A)
        n_ = A.shape[0]
        sq = Mp @ Mp
        bd = (np.array_equal(sq[:n_, n_:], np.zeros((n_, n_), dtype=np.int64)) and
              np.array_equal(sq[n_:, :n_], np.zeros((n_, n_), dtype=np.int64)))
        w(f"{name:10s}: M'^2 block-diagonal = {bd}; "
          f"S' skew-Hadamard from bordered-C3 seed = "
          f"{is_skew_hadamard((lambda S: np.block([[S, S],[S-2*np.eye(4,dtype=np.int64), 2*np.eye(4,dtype=np.int64)-S]]))(border(C3))) if name=='D_skew' else 'n/a'}")
    # explicit: does T[K2] or SC blowup preserve skew-Hadamard? test on bordered C3 cores
    for name, fn in (("T[K2]", D_blowup), ("SC-blowup", D_scblowup)):
        _, Mp = fn(core_tournament(normalize_first_row(S4)) if False else C3)
        Sp = Mp + np.eye(Mp.shape[0], dtype=np.int64)
        w(f"{name:10s}: S'(C_3) S'(C_3)^T == 6I ? {np.array_equal(Sp @ Sp.T, 6 * np.eye(6, dtype=np.int64))}")
    w("")

    # C7: H data on iso classes n=3..5
    w("--- C7: H(T), H(D(T)), H(T[K2]), H(SCblow) on iso classes ---")
    w(f"{'n':>2} {'idx':>3} {'scores':>18} {'H(T)':>6} {'H(D)':>8} {'H(K2)':>8} {'H(SC)':>8}")
    hdata = []
    for n in (3, 4, 5):
        for idx, A in enumerate(iso_classes(n)):
            hT = H_count(A)
            hD = H_count(D_skew(A)[0])
            hK = H_count(D_blowup(A)[0])
            hS = H_count(D_scblowup(A)[0])
            hdata.append((n, idx, scores(A), hT, hD, hK, hS))
            w(f"{n:>2} {idx:>3} {str(scores(A)):>18} {hT:>6} {hD:>8} {hK:>8} {hS:>8}")
    w("")

    # C8: D(T)^op vs D(T^op) iso probe (n=3,4)
    w("--- C8: D(T)^op ~ D(T^op)? (iso probe, n=3,4) ---")
    for n in (3, 4):
        agree = 0; tot = 0
        for A in iso_classes(n):
            Dop = D_skew(A)[0].T  # ^op = transpose of adjacency
            DTop = D_skew(A.T)[0]
            tot += 1
            agree += is_iso(Dop, DTop)
        w(f"n={n}: D(T)^op iso D(T^op) in {agree}/{tot} classes")
    w("")
    w("=== done ===")
    out.close()

if __name__ == "__main__":
    main()
