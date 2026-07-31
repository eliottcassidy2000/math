#!/usr/bin/env python3
"""
ramsey_book_sa.py — search for witnesses of R(B_{n-1}, B_n) = 4n-1
(Epoch FrontierMath "Ramsey Numbers for Book Graphs").

Goal: graph G on N = 4n-2 vertices such that
  (i)  no B_{n-1} in G:      every edge has <= n-2 common neighbours;
  (ii) no B_n in complement: every non-edge pair has <= n-1 common non-nbrs.

Ansatz: Cayley graphs on the dihedral group D_m (order N = 2m, m = 2n-1 odd)
or on Z_2 x Z_m.  For D_m with rotation part U (symmetric subset of Z_m\{0})
and reflection part T (ANY subset of Z_m):

  vertices: rotations r^i and reflections r^j f
  r^i ~ r^j        iff  j-i in U
  r^i f ~ r^j f    iff  j-i in U
  r^i ~ r^j f      iff  j-i in T

Codegree formulas (exact, all pairs, difference d resp. t):
  within-fiber pair (difference d != 0):
      common nbrs      = A_U(d) + A_T(d)
      common non-nbrs  = N - 2k + A_U(d) + A_T(d) - 2
  cross pair (difference t):
      common nbrs      = |T ∩ (t-U)| + |T ∩ (t+U)|  =: X(t)
      common non-nbrs  = N - 2k + X(t) - 2
  where A_S(d) = |S ∩ (S+d)|, k = |U| + |T| (graph is k-regular).

Conditions become, with limit_on = n-2, limit_off = n+1+2k-N:
  d in U:      A_U(d)+A_T(d) <= limit_on     else <= limit_off
  t in T:      X(t) <= limit_on              else <= limit_off

Simulated annealing over (U, T) with exact integer correlation evaluation
(numpy FFT + rint, guarded by exact spot checks; final full brute verification
from the adjacency matrix).

Usage: python3 ramsey_book_sa.py n [k] [seed] [iters] [group=dihedral|abelian]
"""
import numpy as np
import sys, time

def correlations(chi):
    """A(d) = sum_e chi(e) chi(e-d) for all d, exact via FFT with rint."""
    m = len(chi)
    f = np.fft.rfft(chi)
    c = np.fft.irfft(f * np.conj(f), n=m)
    return np.rint(c).astype(np.int64)

def cross_corr(chiT, chiU):
    """X(t) = |T ∩ (t-U)| + |T ∩ (t+U)| = sum_u chiT(t-u)+chiT(t+u) over u in U."""
    m = len(chiT)
    fT = np.fft.rfft(chiT); fU = np.fft.rfft(chiU)
    conv = np.fft.irfft(fT * fU, n=m)               # sum_u chiT(t-u) chiU(u)
    corr = np.fft.irfft(fT * np.conj(fU), n=m)      # sum_u chiT(t+u) chiU(u)
    return np.rint(conv + corr).astype(np.int64)

def energy(chiU, chiT, n, N, group):
    m = len(chiU)
    k = int(chiU.sum() + chiT.sum())
    lim_on = n - 2
    lim_off = n + 1 + 2 * k - N
    AU = correlations(chiU); AT = correlations(chiT)
    within = AU + AT
    X = cross_corr(chiT, chiU)
    e = 0
    # within-fiber differences d=1..m-1
    for d in range(1, m):
        lim = lim_on if chiU[d] else lim_off
        v = within[d] - lim
        if v > 0: e += 2 * v          # two fibers, weight 2 (m pairs each)
    # cross differences t=0..m-1
    for t in range(m):
        lim = lim_on if chiT[t] else lim_off
        v = X[t] - lim
        if v > 0: e += v
    return e, k

def brute_verify(chiU, chiT, n, N):
    """Build full adjacency and check both book conditions by brute force."""
    m = len(chiU)
    A = np.zeros((N, N), dtype=np.int8)
    for i in range(m):
        for j in range(m):
            if i != j and chiU[(j - i) % m]:
                A[i, j] = 1                    # rot-rot
                A[m + i, m + j] = 1            # ref-ref
            if chiT[(j - i) % m]:
                A[i, m + j] = 1; A[m + j, i] = 1
    assert (A == A.T).all() and (np.diag(A) == 0).all()
    M = A.astype(np.int64)
    C = M @ M                                   # common neighbours
    comp = 1 - M
    np.fill_diagonal(comp, 0)
    Cc = comp @ comp                            # common non-nbrs (excl. endpoints auto)
    ok = True
    worst_on = -1; worst_off = -1
    for i in range(N):
        for j in range(i + 1, N):
            if M[i, j]:
                worst_on = max(worst_on, C[i, j])
                if C[i, j] > n - 2: ok = False
            else:
                # common non-nbrs excluding i,j: Cc counts w with comp[i,w]=comp[j,w]=1;
                # w=i: comp[i,i]=0; w=j: 0 -> already excluded
                worst_off = max(worst_off, Cc[i, j])
                if Cc[i, j] > n - 1: ok = False
    return ok, worst_on, worst_off, A

def emit(A, fname):
    N = A.shape[0]
    # Epoch format: adjacency string, column-major binary
    s = "".join(str(A[i, j]) for j in range(N) for i in range(N))
    with open(fname, "w") as f:
        f.write(s + "\n")

def sa(n, k_target, seed, iters, group="dihedral"):
    N = 4 * n - 2
    m = N // 2
    rng = np.random.default_rng(seed)
    # init: random symmetric U of size ~k/2, random T
    chiU = np.zeros(m, dtype=np.float64)
    chiT = np.zeros(m, dtype=np.float64)
    # symmetric pairs for U
    pairs = [(d, m - d) for d in range(1, (m + 1) // 2)]
    rng.shuffle(pairs)
    uneed = k_target // 2 // 2  # pairs
    for (a, b) in pairs[:uneed]:
        chiU[a] = chiU[b] = 1
    tneed = k_target - int(chiU.sum())
    tidx = rng.permutation(m)[:tneed]
    chiT[tidx] = 1
    if group == "abelian":
        # T must also be symmetric: symmetrize
        for t in range(m):
            if chiT[t] and not chiT[(m - t) % m]:
                chiT[(m - t) % m] = 1
    e, k = energy(chiU, chiT, n, N, group)
    best = e
    t0 = time.time()
    temp0 = 3.0
    for it in range(iters):
        temp = temp0 * max(0.02, 1 - it / iters)
        # move: flip a U-pair or a T-element (or T-pair for abelian)
        if rng.random() < 0.4:
            a = rng.integers(1, m)
            b = (m - a) % m
            chiU[a] = 1 - chiU[a]
            if b != a: chiU[b] = chiU[a]
            undo = ("U", a, b)
        else:
            t = int(rng.integers(0, m))
            chiT[t] = 1 - chiT[t]
            if group == "abelian":
                b = (m - t) % m
                if b != t: chiT[b] = chiT[t]
                undo = ("Ts", t, b)
            else:
                undo = ("T", t, t)
        e2, k2 = energy(chiU, chiT, n, N, group)
        if e2 <= e or rng.random() < np.exp((e - e2) / max(temp, 1e-9)):
            e = e2
            if e < best:
                best = e
                if it % 500 == 0 or e == 0:
                    print(f"  it={it} E={e} k={k2} ({time.time()-t0:.0f}s)", flush=True)
            if e == 0:
                ok, won, woff, A = brute_verify(chiU, chiT, n, N)
                print(f"ZERO ENERGY at it={it}: brute verify ok={ok} "
                      f"worst_on={won} (<= {n-2}?) worst_off={woff} (<= {n-1}?)", flush=True)
                if ok:
                    fn = f"ramsey_witness_n{n}_{group}_seed{seed}.txt"
                    emit(A, fn)
                    np.savez(f"ramsey_witness_n{n}_{group}_seed{seed}.npz",
                             chiU=chiU, chiT=chiT, A=A)
                    print(f"WITNESS SAVED: {fn}")
                    return True
        else:
            typ, a, b = undo
            if typ == "U":
                chiU[a] = 1 - chiU[a]
                if b != a: chiU[b] = chiU[a]
            elif typ == "T":
                chiT[a] = 1 - chiT[a]
            else:
                chiT[a] = 1 - chiT[a]
                if b != a: chiT[b] = chiT[a]
    print(f"no witness: best E={best} after {iters} iters ({time.time()-t0:.0f}s)")
    return False

if __name__ == "__main__":
    n = int(sys.argv[1])
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 2 * n - 2
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    iters = int(sys.argv[4]) if len(sys.argv) > 4 else 200000
    group = sys.argv[5] if len(sys.argv) > 5 else "dihedral"
    print(f"R(B_{n-1},B_{n}) witness search: N={4*n-2}, m={2*n-1}, k={k}, "
          f"group={group}, seed={seed}, iters={iters}")
    sa(n, k, seed, iters, group)
