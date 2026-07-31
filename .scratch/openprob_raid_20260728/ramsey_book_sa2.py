#!/usr/bin/env python3
"""
ramsey_book_sa2.py — generalized two-level construction + stronger annealing
for R(B_{n-1}, B_n) >= 4n-1 witnesses on N = 4n-2 = 2m vertices, m = 2n-1.

Model: vertices F_m x {0,1}.
  (x,0) ~ (y,0)  iff  y-x in U0   (U0 symmetric subset of Z_m\{0})
  (x,1) ~ (y,1)  iff  y-x in U1   (U1 symmetric)
  (x,0) ~ (y,1)  iff  y-x in D    (D arbitrary subset of Z_m)
Contains dihedral (U0=U1) and Z_2xZ_m Cayley (U0=U1, D symmetric) models.

Exact codegree conditions (n-book-free both sides):
  within level L, diff d != 0:  W_L(d) = A_{UL}(d) + A_D(d)
       d in UL:  W_L(d) <= n-2      else  W_L(d) <= n+1+2k_L'-N  (see code)
  cross pair, diff t: X(t) = |U0 ∩ (t-D)| + |D ∩ (t+U1)|
       t in D:   X(t) <= n-2        else  X(t) <= ...
Degrees: deg0 = |U0|+|D|, deg1 = |U1|+|D| (not forced equal; limits use both).

For non-regular graphs the off-edge limit depends on the pair's degrees:
common non-nbrs((a,b)) = N - 2 - deg(a) - deg(b) + codeg(a,b)  <= n-1
  => codeg <= n+1 + deg(a)+deg(b) - N.
Annealing with swap moves (size-preserving), flip moves, restarts, tabu-lite.

Usage: python3 ramsey_book_sa2.py n [iters] [seed] [restarts]
"""
import numpy as np
import sys, time

def acorr(chi):
    m = len(chi)
    f = np.fft.rfft(chi)
    return np.rint(np.fft.irfft(f * np.conj(f), n=m)).astype(np.int64)

def xcorr_conv(a, b):   # sum_u a(t-u) b(u)
    m = len(a)
    return np.rint(np.fft.irfft(np.fft.rfft(a) * np.fft.rfft(b), n=m)).astype(np.int64)

def xcorr_corr(a, b):   # sum_u a(t+u) b(u)
    m = len(a)
    return np.rint(np.fft.irfft(np.fft.rfft(a) * np.conj(np.fft.rfft(b)), n=m)).astype(np.int64)

def energy(U0, U1, D, n, N):
    m = len(U0)
    k0 = int(U0.sum() + D.sum()); k1 = int(U1.sum() + D.sum())
    lim_on = n - 2
    W0 = acorr(U0) + acorr(D)
    W1 = acorr(U1) + acorr(D)
    X = xcorr_conv(U0, D)*0  # placeholder replaced below
    # X(t) = |U0 ∩ (t-D)| + |D ∩ (t+U1)| = sum_u U0(t-u)D(u) + sum_u D(t+u)U1(u)
    X = xcorr_conv(U0, D) + xcorr_corr(D, U1)
    e = 0
    off0 = n + 1 + 2*k0 - N
    off1 = n + 1 + 2*k1 - N
    offX = n + 1 + k0 + k1 - N
    d = np.arange(1, m)
    lim0 = np.where(U0[1:] > 0, lim_on, off0)
    lim1 = np.where(U1[1:] > 0, lim_on, off1)
    e += int(np.maximum(W0[1:] - lim0, 0).sum())
    e += int(np.maximum(W1[1:] - lim1, 0).sum())
    limX = np.where(D > 0, lim_on, offX)
    e += int(np.maximum(X - limX, 0).sum())
    return e

def brute_verify(U0, U1, D, n, N):
    m = len(U0)
    A = np.zeros((N, N), dtype=np.int8)
    for i in range(m):
        for j in range(m):
            dd = (j - i) % m
            if i != j:
                if U0[dd]: A[i, j] = 1
                if U1[dd]: A[m+i, m+j] = 1
            if D[dd]:
                A[i, m+j] = 1; A[m+j, i] = 1
    assert (A == A.T).all() and (np.diag(A) == 0).all()
    M = A.astype(np.int64)
    C = M @ M
    comp = 1 - M; np.fill_diagonal(comp, 0)
    Cc = comp @ comp
    okon = True; okoff = True; won = -1; woff = -1
    for i in range(N):
        for j in range(i+1, N):
            if M[i, j]:
                won = max(won, C[i, j])
                if C[i, j] > n - 2: okon = False
            else:
                woff = max(woff, Cc[i, j])
                if Cc[i, j] > n - 1: okoff = False
    return okon and okoff, won, woff, A

def qr_set(m):
    qr = np.zeros(m)
    for x in range(1, m):
        qr[(x*x) % m] = 1
    return qr

def run(n, iters, seed, restarts):
    N = 4*n - 2; m = 2*n - 1
    rng = np.random.default_rng(seed)
    bestE_global = 10**9
    for r in range(restarts):
        U0 = np.zeros(m); U1 = np.zeros(m); D = np.zeros(m)
        # seeding: D = QR (+0?), U0/U1 random symmetric with |U|=m-1-|QR| ~ m/2
        if r % 3 == 0:
            D = qr_set(m).astype(np.float64)
        else:
            idx = rng.permutation(m)[:m//2]
            D[idx] = 1
        for a in range(1, (m+1)//2):
            if rng.random() < 0.5:
                U0[a] = U0[m-a] = 1
            if rng.random() < 0.5:
                U1[a] = U1[m-a] = 1
        e = energy(U0, U1, D, n, N)
        best = e
        T0 = 2.5
        stall = 0
        for it in range(iters):
            T = T0 * max(0.01, 1 - it/iters)
            mv = rng.random()
            if mv < 0.25:
                a = rng.integers(1, m); b = (m - a) % m
                U0[a] = 1 - U0[a]; U0[b] = U0[a]
                undo = lambda: (U0.__setitem__(a, 1-U0[a]), U0.__setitem__(b, U0[a]-0))
                which = ("U0", a, b, None, None)
            elif mv < 0.5:
                a = rng.integers(1, m); b = (m - a) % m
                U1[a] = 1 - U1[a]; U1[b] = U1[a]
                which = ("U1", a, b, None, None)
            elif mv < 0.8:
                t = int(rng.integers(0, m))
                D[t] = 1 - D[t]
                which = ("D", t, None, None, None)
            else:
                # swap inside D (size-preserving): move a 1 to a 0
                ones = np.flatnonzero(D > 0); zers = np.flatnonzero(D == 0)
                if len(ones) and len(zers):
                    t1 = int(ones[rng.integers(len(ones))]); t0 = int(zers[rng.integers(len(zers))])
                    D[t1] = 0; D[t0] = 1
                    which = ("Ds", t1, t0, None, None)
                else:
                    continue
            e2 = energy(U0, U1, D, n, N)
            if e2 <= e or rng.random() < np.exp((e - e2)/max(T, 1e-9)):
                e = e2
                if e < best:
                    best = e; stall = 0
                if e == 0:
                    ok, won, woff, A = brute_verify(U0, U1, D, n, N)
                    print(f"[r{r}] ZERO at it={it}: verify={ok} won={won} woff={woff}", flush=True)
                    if ok:
                        fn = f"ramsey_witness_n{n}_seed{seed}_r{r}"
                        np.savez(fn + ".npz", U0=U0, U1=U1, D=D, A=A)
                        with open(fn + ".txt", "w") as f:
                            f.write("".join(str(A[i, j]) for j in range(N) for i in range(N)) + "\n")
                        print(f"WITNESS SAVED {fn}", flush=True)
                        return True
            else:
                typ, a, b, _, _ = which
                if typ == "U0":
                    U0[a] = 1 - U0[a]; U0[b] = U0[a]
                elif typ == "U1":
                    U1[a] = 1 - U1[a]; U1[b] = U1[a]
                elif typ == "D":
                    D[a] = 1 - D[a]
                else:
                    D[a] = 1; D[b] = 0
            stall += 1
        bestE_global = min(bestE_global, best)
        print(f"[r{r}] done bestE={best}", flush=True)
    print(f"NO WITNESS: best energy {bestE_global} over {restarts} restarts", flush=True)
    return False

if __name__ == "__main__":
    n = int(sys.argv[1])
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 300000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    restarts = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    print(f"n={n} N={4*n-2} m={2*n-1} iters={iters} seed={seed} restarts={restarts}")
    run(n, iters, seed, restarts)
