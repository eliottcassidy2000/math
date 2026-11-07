#!/usr/bin/env python3
# audit J1 (THM-4911 counts): gentourng-independent cross-check for n <= 7: labeled equicentric (resp. strong equicentric) tournaments = sum over equicentric classes of n!/|Aut(T)|.
import itertools, math, subprocess
import numpy as np

GT = "/opt/homebrew/bin/gentourng"
def classes(n):
    m = n * (n - 1) // 2; iu = np.triu_indices(n, 1)
    out = subprocess.run([GT, "-q", str(n)], capture_output=True).stdout
    a = np.frombuffer(out, dtype=np.uint8).reshape(-1, m + 1)[:, :m] - 48
    A = np.zeros((len(a), n, n), dtype=np.int8); A[:, iu[0], iu[1]] = a; A[:, iu[1], iu[0]] = 1 - a
    return A
def landau_strong(A):
    d = np.sort(A.sum(-1), axis=-1).astype(np.int64); cs = np.cumsum(d, axis=-1); n = A.shape[-1]; k = np.arange(1, n)
    return (cs[..., :-1] > k * (k - 1) // 2).all(-1)
for n in range(3, 8):
    m = n * (n - 1) // 2; iu = np.triu_indices(n, 1)
    codes = np.arange(2 ** m, dtype=np.int64)
    bits = ((codes[:, None] >> np.arange(m)[None, :]) & 1).astype(np.int8)
    A = np.zeros((len(codes), n, n), dtype=np.int8); A[:, iu[0], iu[1]] = bits; A[:, iu[1], iu[0]] = 1 - bits
    M = (A - A.transpose(0, 2, 1)).astype(np.int16); x = M.sum(2); y = np.einsum('bij,bj->bi', M, x)
    eq = (y == y[:, :1]).all(1); st = landau_strong(A)
    lab_eq, lab_st = int(eq.sum()), int((eq & st).sum())
    del A, M, bits
    C = classes(n); perms = [np.array(p) for p in itertools.permutations(range(n))]
    orb_eq = orb_st = 0; ncls = 0
    for T in C:
        Mt = (T - T.T).astype(np.int16); xt = Mt.sum(1); yt = Mt @ xt
        if not (yt == yt[0]).all(): continue
        ncls += 1
        aut = sum(1 for p in perms if (T[np.ix_(p, p)] == T).all())
        orb_eq += math.factorial(n) // aut
        if landau_strong(T[None])[0]: orb_st += math.factorial(n) // aut
    print(f"n={n}: labeled equicentric {lab_eq} vs sum n!/|Aut| over the {ncls} classes {orb_eq}: {lab_eq == orb_eq}; strong: {lab_st} vs {orb_st}: {lab_st == orb_st}")
