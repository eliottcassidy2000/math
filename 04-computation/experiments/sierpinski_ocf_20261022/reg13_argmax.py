# Which regular 13-tournament has the most Hamiltonian paths, and does the parameter-free two-loop predictor
# kappa_2 = rho (1 - y_pred/13) find it among all 1,495,297 classes?  Batched numpy (chunks of 50k); the top
# candidates are re-scored at second order (three loops, HYP-9367).  Usage: reg13_argmax.py codes H.
# mac-mini-2026-10-22-sierp.
import sys, math
import numpy as np
sys.path.insert(0, '.')
n = 13
codes = [l.strip() for l in open(sys.argv[1]) if l.strip()]
Hs = np.array([int(l) for l in open(sys.argv[2]) if l.strip()], dtype=np.int64)
assert len(codes) == len(Hs), (len(codes), len(Hs))
iu = np.triu_indices(n, 1)
kap = Hs * 2.0 ** (n - 1) / math.factorial(n)
rho_all = np.zeros(len(codes)); yp_all = np.zeros(len(codes))
CH = 50000
I = np.eye(n)
for st in range(0, len(codes), CH):
    ch = codes[st:st + CH]; N = len(ch)
    bits = np.frombuffer(''.join(ch).encode(), dtype=np.uint8).reshape(N, -1) - ord('0')
    A = np.zeros((N, n, n)); A[:, iu[0], iu[1]] = bits; A[:, iu[1], iu[0]] = 1 - bits
    Sh = (A - np.transpose(A, (0, 2, 1))) / n
    S2 = Sh @ Sh
    G = -S2 @ np.linalg.inv(I + S2); H = S2 @ np.linalg.inv(I - S2); O = Sh @ np.linalg.inv(I - S2)
    gd = np.einsum('kii->ki', G); hd = np.einsum('kii->ki', H); d = hd - gd
    F8 = np.sum(2 * gd * hd - 1.5 * gd ** 2, axis=1)
    TH = np.sum(G * (H * H - O * O - (2 / 3) * G * G), axis=(1, 2))
    DB = -np.einsum('ki,kij,kj->k', d, G, d)
    x = np.linalg.eigvals(Sh).imag
    D1 = np.sum(2 * x ** 2 / (1 - x ** 4), axis=1); D2 = np.sum(4 * x ** 2 * (1 + x ** 4) / (1 - x ** 4) ** 2, axis=1)
    EL = (D2 + D1 ** 2 - D1) / 2
    rho_all[st:st + N] = np.sqrt(np.prod((1 + x ** 2) / (1 - x ** 2), axis=1))
    yp_all[st:st + N] = -n * (F8 + TH + DB) - EL
k2 = rho_all * (1 - yp_all / n)
order_true = np.argsort(-kap); order_k2 = np.argsort(-k2); order_rho = np.argsort(-rho_all)
imax = order_true[0]
rank = lambda order, i: int(np.where(order == i)[0][0]) + 1
print(f"n=13: {len(codes)} regular classes; max H = {Hs[imax]} (kappa {kap[imax]:.6f}); R_13 has H = 3711175")
print(f"  rank of the true maximiser: by rho {rank(order_rho, imax)}, by two-loop kappa_2 {rank(order_k2, imax)}")
print(f"  spearman(kappa, rho) = {np.corrcoef(np.argsort(np.argsort(kap)), np.argsort(np.argsort(rho_all)))[0, 1]:+.4f}, spearman(kappa, kappa_2) = {np.corrcoef(np.argsort(np.argsort(kap)), np.argsort(np.argsort(k2)))[0, 1]:+.4f}")
print("  top 10 by H:  H, rank by kappa_2, rho, y_exact, y_pred")
for i in order_true[:10]:
    print(f"    {Hs[i]}  {rank(order_k2, i):7d}  {rho_all[i]:.6f}  {n * (1 - kap[i] / rho_all[i]):.4f}  {yp_all[i]:.4f}   {codes[i]}")
print("  top 10 by kappa_2: H, true rank")
for i in order_k2[:10]:
    print(f"    {Hs[i]}  {rank(order_true, i):7d}")
# second-order re-scoring of the top 300 by kappa_2
from second_order import predict2
top = order_k2[:300]
k3 = []
for i in top:
    A = np.zeros((n, n)); b = np.array([1 if ch == '1' else 0 for ch in codes[i]]); A[iu[0], iu[1]] = b; A[iu[1], iu[0]] = 1 - b
    rho, first, second, X2 = predict2(A - A.T)
    k3.append(rho * (1 + first + second))
k3 = np.array(k3)
best = top[int(np.argmax(k3))]
print(f"  second-order re-scoring of the top 300 by kappa_2: its argmax has H = {Hs[best]} (true rank {rank(order_true, best)}); "
      f"true maximiser in this top-300: {imax in set(top.tolist())}")
o3 = top[np.argsort(-k3)]
print("  top 10 after second order: H, true rank:", [(int(Hs[i]), rank(order_true, i)) for i in o3[:10]])
