"""
collatz_parity_mixed_levels_20261013.py -- the exact transfer operator on locally constant functions of (n mod 2^j, n mod 3^k)
and the closed form of its characteristic polynomial (THM-4629):
  det(I - M(s)) = (1 - eps_e a/2 - eps_o b/2) * prod_(i=1..k) (1 - (eps_e a/2)^(2 3^(i-1))),   independent of j,
a = e^(-s), b = e^(s(theta-1)), eps_e = zeta^alpha, eps_o = zeta^(alpha+beta).
Even branch: n = r mod 2^j (r even), n = rho mod 3^k  ->  n/2 = r/2 mod 2^(j-1) (two lifts mod 2^j), rho/2 mod 3^k.
Odd branch:  n = r mod 2^j (r odd), n = rho mod 3^k   ->  (3n+1)/2 = (3r+1)/2 mod 2^(j-1) (two lifts), c(rho) mod 3^k,
where c(rho) uses the odd representative of rho mod 2 3^k (the image is one class mod 3^(k+1); a depth-k function
sees only its reduction, which is why W_(j,k) is invariant).
"""
import sympy as sp, itertools
a, b = sp.symbols("a b")
chars = {"tau": (-1, -1), "omega": (1, -1), "C": (-1, 1)}

def M_jk(j, k, ee, eo):
    J = 2 ** j; K = 3 ** k; N = J * K
    inv2 = pow(2, -1, K) if k > 0 else 0
    M = sp.zeros(N, N)
    idx = lambda r, rho: r * K + rho
    for r in range(J):
        for rho in range(K):
            # even/odd is decided by r if j >= 1, else half-half
            if j == 0:
                M[idx(r, rho), idx(0, (rho * inv2) % K if k else 0)] += ee * a / 2
                n = rho if rho % 2 == 1 else rho + K
                M[idx(r, rho), idx(0, ((3 * n + 1) // 2) % K if k else 0)] += eo * b / 2
                continue
            if r % 2 == 0:
                s = (r // 2) % (J // 2) if j > 1 else 0
                tgt = (rho * inv2) % K if k else 0
                for lift in ((s, s + J // 2) if j > 1 else (0, 1)):
                    M[idx(r, rho), idx(lift % J, tgt)] += ee * a / 2
            else:
                s = ((3 * r + 1) // 2) % (J // 2) if j > 1 else 0
                n = rho if rho % 2 == 1 else rho + K      # odd representative of rho mod 2 3^k
                tgt = ((3 * n + 1) // 2) % K if k else 0
                for lift in ((s, s + J // 2) if j > 1 else (0, 1)):
                    M[idx(r, rho), idx(lift % J, tgt)] += eo * b / 2
    return M

def predicted(k, ee, eo):
    p = 1 - ee * a / 2 - eo * b / 2
    for i in range(1, k + 1):
        p *= (1 - (ee * a / 2) ** (2 * 3 ** (i - 1)))
    return sp.expand(p)

for name, (ee, eo) in chars.items():
    for (j, k) in [(0, 1), (1, 0), (2, 0), (3, 0), (1, 1), (2, 1), (1, 2), (3, 1), (2, 2), (0, 2)]:
        M = M_jk(j, k, ee, eo)
        det = sp.expand((sp.eye(M.shape[0]) - M).det())
        ok = sp.simplify(det - predicted(k, ee, eo)) == 0
        print(f"{name:5s} (j,k)=({j},{k}) dim {M.shape[0]:3d}: closed form {'MATCHES' if ok else 'DIFFERS'}" + ("" if ok else f"  det = {sp.factor(det)}"), flush=True)
print("done")
