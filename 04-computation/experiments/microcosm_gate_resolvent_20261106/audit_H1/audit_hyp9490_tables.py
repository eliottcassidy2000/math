# audit_H1: recompute HYP-9490 evidence (large-N deviations vs the K1 tower to N = 1201, all-classes c3 law at N ~ 401, rate table) with U_(d,r) from the exponential walk W_r (transpose route) and cover polynomials from subset convolution; also tests the kappa = 0 case of the Consequence.
import math
from fractions import Fraction
from audit_tournament_lib import classes, transitive, c3, cover_poly

RMAX = 600
# W_r(k): M^r e^y = sum_k W_r(k) e^(ky)
W = {1: 1}
Wlist = [dict(W)]
for r in range(RMAX):
    new = {}
    for k, c in W.items():
        if k == 0:
            continue
        new[k + 1] = new.get(k + 1, 0) + k * c
        new[k] = new.get(k, 0) + k * k * c
        new[k - 1] = new.get(k - 1, 0) - k * c
    W = new
    Wlist.append(W)


def U(d, r):
    return sum(c * (k ** d if d > 0 else 1) for k, c in Wlist[r].items())


def tower_H(co, r):
    n = len(co) - 1
    return sum(co[n - 2 * j] * U(n - 2 * j, r) for j in range(n // 2 + 1))


def circ(n, S):
    return [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]


cores = {
    'C3': circ(3, {1}), 'TT3': transitive(3), 'R5': circ(5, {1, 2}), 'TT5': transitive(5),
    'QR7': circ(7, {1, 2, 4}), 'TT7': transitive(7),
}
cop = {k: cover_poly(A) for k, A in cores.items()}
print("cover polynomials (low -> high):", {k: v for k, v in cop.items()})
print("kappa = c3 - c3(K1 tower of same order):",
      {k: c3(A) - ((len(A) - 1) // 2) ** 2 for k, A in cores.items()})
for N in (101, 201, 401, 801, 1201):
    ref = U(1, (N - 1) // 2)
    row = []
    for name in ('C3', 'TT3', 'R5', 'TT5', 'QR7', 'TT7'):
        n = len(cores[name])
        dev = Fraction(tower_H(cop[name], (N - n) // 2), ref) - 1
        row.append(f"{name} {float(dev):+.3e} [N^2 dev {float(dev) * N * N:+.1f}]")
    print(f"N = {N:4d}: " + " | ".join(row))

# rate table
print("rate: rho_(3,r) = U_(1,r)/U_(3,r), N = 3 + 2r")
for N in (23, 43, 83, 123, 163, 203, 243, 401, 801, 1201):
    r = (N - 3) // 2
    rho = Fraction(U(1, r), U(3, r))
    dev = Fraction(U(3, r), U(3, r) + 2 * U(1, r)) - 1    # TT3 against the K1 tower (= C3 tower)
    print(f"  N = {N:4d}: rho = {float(rho):.4e}, N^2 rho/log^2 N = {float(rho) * N * N / math.log(N) ** 2:.3f}, "
          f"|dev|/2 = {float(-dev) / 2:.4e}, 2 N^2 rho = {2 * float(rho) * N * N:.1f}")

# all classes, orders 4..7
print("c3 first-order law over all classes (reference: TT2 tower for even n, K1 tower for odd n)")
for n, N in ((4, 402), (5, 401), (6, 402), (7, 401)):
    r = (N - n) // 2
    if n % 2:
        refco, c3ref = [0, 1], ((n - 1) // 2) ** 2
        refH = U(1, (N - 1) // 2)
    else:
        refco = [0, 0, 1]
        c3ref = sum(2 + 2 * i for i in range((n - 2) // 2))   # c3(Gamma^m TT2) = sum_{i<m} (2 + 2i)
        refH = U(2, (N - 2) // 2)
    rows = []
    zero_rows = []
    for A in classes(n):
        co = cover_poly(A)
        dev = Fraction(tower_H(co, r), refH) - 1
        kap = c3(A) - c3ref
        rows.append((kap, float(dev)))
        if kap == 0:
            zero_rows.append((float(dev), co))
    unit = [d / k for k, d in rows if k != 0]
    eps = sum(unit) / len(unit)
    resid = max(abs(d - k * eps) for k, d in rows)
    nz = [d for d, co in zero_rows if d != 0]
    print(f"  n = {n} ({len(rows)} classes), N = {N}: c3_ref = {c3ref}; dev/kappa in [{min(unit):.4e}, {max(unit):.4e}], "
          f"eps = {eps:.4e}, max |dev - kappa eps| = {resid:.2e}; kappa = 0 classes: {len(zero_rows)}, "
          f"of which {len(nz)} have dev != 0 (max |dev| = {max([abs(x) for x in nz], default=0):.2e})")
