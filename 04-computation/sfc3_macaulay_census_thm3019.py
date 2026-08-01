"""SFC(3) unbounded: Macaulay surjectivity census beyond the repo box.

f = a z^p + b z^q + c z^r.  L(f^m) = M_m(a,b,c) is the FORM of degree m
   M_m = sum_{i+j+l=m} m!/(i! j! l!) (pi+qj+rl)!  a^i b^j c^l.
SFC(3) at window k says M_{k+1}, M_{k+2}, M_{k+3} have no common zero in P^2.

CERTIFICATE (one-way, sound): if the Macaulay map
   S_{D-d_1} + S_{D-d_2} + S_{D-d_3}  ->  S_D,   D = sum(d_i - 1) + 1,
is SURJECTIVE over a field, then the ideal contains all of S_D, so the
projective variety is empty over the algebraic closure.  Surjectivity mod a
prime implies surjectivity over Q (rank can only drop under reduction).
"""
from math import comb, factorial
from itertools import product

P = (1 << 61) - 1     # prime


def monos(d):
    return [(i, j, d - i - j) for i in range(d + 1) for j in range(d - i + 1)]


def form_M(m, p, q, r, mod=P):
    """coefficients of M_m as dict monomial->coeff mod prime"""
    out = {}
    fm = factorial(m)
    for i in range(m + 1):
        for j in range(m + 1 - i):
            l = m - i - j
            c = (fm // (factorial(i) * factorial(j) * factorial(l))) % mod
            c = c * (factorial(p * i + q * j + r * l) % mod) % mod
            out[(i, j, l)] = c
    return out


def macaulay_surjective(p, q, r, k, mod=P):
    ds = [k + 1, k + 2, k + 3]
    D = sum(d - 1 for d in ds) + 1
    tgt = monos(D)
    tidx = {t: i for i, t in enumerate(tgt)}
    rows = []
    for d in ds:
        F = form_M(d, p, q, r, mod)
        for mu in monos(D - d):
            row = [0] * len(tgt)
            for nu, c in F.items():
                if c:
                    key = (mu[0] + nu[0], mu[1] + nu[1], mu[2] + nu[2])
                    row[tidx[key]] = (row[tidx[key]] + c) % mod
            rows.append(row)
    # rank over F_mod by Gaussian elimination
    ncols = len(tgt)
    rank = 0
    piv_col = 0
    R = rows
    nrows = len(R)
    r_i = 0
    for col in range(ncols):
        piv = None
        for i in range(r_i, nrows):
            if R[i][col]:
                piv = i
                break
        if piv is None:
            continue
        R[r_i], R[piv] = R[piv], R[r_i]
        inv = pow(R[r_i][col], mod - 2, mod)
        R[r_i] = [x * inv % mod for x in R[r_i]]
        for i in range(nrows):
            if i != r_i and R[i][col]:
                f = R[i][col]
                R[i] = [(R[i][t] - f * R[r_i][t]) % mod for t in range(ncols)]
        r_i += 1
        rank += 1
        if rank == ncols:
            break
    return rank == ncols, rank, ncols, D


print("SFC(3) Macaulay census (surjectivity => no common projective zero)")
print("repo box was p<q<r <= 9 with k <= 6 (THM-2836); extending r and k here.\n")
fails, tested = [], 0
for k in range(0, 5):
    bad_k = 0
    for r in range(2, 13):
        for q in range(1, r):
            for p in range(0, q):
                ok, rank, nc, D = macaulay_surjective(p, q, r, k)
                tested += 1
                if not ok:
                    fails.append((p, q, r, k, rank, nc))
                    bad_k += 1
    print(f"  window k={k}: supports p<q<r<=12 done  (cumulative {tested} cells, "
          f"non-surjective {len(fails)})", flush=True)
print(f"\n  TOTAL {tested} cells with 0 <= p < q < r <= 12 and k <= 4")
print(f"  non-surjective (i.e. certificate failed) cells: "
      f"{fails[:6] if fails else 'NONE'}{' ...' if len(fails) > 6 else ''}")
