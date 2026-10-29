#!/usr/bin/env python3
"""Independent audit of THM-4848: run observers M_k of circulant tournaments, their character blocks, the
descent of det B_1 to K_H, norms, Smith data and the converse conjugation.

Own implementation (nothing imported from the session's helpers):
  - P[U][a][b] = number of Hamiltonian paths of T[U] from a to b (single vertex: 1);
  - M_k[s,t] = sum over ordered set partitions (U_1..U_k) of V with s_i, t_i in U_i of prod_i P[U_i][s_i][t_i];
  - representatives R = tuples with first coordinate 0; B_j[r,r'] = sum_s M_k[r, tau^s r'] zeta^(j s);
  - det B_1 in Z[zeta_q] by multi-modular evaluation at primes l = 1 mod q and CRT (power-basis coefficients),
    verified on extra primes.
Outputs: a3_runobs.out (text) and matrices for PARI in audit/mats/.
"""
import os, sys, json, random
from itertools import product, permutations
from math import gcd
from sympy import factorint, isprime, nextprime, resultant, symbols, Poly, cyclotomic_poly, legendre_symbol, jacobi_symbol

HERE = os.path.dirname(os.path.abspath(__file__))
MATS = os.path.join(HERE, "mats")
os.makedirs(MATS, exist_ok=True)

def say(*a):
    print(" ".join(str(x) for x in a), flush=True)

# ------------------------------------------------------------------------------------------------ observers
def path_table(n, adj):
    """P[mask] = dict (a,b) -> number of Hamiltonian paths of T[mask] from a to b"""
    P = [dict() for _ in range(1 << n)]
    for a in range(n):
        dp = [None] * (1 << n)
        dp[1 << a] = {a: 1}
        # iterate masks in increasing order (supersets come later)
        for mask in range(1 << n):
            d = dp[mask]
            if not d:
                continue
            for v, c in d.items():
                P[mask][(a, v)] = P[mask].get((a, v), 0) + c
                for w in range(n):
                    if not (mask >> w) & 1 and adj[v][w]:
                        m2 = mask | (1 << w)
                        if dp[m2] is None:
                            dp[m2] = {}
                        dp[m2][w] = dp[m2].get(w, 0) + c
    return P

def tuples(n, k):
    return list(permutations(range(n), k))

def observer(n, adj, k):
    P = path_table(n, adj)
    tup = tuples(n, k)
    idx = {t: i for i, t in enumerate(tup)}
    N = len(tup)
    M = [[0] * N for _ in range(N)]
    full = (1 << n) - 1
    # ordered set partitions via surjective block assignment
    for f in product(range(k), repeat=n):
        masks = [0] * k
        for v, b in enumerate(f):
            masks[b] |= 1 << v
        if any(m == 0 for m in masks):
            continue
        lists = [list(P[m].items()) for m in masks]
        for combo in product(*lists):
            s = tuple(c[0][0] for c in combo)
            t = tuple(c[0][1] for c in combo)
            val = 1
            for c in combo:
                val *= c[1]
            M[idx[s]][idx[t]] += val
    return M, tup, idx

def observer2_fast(n, adj):
    """k = 2 only, faster: loop over subsets"""
    P = path_table(n, adj)
    tup = tuples(n, 2)
    idx = {t: i for i, t in enumerate(tup)}
    N = len(tup)
    M = [[0] * N for _ in range(N)]
    full = (1 << n) - 1
    for U in range(1, full):
        A = list(P[U].items())
        B = list(P[full ^ U].items())
        for (s1, t1), c1 in A:
            for (s2, t2), c2 in B:
                M[idx[(s1, s2)]][idx[(t1, t2)]] += c1 * c2
    return M, tup, idx


def partial_factor(n, lim=10 ** 7):
    """trial division to lim, then a few rounds of Pollard rho with a step cap; reports a composite cofactor as 'c<digits>'"""
    from sympy import factorint, isprime
    f = factorint(n, limit=lim, use_rho=False, use_pm1=False, use_ecm=False)
    out = {}
    for p, e in f.items():
        if p > lim and not isprime(p):
            key = f"c{len(str(p))}"
            out[key] = out.get(key, 0) + e
        else:
            out[p] = out.get(p, 0) + e
    return out

def circulant(S, n):
    return [[1 if ((j - i) % n) in S else 0 for j in range(n)] for i in range(n)]

# ------------------------------------------------------------------------------------------------ blocks
def blocks_data(M, tup, idx, n):
    R = [t for t in tup if t[0] == 0]
    def shift(t, s):
        return tuple((x + s) % n for x in t)
    # Mshift[r][r'][s] = M[r, tau^s r']
    Ms = [[[M[idx[r]][idx[shift(rp, s)]] for s in range(n)] for rp in R] for r in R]
    return R, Ms

def det_mod(A, l):
    A = [[x % l for x in row] for row in A]
    n = len(A)
    det = 1
    for i in range(n):
        piv = None
        for r in range(i, n):
            if A[r][i]:
                piv = r
                break
        if piv is None:
            return 0
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            det = -det
        det = det * A[i][i] % l
        inv = pow(A[i][i], -1, l)
        rowi = A[i]
        for r in range(i + 1, n):
            f = A[r][i]
            if f:
                f = f * inv % l
                rowr = A[r]
                for c in range(i, n):
                    rowr[c] = (rowr[c] - f * rowi[c]) % l
    return det % l

def bareiss_det(A):
    A = [row[:] for row in A]
    n = len(A)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = None
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    sw = r
                    break
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]

def primitive_root_of_unity(q, l):
    # element of exact order q in F_l (l = 1 mod q)
    while True:
        g = random.randrange(2, l - 1)
        r = pow(g, (l - 1) // q, l)
        if r != 1:
            return r

def solve_mod(A, b, l):
    n = len(A)
    M = [[x % l for x in A[i]] + [b[i] % l] for i in range(n)]
    for i in range(n):
        piv = next(r for r in range(i, n) if M[r][i])
        M[i], M[piv] = M[piv], M[i]
        inv = pow(M[i][i], -1, l)
        M[i] = [x * inv % l for x in M[i]]
        for r in range(n):
            if r != i and M[r][i]:
                f = M[r][i]
                M[r] = [(x - f * y) % l for x, y in zip(M[r], M[i])]
    return [M[i][n] for i in range(n)]

def detB1_coeffs(Ms, q, nprimes_extra=2, seed=5):
    """power-basis coefficients (length q-1) of det B_1 = det(sum_s M[r,tau^s r'] zeta^s) in Z[zeta_q]"""
    random.seed(seed)
    d = q - 1
    l = 2 ** 61
    mods = []
    residues = []
    prev = None
    stable = 0
    while True:
        l = nextprime(l)
        while (l - 1) % q:
            l = nextprime(l)
        rho = primitive_root_of_unity(q, l)
        vals = []
        for j in range(1, q):
            z = pow(rho, j, l)
            pw = [pow(z, s, l) for s in range(q)]
            B = [[sum(x * y for x, y in zip(row_rp, pw)) % l for row_rp in row] for row in Ms]
            vals.append(det_mod(B, l))
        # solve sum_i c_i (rho^j)^i = vals_j
        A = [[pow(rho, i * j, l) for i in range(d)] for j in range(1, q)]
        c = solve_mod(A, vals, l)
        mods.append(l)
        residues.append(c)
        # CRT
        Mtot = 1
        for m in mods:
            Mtot *= m
        coeffs = []
        for i in range(d):
            x = 0
            for m, r in zip(mods, residues):
                Mi = Mtot // m
                x += r[i] * Mi * pow(Mi, -1, m)
            x %= Mtot
            if x > Mtot // 2:
                x -= Mtot
            coeffs.append(x)
        if coeffs == prev:
            stable += 1
            if stable >= nprimes_extra:
                return coeffs, len(mods)
        else:
            stable = 0
        prev = coeffs

def poly_mod_cyclo(c, q):
    """reduce sum c_i x^i (any length) modulo Phi_q = 1 + x + ... + x^(q-1); return length q-1"""
    a = [0] * q
    for i, ci in enumerate(c):
        a[i % q] += ci
    # x^(q-1) = -(1 + ... + x^(q-2))
    top = a[q - 1]
    return [a[i] - top for i in range(q - 1)]

def sigma(c, m, q):
    out = [0] * q
    for i, ci in enumerate(c):
        out[(i * m) % q] += ci
    return poly_mod_cyclo(out, q)

def norm_Qzeta(c, q):
    x = symbols('x')
    f = Poly(list(reversed(c)), x) if any(c) else Poly(0, x)
    return int(resultant(Poly(cyclotomic_poly(q, x), x), f))

def eta_coords(c, q, H):
    """write an element of K_H = Q(sqrt(-q)) (H = QR) as a + b*eta, eta = sum_{h in QR} zeta^h"""
    eta = [0] * q
    for h in H:
        eta[h] += 1
    eta = poly_mod_cyclo(eta, q)
    # c = a*1 + b*eta in power basis; eta has coefficient 1 at each h in QR (if QR contains no q-1) ...
    # solve with two coordinates: use coefficient of x^0 and of x^h0 for some h0 in QR
    # general: c - b*eta must be a constant
    h0 = min(h for h in H if h != q - 1)
    b = c[h0] - (c[0] * 0)
    # c[h0] = a*0 + b*eta[h0] + (a contributes only at index 0 unless reduction) -> handle generally:
    # a*[1,0,..] + b*eta ; solve from indices 0 and h0
    # eta[0] may be nonzero after reduction if q-1 in H
    import fractions
    A = [[1, eta[0]], [0, eta[h0]]]
    rhs = [c[0], c[h0]]
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    a = fractions.Fraction(rhs[0] * A[1][1] - A[0][1] * rhs[1], det)
    bb = fractions.Fraction(A[0][0] * rhs[1] - A[1][0] * rhs[0], det)
    assert a.denominator == 1 and bb.denominator == 1
    a, bb = int(a), int(bb)
    recon = [a * (1 if i == 0 else 0) + bb * eta[i] for i in range(q - 1)]
    assert recon == list(c), "element not in Z[eta]"
    return a, bb

def quad_norm(a, b, q):
    # eta^2 + eta + (q+1)/4 = 0
    return a * a - a * b + ((q + 1) // 4) * b * b

def split_valuations(a, b, p, q, e_tot):
    """(v_P, v_Pbar) of a + b*eta at the two primes above a split p; P = (p, eta - r) with r the smaller root mod p"""
    E = e_tot + 2
    mod = p ** E
    c = (q + 1) // 4
    from sympy.ntheory.residue_ntheory import sqrt_mod
    if p == 2:
        roots = [r for r in (0, 1) if (r * r + r + c) % 2 == 0]
    else:
        s = sqrt_mod((-q) % p, p)
        inv2 = pow(2, -1, p)
        roots = sorted({((-1 + s) * inv2) % p, ((-1 - s) * inv2) % p})
    assert all((r * r + r + c) % p == 0 for r in roots)
    res = []
    for r0 in roots:
        # Hensel lift r0 to a root of X^2 + X + c mod p^E
        r = r0
        for _ in range(E + 1):
            f = (r * r + r + c) % mod
            fp = (2 * r + 1) % mod
            r = (r - f * pow(fp, -1, mod)) % mod
        val = (a + b * r) % mod
        v = 0
        while val % p == 0 and v < E:
            val //= p
            v += 1
        res.append(v)
    return tuple(res)

def kron(D, p):
    return jacobi_symbol(D % p, p) if p != 2 else (0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1))

# ------------------------------------------------------------------------------------------------ local Smith
def local_smith_py(A, p, E):
    """elementary divisor exponents (< E) of integer matrix A over Z_(p), working mod p^E; returns dict e->mult and
    the number of divisors with exponent >= E"""
    mod = p ** E
    A = [[x % mod for x in row] for row in A]
    n, m = len(A), len(A[0])
    res = {}
    rnk = 0
    rows = list(range(n))
    def val(x):
        if x == 0:
            return E
        v = 0
        while x % p == 0:
            x //= p
            v += 1
        return v
    for i in range(min(n, m)):
        best = None
        for r in range(i, n):
            row = A[r]
            for c in range(i, m):
                if row[c]:
                    vv = val(row[c])
                    if best is None or vv < best[0]:
                        best = (vv, r, c)
                        if vv == 0:
                            break
            if best is not None and best[0] == 0:
                break
        if best is None:
            break
        vv, r, c = best
        A[i], A[r] = A[r], A[i]
        for row in A:
            row[i], row[c] = row[c], row[i]
        piv = A[i][i]
        unit = piv // p ** vv
        inv = pow(unit, -1, mod)
        A[i] = [x * inv % mod for x in A[i]]
        # now A[i][i] = p^vv; eliminate column below (entries divisible by p^vv)
        for r2 in range(i + 1, n):
            if A[r2][i]:
                f = A[r2][i] // p ** vv
                A[r2] = [(x - f * y) % mod for x, y in zip(A[r2], A[i])]
        # eliminate row to the right
        for c2 in range(i + 1, m):
            if A[i][c2]:
                f = A[i][c2] // p ** vv
                for r2 in range(i, n):
                    A[r2][c2] = (A[r2][c2] - f * A[r2][i]) % mod
        if vv > 0:
            res[vv] = res.get(vv, 0) + 1
    zero_rank = min(n, m) - sum(1 for _ in range(0))
    return res

def local_smith_np(A, p, E):
    """local Smith form over Z_(p) (exponents < E): phase 1 eliminates with unit pivots mod p^E (numpy, mod < 2^31),
    phase 2 runs the exact routine on the residual block, which is divisible by p."""
    import numpy as np
    mod = p ** E
    assert mod < 2 ** 31
    M = np.array([[x % mod for x in row] for row in A], dtype=np.int64)
    n = M.shape[0]
    i = 0
    while i < n:
        sub = M[i:, i:]
        mask = (sub % p) != 0
        if not mask.any():
            break
        flat = int(np.argmax(mask))
        r, c = divmod(flat, sub.shape[1])
        r += i; c += i
        if r != i:
            M[[i, r], :] = M[[r, i], :]
        if c != i:
            M[:, [i, c]] = M[:, [c, i]]
        inv = pow(int(M[i, i]), -1, mod)
        M[i, i:] = (M[i, i:] * inv) % mod
        col = M[i + 1:, i].copy()
        if col.any():
            M[i + 1:, i:] = (M[i + 1:, i:] - (np.outer(col, M[i, i:]) % mod)) % mod
        M[i, i + 1:] = 0  # column operations clear the pivot row (unit pivot)
        i += 1
    rest = M[i:, i:].tolist()
    if not rest:
        return {}, 0, i
    res = local_smith_py(rest, p, E)
    found = sum(res.values())
    return res, len(rest) - found, i

def rank_mod_np(A, l):
    import numpy as np
    M = np.array([[x % l for x in row] for row in A], dtype=np.int64)
    n, m = M.shape
    r = 0
    for c in range(m):
        piv = None
        nzr = np.nonzero(M[r:, c])[0]
        if len(nzr) == 0:
            continue
        pr = r + nzr[0]
        M[[r, pr], :] = M[[pr, r], :]
        inv = pow(int(M[r, c]), -1, l)
        M[r, :] = (M[r, :] * inv) % l
        col = M[r + 1:, c].copy()
        M[r + 1:, :] = (M[r + 1:, :] - np.outer(col, M[r, :]) % l) % l
        r += 1
        if r == n:
            break
    return r

def det_mod_np(A, l):
    import numpy as np
    M = np.array([[x % l for x in row] for row in A], dtype=np.int64)
    n = M.shape[0]
    det = 1
    for i in range(n):
        nzr = np.nonzero(M[i:, i])[0]
        if len(nzr) == 0:
            return 0
        pr = i + nzr[0]
        if pr != i:
            M[[i, pr], :] = M[[pr, i], :]
            det = -det
        det = det * int(M[i, i]) % l
        inv = pow(int(M[i, i]), -1, l)
        M[i, :] = (M[i, :] * inv) % l
        col = M[i + 1:, i].copy()
        M[i + 1:, :] = (M[i + 1:, :] - np.outer(col, M[i, :]) % l) % l
    return det % l

def write_gp(M, name):
    path = os.path.join(MATS, name + ".gp")
    with open(path, "w") as f:
        f.write("M = [" + ";".join(",".join(str(x) for x in row) for row in M) + "];\n")
    return path

# ------------------------------------------------------------------------------------------------ main
def analyse(label, S, n, k, H, paley=False, save=True, report_norm=True):
    adj = circulant(S, n)
    if k == 2:
        M, tup, idx = observer2_fast(n, adj)
    else:
        M, tup, idx = observer(n, adj, k)
    N = len(tup)
    R, Ms = blocks_data(M, tup, idx, n)
    B0 = [[sum(row_rp) for row_rp in row] for row in Ms]
    say(f"== {label}: q={n}, S={sorted(S)}, k={k}, M_k is {N}x{N}, |R| = {len(R)}, total mass {sum(map(sum, M))}")
    detB0 = bareiss_det(B0)
    say(f"   det B_0 = {detB0} = {factorint(detB0) if detB0 else 0}")
    coeffs, used = detB1_coeffs(Ms, n)
    say(f"   det B_1 power-basis coefficients (CRT over {used} primes l = 1 mod {n}): {coeffs}")
    inv_ok = all(sigma(coeffs, m, n) == coeffs for m in H)
    say(f"   sigma_m(det B_1) = det B_1 for all m in H = {sorted(H)}: {inv_ok}")
    others = [m for m in range(1, n) if m not in H]
    moved = any(sigma(coeffs, m, n) != coeffs for m in others) if others else None
    say(f"   moved by some sigma outside H: {moved}")
    NB1 = norm_Qzeta(coeffs, n)
    say(f"   N(det B_1) = {NB1} = {partial_factor(NB1) if NB1 else 0}")
    res = {"M": M, "B0": B0, "detB0": detB0, "coeffs": coeffs, "N": NB1, "tup": tup, "idx": idx, "R": R, "Ms": Ms}
    if paley:
        a, b = eta_coords(coeffs, n, H)
        nq = quad_norm(a, b, n)
        res["a"], res["b"], res["nq"] = a, b, nq
        say(f"   beta = a + b eta: a = {a}, b = {b}")
        say(f"   quadratic norm a^2 - ab + {(n + 1) // 4} b^2 = {nq} = {factorint(nq)}")
        say(f"   N(det B_1) == quadratic norm^{len(H)}: {NB1 == nq ** len(H)}")
        g = gcd(a, b)
        say(f"   content gcd(a,b) = {g} = {factorint(g)}")
        for p, e in sorted(factorint(nq).items()):
            kr = kron(-n, p)
            typ = {1: "split", -1: "inert", 0: "ramified"}[kr]
            extra = ""
            if kr == 1:
                vv = split_valuations(a, b, p, n, e)
                extra = f" (v_P, v_Pbar) = {vv}, sum {sum(vv)}, {'one-sided' if min(vv) == 0 else 'TWO-SIDED'}"
            say(f"     p = {p}: exponent {e} in quadratic norm, {typ}, p mod q = {p % n}{extra}")
    if save:
        write_gp(M, label.replace(" ", "_") + "_M")
        write_gp(B0, label.replace(" ", "_") + "_B0")
    return res

def main():
    which = sys.argv[1:] or ["P7_2", "NQR7_2", "R7_2", "P7_3", "NQR7_3", "P11_2", "NQR11_2", "R11_2"]
    QR7, NQR7 = {1, 2, 4}, {3, 5, 6}
    QR11 = {1, 3, 4, 5, 9}; NQR11 = {2, 6, 7, 8, 10}
    results = {}
    for w in which:
        if w == "P7_2":
            results[w] = analyse("P7_k2", QR7, 7, 2, QR7, paley=True)
        elif w == "NQR7_2":
            results[w] = analyse("NQR7_k2", NQR7, 7, 2, QR7, paley=True)
        elif w == "R7_2":
            results[w] = analyse("R7_k2", {1, 2, 3}, 7, 2, {1})
        elif w == "P7_3":
            results[w] = analyse("P7_k3", QR7, 7, 3, QR7, paley=True)
        elif w == "NQR7_3":
            results[w] = analyse("NQR7_k3", NQR7, 7, 3, QR7, paley=True)
        elif w == "P7_4":
            results[w] = analyse("P7_k4", QR7, 7, 4, QR7, paley=True)
        elif w == "P11_2":
            results[w] = analyse("P11_k2", QR11, 11, 2, QR11, paley=True)
        elif w == "NQR11_2":
            results[w] = analyse("NQR11_k2", NQR11, 11, 2, QR11, paley=True)
        elif w == "R11_2":
            results[w] = analyse("R11_k2", {1, 2, 3, 4, 5}, 11, 2, {1})
    # converse conjugation checks
    def conj(a, b):
        # eta-bar = -1 - eta
        return (a - b, -b)
    for x, y in (("P7_2", "NQR7_2"), ("P7_3", "NQR7_3"), ("P11_2", "NQR11_2")):
        if x in results and y in results:
            a, b = results[x]["a"], results[x]["b"]
            say(f"converse: beta({y}) = ({results[y]['a']}, {results[y]['b']}); conj(beta({x})) = {conj(a, b)}; "
                f"equal: {(results[y]['a'], results[y]['b']) == conj(a, b)}")
    # det M_k = det B_0 * N(det B_1), checked modulo random primes
    for key, r in results.items():
        ok = True
        for l in (2147483629, 2147483587, 2147483579):
            dm = det_mod_np(r["M"], l)
            rhs = r["detB0"] * r["N"] % l
            ok &= (dm == rhs)
        say(f"{key}: det M_k = det B_0 * N(det B_1) modulo three 31-bit primes: {ok}")
    with open(os.path.join(HERE, "a3_runobs_summary.json"), "w") as f:
        json.dump({k: {"detB0": str(v["detB0"]), "coeffs": [str(c) for c in v["coeffs"]], "N": str(v["N"]),
                       **({"a": str(v["a"]), "b": str(v["b"])} if "a" in v else {})} for k, v in results.items()}, f, indent=1)
    return results

if __name__ == "__main__":
    main()
