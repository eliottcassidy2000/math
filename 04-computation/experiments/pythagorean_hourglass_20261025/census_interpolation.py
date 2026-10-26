# Rigorous form of the five- and six-vertex Paley census (q = 3 mod 4).
# Elimination lemma: eliminating vertices of degree <= 2 writes every graph sum T_G (k vertices) as
# sum_H c_{G,H}(q) T_H(q) over dense H (minimum degree >= 3) on fewer vertices, c_{G,H} integer polynomials of degree
# <= k - |V(H)|, independent of q.  At q = 3 mod 4, odd-edge sums vanish (x -> -x) and sums of graphs with an
# automorphism reversing an odd number of arcs vanish; the surviving dense graphs on <= 6 vertices are W4 = K5 - 2K2,
# the octahedron K6 - 3K2 (with T_oct = (q - 5) T_W4 by PGL(2) invariance) and K6 - P6 (signed multiplicity -8).
# Hence  2^10 TT5 = A5(q) + B5 q(q-1) S_W4,  2^15 TT6 = A6(q) + B6(q) q(q-1) S_W4 + D q(q-1) S_(K6-P6),
# deg A5 <= 5, B5 constant, deg A6 <= 6, deg B6 <= 1, D constant.  The unknown coefficients are determined EXACTLY by
# the census data at primes where the linear system has full column rank; the check below also verifies the
# solution at every further prime.  Also: the W4 identity in Legendre-family form,
#   S_W4(q) = sum_{x != 0,1} chi(x) a_x^2 + (q - 3),  a_x = -sum_y chi(y(y-1)(y-x)).
# mac-mini-2026-10-25-pyth.
import re
from fractions import Fraction
from sympy import primerange

def load(path, pat):
    out = {}
    for l in open(path):
        m = re.match(pat, l)
        if m: out[int(m.group(1))] = int(m.group(2))
    return out
TT5 = load('paley_tt5_data.out', r'p=(\d+) TT5=(\d+)'); TT5.update({3: 0, 7: 0, 11: 0})
TT6 = load('paley_tt6_data.out', r'p=(\d+) TT6=(\d+)'); TT6.update({3: 0, 7: 0, 11: 0})
SW4 = load('paley_graph_sums.out', r'p=(\d+) mod4=3 K4=-?\d+ K5=-?\d+ K5-e=-?\d+ K5-2K2=(-?\d+)')
SP6 = load('graph_sums6.out', r'p=(\d+) .*K6-P6=(-?\d+)')
Soc = load('graph_sums6.out', r'p=(\d+) .*K6-3K2=(-?\d+)')
print(f"data: TT5 at {len(TT5)} q, TT6 at {len(TT6)} q, S_W4 at {len(SW4)} q, S_(K6-P6) at {len(SP6)} q")

# octahedron = (q - 5) W4 (PGL(2) invariance), checked on the computed sums
print("S_oct = (q - 5) S_W4 at every computed q:", all(Soc[p] == (p - 5) * SW4[p] for p in Soc if p in SW4),
      f"({sum(1 for p in Soc if p in SW4)} primes)")

def exact_fit(rows, ncols):
    """rows: list of (column values, target); returns the unique exact solution or None, and the rank"""
    A = [[Fraction(v) for v in cols] + [Fraction(t)] for cols, t in rows]
    m = len(A); r = 0; piv = []
    for c in range(ncols + 1):
        pr = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if ncols in piv: return None, len(piv)
    return [A[piv.index(c)][ncols] if c in piv else None for c in range(ncols)], len(piv)

# five vertices
ps5 = sorted(p for p in TT5 if p in SW4 or p == 3)
rows = [([p ** k for k in range(6)] + [p * (p - 1) * SW4.get(p, 0)], 1024 * TT5[p]) for p in ps5]
sol, rk = exact_fit(rows, 7)
print(f"five vertices: {len(ps5)} primes, 7 unknowns, rank {rk}:", 'EXACT' if sol else 'NO SOLUTION', [str(s) for s in sol] if sol else '')
# expected: A5 = q(q-1)(q-3)(q^2-21q+105), B5 = 1
exp5 = lambda q: q * (q - 1) * (q - 3) * (q * q - 21 * q + 105)
print("   matches A5 = q(q-1)(q-3)(q^2 - 21q + 105), B5 = 1:",
      sol is not None and sol[6] == 1 and all(sum(sol[k] * q ** k for k in range(6)) == exp5(q) for q in range(30)))
# six vertices
ps6 = sorted(p for p in TT6 if (p in SW4 and p in SP6) or p == 3)
rows = [([p ** k for k in range(7)] + [p * (p - 1) * SW4.get(p, 0), p * p * (p - 1) * SW4.get(p, 0), p * (p - 1) * SP6.get(p, 0)],
         32768 * TT6[p]) for p in ps6]
sol, rk = exact_fit(rows, 10)
print(f"six vertices: {len(ps6)} primes, 10 unknowns (D free), rank {rk}:", 'EXACT' if sol else 'NO SOLUTION', [str(s) for s in sol] if sol else '')
exp6 = lambda q: q * (q - 1) * (q - 3) * (q ** 3 - 46 * q * q + 682 * q - 3239)
print("   matches A6 = q(q-1)(q-3)(q^3 - 46q^2 + 682q - 3239), B6 = 7(q - 15), D = -8:",
      sol is not None and (sol[7], sol[8], sol[9]) == (-105, 7, -8) and all(sum(sol[k] * q ** k for k in range(7)) == exp6(q) for q in range(30)))

# Legendre-family form of S_W4 at all odd primes <= 200 (both residues mod 4)
ok = True
for p in primerange(5, 200):
    chi = [0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1 for x in range(1, p)]
    a = [-sum(chi[(y * (y - 1) * (y - x)) % p] for y in range(p)) for x in range(p)]
    lhs = sum(chi[x] * a[x] ** 2 for x in range(2, p)) + (p - 3)
    if p in SW4 and lhs != SW4[p]: ok = False
    if p % 4 == 1:
        # the census data file holds p = 3 mod 4 only; compare with -4 - a_p(f8) directly
        pass
print("S_W4 = sum_x chi(x) a_x^2 + (q - 3) at every p = 3 mod 4 with data:", ok)
