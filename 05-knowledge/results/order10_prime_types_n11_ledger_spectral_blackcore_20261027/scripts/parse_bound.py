"""Explicit exception lower bounds for reparsing a doubly regular (Paley) tournament of order q with a bounded bank.
nu(T) = nuclear norm of the skew-adjacency matrix.  Perturbation: if T and T' differ on e pairs then
|nu(T) - nu(T')| <= ||S_T - S_T'||_* <= 4e (if 2e <= n) and <= sqrt(8 e n) (always).  Hence with Delta = nu(T) - nu(T'):
e >= Delta/4 when Delta <= 2n, and e >= Delta^2/(8n) when Delta > 2n.
Baselines: linear parse TT_k[B_1..B_k], |B_i| <= b:  nu(T') <= nu(TT_n) + n * kappa(b), kappa(b) = max_{m<=b} (nu_max(m) + nu(TT_m))/m,
nu_max(m) = (m-1)sqrt(m) (m odd) or m sqrt(m-1) (m even) [the Cauchy-Schwarz bound].
Bounded-arity substitution trees (any shape): nu(T') <= (2b-3) n ceil(log2 n) [concurrent packet dyadic_module_spectrum T1].
Score bound for linear parses: e >= L(T) - L(T') >= (q^2-1)/8 - q(b-1)/2."""
import math
def nu_tt(n): return sum(abs(1 / math.tan((2 * k - 1) * math.pi / (2 * n))) for k in range(1, n + 1) if (2 * k - 1) != n)
def nu_max(m): return (m - 1) * math.sqrt(m) if m % 2 else m * math.sqrt(m - 1)
def kappa(b): return max((nu_max(m) + nu_tt(m)) / m for m in range(1, b + 1))
def ebound(delta, n):
    if delta <= 0: return 0
    return math.ceil(delta / 4 - 1e-9) if delta <= 2 * n else math.ceil(delta * delta / (8 * n) - 1e-9)
print("q | nu(P_q)=(q-1)sqrt q | nu(TT_q) | linear b=1: spectral e>= | score (q^2-1)/8 | linear b=3: spectral | score | b=3 tree (T1 const): spectral | b=7 tree: spectral")
for q in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 127, 251, 499, 1019, 2039, 4079, 8191]:
    nuP = (q - 1) * math.sqrt(q); nT = nu_tt(q); L = (q * q - 1) // 8
    e1 = ebound(nuP - nT, q)
    e3 = ebound(nuP - nT - q * kappa(3), q); s3 = max(0, L - q * (3 - 1) // 2)
    t3 = ebound(nuP - (2 * 3 - 3) * q * math.ceil(math.log2(q)), q)
    t7 = ebound(nuP - (2 * 7 - 3) * q * math.ceil(math.log2(q)), q)
    print(f"{q:5d} | {nuP:10.2f} | {nT:10.2f} | {e1:7d} | {L:8d} | {e3:7d} | {s3:7d} | {t3:7d} | {t7:7d}")
print("known feedback numbers f(P_q): 7:7 11:20 19:64 23:92 31:180 (THM-4810 / A003141)")
print("asymptotics: spectral e >= q^2/8 - (q^{3/2}/pi) ln q (1+o(1)) for every fixed bank; score e >= q^2/8 - O(bq) for linear parses only")
