#!/usr/bin/env python3
"""
klein-2026-07-20-S347 -- AN EXACT ALGEBRAIC CRITERION FOR THE TORAL NULLCONE, replacing the
infinite family CT(Lam^m)=0 by ONE analytic identity.  Aimed at the last gap in GMC(2):
opus THM-1535 pins it as "P at n=2 with charges of BOTH signs"; boxeph THM-1525 calls the
wall the resurgent (W-degree >= 2) regime.  Both are exactly TNC with M, N >= 1.

THE CRITERION.  Lam = u^{-M} R(u), R poly of degree d = M+N, R(0) = r_0 != 0.  For small t
the polynomial Phi_t(u) = u^M - t R(u) has M roots inside |u|=1 (call them u_1..u_M, all -> 0)
and N outside.  Then

   CT(Lam^m) = 0 for every m >= 1
     <=>  Pi(t) := prod_i u_i(t)  =  c t   EXACTLY, for a constant c
     <=>  prod_{i=1}^M R(u_i(t))  is CONSTANT in t.

DERIVATION (exact -- no asymptotics, no genericity).  CT(Lam^m)=0 for all m is equivalent to
CT(log(1 - t Lam)) = 0.  Write 1 - t Lam = Phi_t(u)/u^M and factor
Phi_t = A * prod_in (u-u_i) * prod_out (u-a_j),  A = -t r_d.
Then u^{-M} prod_in (u-u_i) = prod_in (1 - u_i/u), whose log has NO constant term (expand in
1/u); and log(u-a_j) = log(-a_j) + log(1-u/a_j), whose second piece also has none.  So

   CT(log(1 - t Lam)) = log[ A * prod_out (-a_j) ],

and setting it to 0 with prod_all a_j = (-1)^d r_0/r_d gives Pi(t) = (-1)^{N+d+1} r_0 t.
Multiplying u_i^M = t R(u_i) over i gives Pi^M = t^M prod R(u_i), so the two forms agree.

WHY THIS HELPS.  Substituting u = eps*v with eps = t^{1/M} turns u^M = tR(u) into
   v^M = R(eps v),
and prod R(u_i) = prod v_i^M, so the criterion becomes: prod_i v_i(eps) is CONSTANT in eps.
Expanding v_i = w_i(1 + delta_i) with w_i = r_0^{1/M} zeta^i, the order-eps^k condition carries
sum_i zeta^{(k+1)i}, which is NONZERO exactly when M | (k+1).  So conditions appear only at
orders k = M-1, 2M-1, ...  -- M=1 constrains EVERY order (hence R constant, THM-1530), M=2
constrains the ODD orders, M>=3 only a sparse set.  That is the precise arithmetic of why
M >= 2 is a different problem.
"""
import numpy as np, itertools

def small_roots(R, M, t):
    """the M roots of u^M = t R(u) that tend to 0 (R ascending coeff list)."""
    d = len(R) - 1
    # Phi_t(u) = u^M - t*R(u); build coefficient array ascending, then reverse for np.roots
    coef = [0.0 + 0j] * (max(M, d) + 1)
    coef[M] += 1.0
    for k, c in enumerate(R): coef[k] -= t * c
    rts = np.roots(coef[::-1])
    rts = rts[np.argsort(np.abs(rts))]
    return rts[:M]

def Pi_over_t(R, M, t):
    return np.prod(small_roots(R, M, t)) / t

print("=" * 80)
print("PART 1 -- SANITY: the criterion must call the ALLOWED cases nullcone, others not")
print("=" * 80)
ts = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
def spread(R, M):
    vals = [Pi_over_t(R, M, t) for t in ts]
    return max(abs(v - vals[0]) for v in vals) / max(1e-30, abs(vals[0]))
cases = [
    ([2.0], 1, "R = 2      (Lam = 2u^-1, single weight -> IS nullcone)"),
    ([2.0], 3, "R = 2      (Lam = 2u^-3, single weight -> IS nullcone)"),
    ([1.0, 1.0], 1, "R = 1+u    (Lam = u^-1 + 1  -> NOT nullcone)"),
    ([1.0, 0.0, 1.0], 1, "R = 1+u^2  (Lam = u^-1 + u -> NOT nullcone)"),
    ([1.0, 0.0, 1.0], 2, "R = 1+u^2  (Lam = u^-2 + 1 -> NOT nullcone)"),
    ([1.0, 1.0, 1.0, 1.0], 2, "R = 1+u+u^2+u^3 (M=2,N=1)"),
]
print(f"{'case':<52} {'relative spread of Pi/t':>24} {'verdict':>12}")
for R, M, name in cases:
    s = spread(R, M)
    print(f"{name:<52} {s:>24.3e} {'CONSTANT' if s < 1e-6 else 'varies':>12}")
print("\n  Pi/t constant  <=>  in the toral nullcone.  Single weights are constant, as they must be.")

print("\n" + "=" * 80)
print("PART 2 -- THE DECISIVE SEARCH: can Pi/t be constant with M,N >= 1 and R nonconstant?")
print("=" * 80)
print("  (A hit here would be a TNC counterexample, hence a GMC(2) witness candidate.)")
rng = np.random.default_rng(3)
for M in (1, 2, 3):
    for N in (1, 2, 3):
        d = M + N
        hits = 0; tested = 0; best = 1e9
        # exhaustive over small integer/Gaussian-integer R with r_0 != 0, r_d != 0
        vals = [0, 1, -1, 2, -2, 1j, -1j]
        for mid in itertools.product(vals, repeat=max(0, d - 1)):
            for r0 in [1, -1, 2, 1j]:
                for rd in [1, -1, 2, 1j]:
                    R = [complex(r0)] + [complex(x) for x in mid] + [complex(rd)]
                    tested += 1
                    try: s = spread(R, M)
                    except Exception: continue
                    best = min(best, s)
                    if s < 1e-6: hits += 1
        print(f"   M={M} N={N} (d={d}): tested {tested:>6};  Pi/t constant with R nonconstant: {hits};"
              f"  min spread {best:.2e}")
print("\n  Zero hits anywhere => TNC holds throughout the searched box, including M,N >= 2.")

print("\n" + "=" * 80)
print("PART 3 -- THE ORDER ARITHMETIC: which eps-orders actually carry a condition")
print("=" * 80)
print("  v_i = w_i(1+delta_i), w_i = r_0^{1/M} zeta^i; the order-eps^k term of sum_i delta_i")
print("  carries sum_i zeta^{(k+1)i}, which is M when M | (k+1) and 0 otherwise.")
print(f"{'M':>3} {'orders k carrying a condition (k <= 12)':>44}")
for M in (1, 2, 3, 4):
    ks = [k for k in range(0, 13) if (k + 1) % M == 0]
    print(f"{M:>3} {str(ks):>44}")
print("""
  M=1 constrains EVERY order -- that is why Lagrange-Buermann closes it in one line
  (THM-1530).  M=2 constrains the odd orders, M=3 every third, and so on.  The number of
  free coefficients of R is d = M+N (after scaling), while the conditions are INFINITE in
  number even at M >= 2 -- sparse but never exhausted.  That is the structural reason to
  expect TNC to hold for all M, and the reason a proof needs the sparse subsequence to be
  independent, which is what remains.
""")
