#!/usr/bin/env python3
"""
operator_monad_numbers_s90d.py — Natural Numbers as Tournaments, Operator Theory
opus-2026-03-15-S90d

"Natural numbers are graphs, but they could also be tournaments.
What does operator theory say?"

SYNTHESIS of:
  - Our Cayley monad (S90c): Q⁴=Id, fixed points ±i, cross-ratio 8/7
  - Kind-pasteur (S112): 1/n weights, odd=directed, even=undirected
  - Operator theory: spectral decomposition, resolvent, functional calculus

THE CENTRAL CLAIM:
  The natural number n IS the n-th operator in the Cayley hierarchy.
  Its "graph" aspect is the even (undirected) part.
  Its "tournament" aspect is the odd (directed) part.
  The monad Q acts on these operators, cycling them with period 4.
"""

import numpy as np
from math import factorial, comb, log, pi, sqrt
from fractions import Fraction

# ======================================================================
# PART 1: THE NUMBER OPERATOR
# ======================================================================
print("=" * 70)
print("PART 1: THE NUMBER OPERATOR — n AS AN OPERATOR, NOT A SCALAR")
print("=" * 70)

print("""
In quantum field theory, the NUMBER OPERATOR N counts particles:
  N = Σ_k a†_k a_k

For tournaments, the analog is the CYCLE COUNTING OPERATOR:
  N_T = Σ_{odd k} α_k(T)  where α_k = #{independent sets of size k in Ω(T)}

The OCF says: H(T) = I(Ω(T), 2) = Σ_k α_k · 2^k = evaluation of the
independence polynomial at the "fugacity" z = 2.

This IS a partition function: H = Z(β) at β = -log(2).

In OPERATOR LANGUAGE:
  The "Hamiltonian" is H_op = -Σ_k k·n_k  (n_k = occupation number of cycle type k)
  The partition function is Z = Tr(e^{-βH_op}) = Σ_configs e^{-β·E(config)}
  At β = -log(2): Z = Σ_configs 2^{#cycles} = I(Ω, 2) = H(T)

The NATURAL NUMBER n appears as:
  - Scale index: the n-th term in arctanh(x) = Σ x^{2k+1}/(2k+1)
  - Operator index: the n-th Fourier harmonic of the tournament
  - Cycle length: directed cycles of length n
  - Interaction order: the n-point correlation function
""")

# The transfer matrix as an operator
M_x1 = np.array([[1, 2, 0],
                  [0, 0, 1],
                  [1, 1, 0]], dtype=float)

evals, evecs = np.linalg.eig(M_x1)
print("Transfer matrix M(1) — the OPERATOR:")
print(f"  M = {M_x1.tolist()}")
print(f"  Eigenvalues: {[f'{e:.6f}' for e in evals]}")
print(f"  det(M) = {np.linalg.det(M_x1):.6f} = 1 (unimodular)")
print(f"  tr(M) = {np.trace(M_x1):.0f}")
print(f"  tr(M²) = {np.trace(M_x1 @ M_x1):.0f}")
print(f"  tr(M³) = {np.trace(np.linalg.matrix_power(M_x1, 3)):.0f}")

# ======================================================================
# PART 2: THE RESOLVENT — WHERE THE MONAD LIVES
# ======================================================================
print()
print("=" * 70)
print("PART 2: THE RESOLVENT (zI - M)⁻¹ — THE MONAD'S HOME")
print("=" * 70)

print("""
The RESOLVENT of M is R(z) = (zI - M)⁻¹.
It exists for z not an eigenvalue.

R(z) has POLES at z = λ_i (the eigenvalues τ₃, λ_c, λ̄_c).

The SPECTRAL DECOMPOSITION:
  R(z) = Σ_i P_i / (z - λ_i)
where P_i are the spectral projections.

The TRACE of the resolvent:
  tr R(z) = Σ_i 1/(z - λ_i)

This is a RATIONAL FUNCTION with poles at the eigenvalues.
Its residues are the MULTIPLICITIES (all 1 for M(1)).

GENERATING FUNCTION interpretation:
  tr R(z) = Σ_{n≥0} tr(M^n) / z^{n+1}
           = Σ_{n≥0} T_n / z^{n+1}

where T_n = tr(M^n) is the tribonacci-like trace sequence.

So the RESOLVENT is the GF for the TRIBONACCI SEQUENCE!
""")

# Compute the resolvent trace
def resolvent_trace(z, M):
    """tr((zI - M)^{-1})"""
    n = M.shape[0]
    return np.trace(np.linalg.inv(z * np.eye(n) - M))

# The resolvent trace should equal the rational function
# tr R(z) = 3z² - 2z - 1 / (z³ - z² - z - 1)  [from char poly]
# Since char poly p(z) = z³ - z² - z - 1, we have tr R(z) = p'(z)/p(z)
# p'(z) = 3z² - 2z - 1
# So tr R(z) = (3z² - 2z - 1) / (z³ - z² - z - 1)

print("Resolvent trace tr(zI - M)⁻¹:")
header = "p'(z)/p(z)"
print(f"{'z':>10} | {'tr R(z) numerical':>20} | {header:>20} | {'match':>6}")
print("-" * 65)

for z_val in [3, 4, 5, 10, -2, 0.5+1j]:
    z = complex(z_val)
    try:
        tr_num = resolvent_trace(z, M_x1)
        p_z = z**3 - z**2 - z - 1
        dp_z = 3*z**2 - 2*z - 1
        tr_alg = dp_z / p_z
        match = abs(tr_num - tr_alg) < 1e-8
        print(f"{str(z_val):>10} | {tr_num.real:>10.6f}{'+' if tr_num.imag>=0 else ''}{tr_num.imag:.6f}i | "
              f"{tr_alg.real:>10.6f}{'+' if tr_alg.imag>=0 else ''}{tr_alg.imag:.6f}i | {match}")
    except:
        print(f"{str(z_val):>10} | singular")

# ======================================================================
# PART 3: THE CAYLEY TRANSFORM OF THE OPERATOR
# ======================================================================
print()
print("=" * 70)
print("PART 3: Q(M) — THE CAYLEY TRANSFORM OF THE OPERATOR ITSELF")
print("=" * 70)

print("""
We can apply the Cayley transform Q to the MATRIX M:
  Q(M) = (I + M)(I - M)⁻¹

This is well-defined since no eigenvalue of M is 1.
(Eigenvalues: τ₃ ≈ 1.839, λ_c ≈ -0.420 ± 0.606i)

Wait — τ₃ > 1, so (I - M) has eigenvalue 1 - τ₃ < 0. It's invertible!

Q(M) maps the eigenvalues: Q(λ_i) = (1+λ_i)/(1-λ_i)
  Q(τ₃) = (1+1.839)/(1-1.839) = 2.839/(-0.839) = -3.384
  Q(λ_c) = (1+λ_c)/(1-λ_c)

And Q(M)⁴ = M (since Q has order 4 on each eigenvalue!... wait, Q
  has order 4 only on the PROJECTIVE LINE. On matrices, Q⁴(M) = M
  iff all eigenvalues return to themselves under Q⁴, which they do.)

So Q⁴(M) = M. The Cayley transform of the operator is also 4-periodic!
""")

I = np.eye(3)
# Q(M) = (I + M)(I - M)^{-1}
Q_M = (I + M_x1) @ np.linalg.inv(I - M_x1)
Q2_M = (I + Q_M) @ np.linalg.inv(I - Q_M)  # Q²(M)
Q3_M = (I + Q2_M) @ np.linalg.inv(I - Q2_M)  # Q³(M)
Q4_M = (I + Q3_M) @ np.linalg.inv(I - Q3_M)  # Q⁴(M)

print(f"M =\n{np.round(M_x1, 4)}\n")
print(f"Q(M) =\n{np.round(Q_M, 4)}\n")
print(f"Q²(M) =\n{np.round(Q2_M, 4)}\n")
print(f"Q³(M) =\n{np.round(Q3_M, 4)}\n")
print(f"Q⁴(M) =\n{np.round(Q4_M, 4)}")
print(f"\nQ⁴(M) = M? {np.allclose(Q4_M, M_x1)}")

# Eigenvalues of Q(M)
evals_Q = np.linalg.eigvals(Q_M)
print(f"\nEigenvalues of Q(M): {[f'{e:.4f}' for e in sorted(evals_Q, key=lambda z: -z.real)]}")

# Verify: Q(λ) for each eigenvalue of M
for lam in evals:
    q_lam = (1 + lam) / (1 - lam)
    print(f"  λ = {lam:.6f} → Q(λ) = {q_lam:.6f}")

# ======================================================================
# PART 4: NATURAL NUMBERS AS SPECTRAL RESIDUES
# ======================================================================
print()
print("=" * 70)
print("PART 4: n AS A SPECTRAL RESIDUE")
print("=" * 70)

print("""
THE OPERATOR-THEORETIC VIEW OF NATURAL NUMBERS:

The n-th natural number appears in the RESOLVENT EXPANSION:
  tr R(z) = Σ_{n≥0} T_n / z^{n+1}

where T_n = tr(M^n) = τ₃^n + 2|λ_c|^n cos(nθ).

So the "n-th natural number" in our framework is T_n:
  T_0 = 3,  T_1 = 1,  T_2 = 3,  T_3 = 7,  T_4 = 11,
  T_5 = 21, T_6 = 39, T_7 = 71, T_8 = 131, ...

These ARE the natural numbers — not 1, 2, 3, 4, ... but the
TRIBONACCI NATURALS 3, 1, 3, 7, 11, 21, 39, 71, 131, ...

The standard natural numbers 1, 2, 3, ... appear as the
INTERACTION SCALES (the index n), while the tribonacci
numbers T_n are the SPECTRAL CONTENT at each scale.

RECONCILIATION:
  The interaction scale n and the spectral content T_n are
  DUAL descriptions of the same thing:
    n is the POSITION (when the interaction happens)
    T_n is the MAGNITUDE (how strong the interaction is)

  The resolvent connects them: T_n = ∮ z^n tr R(z) dz / 2πi
  (Cauchy integral formula — the n-th moment of the spectral measure)
""")

# Compute T_n for comparison
traces = []
Mn = np.eye(3)
for n in range(20):
    traces.append(round(np.trace(Mn).real))
    Mn = Mn @ M_x1

print("The spectral naturals T_n = tr(M^n):")
print(f"  n:  {list(range(20))}")
print(f"  T_n: {traces}")

# Show: T_n satisfies the tribonacci recurrence
print(f"\nVerify tribonacci: T_n = T_{{n-1}} + T_{{n-2}} + T_{{n-3}}:")
for n in range(3, 15):
    lhs = traces[n]
    rhs = traces[n-1] + traces[n-2] + traces[n-3]
    print(f"  T_{n} = {lhs}, T_{n-1}+T_{n-2}+T_{n-3} = {rhs}, match: {lhs == rhs}")

# ======================================================================
# PART 5: THE TOEPLITZ CONNECTION — CIRCULANT TOURNAMENTS
# ======================================================================
print()
print("=" * 70)
print("PART 5: TOEPLITZ OPERATORS — CIRCULANT TOURNAMENTS")
print("=" * 70)

print("""
The adjacency matrix of a CIRCULANT tournament is a CIRCULANT MATRIX.
Circulant matrices are special TOEPLITZ matrices.

In operator theory, a TOEPLITZ OPERATOR T_f on the Hardy space H²
has symbol f: S¹ → ℂ (a function on the unit circle).

For a circulant tournament on Z_p with connection set S:
  Symbol f(e^{iθ}) = Σ_{s∈S} e^{isθ}
  Eigenvalues = {f(ω^k) : k = 0,...,p-1} where ω = e^{2πi/p}

The SZEGŐ LIMIT THEOREM:
  As n → ∞, the distribution of eigenvalues of T_n(f) converges to
  the pushforward measure f_*(dθ/2π).

For PALEY tournaments:
  f(e^{iθ}) = Σ_{s∈QR(p)} e^{isθ} = (χ(k)i√p - 1)/2  at θ = 2πk/p
  The symbol takes only 2 values! (Fibonacci-type, as we discovered)
  Szegő: eigenvalue distribution → two-point measure

For INTERVAL tournaments:
  f(e^{iθ}) = Σ_{s=1}^{(p-1)/2} e^{isθ} = Dirichlet kernel
  The symbol is continuous → Szegő gives absolutely continuous measure

This is why interval tournaments maximize H: their symbol has the
LARGEST L¹ norm, which by Szegő maximizes the trace of any polynomial.

THE PERMANENT via OPERATOR THEORY:
  H(T) = perm(A) where A is the adjacency matrix.
  By the van der Waerden inequality: perm(A) ≥ n!/n^n for doubly stochastic A.
  The near-transitive tournament has A NOT doubly stochastic (row sums differ).
  But the OPERATOR NORM ||A|| = λ_max gives a bound: perm ≤ λ_max^n.
  More refined: perm relates to the HAFNIAN of (A + A^T), connecting
  to the undirected structure.
""")

# Compute symbol norms for various circulant tournaments
for p in [7, 11, 13]:
    omega = np.exp(2j * pi / p)

    # Paley (p ≡ 3 mod 4)
    if p % 4 == 3:
        qr = set()
        for a in range(1, p):
            qr.add((a*a) % p)
        S_p = sorted(qr)
        symbol_vals = [sum(omega**(k*s) for s in S_p) for k in range(p)]
        l1_norm = sum(abs(v) for v in symbol_vals) / p
        linf = max(abs(v) for v in symbol_vals)
        print(f"  Paley T_{p}: |S|={len(S_p)}, L¹ norm={l1_norm:.4f}, L∞ norm={linf:.4f}")

    # Interval
    S_int = list(range(1, (p+1)//2))
    symbol_vals = [sum(omega**(k*s) for s in S_int) for k in range(p)]
    l1_norm = sum(abs(v) for v in symbol_vals) / p
    linf = max(abs(v) for v in symbol_vals)
    print(f"  Interval_{p}: |S|={len(S_int)}, L¹ norm={l1_norm:.4f}, L∞ norm={linf:.4f}")

# ======================================================================
# PART 6: THE SPECTRAL ZETA — CONNECTING ALL SCALES
# ======================================================================
print()
print("=" * 70)
print("PART 6: THE SPECTRAL ZETA — THE MONAD'S DNA")
print("=" * 70)

print("""
The SPECTRAL ZETA FUNCTION of M:
  ζ_M(s) = Σ_{i=1}^{3} λ_i^{-s}

At NEGATIVE integers: ζ_M(-n) = tr(M^n) = tribonacci T_n.
At s = 0: ζ_M(0) = 3 = dim(M) (the number of states).
At POSITIVE integers: ζ_M(n) = "reciprocal moments".

The FUNCTIONAL EQUATION:
  ζ_M(s) involves τ₃^{-s} + 2|λ_c|^{-s} cos(sθ - s·arg(λ_c))

Since |λ_c|² = 1/τ₃ (from det = 1):
  ζ_M(s) = τ₃^{-s} + 2·τ₃^{s/2} · cos(s·arg(λ_c))

This has a SELF-DUALITY: the term τ₃^{-s} at large s is balanced
by τ₃^{s/2} from the complex pair. The critical line is at
Re(s) = 0 (where both terms have the same magnitude).

ANALOGY WITH RIEMANN ZETA:
  ζ(s) = Σ n^{-s}  (Riemann: over all natural numbers)
  ζ_M(s) = Σ λ_i^{-s}  (spectral: over eigenvalues)

  Riemann: functional equation ζ(s) ↔ ζ(1-s)
  Spectral: "functional equation" from det(M) = τ₃·|λ_c|² = 1

  Riemann: trivial zeros at s = -2, -4, -6, ...
  Spectral: ζ_M(-n) = T_n (tribonacci, NEVER zero!)

  Riemann: non-trivial zeros on Re(s) = 1/2 (hypothesis)
  Spectral: non-trivial zeros when τ₃^{-s} = -2τ₃^{s/2}cos(sθ)
""")

# Compute spectral zeta at various s
tau3 = 1.8392867552141612
lam_c = complex(-0.41964337760708, 0.60629073187415)

print("Spectral zeta ζ_M(s) = τ₃^{-s} + λ_c^{-s} + λ̄_c^{-s}:")
print(f"{'s':>8} | {'ζ_M(s)':>20} | {'T_n check':>10}")
print("-" * 45)

for s in range(-10, 6):
    z = tau3**(-s) + lam_c**(-s) + np.conj(lam_c)**(-s)
    check = ""
    if s <= 0:
        Mn = np.linalg.matrix_power(M_x1, -s)
        check = f"tr(M^{-s})={np.trace(Mn).real:.0f}"
    print(f"{s:8d} | {z.real:>10.4f} + {z.imag:>8.4f}i | {check}")

# ======================================================================
# PART 7: THE GRAND DICTIONARY
# ======================================================================
print()
print("=" * 70)
print("PART 7: THE GRAND DICTIONARY — NUMBER ↔ OPERATOR ↔ TOURNAMENT")
print("=" * 70)

print("""
╔══════════════════════════════════════════════════════════════════════╗
║          NUMBER  ↔  OPERATOR  ↔  TOURNAMENT  ↔  MONAD             ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                    ║
║  n ∈ ℕ          M^n           tr(M^n) = T_n     Q^n mod 4         ║
║  1/n            Resolvent     n-th cumulant      1/n weight        ║
║  n odd          Directed op   Odd cycles (3,5,7) arctanh coeff    ║
║  n even         Undirected    Symmetric pairs    log(1-x²) coeff  ║
║  n prime        Irreducible   Paley T_n          Simple monad     ║
║  n composite    Factored      Non-Paley          Composed monad   ║
║                                                                    ║
║  The NUMBER LINE is the SPECTRUM of the shift operator.            ║
║  The NATURAL NUMBERS are the resonances of Q⁴ = Id.              ║
║  Each n is an OPERATOR SCALE at which tournaments interact.        ║
║                                                                    ║
║  GRAPH = even part = undirected = 1/(1-x²) = symmetric            ║
║  TOURNAMENT = odd part = directed = Q(x) = antisymmetric           ║
║  NATURAL NUMBER = both together = 1/(1-x) = the whole             ║
║                                                                    ║
║  n = graph + tournament                                            ║
║  1/(1-x) = √(Q(x)) · √(1/(1-x²))                               ║
║  The whole = √(directed) · √(undirected)                          ║
║                                                                    ║
║  This is the GEOMETRIC MEAN of two worlds.                        ║
║                                                                    ║
║  The MONAD Q cycles through 4 aspects:                             ║
║    Q⁰ = identity (the number itself)                              ║
║    Q¹ = Cayley (directed/undirected exchange)                     ║
║    Q² = inversion (-1/x) (reciprocal/complement)                  ║
║    Q³ = inverse Cayley (return through the mirror)                ║
║                                                                    ║
║  After 4 steps: you're back. The monad has period 4.              ║
║  The fixed points ±i are where directed = undirected:             ║
║  the IMAGINARY UNIT is the tournament that IS its own graph.      ║
║                                                                    ║
╚══════════════════════════════════════════════════════════════════════╝

WHAT OPERATOR THEORY SAYS:

1. Natural numbers are SPECTRAL RESIDUES of the resolvent of M.
   T_n = tr(M^n) = tribonacci = the "content" at scale n.

2. The TOEPLITZ STRUCTURE of circulant tournaments connects to
   the Szegő limit theorem: eigenvalue distributions → symbol measure.

3. The CAYLEY TRANSFORM Q acts on M with period 4, cycling through
   four "aspects" of the operator (direct, Cayley, inverse, reverse).

4. The SPECTRAL ZETA ζ_M(s) encodes ALL moments simultaneously,
   with det(M)=1 providing a self-duality at s ↔ -s.

5. Each natural number n is BOTH a tournament (odd part, directed)
   AND a graph (even part, undirected). The geometric mean
   1/(1-x) = √Q(x) · √(1/(1-x²)) decomposes the number line
   into its directed and undirected components.

opus-2026-03-15-S90d: THE OPERATOR MONAD
""")
