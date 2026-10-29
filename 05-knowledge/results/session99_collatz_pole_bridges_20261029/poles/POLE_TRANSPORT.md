# Differential dilation transport of normalized tournament poles

Status: PROVED by the algebra below; exact finite verification is separately scoped.
Inherited input: Session98 `limits/ALL_POLES_CORE_EXPANSION.md`, deriving from
THM-514 (`the-a000568-unit-tail-core-kernel-and-prime-tail-ladder`) and THM-516
(`burnside-core-kernel-phi-reframe`). No priority claim. No Collatz convergence claim.

## 1. Objects and exact single-core insertion

For a core μ with odd parts at least three, let m=|μ|, c=length(μ),
k=(m-c)/2, t_p its multiplicity of p, and h_p(μ)=Σ_{a in μ}gcd(p,a).
Use Session98's E, L=binomial(m,2)-E, C_μ=m!2^(-L)/z_μ, and

    R_μ(z)=C_μ z^m/(1-z/4^k)^(m+1),     D=z d/dz.

The empty core has R_empty=1/(1-z). For every odd p≥3,

    R_(μ∪{p})(z) = [2^(-m+h_p(μ)-(p-1)^2/2)/(p(t_p+1))]
                    z^p ∏_{j=1}^p(D+j) R_μ(z/2^(p-1)).       (1)

All operators act on the function immediately to their right; z^p is multiplied
AFTER applying them. This ordering matters.

Proof. Adding p changes (m,k) to (m+p,k+(p-1)/2), and changes L by
pm+(p-1)^2/2-h_p. The centralizer denominator multiplies by p(t_p+1).
For any q and m≥0, successive Euler operations give

    ∏_{j=1}^p(D+j)[z^m/(1-qz)^(m+1)]
      = ((m+p)!/m!) z^m/(1-qz)^(m+p+1).

Indeed the j-th operation acts on denominator power m+j and multiplies
by m+j while raising that power by one. Scaling z first contributes
2^(-(p-1)m). Combining these factors proves (1).

Thus an actual core-insertion edge transports a pole 4^k to
4^(k+(p-1)/2), and its rational summand's exact order m+1 to m+p+1.
The affine update is additive, not the Collatz update k→3k+1.

## 2. A closed differential-dilation equation for the extremal pole tower

Let C_k=(3k)!/(3^k k!)·2^(-3k²+k) and

    T(z,u)=Σ_{k≥0} u^k C_k z^(3k)/(1-z/4^k)^(3k+1).

Then

    ∂_u T(z,u) = (z³/12)(D+1)(D+2)(D+3) T(z/4,u),
    T(z,0)=1/(1-z).                                      (2)

Proof. For μ=3^k, m=h_3=3k and t_3=k. Formula (1) becomes
(k+1)R_(3^(k+1))=(z³/12)∏_{j=1}^3(D+j)R_(3^k)(z/4).
Multiply by u^k and sum. Conversely, coefficient extraction in u shows
(2) with its initial condition uniquely determines this tower as a formal
u-series of rational functions. It is therefore a generating recurrence,
not just a consequence about pole locations.

The series converges normally for z in any compact avoiding {4^k:k≥0}
and u in any compact subset of C: the Session98 Gaussian-in-m core bound
absorbs the extra exponential |u|^k. Local derivatives may therefore be
interchanged with summation. Equation (2) holds meromorphically in z and
entirely in u. It gives a genuine exact mechanism for pole location ×4 and
order +3; it does not say the full B(z) consists only of this tower.

The companion exponential core series

    F(t)=Σ_{k≥0} 2^(-3k²+k)t^(3k)/(3^k k!)

is entire and satisfies the pantograph equation

    F'(t)=(t²/4)F(t/4),   F(0)=1.                         (3)

This follows directly by coefficient ratios (and uniquely determines all
coefficients, including the vanishing residue classes 1 and 2 mod 3).
Equation (3) keeps core weights but omits their fixed-point tails; equation
(2) restores those tails and therefore carries the pole information.

## 2a. Every odd single-length tower has the same mechanism

Fix an odd integer ℓ≥3, and retain only cores ℓ^k. Put

    q=2^(ℓ-1),
    C_(ℓ,k)=(ℓk)!/(ℓ^k k!) · 2^[-(ℓ-1)(ℓk²-k)/2],
    T_ℓ(z,u)=Σ_{k≥0} u^k C_(ℓ,k) z^(ℓk)/(1-z/q^k)^(ℓk+1).

Then the exact generalization is

    ∂_u T_ℓ = [2^(-(ℓ-1)^2/2)/ℓ] z^ℓ
                 ∏_{j=1}^ℓ(D+j) T_ℓ(z/q,u),
    T_ℓ(z,0)=1/(1-z).                                    (2a)

For u=1 it is meromorphic on C, with poles EXACTLY q^k and exact order
ℓk+1. Its leading amplitude at q^k is

    lim_(z→q^k)(1-z/q^k)^(ℓk+1)T_ℓ(z,1)
      = (ℓk)!/(ℓ^k k!) · 2^[(ℓ-1)(ℓk²+k)/2] > 0.        (2b)

Proof. The loss formula follows by subtracting
E=k(ℓ-1)/2+ℓ binomial(k,2) from binomial(ℓk,2).
In (1), h_ℓ=m=ℓk and t_ℓ=k, yielding (2a). Normal convergence
follows either from the inherited core bound or directly from the negative
quadratic exponent in C_(ℓ,k), uniformly for bounded u. At a given q^k,
exactly one summand has a pole; all other summands form a holomorphic
function locally by normal convergence. This proves exact order and (2b).
Consequently every such restricted T_ℓ(z,1) is non-D-finite, independently
of the arithmetic behavior of the map x↦oddpart(ℓx+1).

This is a hostile control against attributing arithmetic convergence to
meromorphic continuation, geometric pole spacing, or linear pole-order
growth alone. The ℓ=5 case has these same properties with order 5k+1;
its associated odd-part arithmetic need not behave like the 3x+1 case.
The root amplitude/readout note supplies the explicit arithmetic cycle.

## 3. Full core-bank differential hierarchy

Introduce independent cycle markers x_j for odd j≥3 and write

    M(z;x)=Σ_μ R_μ(z)∏_j x_j^(t_j(μ)).

This is coefficientwise well-defined in z: only cores of mass at most n
contribute to z^n. Then for every odd p≥3,

    ∂_(x_p) M(z;x)
      = [2^(-(p-1)^2/2)/p] z^p ∏_{r=1}^p(D+r)
          M(z/2^(p-1); (2^(gcd(p,j)-j)x_j)_j).            (4)

Proof. Remove one p-cycle from each differentiated monomial. Its factor
t_p+1 cancels the same denominator in (1). The remaining response weight
2^(-m+h_p)=∏_j(2^(gcd(p,j)-j))^(t_j) is exactly the marker substitution.
The initial condition M(z;0)=1/(1-z), with all these equations, determines
coefficients by total cycle count. Setting every x_j=1 recovers the full
normalized B(z), coefficientwise and via Session98 normal convergence.
Equation (4) is an infinite-state hierarchy, NOT a scalar finite-order ODE.

The substitution is a contraction on each nonnegative marker: gcd(p,j)≤j.
It fixes markers with j dividing p and contracts the others. For p=3 the
only retained core marker fixed is x_3, explaining the scalar closure (2).
For a general fixed p the singleton p-cycle tower also closes, with z
scaled by 2^(p-1), pole order increasing by p, and prefactor
2^(-(p-1)^2/2)/p. The top tower is distinguished because p=3 gives maximum
order at fixed defect, not because 3 is a mysterious dynamical multiplier.

## 4. Loss ledger and hostile test

Projection μ↦(k,m) preserves pole location and summand order, but loses h_p
and multiplicities. The cores (3,7) and (5,5) both have (k,m)=(4,10), yet
h_3 is 4 versus 2. The marker-substitution factor in (4) is 2^-6 versus
2^-8, so no insertion rule depending ONLY on (k,m) and applying the same
per-core multiplier can be exact. This does not prohibit more sophisticated
aggregate identities: it refutes the proposed insufficient local state.
The needed sidecar is precisely the gcd response, as in inherited THM-516
and Session98 `recurrence/PROJECTIVE_RESPONSE.md`.

The connection to Collatz is structural but sharply bounded: a dilation
and differential operator transports singularity data, while the root's
amplitude analysis extracts a literal 3k+1 factor from consecutive leading
amplitudes. Neither identifies this insertion graph with the Collatz graph.
Here k grows monotonically under insertion and no parity/valuation branch
selects a next core; Collatz has such branches and may decrease its state.
The full hierarchy moreover retains gcd data absent from pole labels.

## Verification

Run `python3 verify_pole_transport.py` in this directory. It uses exact
Fractions and integer coefficient formulas, with no symbolic-analysis
assumption or floating-point tolerance. It checks (1) on every core of mass
at most 18 for insertions p=3,5,7,9 through Taylor degree 35; the tower (2)
through k=10; the general ℓ=3,5,7,9 towers through k=5 and degree 55; (3) through degree 60; and the hostile same-(k,m) example.
These are finite controls, not substitutes for the all-degree proofs above.
