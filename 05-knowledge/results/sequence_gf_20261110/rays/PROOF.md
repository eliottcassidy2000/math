# A positive Mahler compiler for certified path-count rays

Session249, 10 November 2026. **PROVED in this note**, using the explicitly
inherited ray-membership theorem below. No new all-line coverage theorem
or external novelty claim. Chan–Pak remains OPEN.

## 1. Inheritance and the choice of index

[THM-5126](../../../../01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md)
proves that all values1+c 2^s, c odd, 1<=c<=129, s>=1, are Hamiltonian-path
counts except7 and21. It does so with actual rooted ray seeds and finite
patches below their onsets. We use this as an input, not repeat its census.
Its rooted-tail and doubling-vertex constructions increase vertex order
by one for each doubling of H-1 once a seed's construction starts.

This note's sequence is indexed by the numerical path value N, not the
number of tournament vertices n. Those two indices lead to radically
different recurrence classes for the same family.

## 2. Ray source theorem

Let C be any finite set of positive odd integers, and choose thresholds
s_c>=1. Put

    R(z) = sum_(c in C) z^(c 2^s_c),
    B(z) = sum_(c in C) sum_(s>=s_c) z^(c 2^s).

**Proposition1.** In the ring of formal power series,

    B(z) = R(z) + B(z^2),
    B(z) = sum_(k>=0) R(z^(2^k)).                         (1)

The coefficients are0 or1. The solution with constant term0 is unique.
With P(z)=z B(z), the positive exponents mark exactly the ray values
N=1+c2^s. If the input rays have tournament certificates, every positive
coefficient of P carries such a certificate.

Proof. Iterating z->z^2 gives (1), locally finite because R(0)=0. Every
positive integer has a unique odd part, and hence distinct c or distinct
s produce distinct exponents. Comparing the coefficient at m gives
b_m=r_m+1_(2|m)b_(m/2), determining b_m from smaller indices. This also
proves uniqueness. The shift by one changes m=H-1 into H.

The same proof applies to an infinite set C if R is a well-defined formal
series: only finitely many positive seed exponents may occur below each
bound. For counting distinct values, choose at most one onset per odd part.
Multiple construction certificates on the same line require taking the
minimum onset before using the0/1 assertion.

**Weighted construction version.** Suppose the given seed for c starts at
value1+c2^s_c on n_c vertices, and each subsequent step adds one vertex.
Then

    B(z,u)=sum_c sum_(k>=0) u^(n_c+k) z^(c2^(s_c+k))
          = R(z,u) + u B(z^2,u),
    R(z,u)=sum_c u^n_c z^(c2^s_c).                       (2)

This counts one selected witness per ray value, not all isomorphism
classes. At z=1 the formal u-series is sum_c u^n_c/(1-u) for finite C.
Thus witness enumeration by graph order is rational while membership by
path value below is not even P-recursive. The source must retain actual
n_c and construction onset, not the earliest independently known value.

## 3. The exceptional point11 and an exact65-line equation

Let C={1,3,5,...,129}. Define the polynomial with65 positive coefficients

    R_65(z) = sum_(c odd, 1<=c<=129) z^(2c)
              - z^6 - z^10 + z^12 + z^40.              (3)

The subtraction just replaces two seeds: on c=3 use exponent12 instead
of6; on c=5 use40 instead of10. Then the indicator GF of the inherited
certified set (including1) is

    P_65(z) = z + z^11 + z B_65(z),
    B_65(z) = R_65(z) + B_65(z^2).                     (4)

Indeed c=3 starts at s=2, omitting only7. c=5 starts at s=3, omitting11
and21, and11 is restored separately. All other lines start at s=1.
This is an exact generating function for a proved infinite subset of the
spectrum, not an equation for A000568 or for the full spectrum.

The early starts in (3) refer to set membership. They do NOT assert that
THM-5126 has one graph operation beginning at every such start: its actual
root seeds can start later, and its finite patches supply the earlier
points. For construction-size weights, use (2) with the actual onsets and
add those patches with their own witness sizes.

The point11 is structurally necessary here. The full conjectural target
S_*={positive odd integers} minus {7,21} fails closure under doubling H-1:
11 ->21. If P_*(z)=z/(1-z^2)-z^7-z^21 and B_*=(P_*-z)/z, then

    B_*(z)-B_*(z^2)
      = z^2/(1-z^4)-z^6+z^12-z^20+z^40.               (5)

The negative coefficient at20 records the edge11 ->21. Once the
exceptional point11 is removed, the source becomes nonnegative:

    R_*(z)=z^2/(1-z^4)-z^6-z^10+z^12+z^40,
    P_*(z)=z+z^11+z sum_(k>=0) R_*(z^(2^k)).           (6)

Thus the desired support has a very simple positive Mahler presentation.
The missing theorem is a tournament-realization lift of its unbounded
source, not an algebraic simplification of its target GF. Formula(6) alone
proves no new membership and does not assert that all its lowest points
are unit-offset graph seeds. An all-line seed theorem still needs its
finite below-ray points addressed on every line (THM-5126's warning).

## 4. Automaticity, exact counts, and non-P-recursiveness

**Proposition2.** For finite nonempty C, the indicator of the union of
rays is2-automatic and not P-recursive over characteristic zero. Finite
patches or deletions that leave infinitely many points do not change
these conclusions.

Proof of automaticity. In the standard binary expansion of N=1+c2^s,
s>=1, the word is bin(c) followed by s-1 zeros and then1. For each c its
language is the fixed prefix bin(c), at least s_c-1 zeros, and final1.
The finite union is regular. Finite changes remain regular. A complete
finite automaton recognizing the binary expansions is precisely a
2-automatic indicator (choose the convention for leading zeros).

Proof of non-P-recursiveness. At most |C|(1+log_2 X) terms below X are
nonzero, yet there are infinitely many nonzero terms. Therefore there
are arbitrarily long zero runs at arbitrarily large indices. If
sum_(j=0)^r p_j(N)a_(N+j)=0 eventually, with p_r nonzero, then beyond
its finitely many integer roots any r consecutive zeros force the next
zero and inductively the entire tail zero. Contradiction. The r=0 case
also immediately forces eventual zero. No analytic theorem is needed.
Consequently the OGF is neither rational nor D-finite in its value index.

The number of ray values <=X is exactly

    sum_c max(0, floor(log_2((X-1)/c))-s_c+1).          (7)

For C=1,3,...,129, add1 for N=1 and1 for N=11 when below X, and use
s_3=2, s_5=3, all other s_c=1. In particular the count is
65 log_2 X+O(1), so this proved infinite subset has density zero.
The compiler computes the floors using integer bit lengths, not floats.
The quasiperiodic bounded fluctuation in log_2 X is retained by (7).

**Corollary: an index-sensitive recurrence obstruction.** A finite bank
of productive affine ray constructions cannot have a P-recursive
value-indicator sequence, despite its rational graph-order witness series
and its finite automaton. This supplies a concrete stopping rule for
recurrence guessing, not a failure of exact generating functions.

## 5. Reproduction and connection contract

Run `python3 rays/compile_rays.py` from the packet directory. It prints
source support, the first values, integer counts at2^k, and exact finite
controls comparing the Mahler recursion, direct odd-part membership, and
the full-target formal identity. These are checks of the new algebra,
not a repeated tournament census or an independent audit of THM-5126.

Source: certified affine graph families. Target: positive Mahler series
and regular binary language. Map: forget the graph but retain H-1, its odd
part, seed onset and optionally witness order. Preserved: membership and
one selected witness cost. Lost: all class multiplicities, automorphisms,
boundary states and ability to apply a different graph operation. Restored
by retaining seed certificates and rooted transition data. Cheap hostile:
the missing point21 obstructs the naive propagation of11.

This connects the repo's dyadic construction work with automatic-sequence
methods through an explicit arithmetic map, rather than a digit analogy.
The general source compiler, positive residual equation, weighted version,
and index-sensitive obstruction are the session contributions. The65-line
coverage is wholly inherited. General positive-source realizability and
Chan–Pak completeness remain OPEN.
