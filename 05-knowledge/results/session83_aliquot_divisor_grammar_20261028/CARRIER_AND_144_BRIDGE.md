# One guarded composition law, and the exact 2–3–11 / 2–5–7 bridge

Status: PROVED elementary composition identities; FINITE-EXACT complete sigma fiber. This is a synthesis of inherited mechanisms and the present classification, not a claim of a new universal theory or a tournament isomorphism.

## 1. The observer that retains what the next operation needs

For n>=1 retain source n, sum S=σ(n), divisor count d=τ(n), and prime support. For coprime a,b, the exact operation is

\[
(n,S,d)(m,R,e)=(nm,SR,de).
\]

The support guard is essential. The derived observations are aliquot sum s=S-n, proper-divisor term count T=d-1, excess Delta=S-2n, and whether S is square. Consequently

\[
s(ab)=a s(b)+b s(a)+s(a)s(b),\qquad
T(ab)=(T(a)+1)(T(b)+1)-1.
\]

This is an operation-preserving carrier for these stated observations. It is not a minimal encoding theorem: n itself determines its factorization and all the other coordinates, but retaining them makes composition explicit. Removing n and keeping only s fails at27 and35; removing the prime-support guard fails when factors share a prime; removing S's square-class information loses square edges.

The current [Möbius–Liouville carrier work](https://github.com/eliottcassidy2000/math/blob/130bfba5d836e5a7a4d9b2c1b964eb991f81cd9c/05-knowledge/results/liouville_mobius_eleven_20261027/THEORY.md) retains the same coprimality obstruction in μ(ab)=μ(a)μ(b)1_(gcd(a,b)=1), while λ is completely multiplicative. These signs do not encode our class:6 is II and210 is III, yet both have μ=λ=+1. Conversely8 and6 both have three proper divisors but have λ=-1 and+1. These are exact compatibility/loss statements, not an identification of the two carriers.

## 2. Square classes are a parity vector, not a single bit

Write Q(S) for the squarefree part of S, or equivalently the vector of prime-exponent parities of S over F2. Then

\[
Q(\sigma(ab))=Q(\sigma(a))+Q(\sigma(b))
\]

in that vector space, provided gcd(a,b)=1. Thus the product has square sigma iff the two input signatures agree. It is insufficient to know only whether each input sigma is square: two nonsquares can either cancel or fail to cancel.

The smallest relevant cancellation here is

\[
\sigma(2)=3,\quad \sigma(11)=12,\quad \sigma(22)=36.
\]

Both input signatures are the prime3. Meanwhile σ(3)=4 has zero signature. Multiplication by3 preserves the square-sigma predicate whenever3 does not divide the source. This is the inherited coprime square-sigma multiplication mechanism, specialized to a change of the user's source class:

\[
22\longrightarrow14\quad\mapsto\quad66\longrightarrow78.
\]

The sources are22=2*11 (II, three terms) and66=2*3*11 (III, seven terms); their totals are6² and12². The multiplier changes the total by4, rather than multiplying both endpoints by3. Exact endpoint transport is s(3n)=3s(n)+σ(n) when3 does not divide n.

A second prime triple lands on the same total:

\[
\sigma(70)=\sigma(2)\sigma(5)\sigma(7)=3\cdot6\cdot8=144
=\sigma(66)=3\cdot4\cdot12.
\]

So {2,3,11} and {2,5,7} give two distinct seven-term sources with a common square sigma, but different targets78 and74. For70 the square-class signatures3,6,2 also cancel. Keeping full prime-exponent parity explains the cancellation without treating all triples as interchangeable.

## 3. The entire square-total144 fiber is small and exact

\[
\{n:\sigma(n)=144\}=\{66,70,94,115,119\}.
\]

Their aliquot edges are respectively

\[
66\to78,\quad70\to74,\quad94\to50,\quad115\to29,\quad119\to25.
\]

The first two sources are III and the last three are II. There is no class-I source in this sigma fiber. All five are distinct-positive-vertex square-sum edges, with the same total12².

Completeness is certified by a direct finite argument: σ(n)>=n+1 for every n>1, so σ(n)=144 forces2<=n<=143. Independently summing every divisor of all143 possible positive sources therefore certifies the full list. This is a complete finite fiber, not an unbounded search assumption.

This also separates two minima in the target census:66 is the least class-III square-sigma source, whereas74 is the least target of a class-III square-sigma edge, coming from70. The bounded checks through the displayed witnesses establish these minima by checking every smaller source or every smaller target respectively; no monotonicity of s is assumed.

## 4. The seven-term classification closes one intersection exactly

Using [SEVEN_TERM_ABUNDANCE.md](SEVEN_TERM_ABUNDANCE.md), abundant seven-term sources are seven fixed exceptions plus6q for prime q>=5. Only70 among the exceptions has square sigma. For the family,

\[
\sigma(6q)=12(q+1)\text{ is square}
\iff q=12k^2-1\text{ for an integer }k\ge1.
\]

Indeed a square divisible by12 has root6u, giving q+1=3u²; q odd forces u even. Conversely the formula gives σ(6q)=144k². Thus this intersection is exactly

\[
\{70\}\ \cup\ \{6(12k^2-1): k\ge1,\ 12k^2-1\text{ prime}\}.
\]

The seven-term abundant family itself is unconditionally infinite, but infinitude of its square-sigma subfamily is not established. The quadratic-prime guard cannot be dropped. This subfamily is also exactly the multiplication-by3 lift of the even-target class-II square edges, with70 as the exceptional extra seven-term abundant source.

## 5. What this does for the1001 anchor

The verified1001→343→57 chain has two distinct descriptions: a weighted Boolean divisor box evaluates to343, and343 has square sigma400. The earlier1001+1001000=1001² edge instead comes from the universal companion n(n-1). The carrier explains multiplication and signature cancellation; the companion formula in [LATTICE_COMPANIONS.md](LATTICE_COMPANIONS.md) explains that second construction. Neither operation identifies those two edges or makes the aliquot map preserve every square-sum edge.

This is the proposed organizing principle: **factorization supplies an exact guarded composition; evaluation creates fibers; later predicates require the coordinates discarded by evaluation.** The same principle is useful in the tournament port work, but no pairwise orientation or tournament equivalence is claimed for the integers here.

## Checks

Run `python3 check_carrier_144.py`, also under `-O`. The checker exhausts the sigma144 fiber by direct divisor addition, verifies the least III square-sigma source below the witness66, checks coprime product/signature rules on1..100, and verifies the classified seven-term square intersection on all sources through100000. Its finite tests complement, rather than replace, the algebraic proofs and cited low-layer classification.
