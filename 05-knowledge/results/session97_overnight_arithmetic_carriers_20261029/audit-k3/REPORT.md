# Independent audit: K3 Frobenius obstruction and the 948 lattice

Status: **FINITE-EXACT arithmetic; standard geometric implications proved below**.
Audit scope: phase2 claims only. No claim of external novelty, generic rank17
for ICARM273/302, or geometric non-isomorphism of the three record surfaces.
No producer Python module is imported by either checker.

## 1. Good proper models: the missing geometric link is valid

Use the generalized Weierstrass equation over P1 with fundamental line bundle
L=O(2): a_i is a section of L^i. The projective cubic bundle supplies a proper
model, including the identity section. There are two base charts; on the
infinity chart the coefficients are u^(2i) a_i(1/u). Consequently the infinity
fiber uses the coefficients of degrees 2,4,6,8,12, not an affine large value.

For p=233, independently compute the discriminant of the cubic

    F(X)=4X(X+D)(X+E)+(LX+B)^2,   Delta=disc_X(F)/16.

Each Delta has degree24 and is squarefree both over Q and modulo233. Its
leading coefficient modulo233 is nonzero. The checker explicitly verifies
coprimality with c4 modulo233. Thus every geometric bad fiber is nodal I1,
not cuspidal, and the infinity fiber is smooth. The simple vanishing of the
discriminant makes the total Weierstrass surface smooth at each node: the
local smoothing parameter has a simple zero. Away from those nodes smoothness
is immediate; the section at infinity is smooth as well. Hence this is a
smooth proper surface over Z_(233) (checking both geometric fibers suffices
here); there is no hidden exceptional-component contribution to point counts.
Its generic minimal elliptic surface has chi(O)=2, trivial canonical bundle,
and irregularity0, so is K3.

## 2. Independently exhaustive point counts

`check_p233.py` has separately transcribed arrays for L,p,q,D,E. It uses the
cubic resultant discriminant, not the producer's b-invariant formula. For
every one of234 base points, it enumerates every x and every y in F233 in the
original generalized Weierstrass equation, and adds its one point at infinity.
No character shortcut or producer count code is used. Both normal Python and
`python -O` give identical receipts `p233.json`.

| surface | points | H2 Frobenius trace |
|---|---:|---:|
|273|58509|4219|
|302|58404|4114|
|4239|58374|4084|

The traces subtract 1+233^2. A Q-isomorphism induces a G_Q-equivariant
isomorphism on etale cohomology. Smooth proper base change at233 identifies
these representations with those of the good special fibers; the trace
formula then forces equal counts. Their different counts prove **pairwise
non-Q-isomorphism of the underlying K3 surfaces**, not merely different
fibrations. They do **not** exclude Qbar-isomorphism: rational forms can have
different Frobenius actions. No need to claim the Q-isomorphism extends to a
chosen integral model; the cohomological argument bypasses that issue.

Sources of coefficient arrays: ICARM curve pages
https://elliptic-rank.icarm.cloud/curve/273,
https://elliptic-rank.icarm.cloud/curve/302,
https://elliptic-rank.icarm.cloud/curve/4239.

## 3. Classic determinant948 model: independent arithmetic route

Primary source: Noam Elkies, arXiv:2608.25406v1, sections2.1–2.2,
https://arxiv.org/html/2608.25406v1. The source is explicitly separate from the
record surface4239. The auditor directly parsed and compared its saved primary HTML pre blocks.
The essential mathematical arrays are preserved in `elkies948_inputs.json`,
with versioned primary URL, access date, and SHA256 of that source snapshot.
`check_948.py` reads this compact input, not the HTML; it does not load the
producer's mathematical JSON or execute its checker.

The y-polynomials are reconstructed by triangular leading-coefficient
recursion for a polynomial square root, rather than squarefree factorization.
All17 exact curve identities and the separately published y1 match. For each
of136 distinct pairs, compute the reduced x-coordinate of P_i+P_j by the
group law. For reduced N/D, the homogeneous numerator degree is
max(deg N,deg D+4). The source's height lemma gives this degree as the height;
polarization gives the pairing (height(P_i+P_j)-8)/2. This independently
recovers the full published Gram, determinant948, and positive definiteness.
This route does not reuse the producer's coordinate-difference gcd formula.

The15 slopes all verify; the alternative point on the secant is on the curve
and distinct. Replacing any chosen P_i with -P1-P_i is an integral basis
change whose determinant is a product of15 signs. Thus all32768 choices keep
the generated lattice and determinant, while losing the named section choices.
This proves rank17 for the displayed subgroup, not its saturation without an
additional argument. The two choices really are distinguishable sections;
15 bits are required to specify the ordered choice tuple from the fixed slopes.
No relation to the count of supersingular primes follows from this cardinality.

## 4. Self-contained publication manifest and reproduction

Publish only these six files:

- `REPORT.md`
- `check_p233.py`
- `p233.json`
- `check_948.py`
- `948.json`
- `elkies948_inputs.json` (7963 bytes; mathematical facts, no paper prose)

Dependencies: Python3 standard library and **SymPy1.14.0**. All input paths
resolve relative to `__file__`; no network, cached HTML/PDF, producer modules,
or files outside this packet are required. Run both checkers normally and with
`python -O`; the receipts agree byte-for-byte. The primary HTML snapshot itself
is deliberately excluded from publication. Direct primary-source extraction
and comparison were performed during this audit, before packaging.
