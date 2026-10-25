# Session63 independent integral/Fourier audit

**Verdict: PASS for the claims and scopes below.** Prepared 2026-10-25.
The mathematics was reconstructed from the explicit statements; the exact
checker imports no producer code. No repository modifications or pushes.

## 1. Integral signed-orbit theorem

Let a finite group A act on a finite free integral lattice M by signed
permutations of specified basis vectors. For every orbit of basis lines choose
a representative s and let epsilon_s:H_s -> {+1,-1} be its stabilizer action.
Write a,b for the numbers of trivial and nontrivial stabilizer signs.

For every commutative ring R,

    (M tensor R)^A = R^a direct-sum R[2]^b,

up to the coordinate choices implicit in the displayed copies of R. The
natural map M^A tensor R -> (M tensor R)^A is injective and has cokernel
R[2]^b. Proof: an invariant on a transitive signed orbit is determined by one
coefficient r. Transport determines the other coefficients; consistency is
exactly r=epsilon_s(h)r for all h in H_s. A trivial sign imposes no relation;
a nontrivial sign imposes precisely 2r=0. Over Z the latter orbit contributes
zero. No averaging denominator or flatness of R is used.

The independent H1 proof is also valid. Use cocycles
c(gh)=c(g)+g c(h). For S=sum_g c(g), hS=S-|A|c(h); hence v=-S/|A| gives
c(h)=hv-v. It follows that

    H1(A,M) = (M_Q/M)^A / image(M_Q^A).

The orbit-coefficient calculation over Q/Z contributes Q/Z for a trivial
stabilizer sign, and (Q/Z)[2]=Z/2 for a nontrivial sign. Rational invariants
surject onto the first contribution and vanish on the second. Therefore
H1(A,M)=(Z/2)^b. The usual mod-p exact sequence

    0 -> M^A tensor F_p -> (M tensor F_p)^A -> H1(A,M)[p] -> 0

has connecting cocycle (g x-x)/p for an integral lift x. This is a separate
claim from the failure of an integral change-of-basis matrix modulo p.

For M=exterior^r Z[X], the stabilizer sign is the sign of its permutation of
the selected r-element subset. If A has odd order, every homomorphism from
any stabilizer to {+1,-1} is trivial. Thus all exterior invariant modules
commute with every commutative coefficient ring. Modular cancellation depth
is not the dimension of a new invariant space.

## 2. C3 multiplicative hostile

For cyclic permutation of e1,e2,e3, invariant bases in degrees 0,1,2,3 are
1; x=e1+e2+e3; y=e1 wedge e2+e2 wedge e3+e3 wedge e1; z=e1 wedge e2 wedge e3.
Each is a reoriented orbit sum. The only nonzero product of positive degrees
is xy=3z (and yx=3z). Expansion gives exactly three cyclically ordered terms
of sign +1. Thus the degreewise invariants do not change mod3, but that
product vanishes. The x-wedge complex is Z --1--> Z --0--> Z --3--> Z,
with only H^3=Z/3. This is not an invariant-basechange defect.
Its Hilbert polynomial (1+t)(1+t^2) has root depth 1 in characteristic zero
and 3, and depth 3 in characteristic 2; the latter also creates no invariant
basechange defect.

## 3. THM-4531: sign scope and two different integral lattices

Exact inherited source inspected:
`01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md`
and section 1 of
`05-knowledge/results/procgen_tbij_20261001_natural_bijections.md`.

THM-4531 uses endpoint-reversal parity, equivalently sign on the set of
2|E(X)| directed arcs. It is NOT the sign on the undirected edge set, and is
NOT in general the vertex-permutation sign. A single edge with its endpoints
transposed has endpoint sign -1 and edge-permutation sign +1. A single edge
with two isolated vertices transposed has endpoint sign +1 and vertex sign -1.
The primary source confirms the endpoint-reversal definition:
[Royle--Praeger--Glasby--Freedman--Devillers, arXiv:2204.01947v2](https://arxiv.org/abs/2204.01947),
published in J. Algebraic Combin.57 (2023),515--524. This check imports the
definition and bibliographic identity only, not a novelty claim for this audit.

Let P=Z^{T_U} be the honest point-function lattice, and let L be the Z-span
of the Fourier characters phi_X inside P. For d=dim_F2(F2^E/U), its evaluation
matrix is a Walsh matrix up to choices of row/column labels and signs.
For d>0, F^T F=2^d I gives [P:L]=2^(d 2^(d-1)); for d=0 the index is 1.
Modulo 2 this matrix has rank 1. Therefore a characteristic-zero Fourier
basis does not identify these two integral lattices at 2.

L is a signed permutation lattice with b the number of non-even graph orbits
in U-perp; H1(S_n,L)=(Z/2)^b. P is a genuine permutation lattice, so
H1(S_n,P)=0 and its invariants commute with every ring. These facts coexist
without contradiction. For U=0 the independent complete n=2,3,4 controls give
(a,b)=(1,1),(2,2),(4,7), respectively. The mod2 invariant dimensions of L
are 2,4,11, whereas those of P are 1,2,4.

## 4. Fourier lattice, positivity, and operation closure

Set N=2^r and H_{st}=(-1)^(s dot t) for F2^r. The integral Smith factors of H
are 2^j with multiplicity binomial(r,j). An elementary proof begins with

    [[1,0],[1,-1]] [[1,1],[1,-1]] [[1,-1],[0,1]] = diag(1,2).

The outer matrices are integral unimodular. Taking r-fold tensor products
gives diagonal entries 2^(popcount(s)); sorting gives a Smith diagonal since
successive entries divide one another. This simultaneously proves the index
formula and shows that the obstruction is purely 2-primary.

For a real table A, a nonnegative integer histogram n exists with A=Hn iff
HA is entrywise a nonnegative integer multiple of N. Equivalently A belongs
to H Z^N and C_{ab}=A_{a+b} is positive semidefinite. Indeed the character
vector s is an eigenvector of C with eigenvalue (HA)_s=N n_s.
The two conditions are independent: A=(1,0) for r=1 yields PSD identity C but
n=(1/2,1/2); A=H(-1,2)=(1,-3) belongs to the lattice but has a negative
histogram coordinate. This characterizes abstract count tables, NOT tables
realizable by a specified geometric cover, field, or branch arrangement.

If a response is tested against every point context, its full histogram is
determined by delta contexts. No fixed omission of a nonconstant Fourier
coordinate suffices universally, even at fixed mass: compare the nonnegative
integer histograms 1 and 1+chi_s. Their transforms differ only at s. If total
mass is already supplied, the constant Fourier coordinate is of course not
additional missing information. This is a statement about universal exact
context determination, not about minimal data for a restricted family.

For pointwise gluing f*g (here the star denotes ordinary pointwise product),
the Fourier table is normalized convolution:

    H(f*g)_s = (1/N) sum_t (Hf)_t (Hg)_(s+t).

Thus the entire integral response lattice and nonnegative cone are preserved
under this specified gluing operation. The factor 1/N is load-bearing.

## 5. Exact controls and receipts

Reproduction:

    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json

**FINITE-EXACT:** 521 checks PASS; normal and optimized outputs are identical.
The checker includes 112 signed-involution/ring cases (all signed involutions
of ranks 1..3 over Z/2,Z/3,Z/4,Z/6), all graph and tournament orbits n=2..4,
40 Fourier-table product/eigenvalue controls at r=0..4, explicit unimodular
Smith reductions at r=0..4, and exact Legendre counts.

For the branch-free curve y^2=d x(x-1)(x-lambda), with d=1 and the least
positive nonsquare, literal finite-field enumeration gives:

| Field | lambda | ordered twist counts | mean | variance |
|---|---:|---|---:|---:|
| F5 | 2 | (4,0) | 2 | 4 |
| F5 | 3 | (0,4) | 2 | 4 |
| F7 | 2 | (4,4) | 4 | 0 |
| F7 | 3 | (0,8) | 4 | 16 |

The F7 pair has equal branch-incidence data and mean but different variance;
the F5 pair has equal mean and variance but different labeled profiles.

For the regular (C5)^2 subset-orbit polynomial, the independently computed
coefficients at t=-1 have gcd filtration
Delta_1..5=0, Delta_6..13=3000, Delta_14..17=600,
Delta_18..21=120, Delta_22..25=24, Delta_26=1.
This chain was first independently computed using JavaScript BigInt and then
recomputed in the Python checker; complete jets are in normal.json.

Final script SHA256:
`b0fd50eeaa52594151744337c3eec221a0a4648a4982f155c731524ff3b64050`.
Final normal/optimized JSON SHA256:
`7a34f738f9ed0fa76a68d9d88f29bab906d661eec10004aaab0c6dd2cbd8f625`.

This audit does not infer geometric realizability from count positivity,
new invariant dimensions from modular root depth, or integral basis identity
from a characteristic-zero Fourier change of basis.

## 6. Final producer read and replay

The full authoritative `session63-specialization/RESULT.md` and
`session63-root/FOURIER_EFFECTIVITY.md` were read after execution service
recovery. **PASS:** no mathematical correction required. In particular the
specialization note's modulo 5,25,125,625 depths 21,17,13,5 follow from the
independently computed gcd chain, and its C3 product statement explicitly
includes both xy and yx.

The root's `check_fourier_effectivity.py` was replayed separately under
ordinary Python and `-O`: **7552 checks PASS**, outputs byte-identical,
recorded in `root_fourier_replay.json` and
`root_fourier_replay_optimized.json`. This replay is additional to the
independent 521-check script, not a claim of a second independent derivation
of the root's exhaustive table loop.

Primary PDF checks also confirm the producer's exact source locations:
Royle et al., section 4, Lemma4.1 identifies the endpoint sign as a
homomorphism; Brown's
[author-hosted lectures](https://pi.math.cornell.edu/~kbrown/papers/cohomology_hangzhou.pdf),
section3.3, printed page16 gives the coefficient cohomology long exact
sequence, and Proposition3.3, printed pages19--20 gives Eckmann--Shapiro.
The direct proofs above do not depend on these citations.

Read-only repository status at the end of this audit was clean on
`codex/session63-unifying-theories...origin/main`. Publication and final
upstream verification belong to the parent task; this lane changed only its
unique scratch directory.
