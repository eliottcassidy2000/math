# Finite affine representations and the missing positivity under refinement

2026-10-05. **PROVED:** the exact residue-closure obstruction, colored flux
repair, arbitrary finite guarded lifts, opposite-drift aliases, Fourier
orbit decomposition, and the summable-measure refinement criterion below.
**INHERITED PROVED:** the corrected affine group and its rooted lifting
theorem. **FINITE-EXACT:** the explicit implementation universe.
**OPEN:** positivity of the Collatz weight at every source. No Langlands or
modularity theorem is asserted to supply that positivity.

[Program](../../04-computation/experiments/collatz_representation_positivity_20261005.py)
and [saved output](collatz_representation_positivity_20261005.out).

## 1. Inheritance and the objects that must remain distinct

The closest proved mechanism is
[affine guarded lifts, sections 2–4](collatz_affine_guarded_lifts_20261004.md):
the exact finite affine group, its correlated-CRT correction at modulus95,
and actual ancestors realizing every finite group element from a supplied
3-unit hub. The canonical hostile is its formal inverse E(1)=1/3:
being a permutation modulo an odd modulus does not establish integer
legality. The corrected near miss is to promote positive mass in every
finite residue class to a positive atom at every integer. The least-used
sidecars are the exact valuation and the lower bound through arbitrarily
fine refinement.

The current weight is the incoming package's W in
[adaptive mixture flow, P2–P3](collatz_adaptive_mixture_flow_20261005.md).
It is a nonnegative summable sequence on positive odds, is positive
exactly at sources with a finite ROOT certificate, and obeys the incoming
inequality with root deleted. Its improved norm bound and beta formulas
are inherited, not needed to prove new positivity here.

Other nearby mechanisms are [THM-4501, recursive motif families and
frequency](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md),
which supplies dyadic density with no given-source conclusion, and
[finite lookahead obstruction](collatz_finite_lookahead_weight_obstruction_20261005.md),
which excludes its stated finite observer as a globally positive paid
weight. Neither says every finite algorithm is insufficient.

The board is **affine map / valuation guard / representation / incoming
flux / positive residue mass / refinement lower bound**. The anchor is
positivity; the niche is an exact representation retaining the translation
phase; the wildcard is whether a faithful finite representation pays the
source predicate. Its cheapest hostile is a pair of identical finite
matrices with opposite ordinary drift.

Here m always denotes an odd modulus coprime to6, not a sibling depth,
an odd-source index, or a modular-form level. The representation acts on
the m residues of the actual integer n.

## 2. A faithful finite representation, with an exact arithmetic boundary

Fix m>1, gcd(m,6)=1. For a positive valuation a write

    F_a(x)=(3x+1)/2^a,
    f_a(r)=2^(-a)(3r+1) mod m.

The rational affine map F_a is an actual U-step only on

    C_a={n>0 odd: 3n+1 = 2^a mod 2^(a+1)}.             (1)

Let h=ord_m(2). The finite permutation f_a depends on a modulo h.
The inherited group theorem gives

    H_m=<2,3> in (Z/m)^*,
    G_m={(u,v):u in H_m, v in Z/m},
    (u,v)(u',v')=(uu',uv'+v).                         (2)

For completeness, f_1 f_2^(-1) is multiplication by2, and f_1^(-1)
is E(x)=(2x-1)/3. The commutator of multiplication by2 with E is
translation by -1/3, which generates all translations modulo m.
This recovers (2) directly. The positive generators also generate their
inverses as powers, since the group is finite.

Define R(g)e_r=e_(g(r)) on C^(Z/m). This permutation representation is
faithful: the values of g at0 and1 recover v and u. This is faithfulness
to the finite affine map; it does not retain the original integer, the
exact exponent, the dyadic domain, or a route home.

**Exact residue-closure theorem.** The actual value U(n) modulo m is
constant over the positive odd source fibre n=r modulo m if and only if
3r+1=0 modulo m. In that exceptional fibre it is always0. In every other
fibre the next residue is ambiguous.

**Proof.** For every a>=1, the exact cell (1) intersects every residue
class modulo m, in infinitely many positive integers, by coprime CRT.
If 3r+1=0, every branch has target0. Otherwise f_1(r) differs from
f_2(r), since their difference is (3r+1)/4, nonzero modulo m. This
proves both directions, including composite m. Thus no uncolored
deterministic residue transition represents U on all these integers.

**The source-side repair.** For a finitely supported mass p on positive
odds, retain the joint masses

    p_c(r)=sum_(n=r mod m, v2(3n+1)=c mod h) p(n).

Then the exact residue pushforward of the full incoming U-operator is

    pi_m(U_*p)=sum_(c mod h) R(f_c) p_c.              (3)

For c=0 use any positive exponent congruent0, such as h. The ordinary
residue marginal sum_c p_c does not determine (3). The color records
only a finite valuation phase; its next color is not determined by this
state in general. To use a killed-root incoming operator, also retain
the literal source/target ROOT flag and delete the corresponding edges.
Residue1 is not the literal integer1.

**Same matrix, opposite drift.** The letters a=1 and a=1+h have
identical full permutation matrices, while their positive real slopes
are respectively 3/2>1 and 3/2^(1+h)<1. Here h>=2 because m>1.
Both exact source cells contain every residue modulo m; the cells are
disjoint, so this never asserts two different actual valuations at the
same integer. At m=11 the slopes are3/2 and3/2048. Finite matrix
faithfulness therefore cannot detect ordinary contraction after the
exact exponent has been reduced.

## 3. Every finite labeled path lifts, and the omitted sidecar is unbounded

For a positive valuation word w=(a_1,...,a_t), put A=sum a_i and

    F_w(n)=(Pn+B)/Q,  P=3^t, Q=2^A.

Chronological appending of a updates (P,Q,B) to (3P,2^a Q,3B+Q).
The complete exact source guard is the one odd residue

    Pn+B = Q mod 2Q.                                  (4)

Indeed the cleared final identity with an odd endpoint first forces
every earlier division to be integral. At an intermediate even state,
the next numerator would be odd, contradicting the next division; the
last endpoint is odd by (4). Hence every division has its stated exact
valuation. This is the inherited final-oddness argument with its last
bit retained.

For any residue r modulo m, CRT between (4) and n=r modulo m gives a
unique class modulo2mQ. Adding sufficiently many periods makes every
prefix exceed any prescribed positive height H: it suffices to take
n>HQ, since each prefix is at least n/Q. Thus every finite labeled
residue path has infinitely many actual positive lifts, even with all
nodes above H. This theorem constructs a source for the path. It does
not certify a separately supplied source, turn a terminal residue1 into
ROOT, or realize a single positive integer for an infinite path.

The minimal requirements depend on the requested output:

| Source -> target map | Preserved predicate | Lost information | Required sidecar |
|---|---|---|---|
| Affine (P,Q,B) -> G_m | Composition and exact residue action | Ordinary slope, carry size, word, source | Exact coefficients and valuation word |
| G_m -> permutation matrix | Entire finite affine map, faithfully | Nothing further about that map | Still no integer guard |
| Actual source -> residue | Congruences | Actual next residue in m-1 fibres | Valuation phase for one step |
| Guarded word -> finite residue path | All marked finite residue edges | Which integer realizes them | Source plus (4), unbounded A |
| Joint source/valuation mass -> residue flux | Exact one-step incoming mass | Future valuation phases | Refine from actual sources again |
| ROOT certificate -> W atom | Actual rooted positivity | Certificate cannot be replaced by residue | First-hit word and exact source |

These are composition-preserving reductions and labeled lifts. They are
not an equivalence of the finite residue group with the category of
actual certified integer paths. In particular, group inverses need not
be legal integer inverse steps at a supplied source.

## 4. Fourier phases, irreducible blocks, and the two-dimensional boundary

For chi_k(r)=exp(2 pi i kr/m), the Koopman action is exactly

    chi_k o (r -> ur+v) = exp(2 pi i kv/m) chi_(uk).
    chi_k o f_a = exp(2 pi i k2^(-a)/m) chi_(3·2^(-a)k). (5)

Equivalently, the Fourier transform of a pushed mass has the same phase
and frequency formula. Thus the finite representation becomes a
monomial action in the character basis. Omitting the phase in (5)
would discard translations and the ordered affine carry.

**PROVED decomposition.** Its irreducible invariant subspaces over C
are the spans of characters in the H_m-orbits on Z/m. To prove this,
translations diagonalize with distinct characters. Their finite Fourier
averages project onto each character line. Any invariant subspace must
therefore contain the individual character components of its vectors;
the multiplier group transports one component to its entire orbit.
This proves irreducibility on each orbit and the stated direct sum.

At m=11, 2 has order10, so the decomposition has dimensions1+10.
More strongly, **every complex representation of G_11 with dimension
less than10 annihilates the translation subgroup**: a nontrivial
translation eigencharacter has ten different conjugates under H_11,
all of which would have to occur. This is a precise obstruction to
using a two-dimensional complex representation as a faithful carrier
of this finite affine action.

At m=95 the dimensions are1+4+18+36+36. The two36-dimensional blocks
retain the inherited correlated-CRT distinction: |H_95|=36, while
|H_5|·|H_19|=72. Replacing the global multiplier group by that full
local product changes the representation.

This elementary finite representation is not an elliptic-curve Galois
representation or an identification of Collatz letters with Hecke
operators. Already the chronological words12 and21 have different
maps (9n+5)/8 and (9n+7)/8. A multiplicative encoding into a commutative
algebra identifies those two products and cannot retain the carry.
Any external correspondence would need its own lawful map, hypotheses,
and retained source predicate; matching dimensions or the number11
does not provide one. The statements of this section are proved above
without importing a representation-theoretic convergence theorem.

## 5. Positive finite quotients do not supply positive integer coefficients

For any nonnegative summable weight w on the positive odds, define

    mu_m(r)=sum_(n=r mod m) w(n).

For the actual W, **mu_m(r)>0 at every residue of every m coprime to6**.
This follows from the inherited guarded all-element lift: start with
the certified hub1, choose a finite affine translation carrying1 to r,
and obtain a positive odd ROOT ancestor in that residue. Its W atom
is positive. This claim needs only that one certified representative;
it does not say every integer in the residue is rooted.

Let muhat_m(k)=sum_r mu_m(r) exp(2 pi i kr/m). The circulant matrix

    C_(s,t)=muhat_m(s-t)

has eigenvalues m mu_m(r). Therefore every such finite W-matrix is
strictly positive definite. This is a genuine positive result about
the finite Fourier quotient, not the missing pointwise positivity.

**Exact hostile.** Put q(3)=0 and q(n)=2^(-n) at every other positive
odd integer. Its total mass is13/24. Every residue of every odd modulus
has positive q-mass, hence every associated Fourier circulant is
strictly positive definite. Nevertheless the coefficient at3 is zero.
For odd m>3 its residue mass at3 is exactly

    mu_m(3)=1/[8(2^(2m)-1)] > 0,                       (6)

which tends to zero. This is not a proposed Collatz flow or a Collatz
counterexample. It refutes the logical implication from all finite
residue positivity, even strict matrix positivity, to atom positivity.

**Lossless numerical repair.** The complete unbounded family of
numerical pushforwards does determine a summable w. Extend w by zero
at positive even integers. For a fixed n and every m>n,

    0 <= mu_m(n)-w(n) <= sum_(j>=m) w(j).              (7)

Every other positive integer congruent to n is at least m, which proves
the bound. The right side tends to zero by summability. One nested
tower, for example m=5^j, already suffices. Thus losslessness of the
full tower is compatible with failure of every finite positivity test.

For fixed n, the exact missing condition is

    w(n)>0  iff there exist rational epsilon>0 and j0
             such that mu_(5^j)(n mod 5^j)>=epsilon
             for every j>=j0.                         (8)

The forward direction uses w(n) itself as a lower bound; the reverse
direction takes the limit in (7). An effective summable tail makes (7)
a quantitative reconstruction algorithm. It supplies an upper error
bound, not the positive epsilon in (8). The already proved computable
weight constructions can therefore be fully compatible with a still
unresolved nonvanishing question.

This is the useful categorical boundary: faithful finite affine maps,
exact character phases, and consistent numerical refinement can all
be retained, while a strict positivity predicate is lost on passage to
the limit. No unproved orbit is supplied by changing representations.

## 6. Exact computation and the next admissible obligation

Run:

```text
python -B 04-computation/experiments/collatz_representation_positivity_20261005.py
python -B -O 04-computation/experiments/collatz_representation_positivity_20261005.py
```

The independent finite universe is m in {5,7,11,19,35,95}; complete
affine-group closure from the two accelerated generators; complete
multiplier closure from2 and3; permutation faithfulness; all residue
fibres and both opposite-drift exact cells; Fourier exponents at every
frequency for three distinct source residues; and the full frequency
orbit decomposition. Eight declared valuation words, every residue,
and three lifts each give4,128 literal prefix controls above height17.
The finite-mass color identity is checked independently on all128 odd
sources below256. There are172 explicitly replayed ROOT representatives,
one per residue for each modulus. This finite list does not prove the
all-modulus rooted lifting theorem, which is inherited with its proof.

The hostile measure is summed by exact geometric series at seven
declared moduli, not estimated by a cutoff. Thirteen exact-input hostiles
check bool/float-equivalent and malformed domain boundaries. The output
records26,980 checks, identical normally and under optimized Python.

The next positivity obligation must therefore retain an actual source
and produce a lower bound stable through refinement, or transport an
independently verified ROOT certificate to it. A positive finite
eigenvector, a solvable residue group, or a full-support residue marginal
alone does not discharge either obligation.
