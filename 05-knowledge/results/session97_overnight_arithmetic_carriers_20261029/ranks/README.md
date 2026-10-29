# Session97: high-rank curves, exact local certificates, and the meaning of 21

Status: **CITED / FINITE-EXACT / PROVED elementary certificate lemma / HEURISTIC scope separated**.
Audit date: 2026-10-29. Repository inheritance pin supplied by parent:
`5d6eb91dd76eee86b91e04d856daa32c32b3d61d`.

## Main result of this pass

We independently verified the rank lower bounds 30, 31, and 32 for the three
public record curves 273, 302, and 4239 using exact rational point membership
and full-rank binary matrices of local quadratic characters. This is not merely
a floating regulator or a check that the listed points lie on the curves.

| ICARM curve | independent points | largest certificate prime | prime excluding rational 2-torsion |
|---|---:|---:|---:|
| [273](https://elliptic-rank.icarm.cloud/curve/273) | 30 | 269 | 59 |
| [302](https://elliptic-rank.icarm.cloud/curve/302) | 31 | 347 | 59 |
| [4239](https://elliptic-rank.icarm.cloud/curve/4239) | 32 | 269 | 7 |

These are **lower bounds**, not unconditional exact ranks. The curve pages
report exactness conditional on GRH and BSD. Curve4239 is a specialization of
a generic-rank17 fibration on the K3 surface of discriminant1152, credited to
Elkies, Klagsbrun, and Li. Its posted explicit-formula certificate gives an
analytic upper bound below34 and root number+1. We have not rerun that analytic
certificate. [Record4239](https://elliptic-rank.icarm.cloud/curve/4239)
[Dujella's primary record history](https://dujella.github.io/tors/rankhist.html)
separately lists the highest unconditionally exact known rank as22 (Li,
September2026). Do not conflate that category with the record lower bound32.

The website's current API also describes exact local-character certification;
our implementation and verifier were written independently from its code.
[ICARM API](https://elliptic-rank.icarm.cloud/api)

## Why the finite certificates prove independence

Let E/Q have integral Weierstrass coefficients a1,a2,a3,a4,a6. Put
b2=a1²+4a2, b4=2a4+a1a3, b6=a3²+4a6 and
f(x)=4x³+b2x²+2b4x+b6. At a good odd prime p, a root r of f gives
an Fp-rational point of order2 on the completed-square model.
The standard quadratic-character map attached to that point sends an ordinary
point P to the Legendre symbol of x(P)−r. In additive F2 notation it is a group
homomorphism. The value at the special order2 point uses the derivative of the
monic completed-square cubic, and the identity maps to0.

For completeness, the homomorphism follows from the chord rule: substituting
a line into the monic cubic and evaluating at r shows that the product of the
three x-coordinate differences from r is a square. Tangencies and the special
point follow by the corresponding limiting algebraic identities. Composition
with good reduction gives E(Q)→F2. The verifier avoids special values entirely
in its chosen columns: all denominators are p-units and x(P)≠r.

If f has no root at another good odd prime, E(Q) has no rational2-torsion:
prime-to-p torsion injects under good reduction. Its rational torsion therefore
has odd order, so every character kills it. A full-column-rank character matrix
implies any relation modulo torsion has all coefficients even. Divide by2;
because torsion has odd order the resulting half-relation is still torsion.
Repeat. The integer coefficients must be divisible by every power of2, hence
all vanish. Thus the listed points are independent modulo torsion.

This proof needs neither BSD nor GRH. It does not prove saturation at odd
primes, that the listed points are a full Mordell–Weil basis, or a rank upper
bound. Mod2 observation loses sign: replacing a point by its inverse correctly
preserves every character. This is an exact instance of a small certificate
being sufficient for independence but insufficient for faithful group-state
reconstruction.

## Source-preserving family transport

`check_specializations.py` evaluates each page's displayed fibration at its
stated rational parameter, proves c4/c6/discriminant scaling, and transports
all93 named points through the resulting rational Weierstrass isomorphism.
It retains parameter, scale, x-translation, y-shear, and y-translation.
The explicit transports are in `specialization_certificates.json`.

All three fit

    y²+(Lx+B)y=x(x+D)(x+E),
    B=pE+qD+pq(L−p−q),

with quadratic L,p,q and quartic D,E. Curve302's minus-sign presentation must
be normalized by replacing L with −L, **not** also replacing p,q. This was
caught by the exact j/scaling check during implementation and repaired before
filing. This common construction architecture is a genuine relationship;
it does not identify these different fibrations or prove their generic ranks.
The three exact parameter values and scaling factors are:

| curve | T | u in x=u²X+r |
|---|---|---|
|273|12473/32946|1/542719458|
|302|164518/924945|1/285174417675|
|4239|−2868768/362209|4/131195359681|

The public snapshot has6127 entries and66 with stated lower bound≥30.
Those66 have66 distinct rational j-invariants: no two are quadratic twists
of one another. Distinct j does **not** exclude isogeny or membership in a
common parameter family. All66 have empty built-in family tags, but these
only test special j=0/1728 families; empty tags do not mean no construction
family. The leaderboard is selected by rank/height/conductor/discriminant
records, not an unbiased sample. No density regression is justified.
[Dataset documentation](https://elliptic-rank.icarm.cloud/api)

## What rank21 does and does not mean

PPVW predicts finitely many Q-isomorphism classes of rank>21; it does not
predict none. The proof is for their random-matrix model, not actual elliptic
curves. Its threshold balances H^(5/6) curves with rank-tail probability
H^(-(r−1)/24+o(1)), giving21=1+24·5/6. Rank21 is a boundary case.
The separate minimalist prediction is asymptotically half rank0, half rank1.
The paper explicitly warns about high-generic-rank parametric families and
says that no natural alternating matrix attached to an elliptic curve is
known in this model. These qualifications rule out treating tournament skew
matrices as an established arithmetic realization of Sha or rank.
[PPVW, Theorem1.1.1, Remarks1.1.2–1.1.4, Sections7.3,8](https://math.mit.edu/~poonen/papers/bounded-ranks.pdf)

### A useful proposed family diagnostic, not a theorem about actual ranks

Suppose a family has generic rank r0, N(H)=H^(α+o(1)) distinct specializations
up to Q-isomorphism, and one **assumes** the shifted tail law
P(rank≥r0+s)=H^(-(s−1)/24+o(1)). Its expected count has exponent
α−(s−1)/24. Thus the family-specific critical total rank is

    rcritical = r0 + 1 + 24α.

This elementary deduction makes the required missing data explicit: count
DISTINCT curves by a declared arithmetic height, establish generic rank,
control exceptional parameter maps, and justify a family tail law. None is
supplied by matching the numeral21 with tournament exclusions{7,21}.

The honest bridge to Session96 is alternating forms, symplectic cokernel
pairings, parity and retained source maps—not equality of rank thresholds
with Hamilton-path forbidden values. A candidate stronger connection must
produce a map preserving the relevant predicate and survive independent
rank/height and path-count controls.

## Reproduction and controls

Run from this directory:

    python check_characters.py
    python verify_certificates.py
    python check_specializations.py

Normal and `python -O` runs agree; checks use explicit exceptions rather than
optimization-disabled assertions. The independent verifier performs13851
checks, rejects12 hostile mutations (duplicate point, off-curve coordinate,
bit flip, wrong local root across three curves), and correctly accepts three
point-sign reversals. Output: `verification.json`.

`analyze_families.py` requires the public `database.json` snapshot, downloadable
from https://elliptic-rank.icarm.cloud/database.json. Its SHA256 and attribution
are retained in `highrank_family_audit.json`. Do not publish the downloaded
PPVW PDF or the large raw database as research outputs; publish the scripts,
three small factual curve snapshots, certificates, and analysis instead.

Data acknowledgement: NSF Institute for Computer-Aided Reasoning in
Mathematics (ICARM), NSF Grant DMS2425401. Source URLs remain attached to the
research and the family-audit metadata.

## Next decisive work

1. Independently verify the family-level rank17 claim using a named section
   height matrix/intersection pairing, rather than inferring it from one
   specialization. This would quantify the actual rank-jump mechanism.
2. Count distinct j-specializations of the same named fibration with a fixed
   rational parameter-height budget, retaining duplicate/isomorphism guards;
   do not substitute leaderboard entry counts.
3. Test whether local-character certificates transport under the named
   rational isomorphisms and under genuine isogenies, retaining kernels and
   bad-prime exclusions. This links the repo's source/port discipline to a
   concrete arithmetic certificate, without a numerological rank21 claim.

## Second checkpoint: certificates survive the actual coordinate change

`check_character_transport.py` transports every named local root by
`r_family=u² r_minimal+r_translation` at primes where all denominators and u
are units. Since `x_family−r_family=u²(x_minimal−r_minimal)`, the quadratic
characters are unchanged. All93 certificate rows survive with no excluded
primes in these examples. Forgetting the root translation gives465,463,494
wrong observation bits respectively. This is a concrete, source-preserving
arithmetic analogue of retaining shared ports: the curve isomorphism alone
is not enough unless the observed local root is transported with it.

**Additional Sha guard:** do not identify symplectic Lagrangians with maximal
isotropics for a chosen quadratic refinement. In dimension4 over F2 the
former count is15 while the split-quadratic latter count is6. The incoming
repository `F1_QUADRATIC_BOUNDARY.md` at commit1151efc0b supplies this distinction;
BKLPR's Selmer model uses split quadratic data. Thus the symplectic15 count is
not automatically the appropriate finite sample space for a rank/Selmer
claim. The local-character certificates here do not require that identification.
