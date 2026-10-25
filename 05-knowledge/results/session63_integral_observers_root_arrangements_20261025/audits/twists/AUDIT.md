# Session 63 independent audit: twists, quotients, and effectivity

2026-10-25. **Verdict: PASS for the stated scope.** No mathematical failure
found in the audited proofs. Finite tests are controls, not proofs for all
fields or extension degrees. No repository edits or publication were made
by this audit lane.

## 1. Audited artifacts and boundary

The producer packet is `session63-twists/{RESULT.md,ROOT_STACK.md,
ZETA_AND_EFFECTIVITY.md,verify.py}`. The independent recovery packet is
`session63-recovery/{check_twist_histograms.py,twist_finite_normal.json,
twist_finite_optimized.json}`. The inherited complete-intersection/A1 and
nonrationality arguments were read in Session59's
`toric-shared-roots/TORIFICATION_OBSTRUCTION.md`; the complementary
six-element canonical characters were checked against
`six_seven_gale_cube_20261025/CHARACTERS.md`.

This review includes the producer's general character-norm and complex
line-complement classification addenda. It does **not** certify the root
agent's separate integral Fourier-effectivity theorem or any result not
listed here.

The three relevant carriers must remain distinct:

| Source | Quotient | Fresh mean, characteristic not 2 or 3 | Formal value at q=1 |
|---|---|---|---|
| Branch-free root cover V | U = P² minus the seven lines | q²−6q+12 | 7 |
| Whole singular root surface X | P² | q²+q+1 | 3 |
| Minimal resolution Y | Blowup of P² at the three triple points | q²+4q+1 | 6 |

In characteristic 3, the open mean is q²−6q+10 and the resolved mean is
q²+6q+1; the singular mean is unchanged. Characteristic 2 is excluded.

## 2. Proof audit

### Torsor, normalization, and nonfree mass

The six ratios L_i/E have independent square classes: a selected numerator
line has odd valuation in any nonempty product. Thus the degree is exactly
64, not merely bounded by 64. On U every ratio is a unit and adjoining its
square root is étale in odd characteristic. The projective sign group is
constant of order 64 there. The twists are *marked G-twists*, so accidental
isomorphisms between their underlying varieties do not remove multiplicity.

The local square-root charts are smooth over the zero-, one-, and two-line
strata. At a triple intersection, the remaining equation is the nondegenerate
quadratic node w²=a u²+b v², with a,b nonzero. Consequently X is already
normal; no unidentified normalization is being substituted in the counts.

The nonfree mass formula is valid. Each Frobenius-stable geometric G-orbit
contributes exactly |G| to the sum over g of fixed points of gF, because
orbit size times stabilizer size is |G|. This counts coarse quotient points
without requiring rational lifts in the original W. The stabilizer orders
on the four strata are 1, 2, 4, 8, respectively. Each base point has 2^r
compatible twists and 2^(6−r) rational points in each compatible reduced
fiber. Nonreduced scheme fibers are not mistaken for étale fibers.

All branch intersections are rational over the given odd base field: they
are intersections of its own distinct rational lines. Deconing proves
(t₂,t₃)=(12,3), except in characteristic 3 where the two additional affine
triples produce (6,5). There are no fourfold points in either regime.

### Resolution quotient and rationality of exceptional curves

The node blowup is the minimal resolution and is equivariant. Its u-chart
has ring k[u,V,W]/(W²−a−bV²), where V=v/u and W=w/u. The three transformed
signs are independent, so the invariant ring is exactly k[u²,V²]; W² is
eliminated by the equation. The map (x,t)=(u²,V²) is the x-chart of the
blowup of the base origin. The analogous v-chart has the usual overlap.
These charts cover because u=v=0 forces w=0 in the exceptional homogeneous
conic. The constructions descend from the étale root charts and agree away
from the triple locus. Hence Y/G is the stated blowup of P².

The independent point-count proof also passes: each rational node is
replaced by a smooth finite-field conic with q+1 points, adding q; a point
over a nonrational node cannot be rational over the smaller field. The
argument does not incorrectly declare every conic split over an arbitrary
base field. The generalized correction q·#T(F_q) is correct.

Over C, direct stratum arithmetic gives

    χ_c(V)=64·7=448,
    χ_c(X)=64·7+32·(−19)+16·12+8·3=56,
    χ_c(Y)=56+24=80.

The quotient Euler characteristics are 7, 3, 6, not those source values.

### Fixed twists, zeta series, and ordinary motives

For constant abelian G, (gF)^m=g^mF^m. The power map G→G^m has uniform
fibers, proving that the fixed-family degree-m mean is #(W/G^m)(F_(q^m)).
This is not generally #(W/G)(F_(q^m)). Constancy and abelianness are used
explicitly; neither is silently omitted.

For exponent two the quotient count survives at odd degrees, while at even
degrees all fixed twists become the original W. This proves the zeta
product with exponents 64 and 32. The factor 32, rather than 64, comes from
the logarithmic denominator 2r at degree 2r. The general character identity

    ∏_g det(1−χ(g)tF) = det(1−t^d F^d)^(|G|/d)

is valid for a character of order d; it requires no diagonalizability of F.
Its cohomological use has alternating exponents and splitting coefficients.
At extension m, character aliasing χ↦χ^m retains every character trivial
on G^m, not merely the globally trivial character.

The revised motive-hostile paragraph is precise: over F5→F25 all fixed
forms become isomorphic to V, but V has 832 points and U has 487. Point
counting therefore already disproves the proposed ordinary averaged class
identity over finite fields. The complex Euler counterexample is separately
identified as characteristic zero, not as an impossible F_q→C base change.

### Fourier information and torification scope

The exact bin identity N_ε=64 n_ε, Walsh inversion, variance, and all shift
covariances are correct. The mean removes 63 nontrivial character sums;
covariances recover their squares, not their signs. The support inequality
follows directly from Cauchy–Schwarz when the denominator is nonzero (in
the hourglass finite-field cases the complement is nonempty).

The coefficient −4 in the (L−1)-basis obstructs finite split-torus
stratifications of U. For all 64 *fixed* finite-field twists, odd extension
degrees suffice to obtain a polynomial identity and hence a contradiction
to simultaneous torifiability. That average alone does not prove every
individual twist non-torifiable. The separately inherited complex p_g=7
argument does prove individual complex non-torifiability by birational
nonrationality, including every nonempty open.

The positive resolved quotient is correctly handled constructively, not by
positivity alone: since the three good-characteristic centers lie on one
line, its blowup decomposes as A² disjoint union P¹ disjoint union three A¹.
This refines to one two-dimensional torus, six one-dimensional tori and six
points. It says nothing similar about the source Y.

The general n-line classification is correct: the torus-basis coefficient
of degree one is 3−n; three concurrent lines instead have E=T²−1.
The only split-torifiable complements are n≤2 and a nonconcurrent triple.
The constructive positive cases are P², A², A¹×G_m, and G_m². The two
hourglass occurrences of seven are not a universal equality: seven general
lines have χ_c(U)=10 while the root resolution still has p_g=7.

### Root stack and retained algebra

The product of the seven divisor-root stacks agrees with [X/G]. On L_j≠0,
the jth root object is trivial and the other six roots have precisely the
normalized root-cover quotient chart. Transitions r_i→r_i/r_k agree with
the line-bundle data and projective signs. Stabilizers are (μ₂)^r.
The triple-point atlas is singular, so the stack is not spuriously declared
smooth. Roots without sections are correctly distinguished.

The square-free monomial multiplication retains the branch functions;
replacing it by only symmetric-difference labels loses the algebra. The
global canonical calculation gives seven complementary six-element
characters, not singleton characters. No F1 descent, preferred F1 point
functor, integer-square solution, or RH implication follows.

## 3. Executed independent controls

`independent_checks.py` imports no producer code. It enumerates all
normalized projective coordinate vectors in P⁶(F3) and P⁶(F5), evaluates
the four defining quadratic equations for every one of the 64 twists, and
independently tests the four-row Jacobian rank for singular points.
Every labeled compact/open/node bin agrees with the recovery lane's
root-table implementation, not just their means.

| Field | Ambient P⁶ points enumerated | Sum of 64 compact counts | Sum of open counts | Sum of Jacobian-singular counts |
|---|---:|---:|---:|---:|
| F3 | 1,093 | 832 | 64 | 320 |
| F5 | 19,531 | 1,984 | 448 | 192 |

A separate explicit finite G-set of order-64 deck group, with stabilizers
of orders 1, 4 and 64 and commuting Frobenius, checks fixed-family odd/even
counts and the zeta-log identity through degree 16. The example has a
genuine nonzero correction, so a missing factor cannot pass vacuously.

Both normal and `python -O` runs passed and produced byte-identical outputs.
The test uses explicit `require` guards, not optimization-disabled asserts.
To reproduce from the scratch audit directory:

    python independent_checks.py --reference ../session63-recovery/twist_finite_normal.json > normal.json
    python -O independent_checks.py --reference ../session63-recovery/twist_finite_normal.json > optimized.json
    cmp normal.json optimized.json

The reference argument permits relocation of the published packets.

The producer's own `verify.py` was also independently rerun in both modes.
Both outputs match its saved SHA256
`e7583686dd525fb65b65c5e731eedccebe3304ecd586f469f0bec14b651685ed`.
This covers its five prime-field compact/conic checks and 4,080 finite
power-map orbit cases. The recovery script and all nine saved field results
were inspected; this lane did not redundantly rerun that whole script.

## 4. Primary-source verification

- [Stacks Project, Tag 04UV](https://stacks.math.columbia.edu/tag/04UV)
  was opened: its torsor-with-equivariant-map description supports the
  quotient-stack interpretation.
- [Stacks Project, Tag 07S7](https://stacks.math.columbia.edu/tag/07S7)
  was opened: its finite-locally-free free-action hypotheses match the
  torsor use here, not the nonfree branch argument.
- [Cadman, Using Stacks to Impose Tangency Conditions on Curves](https://www.charlescadman.com/raw/pdf/stacks.pdf)
  was retrieved and inspected at Definition 2.2.4, Remark 2.2.5 and
  Examples 2.4.1–2.4.2. These directly support simultaneous roots, the
  affine root quotient, and triviality for a nonvanishing section. The
  producer replaced its former partially retrieved citation accordingly.
- [López Peña–Lorscheid, Torified varieties and their geometries over F1](https://arxiv.org/html/0903.2173v3)
  was opened: the finite-immersion decomposition convention matches the
  torification scope of the audit. The coefficient and rationality
  obstructions above are independently proved, not attributed as new
  theorems of that source.

## 5. Remaining boundaries

No unresolved mathematical blocker remains in this packet. No novelty
certification is made. Finite-field tests do not substitute for the general
geometric proofs; characteristic 2 and arbitrary nonconstant/nonabelian
twisting are outside their claims. This audit does not certify a particular
Diophantine twist's rational points merely from a positive average.
