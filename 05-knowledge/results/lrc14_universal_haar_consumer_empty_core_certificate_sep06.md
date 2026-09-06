# The universal triple bound gives a body-specific scale-three completion

**Status: PROVED COROLLARY RELATIVE TO THE UNIVERSAL LOCAL BOUND.**
The compact/open transfer is independently inherited from THM-4150 and
THM-4373. No arbitrary-body Haar floor, parity entry, synchronization, or
LRC(14) conclusion is asserted.

The actual new consumer of the universal `6/77` theorem is a complete
scale-three analogue of the dyadic Haar transfer. It upgrades the resonant
body-specific corollary of
[THM-4373 — signed-(1,2,1) triple-comb bound](../../01-canon/theorems/THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound.md)
to every distinct positive **odd** tail triple with no member divisible by
three. The exact owner criterion was already proved in the
[scale-three Haar scout](../../07-reflections/lrc14-scale-three-triple-owner-haar-root-20260903.md),
Sections 2-3; that historical note is used for its displayed proof, not for
its now-superseded frontier statement.

## 1. Exact statement and equality-safe proof

Let `C` be a finite set of positive integers, and put

```text
G_C={y in R/Z: ||cy||>=1/14 for every c in C}.
```

Suppose `G_C` is nonempty. Let `T={a,b,c}` consist of three distinct positive
odd integers prime to three. Then

```text
mu(G_C)>=6/77
  ==> G_(3C union T) is nonempty.                     (1)
```

In particular, every ten-speed body satisfying this measure condition
accepts every such triple of tails, giving thirteen nonzero speeds. Body
nonemptiness is automatic for ten speeds by cited `LRCUpTo13`; the new
theorem supplies the triple-comb ceiling, not that cited result.

For each tail and sheet define

```text
D_(w,j)={x: ||w(x+j/3)||<1/14},
F_T=union_(pi in S_3) intersection_(i=1)^3 D_(T_i,pi(i)).
```

The three danger sheets of any tail are disjoint because it is a unit modulo
three. Thus `x in F_T` holds exactly when all three translates `x+j/3` are
spoiled by the three tails. Put `P=m_3^(-1)(G_C)`. All those translates have
the same body clearance, so failure of the full completion is equivalent to

```text
P subset F_T.                                         (2)
```

The set `P` is nonempty and compact; `F_T` is open and proper (`0 notin F_T`).
Multiplication by three preserves normalized Haar measure, hence
`mu(P)=mu(G_C)`. The universal local theorem gives `mu(F_T)<=6/77` after
the primitive tail normalization described below. Equation (2) is impossible
if `mu(P)>mu(F_T)`.

It is also impossible at equality. If a nonempty compact `P` lies in a
proper open `F_T` with equal Haar measures, the open difference `F_T\P` has
measure zero. Every nonempty open subset of the circle has positive Haar
measure, so that difference is empty. Then `P=F_T` is nonempty, proper,
open and closed, contradicting connectedness of the circle. Thus (1)
includes the equality threshold. Endpoints at clearance `1/14` remain safe.

If the tails have a common divisor `d`, it is odd and prime to three.
Writing `T=dT_0`, multiplication by `d` permutes the sheet labels modulo
three and gives

```text
F_T=m_d^(-1)(F_(T_0)),  mu(F_T)=mu(F_(T_0)).
```

Hence primitive normalization loses neither the failure predicate nor its
Haar mass. No primitivity condition on `C` is needed.

There is also a stronger tail-adaptive sufficient criterion. Let
`kappa(T_0)=min_i E_i(T_0)` be the primitive triple's exact network bound.
The same proof gives

```text
mu(G_C)>=kappa(T_0)  ==> G_(3C union T) nonempty.       (3)
```

This can certify a body below the uniform `6/77` threshold without any
additional theorem about its safe-set measure.

## 2. What strict failure would now require

Within this exact typed domain, if `G_C` is nonempty but the thirteen-speed
completion fails, compact/open containment forces

```text
mu(G_C)<mu(F_T)<=kappa(T_0)<=6/77.                    (4)
```

This discharges the universal nonresonant triple-comb obligation in the
Haar scout and THM-4373, and the universal degree-zero projection ceiling
left by
[THM-4409 — third-sheet component network](../../01-canon/theorems/THM-4409-lrc14-third-sheet-component-network-certificate.md)
and
[THM-4414 — contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md).
The incoming [owner-line closure](overnight2_20260906_lrc_owner_lines.md)
remains a valid stronger *count* certificate on its subclass. Its residual
two-dimensional owner polygons no longer leave the selected `6/77` network
target open. They can still matter for sharper counts, capacities, or weights.

No assertion makes arbitrary weighted THM-4409 flows exact: degree zero has
constant length capacities, while separate nonconstant flows still lose
their edge-integral coupling. Universal degree-zero sufficiency removes
their necessity for the present unweighted ceiling.

## 3. The map and its two outstanding inputs

```text
source:    the complete safe body G_C and all three labelled tail owners;
target:    one safe phase for the literal row 3C union T;
map:       y -> its three lifts (y+j)/3, with common j retained;
preserved: weak body safety, strict tail danger, owner permutation,
           common dilation, and Haar measure;
lost:      individual body-component addresses after scalar measure comparison;
sidecar:   actual G_C and either its Haar mass or its intersection with F_T^c;
test:      mu(G_C)>=min_i E_i(T_0), or directly G_C\m_3(F_T) nonempty.
```

Two separate inputs remain. First, an arbitrary prescribed-runner row need
not have ten speeds divisible by three and exactly three odd nonmultiples.
Second, lower-dimensional LRC supplies body nonemptiness, not the quantitative
floor `mu(G_C)>=6/77`. These are independent mathematical obligations.

The prime-three boundary of
[THM-4004 — three-detuned divisor combs](../../01-canon/theorems/THM-4004-lrc14-three-detuned-divisor-comb-profile.md)
does supply a ten-speed divided pack, but its three exceptional speeds can
have mixed parity and its pack has no inherited `6/77` Haar floor. Therefore
the entire prime-three incidence boundary is not closed by (1).

Likewise
[THM-4370 — septimal wall quadrature](../../01-canon/theorems/THM-4370-lrc14-septimal-wall-quadrature-and-valuation-reanchor.md)
and
[THM-4372 — current layer cake and critical-fibre transport](../../01-canon/theorems/THM-4372-lrc14-sharp-current-layer-cake-and-critical-fibre-transport.md)
retain their exact obligations. They concern an even anchor `{2h}` and
twelve odd tails on seven-point fibres, with the anchor-safe weight `G_h`.
Their lower-shell one-sheet term can remain negative, and their deeper
current can vanish. An unweighted three-sheet bound supplies no inequality
for those weighted seven-fibre quantities without an actual common carrier.
The residue primes, number of sheets, retained body, and weight cannot be
identified by analogy.

## 4. The cheapest next consumer probes

The decisive body-side question is now a quantitative safe-set problem:
can one prove `mu(G_C)>=6/77` for the **actual ten-speed bodies produced by
an entry route**, or exhibit their first exact body below that floor?
The stronger universal ten-body floor is **REFUTED** by the
[recovered exact bodies](lrc14_haar_body_empty_core_sep06.md):
`C={1,2,3,5,7,8,9,11,12,13}` has `mu(G_C)=14249/252252<6/77`,
already in [THM-530 / G_P-intersection global witness floor, Section A](../../01-canon/theorems/THM-530-lrc-gp-intersection-global-witness-floor.md).
Even the necessary small-clock divisor sieve does not imply this floor.
The [independent consumer reconstruction](overnight3_20260906_lrc_consumers.md)
retains all its components and six isolated safe points. Actual entry-produced
body components or joint safe-phase geometry are needed. The corollary with
its explicit measure hypothesis remains valid.

The canonical synchronization hostile is already in
[THM-4032 — d=3 affine defect boundary, Section 6](../../01-canon/theorems/THM-4032-lrc14-d3-affine-defect-lattice-boundary.md):
the divided pack `{1,...,10}` is safe at `y=2/11`, yet tails `(1,5,11)`
spoil all three lifts there. The same full row is safe at `x=1/14`.
Any proposed body-to-tail map must recover a different body phase rather
than insist that an arbitrary cited witness lifts. The exact smallest next
intersection test is to reconstruct `G_C\m_3(F_T)` on this already named
row, with its two specified phases as hostile and positive controls.

For a proposed new family, first apply the inherited small-denominator
and divisor sieves. The warning is concrete:
[THM-4154 — Haar-pool inheritance correction](../../01-canon/theorems/THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction.md)
showed that three large dyadic Haar pools were already safe at the fixed
phase `1/12`. Their abstract transfer remained valid, but counting those
families did not establish a new frontier. The same novelty audit is
necessary for scale-three examples.

Finally, Haar sufficiency must not be promoted to equivalence. In the
nearby dyadic route,
[THM-4330 — anchored entry, Section 6.1](../../01-canon/theorems/THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve.md)
has `H={1,...,11}`, tails `(1,9)`, and
`mu(G_H)=10931/194040<4/63`, while `x=5/24` has clearance `1/12`.
This is an actual failure of that Haar gate on a safe row, not a counterexample
to the scale-three body-floor candidate. It motivates retaining exact
body-component addresses if the new measure gate stops.

## 5. Recovered component-width escape and the joint residual

Recover [THM-4052 — affine-component width escape cones](../../01-canon/theorems/THM-4052-lrc14-affine-component-width-escape-cones.md),
Section 3, before trying to derive another large-tail theorem. For a ten-core
`C` with `M=max C`, cited eleven-runner loneliness supplies clearance `1/11`.
The Lipschitz margin contains a closed body-safe interval of length
`3/(77M)`. On the body-phase circle let `Sigma_T=m_3(F_T)`. For three
ternary-unit tails of any parity, its distinct owner-assignment strata are
open and disjoint. Each connected component stays in one danger tooth of
the fastest tail `c=max T`, of length at most `3/(7c)`. The closed/open
comparison therefore proves escape when `c>=11M`, including equality.
This is exactly the inherited theorem, not a new height cone.

More precisely, if `L_C` is the largest closed safe-component length and
`W_T` the largest connected component of `Sigma_T`, then `L_C>=W_T`
precludes full failure. The exact phase condition remains
`G_C` not contained in `Sigma_T`; the two scalar gates are sufficient only.
For odd ternary-unit tails, combining THM-4052 and THM-4434 leaves any
hypothetical failure in the necessary region

```text
c<11M,
mu(G_C)<mu(F_T)<=6/77,
L_C<W_T<=3/(7c),
G_C subset Sigma_T.
```

The recovered clock-filtered body
`C={1,3,4,10,11,13,14,16,17,18}` with `T=(1,5,11)` defeats both scalar
gates: its measure is below `6/77`, and its exact largest component has
length `13/1568<3/77=W_T`. The first maximum interval is
`[29/224,27/196]`, whose length is `13/1568`; the complete component list
is in the [frozen body transcript](lrc14_haar_body_empty_core_sep06.out).
Yet `x=9/19` is safe with clearance `2/19`. Thus even combining full Haar
mass and largest-component width does not remove the need for component
location. This is an application of a recovered safe row, not a new LRC
family or an actual counterexample entry claim.

Reproduce the extra width comparison from the already audited closed sets:

```python
import runpy
m = runpy.run_path("04-computation/lrc14_haar_body_empty_core_sep06.py")
C = (1,3,4,10,11,13,14,16,17,18)
assert m["safe"](C) == m["wall_safe"](C)
print(max(b-a for a,b in m["safe"](C)))  # 13/1568
```

The useful next object is the pair of complete rational component lists
inside the strict relative-height cone, with the actual common lift label.

## 6. Incoming parity-free generic consumer

Incoming `eb50ee68df` promotes
[THM-4437 — all-parity reduction to three low circuits](../../01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md).
After primitive reduction, every positive distinct ternary-unit tail triple
with no signed relation of magnitude pattern `(1,1,1)`, `(1,1,2)` or
`(1,2,2)` has `mu(F_T)<=min E_i<6/77`, with no parity hypothesis.
The same actual-body consumer therefore extends to that complete generic
domain. Its three individual projection equalities do not change the strict
minimum. The earlier mixed-parity hostiles lie in the explicitly retained
low circuits and are preserved.

For the additive family, the separately
[proved sharp physical bound](lrc14_additive_parity_empty_core_sep06.md)
instead gives the sufficient gate `mu(G_C)>=6/55`, including equality by
the same compact/open argument in Section 2. These are tail-specific
sufficient thresholds. Neither proves that a particular entry supplies the
required body measure, and both leave the exact phase intersection as the
consumer when the scalar gate fails. THM-4052's width escape applies to
all these ternary-unit tails independently of parity.
