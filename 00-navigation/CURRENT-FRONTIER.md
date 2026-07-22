# Current Mathematical Frontier

**Rolling state — refreshed 2026-07-22.** This file supersedes dated frontier
snapshots as a statement of present status. A linked theorem is the proof
source; this document records how the pieces compose and what remains.
Status vocabulary: **PROVED** means an in-repo proof; **CITED** a scoped external
import; **FINITE-EXACT** an exhaustive stated finite universe; **VERIFIED**
reproducible evidence rather than a general proof; **CONDITIONAL**, **OPEN**,
**REFUTED**, and **SUPERSEDED** have their literal meanings.

## LRC(14)

### Headline

**OPEN.** In the standard reduction there are 13 distinct nonzero relative
speeds, corresponding to 14 total runners. External work settles at most 13
total runners. The residual is not a missing routine finite run and not a
uniform small-period lemma.

The current proof graph is:

```text
counterexample
  -> genuine support-3..5 relation, height <= 2^20       [THM-2051]
  -> bounded sparse relation code of rank >= 11          [THM-2052]
       |-- rank 12 -> finite maximal-minor box
       `-- rank 11 -> finite two-anchor star atlas
            -> deck + sufficient gate / indexed open 26-disk failure carrier [THM-2053]
            -> signed-hull owner fan / Kelvin-Farey certificate [THM-2055/2056]
            -> primitive packets + one coprime owner interval [THM-2058]
            -> exact core/tail CRT compatibility on each clock [THM-2059]
            -> clock, pair-sum, Fejer, Euler, or rank discharge [OPEN]
```

The rank-eleven residual is therefore a finite but enormous labelled
intersection, not merely a plane or a quadratic form; rank twelve is the
separate finite-box branch and still needs exact decision. The rank-eleven
object's indispensable coordinates are the primitive parameter, integral column data, transverse
deck and bad modulus, signed-hull owner cone, pair-sum clock, and
phase/endpoint sidecar. Failure of any sufficient gate means **uncertified**,
not unsafe. The reusable certificate architecture is
`seed + selector + preserver + pointwise exit`.

### Proved reductions and terminals

- **CITED:** LRC through 13 total runners. The newest step is
  [Sungkawichai–Trakulthongchai](../05-knowledge/reference/CORE-PAPERS.md#lonely-runner-conjecture),
  a computer-assisted proof for 10–12 nonzero speeds.
- **PROVED finite-circuit alternative:**
  [THM-2051](../01-canon/theorems/THM-2051-fejer-bv-small-relation-alternative-for-lrc14.md)
  closes the coarse analytic
  “middle” sought by THM-935/946. THM-965 gives the sharp universal pair-
  covariance floor `delta_(a,b)>=-6/637`, so pairs can be paid exactly. If a
  thirteen-speed row has no exact support-three-through-five relation with all
  coefficient magnitudes at most `2^20`, whole-product Fejer--BV approximation
  gives a positive-measure strict lonely set. THM-2074 counts the exact finite
  hyperplane ledger, proves only `O(B^12)` exceptions (hence density one), and
  gives certified packets on every fixed prime-power tower. The residual is
  templates, but descent or classification within those hyperplanes—including
  the relation-rich AP boundary—remains open.
- **PROVED rank descent and explicit residual atlas:** THM-2052 combines
  THM-763's finite primitive-height bound with a pigeonhole/code argument.
  Every hypothetical counterexample has eleven independent support-at-most-three
  relations of height `91^6`, hence is rank-twelve finite-box or lies in one of
  finitely many two-anchor stars `a_i v_p+b_i v_q+c_i v_i=0`, `c_i!=0`.
  Modulo dilation each star has one projective parameter; pair clusters and
  cross-cluster triples give the same split. THM-2053 proves
  `max_i|a z_i-b u_i|<=(a^2+b^2)/91`, removing `||(a,b)||>=91L`, while adjacent
  normalized columns raise the parameter torus elementarily to margin `1/13`.
  Every row obeys `M(v)>=1/13-R/(2N)`, `N|(v_i+v_j)`; the exact deck `D_N(m)`
  removes longitudinal fibers, with a divisibility down-set and AP cutoff `156`.
  THM-2055 restricts failure to hull-owner sectors while non-hull deck sidecars
  remain essential. THM-2056 Kelvin-inverts them to `(1/91)R^{-1}K^o` and proves
  `2u.v>=A_p(u)+A_p(v)` for whole acute unimodular cones. THM-2057 closes two
  entire one-tail planes: `{a,2a,...,11a,13a,w}` by a `12a` / `14a` / `84a`
  certificate, and `{a,2a,...,12a,w}` by a `13a` / `14a` / `182a` certificate.
  More generally, a missing clock `N<=14` forces `Na|w` over `aC`. Despite
  `640690` primitive determinant failures, its arithmetic carrier has three
  leaves. THM-2059 makes arbitrary-`N` core/tail composition an exact CRT
  histogram dot product; its zero-mode split isolates disjoint-support zeros.
  THM-2060/2064 isolate the odd dyadic seam; THM-2061 folds it; THM-2066 closes
  `59,880` cores through `24`; THM-2068 proves the minimum bank; THM-2072 rules
  out a fixed uniform bank; THM-2073/2076 give a safe-child
  tower. THM-2075 preserves nonempty-core component/Euler/owner data; THM-2078
  closes terminal maximum `<=24`. MISTAKE-238 forbids descending the empty full
  row across zero-child tails; THM-2079 pairs addresses `a,2^r-1-a` and flips
  owner bits. THM-2080 gives the exact mixed `1/14`--`1/7` overlap fold. For an
  odd guard every danger comb spends at least `1/42` of its mass inside the
  guard, with equality only at speed `6h`; six distinct combs therefore cannot
  cover the `5/7` guard complement. This eliminates terminal rank six at every
  height and sharpens the tower to depth at most four. THM-2081 localizes
  Hunter to the guard complement at terminal rank seven: containment forces
  `tau_h(Q)<=2/7-sum_q measure(D_q intersect E_h)`, where `tau_h` is the
  maximum spanning-tree weight of the outside-guard pair intersections. Its
  exact three-frequency replay closes all `4,120` rank-seven pairs through
  terminal height `24`; proving the same strict inequality at all heights is
  the live depth-four target. At rank seven its scalar deficit is exactly the
  negative signed-fold sum divided by `196`, hence equals outside-cover
  multiplicity under containment. THM-2085/2087 force a height-57 complete
  relation cut: either `q/h=r/s`, `r,s<=57`, or a height-6498 guard star.
  THM-2088 makes cut rank seven finite
  (`max(h,max Q)<=91421508108581`); rank six is persistent, with THM-2089 flat
  form `q_i=u_i(z+v_i h)`. THM-2082 retains residue incidence; THM-2086 closes
  `7|h`, five `7|q`, and the lacunary cone, while THM-2094's conditional
  fourfold-moment certificate removes four `7|q`. The live modular profile is
  `7 not|h`, one to three `7|q`, in nonlacunary templates. THM-2090 splices
  the persistent cut to the
  global rank-eleven code: the full row is finite, all thirteen speeds form a
  height-`91^6` last-guard/terminal-anchor star, or `(h,Q)` is literally frozen
  and only the three earlier guards plus two original tails move on one affine
  lattice line. THM-2092 closes that apparent unbounded branch: the depth-four
  needle gives `max(S)<=128 max(Q)/3`, so every frozen block has only finitely
  many full rows. It also lifts THM-2088's bounded-terminal branch to
  `max(S)<=3900651012632789`. THM-2093 closes the remaining global star:
  its actual parameter has the exact support flag of weights `2,3,4,5,6`
  modulo `2,4,8,16,32`; the five outer private coefficients pay dyadic factors
  `16,16,8,4,2`; the six cut rows bound the common guard/anchor gcd by
  `lcm(1,...,57)^6`; and a denominator-cleared THM-2053 gate gives
  `max(S)<2912*lcm(1,...,57)^6*(91^6)^13`. Thus every no-pair rank-seven
  branch is now finite. THM-2097's mixed-threshold two-torus escape also makes
  the bounded guard-ratio branch `q/h=r/s`, `r,s<=57`, finite template-by-
  template; all boxes still need discharge. THM-2098 splits ranks 8--11 between seven transverse bands (budget `5(n-7)/49`) and seven guard-proportional bands. The six persistent height-114 cut
  rows make THM-2065 alone vacuous; location and phase sidecars are essential.
- **PROVED relative decorrelation:** [THM-2054](../01-canon/theorems/THM-2054-relative-fejer-whole-product-decorrelation.md)
  lifts bounded scalar relations; `H=2^19` clears recorded margins.
  MISTAKE-080/082 still require shape-specific torus models.
- **Scope separation:** [THM-1149](../01-canon/theorems/THM-1149-compact-essential-crown-and-farey-blocker.md)
  proves that tight deletion and an all-loose essential crown are different
  branches. Equality classification after a tight deletion cannot be applied
  before extracting that deletion.
- **Local-comb ceiling:** THM-1252 through THM-1274 close or sharply saturate
  most purely local six-comb return arguments. The live residue is global
  endpoint/child transport or a phase-located turn tax, not another unlocated
  local-return charge.
### Exact live obligation

Rank-twelve cells already lie in finite maximal-minor boxes but still need
exact finite decision. For each reduced rank-eleven two-anchor star, intersect:

```text
bad transverse deck D_N(m)<1/14
  x primitive positive parameter in the strict open-disk failure union
  x signed-hull owner/Farey ray
  x pair-sum clock and endpoint-owner word.
```

Discharge every resulting labelled cell by at least one of:

1. an exact `D_N`/clock/binding phase, generalizing THM-2057;
2. a THM-2059 CRT packet overlap or a pair-sum/HYP-2108 endpoint certificate;
3. THM-2054 on a resonance-lift/nonaliasing cell, with the missing
   model-specific plateau proved;
4. an active relation outside the rank-eleven code; or
5. owner-labelled Euler survival in THM-2047's phase-height complex.

[HYP-8871](../05-knowledge/hypotheses/HYP-8871-lrc14-owner-sector-klein-sail-automaton.md)
is the open finite-state/Farey program for this obligation.
[THM-2058](../01-canon/theorems/THM-2058-primitive-phase-packets-and-deck-fan-intervals.md)
is the proved carrier that reduces each fixed bad-denominator/owner fiber to
labelled primitive packets and one coprime interval. It makes the terminal
exact and finite; it does not prove that every surviving interval is empty.

On the dyadic branch, `r=1..4`, `|Q_r|=10..7`, and `max(Q_r)>=25`;
`r=0` is the hereditary eleven-core. At rank seven, THM-2081's tree gate and
THM-2083's bounded-support relation alternative are now the highest-leverage
target: make the bound effective or discharge every template at unbounded
scale. THM-2085/2086 are reserved only. All ranks must retain guard containment,
component addresses, endpoint owners, and both original odd-tail owner words.

### Independent routes that remain live

- **AP-core / tight-deletion supplier.**
   [THM-1017](../01-canon/theorems/THM-1017-ap-core-bridge-reduction.md) proves
   `AP core -> far element -> LRC(14)`. The extraction of that AP core from the
   compact Cover14 residual is open (the current
   [HYP-6820 audit](../05-knowledge/hypotheses/HYP-6820-q25-and-n12-uniformity-audit.md)).
   Uniform emptiness of non-AP/deep multi-defect twelve-speed branches remains
   open.
- **Euler/global phase route.** THM-2047 proves the corresponding labelled phase-height carrier:
  it retains sign, owner, side, height, and exact deletion;
  `chi(G_delta)` detects isolated tight points. THM-2050 shows complete local
  period-14 germs can agree while
  global maxima differ, so first-exit magnitude or gluing is mandatory.
- **Peel and comb routes.** Every zero-measure covering packet obeys THM-731's
  necessary peel inequality `disc_v>=6|G'_{~v}|^2`; THM-2048's integer fiber tax
  strictly sharpens it but is still a pruning gain, not a classifier. Local
  six-comb return machinery is saturated; pursue endpoint/child transport or a
  phase-located turn tax.
- **Effective spectrum.** THM-1289 imports an ineffective one-sided isolated
  gap. THM-1290 is exact through maximum speed 55 (and empties
  `(1/14,3/41)` through 64). Extend computation only with a structural filter
  or make the gap effective.
- **Exact per-row pair-sum evaluator.** THM-1002/2047 reduce a fixed row's
  maximum to rational pair-sum vertices; reflection halves the phase interval.
  HYP-8900 replays several rows exactly, including value `14/183` for the deep
  well. This is not a uniform finite family or a proof of its Wall-A restatement.
- **Frobenius wildcard.** THM-2041 preserves exact packets, but LRC still needs
  a nonzero safe seed and pointwise exit. Corrected HYP-8840's mixed moments
  preserve augmentation; neither supplies an LRC implication by itself.

### Mandatory controls and perspective prompts

- `26*{1,...,12} union {339}` has first good period `q=27`: uniform `q<=25` is
  false. `{1,...,12,5460}` defeats the old `q<=1200` scan.
- THM-2043 proves parity-Hasse completeness modulo 14, but THM-2058's GW
  family preserves that packet, blockedness through `q=13`, threshold `14`,
  every prescribed finite lift depth, and one-tail owner sector while exiting
  strictly at `47/113`. Raw jets and any fixed lift depth are not certificates.
- AP13 and the `12->26` lift have identical complete local unit-phase germs but
  maxima `1/14` and `1/12`.
- A disk/Farey bad ray is only uncertified. THM-2057's AP one-tail planes
  contain many such rays yet are wholly safe by clocks and binding phases.
- The determinant polygon is basis-dependent Euclidean data, not a Heegner
  discriminant. Paley, modular-curve, entropy, and cusp-form coincidences are
  inspiration until an explicit LRC-predicate-preserving map is proved.
- THM-2047's selected inequality slice retains sign, owner, side, height, and
  exact paired deletion; it is not an ordinary toric complement. HYP-8840's
  mixed moments retain augmentation, but no Bessel-to-sinc swap proves LRC.
- Do not default to runners as vertices. Compare columns, gaps, clock orbits,
  Farey cones, residues, wall events, endpoints, circuits, and proof
  obligations; record what each quotient forgets and which sidecar repairs it.

## NC2 and Gaussian moments

### Headline

- **PROVED:** [THM-2022](../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md)
  proves NC2 and hence unrestricted GMC(2).
- **REFUTED for higher dimensions:** explicit GMC(3) counterexamples now exist;
  GMC is false for every dimension at least 3. See the current external preprint
  entry in [`CORE-PAPERS.md`](../05-knowledge/reference/CORE-PAPERS.md).
- **FORMALIZATION PARTIAL:** normalized residue, face transport, finite-field
  contradiction, and `GMC2NC2.heightWitnessSupplier_holds` are kernel-checked
  and root-imported. The clean endpoints `nc2_of_dvdK1` and `gmc2_of_dvdK1`
  have general complex `DvdK1` as their sole explicit premise. On paper,
  [THM-2067](../01-canon/theorems/THM-2067-galois-orbit-product-closes-one-variable-dvdk.md) proves the required one-variable statement internally.

### Why THM-2022 works

For fixed support, torus descent reduces a complex null point to a number
field. THM-2067 makes the lowest balanced Wick face nonzero. At a good prime,
Kummer kills non-dilated channels, strict height kills dilated off-face terms,
and Lucas--Frobenius leaves `Q^p`: **whole-layer preservation**, not atoms.

### Boundaries and repaired bridges

- [THM-2040](../01-canon/theorems/THM-2040-the-de-factorialization-principle.md)
  is now a retired pointer to THM-2022/2033, not an independent theorem. The
  exact symmetric-wall calculation survives in the historical proof record;
  the general principle is prime-local initial form, not a global
  common-factorial or Vandermonde factorization.
- [THM-2033](../01-canon/theorems/THM-2033-the-nc2-wall-is-the-confluent-transitivity-vandermonde.md)
  is a valid determinant/Vandermonde identity for a special moment matrix, not
  for the general scalar moment.
- [THM-2023](../01-canon/theorems/THM-2023-gmc2-hyperbessel-boundary-laguerre-polya.md)
  proves the boundary hyper-Bessel function is Laguerre–Pólya; it is now a
  refinement rather than the missing NC2 step.
- MISTAKE-211/212/214/215 bar four recurring overclaims: scalar moments do not
  separate atoms; a tournament source is sufficient, not iff; Vandermonde node
  repetition is not tournament score repetition; and no global
  de-factorialization exists.
- [THM-2070](../01-canon/theorems/THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation.md)
  embeds every Laurent polynomial as a horizontal face; its aperiodic dihedral
  control has cofinite returns but infinitely many zero constant terms, killing
  the S222/S223 saddle and semigroup shortcuts.

### Live work

Package THM-2022 for publication and formalize THM-2067 to discharge general
complex `DvdK1`. Root-imported modules prove the positive-coefficient,
two-charge, and HYP-8930 fixed-support unique-channel leaves. The last mechanism is
substantive: one unique balanced composition leaves one nonzero monomial.
But HYP-8931's `LowestFaceUniqueChannel P` class is inconsistent by the empty
level-set witness `lambda=0, delta=-1, F=empty` (MISTAKE-240), so its
kernel-checked implication is vacuous and not wired into NC2. Its `98/116`
figure is a bounded census, not a theorem. THM-2070 refutes HYP-8890's saddle/Watson bypass and
HYP-8895's semigroup bypass; only their explicitly positive or two-charge
subcases survive. Continue to transfer the `seed/selector/preserver/exit`
design without pretending a preserver supplies the seed.

## Tournaments

### Reliable structural toolkit

- [THM-1805](../01-canon/theorems/THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels.md):
  Vandermonde is a signed tournament sum; transitive score permutations survive
  and directed 3-cycles are the fundamental cancellation atom.
- [THM-1862](../01-canon/theorems/THM-1862-order-join-reduction-principle.md):
  under order-join, Hamiltonian-path count `H` is multiplicative and directed
  triangle count `c3` is additive.
- [THM-1936](../01-canon/theorems/THM-1936-signed-redei-join-multiplicative.md):
  signed Rédei data is join-multiplicative.
- [THM-1880](../01-canon/theorems/THM-1880-the-a-b-functional-frame-chebyshev-pell-companions.md)
  and [THM-1885](../01-canon/theorems/THM-1885-the-half-and-plus-one-monoid-is-baumslag-solitar.md):
  the transitive skew-characteristic recurrence is Chebyshev–Pell-like under
  `x -> x+1` and `x -> x/2`, producing `BS(1,2)^+`.
- [THM-1926](../01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md):
  the Bowen–Lanford zeta factors over the strong core.
- [THM-1940](../01-canon/theorems/THM-1940-var-lambda2-is-a-4-subtournament-census-invariant.md):
  `var(lambda^2)` is a census of four-vertex subtournaments.
- [THM-1965](../01-canon/theorems/THM-1965-the-tournament-invariant-lattice-mapped.md):
  the mapped invariant lattice is exact through `n <= 6` only.
- [THM-1966](../01-canon/theorems/THM-1966-signed-redei-count-independent-invariant-n7.md):
  at `n=7`, signed Rédei magnitude is independent from spectrum plus `H`.
- [THM-2013](../01-canon/theorems/THM-2013-coordinates-for-the-continuum-cyclic-temperature.md)
  and [THM-2016](../01-canon/theorems/THM-2016-the-deep-continuum-and-the-reducibility-ceiling.md):
  cyclic temperature, the reducibility ceiling, and the invariant-resistant hot
  center are the best current decomposition. `H` is an empirical thermometer,
  not a proved coordinate of every shell.

### Live work and limits

- [THM-1950](../01-canon/theorems/THM-1950-h-ge-disc-reduced-to-strongly-connected.md)
  reduces `H >= disc` to strongly connected tournaments and verifies finite
  cases; the global inequality is open.
- Determine operation-response laws before proposing new scalar invariants.
  Work on substitution, join factorization, SCC cores, switching, duality, and
  minimal independent coordinates at `n>=7`.
- A literal Paley bridge exists only in THM-640's prime quadratic-residue scope.
  Composite modulus 14 and repeated score/node heuristics do not inherit it.
- MISTAKE-227–229 correct the newest false syntheses: the AP chain frame is
  index `11!`, Paley adjacency roots are shifted half-scaled Gauss sums, and
  THM-2053's disk union has no discriminant `-7` or Heegner classification.
- For applications outside tournament theory, first prove that the pairwise
  relation is intrinsic and target-preserving. A forced total orientation can
  destroy exactly the ties or magnitudes the original problem needs.
- Corrected HYP-8835 isolates a useful game/dynamics coordinate without
  overclaiming: tournament optimal support is uniquely odd (skew singularity
  plus the tournament block modulo two), while a positive kernel vector `Mp=0`
  yields the replicator integral `product x_i^(p_i)`. Pure optimum means
  Condorcet winner, not transitivity; an intransitive source-over-3-cycle is the
  minimal hostile control. RPS levels are circles in a simplex, not `T^2`.

## Integer sequences

### Headline

[THM-2000](../01-canon/theorems/THM-2000-support-harmonic-abel-dini-figurate-surface.md)
and [THM-2005](../01-canon/theorems/THM-2005-support-dirichlet-automatic-tournament-atlas.md)
replace “the reciprocal sum of a sequence” with two explicit objects:

```text
support profile       D_A(z) = sum_{a in support} a^{-z}
indexed multiplicity  = D_A(z) + collision tax.
```

At `z=1`, Abel–Stieltjes and logarithmic block occupancy give the exact
convergence criterion. Density zero is insufficient; the critical family
continues through every Bertrand/iterated-log scale. Egyptian splitting
conserves mass exactly at `z=1`.

### Established landmarks

- Ordinary and centered polygonal supports have digamma clocks; square and
  centered-octagonal confluences lead to trigamma forms.
- The maximum-cyclic-triangle support mass is exactly
  `75/4 - 24*log(2)`.
- The condensation-hazard support is `{5,6,...}` with profile
  `zeta(s)-H_4^(s)` and prefix-product parity shuffle tax
  `67/4 - 24*log(2)`.
- The A000568 support mass has a certified bracket narrower than
  `3.11e-44`.
- [THM-2010](../01-canon/theorems/THM-2010-new-tournament-invariant-sequences.md)
  records four-term candidates and a no-match search. That is not proof of
  sequence novelty.
- **EXACT identities / heuristic bridge:** corrected HYP-8820 places caterer
  and cake as Pascal-prefix sums, Moser and bagel as fixed binomial row
  functionals, and Fibonacci as a shallow diagonal. Algebraically,
  for `n>=1`, `bagel(n)=cake(n+1)-2`, hence
  `bagel(n)-cake(n)=T_n-1`. Klein-S313's
  full-rank gap diagonals exactly realize the g-bonacci kernels; finite-rank
  shadows eventually differ. MISTAKE-222 blocks the stronger claim that a
  shared array or matching `-1` identifies the torus and shadow boundaries or
  transfers an LRC/JC predicate. The live bridge test is an explicit common
  boundary-cell/Euler-characteristic valuation, not another prefix match.
- **EXACT arrangement shadow / blocked topology transfer:** corrected HYP-8825
  identifies the Vandermonde as the braid-arrangement defining polynomial,
  its real chambers with labelled transitive tournaments, its leading
  coalescence product, and the gap recurrence with a companion determinant.
  MISTAKE-223 blocks the jumps from these facts to general NC2
  noncancellation, hyper-Bessel factorization, or one bagel/shadow Euler
  characteristic. THM-2023 already proves the `Phi_(p,q)` zero theorem by a
  different route.

### Live work

Classify profiles under support operations rather than compare only their value
at one; study analytic continuation/abscissae and automatic/Mahler structure;
track collision taxes for census sequences; and feed support identities back
into tournament operation laws and LRC residue packets.

## Other active portfolio

- **Jacobian/Dixmier:**
  [THM-1300](../01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md)
  verifies the three-dimensional Keller map and collision; THM-1315 analyzes
  its geometry. Public provenance remains unsettled (MISTAKE-205), and the
  exact `A_3` cotangent certificate does not decide `DC(2)`.
- **Rank-two Poisson / planar and Weyl boundary:**
  [THM-2044](../01-canon/theorems/THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension.md)
  disproves the two-pair Poisson conjecture by symplectic suspension.
  [THM-2045](../01-canon/theorems/THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate.md)
  blocks planar mates for `R=x(a-b*x^r*q^s)`, while
  [THM-2046](../01-canon/theorems/THM-2046-first-order-cotangent-pullbacks-cannot-cross-the-DC2-wall.md)
  proves that this witness cannot descend through multiplication positions and
  first-order dual momenta. [HYP-8803](../05-knowledge/hypotheses/HYP-8803-A3-pullback-versus-A2-quantum-descent.md)
  localizes the remaining nonfiltered gate at extension across `x=0` in an
  exact Ore-Weyl chart. [THM-2049](../01-canon/theorems/THM-2049-the-DC2-Ore-boundary-correction-complex-is-acyclic.md)
  proves its simultaneous `(S,T)` associated-graded correction map is onto in
  every relevant beta grade, kills the grade-six obstruction guess, and gives
  a formal beta-adic lift. Finite polynomial termination and coupling to `D`
  are now decisive. THM-2063 tames maps with a pencil member affine on a source
  fiber. [THM-2071](../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
  closes every quadratic member in every output-pencil/source-foliation
  direction. [THM-2084](../01-canon/theorems/THM-2084-cubic-fiber-low-complement-gauss-manin-gate.md)
  identifies the cubic coefficient ladder with the polynomial parts of
  `(z^3+pz+q)^(n/3)` and closes reduced complementary degrees through `11`:
  a non-tame pair with a cubic pencil member must have reduced complement
  degree at least `13`. The live cubic object is the all-degree Faber
  boundary/flux/primitive triple `(E_n,Phi_n,R_n)`, beginning next at `n=13`.
  These are source-fiber degree gates, not generic cover-degree or Jelonek
  classifications.
  MISTAKE-228/229 block the old atlas and NC2/GMC chain; DC(2)/JC(2) remain open.
- **Gaussian higher dimensions:** THM-1490 is one verified higher-dimensional
  construction; newer three-real-Gaussian examples supersede any claim that
  dimension four is sharp.
- **Full portfolio:** use [`PROBLEM-LEDGER.md`](PROBLEM-LEDGER.md). Every session
  should keep at least one niche or wildcard item outside its anchor, especially
  when an operation, obstruction, or sidecar transfers across domains.

## Cross-domain connection discipline

The most reusable current bridges are not literal object identifications:

| Mechanism | Proven source | Legitimate transfer question |
|---|---|---|
| Whole-layer Frobenius / orbit norm | NC2 balanced face and TNC monodromy orbit | Exact finite-abelian packet preservation and uniform-incidence norms transfer; LRC still lacks the seed and exit, and `-1` in the packet stabilizer forces a signed orientation sidecar. |
| Quotient + sidecar | spectrum/SCC/support quotients | Which lost coordinate restores the target predicate? |
| Bulk / boundary / null | THM-2058 phase height | Strict templates occur on all large prime grids; tight support lies on finitely many level-14 pair-sum clocks; subthreshold packets vanish. |
| Operation-response | tournament joins/support unions | Which observables add, multiply, localize, or collide? |
| One-sided accumulation | LRC spectrum | Can qualitative isolation be made effective by finite state/pinning data? |
| Support Dirichlet profile | integer sequences | Do tournament/LRC census supports have functional or automatic structure? |
| Weighted-sector obstruction | THM-2045 factorized `R` family | Does the unique constant-producing sector fail on a broader Lee–Li inner edge? |

HYP-8810's “JC(2) and LRC(14) share AP-rigidity” is a **wildcard frame, not a
proved reduction**. LRC has the one-way THM-1017 supplier; the planar-JC wall
still needs an exact map and preserved predicate. THM-2045 supplies one exact
JC-side sector obstruction to compare, not an AP reduction.

A bridge must name its map, preserved predicate, loss, sidecar, and hostile control.
