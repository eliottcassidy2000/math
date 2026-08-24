---
id: THM-3878
title: "LRC(14) eleven-plus-two seams collapse by harmonic absorption"
status: >
  PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED.  In THM-3818's
  exact rank-eleven two-component 11+2 branch, THM-668 harmonic absorption
  closes every seam of scale at least three and every scale-two seam with an
  even pair coordinate.  Exact two-lift geometry closes (1,3,2).  The
  unconditional necessary ledger falls from 52,692 to 7,505.  In the extra
  relative-scale slice t>=U, cyclic slack and one auxiliary LRC(13) witness
  leave only 58 certificate survivors.  An independent safe-mass/component
  audit closes none of them.  Subsequent THM-3910 closes 41 scale-one types,
  leaving 16 plus the scale-two `(1,9)` type; THM-3976 independently closes a
  fixed 1,365-body family inside that type.  The remaining 17, the t<U slice,
  and LRC(14) are open.
source: root / THM-3818 cyclic seam and THM-668 dispatch join, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS.  A separate 93,524-gate checker rebuilds
  the atlas through Euler-totient blocks, counts seams by divisor formulas,
  constructs every obstructed odd-pair phase from the cross-centre gap
  1/(2pq), and computes the control row by exact pair-sum events.  It confirms
  45,186 absorption closures, the unique universal scale-two atlas pair
  (1,3), and the final 7,505 necessary triples.  It also proves the (3,7,2)
  control has maximum 1/10 and connected decoder graph, so it is only a
  pair-selector hostile.  Normal, optimized and frozen streams byte-match.
  Two further exact implementations compute cyclic obstruction-component
  widths and the conditional t>=U deletion ledger by wall-cell topology and
  by periodic real-line open unions.  They agree on all 5,855 pairs, the
  5,445+1,649 base closures, 353 auxiliary closures, the final ordered list of
  58 certificate survivors, and the (1,9,2) multiplier hostile.
  A fourth independent audit tests the 58 rows against scalar safe measure,
  inversion symmetry, and component multiplicity in 294 active checks.  It
  finds zero further closures, retains the corrected THM-1042 endpoint
  convention, and isolates the missing mixed moments M1,M2.  Normal,
  optimized and frozen streams byte-match.  THM-3910 and THM-3976 are the
  separately audited responses to that invoice.
depends_on:
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-668-detuned-harmonic-dispatch
  - LRCUpTo13
related:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-3976-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family
script: 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_thm3878.py
output: 05-knowledge/results/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_thm3878.out
script_sha256: 246dcb77753616aa399300daad62adaedfa838a148ea1b63edf5f75e4f4eae69
output_sha256: 2c4ad16022fddbb20bcd3a407bdd32cececc3e8e52ca6201f132a92d05767500
semantic_sha256: 58a950770c4984d4a1a3f4a4031a360042ca471ed406d529dc28aad7812c1de6
independent_audit_script: 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_independent_audit_thm3878.py
independent_audit_output: 05-knowledge/results/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_independent_audit_thm3878.out
independent_audit_script_sha256: b9760e34de5779e5ccf328b12058feb1769965c0e7d614b5fd89f779141a9143
independent_audit_output_sha256: 6dfca0670072ca39c6cba48e895cc79abcb4548b17223989518bbd7c657450a3
independent_audit_semantic_sha256: 22bffbd1a61a3e19596408c6c73bc1db6562896c3e314e327b285e35dfe73b86
cyclic_slack_script: 04-computation/lrc14_eleven_plus_two_cyclic_slack_thm3878.py
cyclic_slack_output: 05-knowledge/results/lrc14_eleven_plus_two_cyclic_slack_thm3878.out
cyclic_slack_script_sha256: 28566aa72801afef9486fcd5731d4e16e44275f26b63a37a36e6a553712fc04b
cyclic_slack_output_sha256: 64b6159df9e3ba159b0e6b3ad5558aeb6627f9016b8991b18c1f8b79b8e62806
cyclic_slack_semantic_sha256: 2466ab9a4400ecd20d3cb9c3bc01c18e01d9cc27ca99185387f9d88ab3ecd2aa
cyclic_slack_independent_script: 04-computation/lrc14_eleven_plus_two_cyclic_slack_independent_audit_thm3878.py
cyclic_slack_independent_output: 05-knowledge/results/lrc14_eleven_plus_two_cyclic_slack_independent_audit_thm3878.out
cyclic_slack_independent_script_sha256: db4fec19bc886df766d108c9f30da4085208b8c2a3f58329484da0bedb181b94
cyclic_slack_independent_output_sha256: bc423d3d6c161398f420ae5c14169788c6dff07d4ba586b2e096ae70e0ef668b
cyclic_slack_independent_semantic_sha256: 8c0f573d0f2ac275c3877fe32d4d74176b0058c9417ce7322c4428d73c1e02ad
scale1_deletion_script: 04-computation/lrc14_eleven_plus_two_scale1_deletion_thm3878.py
scale1_deletion_output: 05-knowledge/results/lrc14_eleven_plus_two_scale1_deletion_thm3878.out
scale1_deletion_script_sha256: dd1f2de9b8f77ac657be98f10c523d46e2c1af4c6b992215aa7c25bc3f42e149
scale1_deletion_output_sha256: a2c2fec59aa1f72f7247dce22e28a1d3776c57f8151d94724d1da6a66ab09502
scale1_deletion_semantic_sha256: 7401b4cbd68d3c2263e4f07f1921e1ec6e8b6e8b1b681cc0cc6315c262745a74
full_slack_independent_script: 04-computation/lrc14_eleven_plus_two_cyclic_slack_full_independent_audit_thm3878.py
full_slack_independent_output: 05-knowledge/results/lrc14_eleven_plus_two_cyclic_slack_full_independent_audit_thm3878.out
full_slack_independent_script_sha256: 408f539b3b0b71b13d0ffa31f06314fe5dcfd303a8cc3d4da1d077a2fcec878b
full_slack_independent_output_sha256: f20c8bf9af7b8376e73aade595348693b40b818efa0131b97dcb17fa48b95726
full_slack_independent_semantic_sha256: 70647b52896622718cd39ba148c4afd80bf6193c847488d5487b4fa2de4014a1
scale2_auxiliary_hostile_script: 04-computation/lrc14_eleven_plus_two_scale2_auxiliary_hostile_thm3878.py
scale2_auxiliary_hostile_output: 05-knowledge/results/lrc14_eleven_plus_two_scale2_auxiliary_hostile_thm3878.out
scale2_auxiliary_hostile_script_sha256: 27dc43dd5400de8e042f591f2972e53974791299b053eab43077536e9260a6d2
scale2_auxiliary_hostile_output_sha256: 4cecdda5661618801d2f82a6cdec6264331b094296dd6b8fb73446bd8f11e21b
safe_mass_component_audit_script: 04-computation/lrc14_eleven_plus_two_safe_mass_component_audit_thm3878.py
safe_mass_component_audit_output: 05-knowledge/results/lrc14_eleven_plus_two_safe_mass_component_audit_thm3878.out
safe_mass_component_audit_script_sha256: 2021e793e0ee734a432a9a95071a5c537d4d235250a1f57eaa91544cdf262893
safe_mass_component_audit_output_sha256: 9ce9c9fd4a0ba2179ac0aa1cef86f93c6e0e06b50d3443013729e6ae6a491b2d
safe_mass_component_audit_semantic_sha256: 843aa403298c0c884f0e2a60b552448a31bd1158185be276e2ef60a1b73478a0
isolated_endpoint_audit_script: 04-computation/lrc14_eleven_plus_two_isolated_endpoint_audit_thm3878.py
isolated_endpoint_audit_output: 05-knowledge/results/lrc14_eleven_plus_two_isolated_endpoint_audit_thm3878.out
isolated_endpoint_audit_script_sha256: f1d665ef39dc327a59d3ea3a6a911b99e4460d6325d115be30f35ecbcb51bb1d
isolated_endpoint_audit_output_sha256: 2608e180974577cfed10960da28e324bad30d9d51777d3a9914bdadb06ecf2a9
isolated_endpoint_audit_semantic_sha256: 3a4a86aac08912db3e97f8e2625d25ffd55cc7e2a1f9f8cffce38e27d1a5a66c
hash_basis: raw LF bytes
---

# THM-3878 -- harmonic absorption collapses the 11+2 seam ledger

**PROVED RELATIVE TO CITED LRCUpTo13 + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  This is a
bounded theorem inside THM-3818's exact oriented rank-eleven two-component
branch.  It proves no statement about other component shapes, physical entry,
owner/arrival, or arbitrary 14-runner rows.  In particular, **LRC(14) remains
OPEN.**

## Bounded theorem

Assume THM-3818's exact `11+2` equality branch

```text
n = s u direct-sum t(p,q),
s,t>0, gcd(s,t)=1,
p<q, gcd(p,q)=1,
```

where `u` has eleven coordinates and `(p,q)` belongs to the `5,855`-ratio
cube-class atlas.  Then a hypothetical counterexample must satisfy exactly
one of the following necessary alternatives:

```text
s=1; or
s=2, p and q odd, and (p,q)!=(1,3).
```

Equivalently, all `s>=3` seams, all even-coordinate `s=2` seams, and the
odd cell `(1,3,2)` are lonely.

### Harmonic absorption

Suppose `s|p`.  Rewrite the row as

```text
n = s( u direct-sum (tp/s) ) direct-sum {tq}.          (1)
```

The parenthesized pack contains twelve nonzero speeds.  Since `gcd(p,q)=1`,
`s|p` implies `gcd(s,q)=1`; together with `gcd(s,t)=1`, this gives
`s` not dividing `tq`.  THM-668 applies to (1) and supplies a common time of
clearance at least `1/13`, strictly above `1/14`.  The case `s|q` is
symmetric.

THM-3818 already proves that every surviving tariff seam with `s>=3` divides
one of `p,q`, so harmonic absorption eliminates all `40,982` such triples.
At `s=2`, exactly one coordinate is even unless both are odd; absorption
eliminates another `4,204` triples.

### Exact odd scale-two geometry

It remains to understand `s=2` with `p,q` odd.  Choose any `y` at which the
eleven-speed component `u` is `1/14`-safe.  Since `t` is odd, its two cyclic
lifts present the pair with phases

```text
z=t y/2,             z+1/2.                            (2)
```

Let

```text
D_w={z in R/Z: ||wz||<1/14}.
```

Both lifts fail precisely when

```text
z in (D_p union D_q) intersect
     ((D_p union D_q)-1/2).                            (3)
```

For odd `w`, `D_w` and `D_w-1/2` are disjoint: their centres are at least
`1/(2w)` apart, while the two open radii total only `1/(7w)`.  A cross term,
say `D_p intersect (D_q-1/2)`, compares centres

```text
a/p,                  (2b+1)/(2q).
```

Because `p,q` are coprime and odd, the numerator

```text
2qa-p(2b+1)
```

is always odd and can equal `1`: Bezout solves
`2qa-2pb=p+1`.  Hence the minimum circular centre gap is exactly

```text
1/(2pq).
```

The two open radii sum to `(p+q)/(14pq)`.  Therefore (3) is empty exactly
when

```text
p+q<=7.                                                (4)
```

Among the THM-3818 atlas's odd pairs, only `(1,3)` satisfies (4).  This proves
that every phase has a pair-safe lift for `(1,3,2)`, while retaining all
eleven `u`-clearances from the chosen `y`.

## Exact census

Two independent atlas builders—trial factorization with a pair loop and an
SPF sieve with a sum/totient loop—give identical objects.  Two independent
seam builders—direct scale scanning and divisor-set union—also agree.

```text
primitive cube-class ratios                         5,855
original necessary seams with s>=2                46,837
  s>=3                                             40,982
  s=2                                               5,855

closed by harmonic absorption                     45,186
  all s>=3                                         40,982
  s=2 with one even coordinate                     4,204

odd-coordinate s=2 residual                        1,651
closed by universal two-branch geometry: (1,3)         1
unresolved s=2                                     1,650
unresolved s=1                                     5,855
combined exact residual                            7,505
```

The census is **FINITE-EXACT**; the absorption implication and odd-pair
criterion are proved for all tuples in their stated symbolic scopes.

## Relative-scale cyclic-slack theorem

The residual ledger admits a much sharper filtration when the relative scale
is retained.  Put `U=max_i u_i` and

```text
A_pq=D_p union D_q,             D_w={z in R/Z: ||wz||<1/14}.
```

Let `beta_1(p,q)` be the maximum circular component length of `A_pq`.  For
odd `p,q`, let `C_pq` be the image under `w=2z` of

```text
A_pq intersect (A_pq-1/2),
```

which is half-periodic, and let `beta_2(p,q)` be its maximum component
length.  Every hypothetical counterexample in the branch of the bounded
theorem obeys

```text
s=1:                 t < 42 beta_1(p,q) U,
s=2, p,q odd:        t < 42 beta_2(p,q) U.               (5)
```

Indeed, cited `LRCUpTo13` applied to the eleven speeds `u` supplies `y_0`
with clearance at least `1/12`.  The clearance is `U`-Lipschitz, so the
closed interval of radius

```text
(1/12-1/14)/U=1/(84U)
```

about `y_0` keeps every `u_i` at clearance at least `1/14`.  Its length is
`1/(42U)`.  If `s=1` and the full row were a counterexample, its image under
`w=ty` would be a compact connected subset of the open set `A_pq`, hence
would fit strictly inside one component.  This gives the first inequality
in (5).  If `s=2`, coprimality makes `t` odd.  The two actual times `y/2`
and `(y+1)/2` present pair phases `z=ty/2` and `z+1/2`; both must be bad.
Passing to `w=2z=ty` gives the second inequality by the same argument.

Two exact constructions of these open sets agree on every atlas pair.  The
first uses rational wall cells and joins adjacent active cells only across an
active wall.  The second constructs periodic open intervals on the real
line, merges only strict overlaps, intersects the shifted copy, and then
takes the half-period quotient.  They give

```text
max beta_1 = 29/196, uniquely at (1,28);
42 beta_1<=1 for 5,445 of 5,855 scale-one pairs;

max beta_2 = 2/63, uniquely at (1,9);
42 beta_2<=1 for 1,649 of 1,650 residual scale-two pairs. (6)
```

Consequently, in the additional slice

```text
t>=U,                                                        (7)
```

(5)--(6) close `5,445` scale-one types and all residual
scale-two types except `(s,p,q)=(2,1,9)`.

### One-auxiliary deletion certificate

For any positive integer `a`, apply cited `LRCUpTo13` to the twelve-speed
pack

```text
u union {at}.
```

Under (7), its largest speed is `at`; a `1/13` witness therefore has a
closed `1/14`-safe interval of length `1/(91at)`.  Under `w=ty`, the image
has length `1/(91a)`.  A counterexample would force that image into a
component of

```text
A_pq intersect {w: ||aw||>=1/14}.                         (8)
```

For `a=p` or `a=q`, this is a genuine deletion witness: the chosen pair
runner joins the twelve-speed pack and the other pair runner must cover the
whole safe interval.  Exact component lengths in (8) close another `353` of
the `410` scale-one pairs left by (6).

For completeness, every positive `a` was tested on the final 57 pairs.  The
search is finite symbolically: if a pair-danger component has length `beta`
and `a>13/(7 beta)`, it contains a complete auxiliary-safe cell of length
`6/(7a)>1/(91a)`, so that `a` cannot satisfy the strict component cut.  The
largest necessary finite cutoff is `77`, attained at `(6,19)`, `(8,21)`, and
`(6,47)`.  The frozen companion's weaker phrase `at most 78` is still valid.
No arbitrary one-auxiliary choice
closes any of the final 57 types.

For the scale-two exception `(1,9)`, the quotient obstruction is exactly

```text
(2/21,8/63) union (55/63,19/21),
```

with both components of length `2/63`.  Thus (5) leaves only

```text
U<=t<4U/3.                                                 (9)
```

Every auxiliary multiplier `1<=a<=58` fails the component cut.  For
`a>=59`, each obstruction arc contains a complete auxiliary-safe cell of
length `6/(7a)`, proving that no larger multiplier can help either.

The exact conditional ledger is therefore

```text
unconditional THM-3878 residual                       7,505
t>=U base cyclic-slack closures                       7,094
additional scale-one deletion closures                  353
conditional one-auxiliary certificate survivors           58
  scale one                                                57
  scale two                                      (2,1,9)     1
```

The ordered survivor-list SHA-256 is

```text
34ffd609ed76d287cf7379c697cc060633345459ea321a896cd0f4c30a6255ec.
```

The word *survivor* refers to this certificate family, not to LRC.  Outside
the extra slice (7), the unconditional ledger remains `7,505`; the region
`t<U` is untouched except for the strict inequalities (5).

### Flatness sidecar saturation

THM-3743 supplies no further scale cut here.  Every atlas pair already has
the internal primitive relation `(q,-p)` of norm `p+q<=356`.  In the direct
sum decoder space, a relation supported on both components conformally
decomposes into its two component restrictions, so every Graver element is
component-supported and carries no information about `s:t`.  The first
failed implication is therefore

```text
one short internal relation  -/->  a crossing relation or a scale bound.
```

A separately forced crossing relation could still rank-raise, but the bare
THM-3743 sidecar deletes zero seams from this atlas.

### Safe-mass saturation and the routed mixed response

An independent route tests the final 58 conditional rows against the exact
safe mass of the eleven-pack, inversion symmetry, and component counts.  It
closes no row.  The scale-one pair-obstruction measures are all at least
`5/21`, while the exact AP11 safe mass is only `10931/194040`; the scale-two
exception has obstruction measure `4/63`.  Thus even the best possible
universal scalar safe-mass floor cannot force noncontainment by a Frechet
comparison.  Both sides are inversion-symmetric, and raw component counts
carry no relative phase.

The exact missing quantities are visible without approximation.  If `G_u`
is the eleven-pack safe set and `D_p,D_q` are the two labelled danger combs,
then

```text
mu(G_u \ (D_p union D_q)) = M_0-M_1+M_2,
M_0=mu(G_u),
M_1=mu(G_u intersect D_p)+mu(G_u intersect D_q),
M_2=mu(G_u intersect D_p intersect D_q).                (10)
```

Current pair geometry supplies unconditional danger moments, not the
`G_u`-conditioned incidences `M_1,M_2`.  The needed sidecar is therefore an
owner/phase-labelled mixed-incidence theorem, equivalently suitable Fourier
correlations with both pullback combs.

THM-3910 supplies fixed-radius auxiliary-center, body-component, and integer
t-sheet responses.  Its first response closes 41 of the 57 scale-one types,
leaving 16 plus the scale-two type.  THM-3976 gives the exact signed-endpoint
cross-phase form, re-proves the known 57-row AP11 scope, and closes the fixed
1,365-body scale-two family
`2E union {t,9t}` for `E subset {1,...,15}`, `|E|=11`, odd `t>=max E`.
Its Graver-fibre and half-translate hostiles prove that relation length and
full autocorrelation still forget the needed relative phase.  Crucially,
THM-3976 closes no arbitrary-body type; the general conditional ledger after
THM-3910 has 17 types.

The failed closure route nevertheless gives a rigorous packet corollary.
For every eleven distinct positive speeds `u=(u_1,...,u_11)`, cited
`LRCUpTo13` supplies a `1/12`-deep time `t_0`.  Partition their joint phases
into `84^11` cells, take a preimage `A` of measure at least `84^-11`, and put
`D=A-A`.  Then

```text
D=-D,  0 in D,  mu(D)>=2*84^(-11),
||u_i d||<1/84 for every i and d in D,
t_0+D subset {t:min_i ||u_i t||>1/14}.                  (11)
```

Same-cell subtraction proves the return bound.  Since `A` lies inside one
coordinate-cell preimage of measure `1/84<1/2`, the circle
Kneser--Macbeath inequality gives the doubled mass.  The strict safe set has
at most `sum_i u_i` positive arcs, so one has length at least

```text
2*84^(-11)/(sum_i u_i) >= 2*84^(-11)/(11 max_i u_i).    (12)
```

This strengthens the scalar floor but does not locate the translate relative
to `D_p,D_q`, so it closes none of the 58 rows and makes no LRC(14) claim.

Isolated closed-safe walls form a separate endpoint channel.  Being finite,
they change none of `M_0,M_1,M_2`, their Fourier coefficients, Frechet bounds,
or positive widths; they cannot be charged as positive component slots.  For
the AP11 control the four walls are `3/14,5/14,9/14,11/14`.  Across all 57
scale-one rows they are simultaneously pair-safe exactly when

```text
14 does not divide tp and 14 does not divide tq.          (13)
```

The exact helpful-residue-count histogram is
`0:2, 6:7, 7:4, 12:33, 13:11`; `(3,14)` and `(5,42)` never benefit.
Moreover the legal hostile `U=11,t=14` maps all four walls to pair phase zero
for every row, so there is no universal scale-one endpoint closure.

The failure is not repaired by adding pair sum or positive multiplicity.  At
the same legal scale `t=15`, pairs `(3,14)` and `(7,10)` both have sum `17`,
danger measure `13/49`, inversion symmetry, 14 positive base components and
210 positive pullback components.  Yet the first kills all four AP11 isolated
walls and the second preserves all four.  Atomic endpoint incidence is not a
function of those scalar statistics.

For the scale-two exception `(2,1,9)`, every allowed odd `t` gives a safe
coherent lift at each of the four AP11 walls (56 residue/wall cases: 32 with
one safe lift and 24 with both).  This is a special-body positive control, not
an arbitrary-pack theorem.  A general isolated wall has address
`(14k+epsilon)/(14u_j)`; its owner `u_j`, sign, numerator, and the residue of
`t` modulo its denominator are the missing atomic sidecar.

## Minimal hostile and controls

The first atlas pair for which the universal two-branch selector fails is

```text
(p,q,s)=(3,7,2),             z=13/20.
```

Indeed,

```text
min(||3z||,||7z||)=1/20,
min(||3(z+1/2)||,||7(z+1/2)||)=1/20,
```

so a different pair runner kills each lift.  This is a hostile to the
**universal pair-only selector**, not a counterexample to LRC(14).

The phase-uniform selector obstruction can be realized at an actual
eleven-speed good point.  This is a method control, not an example in the
rank-eleven two-component equality branch: its full THM-3818 decoder graph
is connected, with 33 decoder edges.  Take

```text
u=(1,2,3,4,5,6,7,8,9,11,12),       y=3/10,
n=2u direct-sum (3,7).
```

Every coordinate of `u` is safe at `y`; the selected lifts `3/20,13/20`
both have row clearance `1/20`.  But the same full row is safe at `1/20`
with clearance `1/10`.  This explicitly prevents the selector hostile from
being misread as an LRC hostile or as a realization of an unresolved
two-component seam.

The off-atlas pair `(1,5)` is a positive geometry control: it also satisfies
`p+q<=7` and is universally two-branch safe, but it is absent from the cube
atlas because its sum has a factor `3`.  This keeps geometric classification
separate from atlas membership.

## Typed connection ledger

```text
source:      THM-3818's oriented 11+2 quotient and cyclic tariff
target:      THM-668's scaled twelve-pack plus one detuned runner
map:         s|p sends tp into the pack as s*(tp/s), symmetrically for q
preserved:   all eleven u-clearances and the absorbed pair-runner clearance
destroyed:   the chosen LRC(<=13) time, winning branch, owner, and arrival
sidecar:     gcd(s,t)=gcd(p,q)=1, ensuring the other runner is detuned and
             its branch orbit is complete
decisive test: exact 46,837-seam census plus rational scale-two wall cells
```

There is a real but limited connection to the factorial local-rule audit.  In
both cases, a local divisor condition is a **selector for the next proof
mechanism**, not a stable global barcode.  At factorial `d=9996`, divisor
places leave degree `3998` and an adaptive nondivisor `p=19` sidecar is
needed.  Here `s|p` does something stronger than filtering: it changes the
decomposition itself by absorbing one runner into a coherent harmonic pack,
after which THM-668 is proof-forcing.  There is no claimed arithmetic map
between factorial primes and LRC scales; the shared content is only this
typed selector/sidecar architecture.

## Remaining frontier

Unconditionally, the bounded result leaves exactly:

- all `5,855` scale-one triples, where no cyclic branch orbit exists; and
- `1,650` odd-coordinate scale-two triples beginning with `(3,7,2)`.

Inside `t>=U`, this theorem leaves the 58 certificate survivors above.  The
exceptional `(2,1,9)` row survives every scalar auxiliary multiplier, while
the 57 scale-one rows survive every such single-auxiliary component cut.
THM-3910 later closes 41 scale-one types.  THM-3976 closes its fixed
1,365-body family inside the exceptional type, but neither that family nor
the AP11 control deletes another arbitrary-body type.  Thus 17 remain.

The hostile shows that the latter require information about which points of
`G(u)` are available, or an owner/arrival sidecar; pair geometry alone cannot
choose a safe branch.  Component shapes `3+10` through `6+7`, the rank-twelve
terminal, physical entry, and LRC(14) are untouched.

## Exact replay

Run from the repository root:

```bash
python3 -B 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_independent_audit_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_independent_audit_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_cyclic_slack_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_cyclic_slack_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_cyclic_slack_independent_audit_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_cyclic_slack_independent_audit_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_scale1_deletion_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_scale1_deletion_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_cyclic_slack_full_independent_audit_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_cyclic_slack_full_independent_audit_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_scale2_auxiliary_hostile_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_scale2_auxiliary_hostile_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_safe_mass_component_audit_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_safe_mass_component_audit_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_isolated_endpoint_audit_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_isolated_endpoint_audit_thm3878.py
```

The 96,102-gate primary companion independently builds the atlas and seam
union twice and checks all 1,651 odd pairs by exact rational wall cells.  The
93,524-gate audit uses different atlas, divisor-count, modular-inverse and
pair-sum-event routes.  Normal and optimized executions byte-match their
raw-LF frozen outputs.  The finite census controls the ledger; the absorption
and odd-pair arguments above prove their stated symbolic scopes.  The two
cyclic-slack implementations and their two deletion ledgers are
structurally independent in interval representation; all normal, optimized,
and frozen streams also byte-match.  Their finite counts support only the
explicit conditional implications above.  The safe-mass audit is a stopping
  boundary plus the proved packet corollary `(11)`--`(12)`.  THM-3910 and
  THM-3976 are subsequent response theorems; only THM-3910 changes the general
  certificate-survivor count, from 58 to 17.  **QED.**
