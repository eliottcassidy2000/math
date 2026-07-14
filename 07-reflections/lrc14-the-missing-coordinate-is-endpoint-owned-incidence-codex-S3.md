# LRC(14): the missing coordinate is endpoint-owned incidence

**Session:** codex-2026-07-14-S3
**Status:** frontier synthesis after THM-761--766; LRC(14) and the `n=12`
sporadic tight branch both remain open at their stated residuals.

The broad atlas already identified the live object as a scale-residue packet
rather than a box of speeds.  The present audit sharpens the fiber carried by
that packet.  Two proposed scalar terminals failed for the same reason:

- `q<=25` remembered a denominator but forgot multiplier and blocker ownership;
- safe measure and runner tournaments remembered mass/order but forgot which
  endpoint component must fit which danger tooth simultaneously.

The missing coordinate is **endpoint-owned incidence**.  At small rational
periods it is a zero-owner/signed-pair blocker deck.  At tight equality it is a
component-tooth/splice hypergraph with exact endpoint phases.  Both sit over a
scale-normal residue packet, and both become false when projected to one score
per runner or modulus.

## 1. What the historical viewpoints contributed

The repository's viewpoints are not discarded analogies; each exposes one
coordinate of the same total object.

| viewpoint | durable information | unsafe forgetting |
|---|---|---|
| torus line / central box | exact LRC predicate and closed boundary | arithmetic reason a wall is present |
| danger arcs / interval complex | component endpoints, owners, open-vs-closed equality | scale-normal recurrence |
| divisor covering / CRT | zero phases and rational witnesses | nonzero signed blockers |
| Farey / Ostrowski / pair-sum rulers | all exact maximizer times | simultaneous obstruction across rulers |
| Bonferroni / coverage depth | exact finite small-body certificates | signed cancellation and coherent scale |
| Fourier / relation lattice | global resonance | location if absolute values are taken |
| autocorrelation / Bernoulli `B2` | signed endpoint interaction for one peel | global classification of cores |
| safe peel | exact recursive descent when a runner is nonbinding | confluence and irreducible normal form |
| pack / cluster / affine suspension | coherent clocks and scale fibers | exceptional residues if over-compressed |
| zonotope / geometry of numbers | genuine global finite height | equality classification inside the bound |
| tournaments / observer switches | ordering, cycles, gauge sensitivity | metric clearance and joint hyperedges |
| Cech / sheaf / normal fan | compatibility and controlled-forgetting tests | the numerical predicate without labels |
| endpoint splice graph | exact equality defect and blocker ownership | none at fixed threshold; size grows |

The recurring historical error was not using the wrong viewpoint.  It was
treating a useful projection as injective.

## 2. Small periods: from denominator ceiling to blocker deck

THM-762/764 prove that for a covering family and `15<=q<=28`, a rational
`q`-witness exists exactly when there is no `q`-divisible speed and the signed
unit-pair deck misses a class.  The counterexamples

```
26*{1,...,12} union {339},
{81,91,131,151,157,196,258,274,313,328,330,339,348}
```

first witness at denominators `27` and `26`.  The latter is loose and
gcd-incoherent.  Exact S105 replay finds `91/8260` least denominators above 25
and a bank maximum of 38.  Hence `q<=25` is false for both coherent and
incoherent reasons.

The correct finite carrier is not a denominator tournament.  Its vertices are
the obligations `(q,a)` or, after the proved `q<=28` compression, signed unit
pairs plus a zero-owner bit.  A modulus tournament may rank how compressed its
deck is, but it cannot decide whether a card is missing.  In the exact audit,
changing the gauge flips `47/55` modulus edges while the blocker verdict is
unchanged.

This redirects the LRC(14) closure problem.  A rational-period route must find
an **adaptive** incomplete deck, or prove that simultaneous deck completeness
forces a named coherent/cluster structure.  THM-566 forbids a universal raw
finite ladder; THM-761 closes many small-exception sheet packets by a different
adaptive mechanism.

## 3. Tight twelve-speed equality: the new reduction stack

For a primitive tight twelve-set `A={a_1<...<a_12}` and its max-peel
`P=A\{a_12}`, the sporadic branch is `M(P)>1/12`.  Four uniform facts now
compose:

1. **Finite height (THM-763).** `sum A<=78^11`.  The branch is genuinely
   finite, but the bound is roughly `6.5e20`, so finite does not mean checked.
2. **Hereditary primitivity (THM-765).** Every leave-one-out core has gcd one.
   Translation sheets eliminate every imprimitive peeled core at all heights.
3. **Ratio and first-window geometry (THM-759/766).** With
   `m=a_1`, `b=a_11`, `w=a_12`, every tight set has
   `w>=12m` and `b/m>=72/7`.  If `b<12m`, it lies in one of eleven exact
   projective tooth cones and satisfies a rational center-alignment band.
4. **Residue pinning.** A tight pair-sum maximum has denominator divisible by
   13.  Inside the separately defined full-residue branch, the twelve speeds
   are labelled by all nonzero classes modulo 13; tightness alone does not
   force that branch, so the endpoint argument must either reach it by a
   proved reduction or treat the other residue multiplicities as well.

The exact max-peel predicate is

```
G_P={t:min_(p in P)||pt||>1/13} subseteq
D_w={t:||wt||<=1/13}.
```

This is simultaneous component-tooth containment.  Measure cannot replace it:
for `P={1,...,10,12}`, `M(P)=1/11` but
`|G_P|=461/8190<2/13=|D_w|`.

## 4. Endpoint Euler defect

For target `1/q`, split each open danger comb into individual teeth

```
I(w,m)={t in R/Z: |wt-m|<1/q},   m in Z/wZ,
```

where the absolute-value expression is read on the compatible real lift
(equivalently `|t-m/w|<1/(qw)`).  The index `m` is part of the tooth data;
the unindexed condition `||wt||<1/q` denotes the whole comb.

Assume the open danger union is not the whole circle, and let `kappa_q(W)` be
the number of connected components of its tooth-intersection graph
(equivalently, of the union).  A **protected splice** is a common end/start
endpoint at which no third open tooth is active; let `P_q(W)` count these
splices.  An exact endpoint sweep gives

```
chi_q(W)=kappa_q(W)-P_q(W)
        = number of open q-safe components.
```

Thus `chi_q=0` is the exact no-strict-witness condition at that threshold.
It remembers both overlap rank and boundary closure.

The non-full-union condition is automatic for twelve speeds at `q=13` by the
settled LRC lower bound.  In the excluded full-circle case, the sweep has no
component start and reports zero while the union has one topological
component; this is why the hypothesis is stated explicitly.

For two speeds `u,v`, the strict tooth-overlap condition is the integer
determinant test

```
q*|m v-n u|_(mod uv) < u+v,
```

and the number of pair intersections is

```
e_q(u,v)=gcd(u,v)
          * (1+2 floor((u+v-1)/(q gcd(u,v)))).
```

The stored audit checks this formula on 360 pairs with zero failures.

For a complete nonzero residue packet modulo 13, label the speed in residue
`r` by `w_r`.  Nominal equality splices come from complementary residue pairs,
with raw capacity

```
P* = 2 sum_(r=1)^6 gcd(w_r,w_(13-r)).
```

The defect decomposes exactly as

```
chi_13 = (kappa_13-P*) + (P*-P_13).
```

The first term is overlap-rank shortage; the second is third-runner blocker
debt.  This decomposition joins the old Cech nerve, pair-sum obligations, gcd
decks, and endpoint-owner work in one integer invariant.

## 5. Exact evidence and the boundary of THM-766

The height-one complete-residue cube `w_r=r+13k_r`, `k_r in {0,1}`, has 4095
nonzero lift patterns.  Exact enumeration gives:

```
4085  kappa>P*                         (rank shortage),
   9  capacity survives but blockers leave chi>0,
   1  chi=0: 2*{1,...,12}, gcd 2.
```

All 4094 primitive rows have positive defect; the minimum is 2.  The named
defects are AP `0`, doubled AP `0`, one-coordinate lifts `2,2,4,12`, and the
all-odd lift `17`.

This also identifies where the first-window theorem stops.  Some failed lifts
and the Goddyn--Wong tight packet lie on the boundary `b=nm`, outside the open
first-window split.  Endpoint defect still distinguishes the failed `n=12`
lifts there.  A uniform proof therefore needs both the projective cone and the
endpoint splice invariant.

The bounded historical 77-core bank is consistent: the exact THM-759 cap
leaves 790 primitive completions; THM-766 eliminates all 40 narrow-core rows,
and exact `M` evaluation eliminates the 750 wide-core rows.  Its minimum is
`1/12`, but this remains a bounded bank.

## 6. Tournament Analysis: useful liar, exact sidecar

For the endpoint audit, use the six complementary-pair splice obligations as
vertices.  Pairwise blocker counts define an orientation; reversing the count
comparison is the gauge, and `1->2->...->6` is the tie Hamiltonian path.  The
script reports scores, directed triangles, SCCs, and edge flips.

The guardrail is decisive: the failed lifts replacing `2` by `15` and `7` by
`98` have the same transitive tournament fingerprint as the AP, but endpoint
defects 2 and 12.  The tournament preserves blocker direction and destroys
unspliced tooth components.  It is telemetry; the endpoint graph plus
protected-splice sidecar is theorem-facing.

The challenged default vertex set is therefore runners.  Better vertices are,
depending on the stage:

- signed residue obligations for small-period closure;
- deck lifts for gcd-sheet obstruction;
- individual teeth and endpoint splices for equality rigidity;
- proof routes only for scheduling, never for deciding clearance.

## 7. The precise next lemmas

### Post-pull bridge: metagraph address versus endpoint stalk

The concurrent HYP-6825 atlas now gives an exact finite fact that clarifies,
rather than replaces, this endpoint picture.  Rooted weighted blue/black
line-WL separates all 272 converse-merged tournament nodes at `n=7`, while
raw line incidence gives only 159 cells.  This is a complete address for the
known finite **base**.  It is not a clearance invariant: the arithmetic
progression and positive-defect lifts above can share the same tournament
fingerprint.  In the constructible-atlas language, `chi_q`, protected-splice
owners, rational endpoint positions, and the closed-threshold flag form a
concrete metric/owner **stalk** over the metagraph address.

This also makes two live pull cards precise.  MPA-19 should use endogenous
pair-sum or splice obligations as vertices and retain their owners; MPA-20
should test whether `(metagraph address, divisor mask, cap ratio,
endpoint-splice defect, peak witness)` is predicate-pure on named LRC banks.
At THM-761's seven-exception wall, the analogous base is a cyclic tiling of
the sheet fibre `Z_c`, not the staircase tiling.  A promising transfer is to
define sheet components and protected owner-splices there and classify its
zero-defect labelled tilings.  Any map between the two tiling systems must
preserve free-sheet nonemptiness and owner transport; shared vocabulary is
not yet a functor.

The incoming four-far cone audit supplies the outer chart.  Exactly four far
speeds give a rank-four semilinear cone over each nine-speed small core, while
larger far counts have higher-dimensional cones.  The endpoint stalk varies
over those charts even at fixed residue address.  Thus the proposed object is
now layered without a dimension sleight of hand:

```text
far-cone / affine-slope chart
    -> metagraph or cyclic-tiling base address
    -> endpoint-owner metric stalk
    -> wall, sheet, peel, and scale transport
    -> witness or named residual.
```

### Equality-rigidity lemma (n=12)

For every primitive full-residue packet
`W={w_r: w_r=r (mod 13)}` satisfying THM-763's height bound and THM-766's
cone/alignment constraints,

> `chi_13(W)=0` implies `W=c*{1,...,12}`; hence primitivity gives
> `W={1,...,12}`.

An equivalent quantitative target is: for every nondilate, either
`kappa>P*`, or at least `P*-kappa+1` nominal complementary-pair splices are
blocked by third runners.  This is the **splice-lattice coherence lemma**.
The `q=5` beater and Goddyn--Wong at `q=14` show that it must be stated for the
full-residue `q=13`, twelve-speed setting, not as a universal tight-set slogan.

### LRC(14) closure lemma

For the scale-normal `f>=4` covering residual left after THM-760/761, prove a
well-founded alternative:

1. an adaptive modulus has an incomplete signed blocker deck;
2. an endogenous pair-sum obligation is unblocked;
3. a safe peel descends;
4. the packet has a coherent/cluster sheet certificate; or
5. a normalized incidence complexity strictly decreases.

The measure of complexity must include exceptional residue words and endpoint
ownership.  Raw height, far count, tournament class, and denominator ceiling
have all been proved insufficient.

## 8. Research verdict

The repository no longer lacks viewpoints.  It lacks an inverse theorem that
reconstructs arithmetic coherence from **complete endpoint-owned obstruction**.
That is the common frontier behind the open scale-normal LRC(14) branch and the
open `n=12` equality-rigidity branch.  The right object is a recursive,
proof-carrying incidence sheaf over scale-residue packets—not a larger speed
box and not a more elaborate unlabelled tournament.

*Artifacts:* `lrc14_q25_uniformity_refutation_codex_S3.py`,
`lrc14_s105_q25_bank_audit_codex_S3.py`,
`lrc13_sporadic_exact_cap_audit_codex_S3.py`, and
`lrc13_endpoint_splice_defect_codex_S3.py`, with stored exact outputs.
