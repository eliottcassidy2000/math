# The `n=12` sporadic branch is finite and highly stratified, not uniformly empty

*codex-2026-07-17-S64 subaudit.  Scope: the primitive tight twelve-speed
rigidity branch used by HYP-6800/HYP-6820.  This is a quantifier and
formalization audit, not a claim to close LRC(14).*

## 1. Exact predicate and verdict

For a finite positive speed set put

```text
M(S)=max_t min_(s in S) ||st||.
```

The branch audited here consists of ordered primitive sets

```text
A={a_1<...<a_12},       M(A)=1/13,
P=A\{a_12},             M(P)>1/12.
```

The last inequality defines the **sporadic max-peel branch**.  It is not the
same problem as the twelve-speed first spectral gap
`M(S) in (1/13,2/25)`: the sporadic set `A` lies at the lower endpoint, while
its nonextremal peel `P` has only eleven speeds.  Tight uniqueness and
first-gap emptiness do not imply one another.  Lean's
`HcompSurface.TightLooseDichotomy` packages both kinds of information into one
stronger open predicate, which is why the historical language can make them
look interchangeable.

The current verdict is:

- **not proved empty;** no theorem in canon or Lean closes the displayed
  quantifiers;
- **not known false;** no sporadic tight twelve-set is known, although the
  analogous thirteen-speed assertion is false by the Goddyn--Wong set
  `{1,...,11,13,24}`;
- **uniformly finite, conditional on the cited lower-dimensional LRC input;**
  THM-763 gives `sum A<=78^11`, and THM-668 makes each candidate an exact
  finite integer calculation;
- **not computationally exhausted;** the bound is about `6.5*10^20`, whereas
  every executed census is a structural face or much smaller box.

The lower-dimensional citation is now concrete rather than folklore:
Sungkawichai--Trakulthongchai, *Eleven, twelve, and thirteen lonely runners*,
arXiv:2604.23906, gives a computer-assisted proof for `k=10,11,12`.  THM-763
still remains conditional at the repository proof-object level: neither that
external computation nor the MSS zonotope descent is formalized in
TournamentH7.

## 2. What is actually closed

The canonical reductions and computations remove large, precisely quantified
faces:

1. THM-759 gives the max-peel completion inequality; THM-763 gives finite
   height; THM-765 gives hereditary primitivity and exact component-tooth
   containment; THM-766 gives the eleven projective tooth cones in the narrow
   span.
2. THM-768 removes the branch in which the maximum is the unique
   `13`-multiple.  `LRCResiduePinning.lean` proves full nonzero residue pinning
   only under the no-`13`-multiple and global tight-from-above hypotheses.
   Independently, the published Goddyn--Wong one-runner acceleration criterion
   closes the **entire** single-acceleration AP family at `n=12`, not merely a
   sampled list: replacing `r in [12]` by `m*r`, `m>=2`, would require every
   integer beginning with `13-r` in its obstruction interval to share a
   nontrivial factor with `r`, but
   `gcd(r,13-r)=gcd(r,13)=1`.  This prime-apex observation explains why the
   known acceleration mechanism cannot populate the branch; it does not
   classify arbitrary multi-coordinate or non-AP-centred sporadics.
3. THM-769 splits shallow full-residue packets from deep binding-scale sheets.
   THM-770 closes only a bounded shallow box.  THM-795 and the later Hamming
   theorems close AP-centred stars at stated radii/scales; they are not a
   completeness theorem saying that every sporadic packet is AP-centred.
4. On the AP-centred proper Hamming-six common-sheet face, the exact recursions
   and owner obstructions close scales `c=1,...,12,14,15`
   (THM-857/861/862, THM-957/958/960/962/963/969/970/974/976/977/978),
   while THM-860 excludes `c=13` by primitivity.  Scale twelve dies by complete
   owner orthogonality; scale fourteen dies at every owner-local gate; and
   scale fifteen has at most four feasible owners in any scalar row.  This
   still leaves `c>=16`, the finite ramified H5
   metric bank, radius-seven and higher correlated component languages,
   non-AP-centred/deep two-sheet packets, shell-five transport, dyadic/collar
   residuals, and higher sheets.
5. `lrc13_sporadic_exact_cap_audit_codex_S3.py` exhausts 790 completions of 77
   nonextremal cores contained in `{1,...,13}` and finds no tight completion.
   The height-twelve owner CSP and the various Hamming recursions are exact at
   their own quantifiers.  None enumerates every candidate below `78^11`.

The sharp missing bridge is therefore a **coverage/completeness theorem**:
every primitive tight twelve-set must enter one of the already closed labelled
packet languages, or a global invariant must decrease until it does.  Finite
ramification of several AP-centred faces is not such a bridge.  Nor is
“finite-decidable”: it means that a terminating exact algorithm exists after a
candidate/packet is fixed, not that the complete finite universe has been run.

## 3. Formalization audit

TournamentH7 records the open boundary honestly when one follows definitions
rather than theorem names:

- `LRC14Grand.PrimitiveTightRigidity` in `LRCTightRigidity.lean` is a `Prop`,
  not a theorem with an inhabitant.
- `HcompSurface.TightLooseDichotomy` in `LRCHcompSurface.lean` is likewise a
  named open input to an assembly theorem.
- `LRCResiduePinning.lean` formalizes pinning, not lift rigidity or deep-sheet
  exclusion.
- the Hamming/scale Lean modules formalize terminal combinatorial or arithmetic
  bricks for their named faces; no root theorem composes them into global
  sporadic emptiness.
- THM-763's zonotope finite-height proof and THM-668's general pair-sum-ruler
  theorem do not yet have end-to-end Lean counterparts.

One historical executable also blurred these levels.
`lrc13_tightness_rigidity_macmini_S107.py` inferred that twelve nonzero
residues at one `q=13` maximizing point were automatically distinct.  They are
not: completeness comes from the separate global tight-from-above pinning
theorem.  Its source/output wording is corrected in this subaudit.  Its actual
finite census and AP control remain valid.

## 4. A new uniform arithmetic brick

There is a small but genuine improvement available before any sheet split.
Let `b=max(P)` and `mu=M(P)>1/n`.  THM-668 supplies a maximizing pair-sum ruler
`q=u+v<=2b` and an integer clearance numerator `m` with `mu=m/q`.  The fraction
need not be reduced.  Strict nonextremality gives

```text
q<n m  =>  n m-q>=1
        =>  mu>=1/n+1/(nq)
        >=1/n+1/(2nb).
```

THM-759 gives, for a tight top completion `w`,

```text
w*((n+1)mu-1)<=b.
```

Combining the two inclusive inequalities yields

```text
w <= 2 n b^2/(2b+n+1).
```

At twelve speeds,

```text
w <= 24b^2/(2b+13) < 12b.
```

`LRCSporadicDiscreteCap.lean` formalizes the denominator-gap arithmetic, the
abstract completion step, their natural-number ruler composition, and the
`n=12` specialization without `sorry` or `native_decide`; the registered
module audits to `propext`, `Classical.choice`, and `Quot.sound`.  It deliberately
takes the THM-668 ruler data and THM-759 completion inequality as hypotheses;
those semantic bridges remain outside the module.  In the 77-core regression
bank, the coarse `w<12b` cap admits 10,813 primitive completions, while this
uniform denominator-quantized cap admits 6,897; using each core's exact `mu`
then leaves the already audited 790.  This is compression, not emptiness.

Endpoint audit: `q<=2b` is inclusive; `mu>1/n` is strict and becomes the
integer step `nm-q>=1`; the gap and final rational cap are inclusive.  For an
integer search the last display is replaced by its floor.

## 5. Tournament and alternate-vertex audit

The branch has been viewed through several binary relations, but none of the
lossy pair shadows is the underlying object.

| vertices | gauge / pair observable | preserves | destroys |
|---|---|---|---|
| runners | speed, lift height, or first boundary arrival | an ordering and useful priority queue | literal safe components, simultaneous cover, and owner multiplicity |
| nonzero residues mod 13 | cyclic order after a unit multiplier; tie path by residue label | shallow pinning and missing-class perturbations | lift height, zero owners, and deep binding scale |
| AP labels / lifted coordinates | Hamming distance and common scale | a specified AP-centred packet language | any proof that an arbitrary packet has such a centre |
| sheet masks or providers | who covers which CRT lift; tie path by label | common-sheet feasibility and affine owner equations | the continuous metric fibre; isomorphic masks can have different exact maxima |
| pair-sum rulers | candidate denominator and active-pair event | exact maximum evaluation when every labelled ruler and residue is retained | a tournament ranking of rulers forgets the maximizing numerator and all-runner band predicate |
| safe components and danger combs | component--comb incidence / overlap debt | the exact completion predicate `G_P subset D_w` and labelled future action | pairwise incidence alone loses higher overlap and continuation state |
| proof obligations / prime-power cuts | which relative-lcm cut is tighter | finite ramification and hereditary capacity failures | literal endpoint geometry and metric erosion |
| projective residue classes | zero intersection of owner-obligation sets after a sign-section gauge | the scale-ten `C6` of incompatibilities and triangular-prism nerve | raw runner labels, literal witness identities, and higher intersection multiplicities |
| owner-feasibility subsets / maximum-union vectors | rank owners by local feasibility and sheet-union size; tie path by owner label | the pre-nerve scale-fourteen/fifteen obstruction, including absolute missing-sheet debt | the resulting transitive tournament forgets the coverage threshold and how many owners are impossible |
| boundary events / atomic cells | first wall crossed, owner, slack, and remaining rays | the full finite recursive predicate when all labels are retained | quotienting to an unlabelled tournament loses continuation equivalence |

For any chosen pair gauge, deterministic tie-breaking supplies a Hamiltonian
path.  Existing runner, ray, and order tournaments are often transitive with a
single path even when exact cover verdicts differ.  Directed cycles can encode
signed provider constraints, but the scale-two census shows the limitation
cleanly: all 64 signed cycles are arithmetically viable, while literal metric
recursion kills 63 and identifies the remaining cover as the ordinary AP.

The most faithful current “tournament replacement” is therefore not a
tournament: it is the labelled bipartite component--comb incidence
hypergraph, decorated by rational endpoints, owner heights, sheet masks, and
remaining legal operations.  At scales nine through twelve a smaller faithful
carrier emerges: the nerve of owner-obligation sets, sometimes after
quotienting by sign to projective residue classes.  At scales fourteen and
fifteen the contradiction occurs one layer earlier, in local sheet-union
feasibility.  Pair orientations are telemetry on these objects; they become
proofs only when the edge meaning—or the absolute coverage threshold—is
retained.

## 6. Next proof obligation

The discrete cap suggests ranking a state not by raw maximum speed but by the
integer **denominator debt** `n*m-q` on each active pair-sum ruler, coupled to
component overlap debt.  The former is positive and quantized on the sporadic
branch; the latter detects wasted comb mass.  A viable recursive theorem would
show that every legal insertion either increases denominator debt, decreases
the labelled component obligation, or enters an already closed sheet packet.
Neither coordinate is monotone alone, so this is a proposed carrier, not a
claimed proof.

The immediate mechanical priorities remain: formalize the semantic
THM-668-to-ruler-data bridge, formalize THM-759's completion inequality, and
then use `LRCSporadicDiscreteCap` as their terminal arithmetic consumer.  That
would make the new cap an end-to-end Lean theorem while leaving the global
classification boundary explicit.
