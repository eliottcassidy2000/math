# The n=12 sporadic branch after the c=35 closure, and what LRC(14) still needs

*codex-2026-07-19-S78 exact-status audit.  This separates paper theorems,
external finite certificates, Lean consumers, conjectural bridges, and direct
LRC(14) obligations.  It records THM-1249 but does not claim uniform n=12
rigidity or LRC(14).*

## 0. Verdict

The primitive twelve-speed sporadic branch is still open.  It is rigorously
finite and sharply stratified, but the finite universe has not been exhausted.
The exact global bound is

```text
A={a_1<...<a_12}, gcd(A)=1, M(A)=1/13
       ==> sum A <= 78^11 = 650190514836423555072.          (1)
```

The executed computations cover named faces far below or structurally inside
(1); they are not a run over all candidates in (1).

This audit does close one previously named finite face.  THM-1249 proves that
the primitive proper AP-centred common-scale-35 Hamming-six face is empty:

```text
3,002,076 scalar contexts
       -> 216 scalar survivors on 24 supports
       -> 1,296 owner obligations with U_5<=31<35
       -> zero covers.                                      (2)
```

The AP-centred common-scale frontier is therefore `c=36`, not `c=35`.
Equation (2) is one slice of the shallow problem; it does not extract an AP
centre from an arbitrary packet and says nothing about deep sheets.

The direct LRC(14) frontier has meanwhile moved away from “find another global
margin.”  THM-1221/1237/1244 give positioned mass, a rank-two handoff forest,
two located overlap quanta, and a private owner region.  THM-1240 forces a
directed blocker cycle.  THM-1248 is presently reserved, not proved; its target
is a finite relative-address word and an exported fourth-support wall event.
The highest-leverage direct theorem is therefore phase-bearing owner reuse or
alternate-gap descent, not further unpositioned incidence averaging.

## 1. Keep three statements separate

For a finite speed set `S`, write

```text
M(S)=max_t min_(s in S) ||st||.
```

Three claims have repeatedly been blended in historical prose.

### 1.1 Sporadic max-peel emptiness

For primitive tight

```text
A={a_1<...<a_12}, M(A)=1/13,
P=A\{a_12},       M(P)>1/12,                              (3)
```

prove that no such `A` exists.  The strict inequality in (3) is the sporadic
max-peel branch.  The complementary extremal-core branch is finite and closes
to the ordinary AP using the lower-dimensional rigidity input and exact
completion check.

### 1.2 Primitive tight uniqueness

Prove

```text
gcd(A)=1 and M(A)=1/13  ==>  A={1,...,12}.                  (4)
```

After the extremal-core dispatch, (3) is the live branch of (4).  In Lean,
`LRC14Grand.PrimitiveTightRigidity` is still a `Prop`, not an inhabited theorem.

### 1.3 The first spectral gap

Prove that no twelve-set has

```text
1/13 < M(A) < 2/25.                                       (5)
```

This is stronger stability information away from the tight locus and is not
equivalent to (3) or (4).  THM-1002 proves (5) only for `max(A)<=18`; the
AP-centred Hamming-one, Hamming-two, and Hamming-three banks are closed by
THM-1004, THM-1005, and the currently duplicate-ID file
`THM-1006-hamming-three-rigidity-of-the-AP.md`, but the global gap remains
open.  The duplicate theorem ID is catalogue debt, not a mathematical bridge.
A proof of (4) alone does not provide a uniform positive gap, and a
putative counterexample to (3) lies at the lower endpoint rather than inside
(5).

## 2. What is genuinely proved globally

| Layer | Exact proved content | What it does not prove |
|---|---|---|
| Finite height | THM-763 proves (1), conditional only on the settled lower-dimensional LRC citation; THM-668/1002 put every candidate maximum on a finite pair-sum ruler. | No feasible enumeration of all sets below `78^11`; no end-to-end Lean formalization of the zonotope descent or general ruler theorem. |
| Top completion | THM-759 gives `a_12<=a_11/(13M(P)-1)<12a_11` on (3); `LRCSporadicDiscreteCap` sharpens this to `a_12<=24a_11^2/(2a_11+13)` using quantized denominator debt. | These are compression bounds, not a contradiction.  The Lean cap takes its semantic ruler/completion inputs as hypotheses. |
| Hereditary structure | THM-765 proves every leave-one-out core of a primitive tight twelve-set is primitive and gives exact component-tooth containment. THM-766 gives the projective first-window cones. | Primitivity plus scalar capacity is compatible with every small sheet number tested. |
| Binding scale | THM-769 writes every tight maximizer as `p/(13s)`, identifies the shallow branch `s=1`, and gives the off-sheet capacity law for `s>=2`. THM-1006 identifies `s=val`, proves `val>=gcd(A)`, and reduces deep emptiness to the reverse inequality `val<=gcd(A)`. | The reverse inequality is precisely the open content law, not a proved part of THM-1006. |
| AP locus | THM-1171 proves a twelve-term AP `{a,a+d,...,a+11d}` is tight iff `a=d`. | It neither extracts an AP from an arbitrary tight packet nor supplies a tight deletion inside an LRC(14) crown. |
| All-height shallow carrier | THM-1143 proves that full-residue tightness is a labelled `A_12` mechanical ballot walk on the thirteen-sheet deficit simplex. | The arbitrary-height regeneration/pumping theorem for arithmetic-realizable wall words remains open. |

The exact finite decision formula and the finite-height theorem make (4)
decidable in principle.  “Decidable” here is a quantifier statement, not a
completed certificate.

## 3. Shallow branch: what the Hamming and ramification programme has removed

At `s=1`, tightness forces a complete nonzero residue system modulo 13.  Write

```text
w_r = r+13h_r,       r in F_13^*.                            (6)
```

The following are closed.

- THM-770: the full finite lift box `0<=h_r<=12` has only the AP/dilation
  solutions after gcd descent.
- THM-1001 and THM-795: a single coordinate cannot wind at any height.
- THM-800 and THM-804/806: the AP-centred Hamming-two and Hamming-three stars
  close uniformly in height and scale.
- THM-810/815/816: the Hamming-four common-scale and order-three coset
  alternatives close.
- THM-845: the proper scale-one Hamming-five chart closes at arbitrary lift
  height.  THM-844 and THM-847 close all order-three and mixed order-one/order-
  three contexts.
- THM-858: the remaining common-sheet Hamming-five order bank is finite,
  `{2,3,7}`-smooth, has maximum effective order at most `10,584`, and obeys
  every relative-complement lcm cut.
- THM-857/859/860: scale-one Hamming six closes, common-dilation conjugacy is
  exact, and every primitive proper AP-centred H6 ramification language has
  common scale at most 840 and a finite literal metric recursion.
- THM-861 through THM-1125 close all legal common scales through 34; primes at
  least 19 are uniformly scalar-impossible and multiples of 13 are primitive-
  impossible.
- THM-1249 now closes `c=35` by the exact `Z/5` owner bound (2).

These facts leave several different shallow obligations.

1. **Finite smooth H5 bank.**  Enumerate or structurally exclude the remaining
   `{2,3,7}`-smooth common-sheet languages through order 10,584, retaining
   complement-lcm fibres and literal endpoints.
2. **AP-centred H6 scale line.**  `c=36` is the first untreated common scale;
   composite non-13 scales continue through 840.  Closing scales one by one is
   valid but low leverage.  A prime-power flag theorem uniform in the scale
   factorization would compress hundreds of finite faces at once.
3. **No assumed centre.**  Prove that an arbitrary shallow tight packet enters
   one of the AP-centred languages.  The current Hamming charts do not provide
   that extraction.
4. **Mechanical-word rigidity.**  Equivalently, prove that every primitive
   arithmetic-realizable ballot walk in THM-1143 regenerates the Farey/AP word.
   The finite state (`1,352,078` deficit vectors) is not enough by itself because
   the labelled Beatty wall chronology is unbounded.
5. **Radius seven and beyond.**  The scalar discrepancy coefficient changes
   sign.  The live state must retain endpoint owners, third moments, correlated
   AP windows, and component--comb overlap debt.

The c=35 proof illustrates the right finite carrier: owner-coloured quotient
fibres plus an absolute union threshold.  Its owner tournaments are all
transitive and hence are only telemetry.

## 4. Deep branch: the content law is metric, not scalar

For sheet number `s>=2`, put

```text
E={v:s divides v}=sU,       F=A\E,
D_w=s/gcd(w,s) for w in F.                                  (7)
```

THM-769 proves the necessary capacity

```text
sum_(w in F) (floor(2D_w/13)+1)/D_w >= 1.                   (8)
```

It follows that `|F|>=2`.  Exactly two exceptions force `s=2`, ten even
speeds, and two odd speeds.  Three exceptions either contain a half-sheet
tightener or force the `s=3` equality packet.  THM-1006 shows that the whole
deep branch would vanish if one proved

```text
tight A  ==> val(A)<=gcd(A).                                (9)
```

For primitive `A`, (9) says `s=1`.  Capacity plus primitivity cannot prove
(9): the resulting gcd/capacity integer system has explicit solutions at
every tested sheet number.  The missing information is which exact lifts and
components are covered simultaneously.

The most developed deep packet is `s=2`:

```text
A=2U union {x,y},       x,y odd.                            (10)
```

Tightness is equivalent to both odd runners being dangerous on every
`1/13`-safe component of `U` and owning opposite nearest-integer parities at
every point.  The following are rigorous partial closures:

- THM-775 forces imprimitive deletion routes into a finite dyadic quotient
  chain with binary safe-child fibres.
- THM-776 exhausts the ten-even/two-odd locus through speed height 100.
- THM-782/789 give anchored return/erosion packets and exact global budgets;
  their liars prove that one fixed anchor or raw component order is too lossy.
- THM-797/803 add divisor-grid and parity support, but exact counterexamples
  escape on nonmaximal components.
- THM-824/831/836 reduce a principal fixed-ratio shell to exact radius and
  owner conditions.  The uniform shell-five classes `15,37,41 mod 52` remain
  open; a `U`-independent single endpoint-grid column is proved insufficient.

Thus the highest-leverage deep theorem is not another cardinality inequality.
It is a persistent-colour/component incompatibility, first for (10), that
retains actual lift set, several numerator columns or a nonendpoint column,
and owner ancestry.  A proof of (9) could be uniform in runner number, but the
shallow regeneration theorem cannot: shallow sporadics already occur at other
values of `n`.

## 5. What is formalized, what is externally certified, and what is only named

The formal boundary is clean when read declaration-by-declaration.

- `LRCTightRigidity.lean` defines `PrimitiveTightRigidity`; it does not prove
  it.
- `LRCHcompSurface.lean` defines `TightLooseDichotomy` and uses it as an input
  to conditional LRC(14) wiring.
- `LRCA12Chipwalk.lean` checks root transport and finite-state invariants, not
  arbitrary-height mechanical-word regeneration.
- `LRCNestedFibreRelaxation.lean` kernel-checks the sound anchor/remainder
  upper relaxation.  The scale-specific row enumerations remain external
  exact certificates.
- `LRCScaleThirtyFiveFibre.lean` checks the final implication
  `relaxed card<=31 ==> not all 35 sheets`; it does not replay 3,002,076 rows
  in the kernel.
- The THM-770 and scale-by-scale Hamming computations are deterministic exact
  finite certificates, often with independent replays, but they are not a
  global formal derivation of (4).
- THM-763's MSS zonotope argument and lower-dimensional LRC citation are paper/
  external inputs, not end-to-end TournamentH7 theorems.

No import at the root turns any of the open Props into a theorem.  An imported
conditional assembly is not a completed proof.

## 6. Relation to the remaining LRC(14) mathematics

The n=12 equality theorem is a supplier on one inverse route, not the only
possible proof of LRC(14).  The conditional compact-crown architecture still
needs both

```text
tight-deletion extraction / all-loose crown collapse,
primitive n=12 equality classification.                              (11)
```

AP-internal rigidity does not supply either arrow.

The newer direct seven-wall route has a different, more local residual.
Combining the latest proved results, a hard row has

- a compact six-ratio packet and one of five macroscopic cut indices
  (THM-1233/1241/1236);
- a protected slowest-spoke component with macroscopic private mass;
- an irredundant handoff chain, a rank-two forest, and two located overlap
  quanta (THM-1244);
- a phase-bearing directed blocker cycle with positive canonical-lift
  holonomy (THM-1240 plus THM-1226);
- sharp guardrails showing that one chosen gap, one quotient period, a bare
  Fano contraction, or all sampled pair blockers do not force a contradiction
  (THM-1239/1242/1243/1246/1247).

THM-1248 is still proof-in-progress.  Its proposed finite relative address and
fourth-support wall export are exactly the right next bridge, but must not be
counted as established until its referee, Lean core, and composition are
frozen.

The direct highest-leverage target is:

> **Owner-reuse/off-grid continuation compression.**  Combine the two located
> THM-1244 seams, its macroscopic private owner, the positive-holonomy blocker
> cycle, and a finite relative-address word.  Force either reuse of one label
> at incompatible oriented germs, a second independent debt edge, or descent
> to a strictly better gap/address cell.  The descent must preserve the exact
> gcd sheet and terminate at an already closed clock packet or a lonely phase.

That theorem would consume the new structure.  More global pair mass, raw
Fano incidence, a scale-free blocker-degree bound, or a low-height inverse
would ignore coordinates already proved essential.

## 7. Tournament and alternate-carrier audit

| Candidate vertices | Preserves | Destroys | Verdict |
|---|---|---|---|
| runners | order, ratios | sheet coverage, wall phase, owner chronology | usually transitive telemetry |
| residues mod 13 | shallow pinning and AP labels | lift heights and deep sheet number | useful shallow quotient only |
| sheet cuts / deficit coordinates | ballot cover predicate | arithmetic realizability if wall labels are removed | correct finite state with an unbounded labelled word |
| CRT owner obligations | absolute sheet-union threshold and ramification cuts | continuous component geometry | proof-bearing for finite H5/H6 faces, including c=35 |
| safe components | exact top-completion predicate | arithmetic future action unless comb progressions are retained | part of the faithful n=12 state |
| wall-crossing events | chronology, oriented germs, integer addresses | global mass if contracted too early | natural direct-LRC14 vertices |
| blocker labels | nontransitive functional cycle | gap address and exact tooth lift | THM-1240 is real progress but not a contradiction |
| located handoff seams | overlap debt and owner reuse | quotient-period summary | closest present proof carrier |

For c=35, owner vertices under either `Z/5` or `Z/7` proof-cost gauge give
score word `(0,1,2,3,4,5)`, no cycles, six singleton SCCs, and one Hamiltonian
path, even though the `Z/5` threshold closes every row and `Z/7` leaves 24
rows live at all six owners.  The absolute threshold, not orientation, is the
proof.

For the direct LRC(14) route, the blocker tournament genuinely has a directed
cycle.  THM-1247 shows that even a blocker-complete q=15/Fano packet can still
have an off-grid lonely time.  The cycle becomes useful only after decorating
its edges by relative determinant, canonical gap address, oriented wall germ,
owner, exact gcd sheet, and position in the THM-1244 handoff forest.

## 8. Rerouted priority list

Several attractive targets should no longer receive primary effort.  A
uniform raw `q<=25` good-period theorem is false; the `j=4` continuum ceiling
is already closed by THM-1203; finite low-height phase sampling is refuted by
the scaled `(64K,75K)` exact-relation family; one preselected crack/gap can be
erased by THM-1239; and THM-1247 shows that raw Fano incidence plus every
sampled blocker can coexist with an off-grid witness.  These are guardrails,
not unfinished proof obligations.

1. **Direct LRC(14):** finish and independently replay THM-1248, then prove
   owner-reuse/off-grid continuation or a well-founded alternate-gap descent
   with THM-1244's two located seams.
2. **Uniform shallow n=12:** prove arithmetic-realizable `A_12` ballot-word
   regeneration.  This would dominate continued AP-centred scale hopping.
3. **Uniform deep n=12:** close the `s=2` persistent opposite-colour cover, or
   prove the content inequality (9) from metric component data.
4. **Broad finite compression:** classify the remaining smooth H5 bank or
   prove a factorization-uniform H6 prime-power flag theorem.  If working one
   scale at a time, `c=36` is now the correct next face.
5. **Formalization:** connect the semantic ruler and completion theorems to
   `LRCSporadicDiscreteCap`, and retain open `PrimitiveTightRigidity`/
   `TightLooseDichotomy` as explicit hypotheses until the actual suppliers
   exist.

The audit's central correction is simple: the project has not missed one more
scalar inequality.  Both the n=12 branch and direct LRC(14) now fail only after
pairwise and quotient summaries have thrown away the same kind of datum—an
owner-labelled, oriented component chronology.  The finite c=35 closure is a
useful exact win; the structural frontier is the transport of that ownership
through wall events.
