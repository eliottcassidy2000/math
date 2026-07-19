# The compact bridge lands in an essential crown, not automatically in the shallow ballot cone

*codex-2026-07-18. Scope: exact structural advance and exact counterexamples;
LRC(14) and the `n=12` tight classification remain open. This note integrates
THM-1125's essential-region identity, the corrected `INVcov` scope of
THM-1099, and the new shallow `A_12` chip-walk/Farey-comb carrier.*

## Executive result

Put `lambda = 1/13` and

```text
B_v = {t in R/Z : ||v t|| < lambda}.
```

If a thirteen-speed family `V` has `M(V)<lambda`, deleting a speed does
**not** automatically produce a tight twelve-speed family. There is an exact
dichotomy:

1. some deletion is tight at `1/13`; or
2. all thirteen deletions are loose, and `V` carries thirteen disjoint,
   positive-width **private essential regions**, one inside a tooth of each
   deleted speed. Call this the **all-loose essential crown**.

The second branch is real. The exact row

```text
V0 = {1,2,3,5,7,8,9,10,11,12,17,19,104}
```

is primitive and compact, has `M=8/105<1/13`, and every one of its thirteen
deletion cores has `M>1/13` (the minimum is `2/25`). It misses only the
Cover14 side of the intended inverse-theorem hypotheses: it covers `2,...,13`
but not `14`.

The already-known lift `7 -> 112=7+105` is therefore exactly the right
stress test. It preserves the complete residue/gap germ at `8/105`, supplies
Cover14, and changes the global value to `3/20`. Thus a one-chart Farey word,
an essential-region picture, or a tournament shadow cannot by itself use
Cover14. The missing bridge must glue the crown to the integer-lift atlas.

## 1. Essential-crown theorem

For `v in V`, write

```text
W_v     = V \ {v},
E_v     = {t : ||w t|| >= 1/13 for every w in W_v},
delta_v = M(W_v)-1/13,
m_v     = max(W_v).
```

The set `E_v` is THM-1125's essential region, now at threshold `1/13` and
with the strict/open convention needed by `M(V)<1/13`.

> **Theorem A (essential crown).** Let `V` be thirteen distinct positive
> integer speeds with `M(V)<1/13`. Then:
>
> 1. every `E_v` is nonempty;
> 2. `E_v subset B_v`, and the thirteen `E_v` are pairwise disjoint;
> 3. if `delta_v>0`, then `E_v` contains a circular interval of length at
>    least `2 delta_v/m_v`, while every connected component of `E_v` has
>    length strictly less than `2/(13v)`. Consequently
>
>    ```text
>    0 < delta_v < m_v/(13v);                              (CROWN)
>    ```
>
> 4. if `N(t)=#{v in V:t in B_v}` and
>    `mu_k=measure{t:N(t)=k}`, then
>
>    ```text
>    mu_1 = sum_(k>=3) (k-2) mu_k.                         (BALANCE)
>    ```
>
>    Up to endpoints, `{N=1}` is the disjoint union of the private regions
>    `E_v`.

**Proof.** LRC(13), already known, gives `M(W_v)>=1/13`, so `E_v` is
nonempty. Since the full family has `M(V)<1/13`, at every point of `E_v` the
deleted speed is strictly dangerous; hence `E_v subset B_v`. If a point lay
in `E_v` and `E_w` for two different owners, each of the thirteen speeds
would be safe there, contradicting `M(V)<1/13`.

If `delta_v>0`, take a maximizer `t_v` of the core. The function
`min_(w in W_v)||wt||` is `m_v`-Lipschitz, so the circular interval of radius
`delta_v/m_v` around `t_v` lies in `E_v`. A connected subset of `B_v` lies
inside one of its separated teeth, each of length `2/(13v)`. Because `E_v`
is closed and that tooth is open, the upper inequality is strict. This proves
`(CROWN)`.

Finally every `B_v` has measure `2/13`, so `integral N=13(2/13)=2`.
Strict coverage gives `N>=1`, hence `sum mu_k=1` and `sum k mu_k=2`.
Subtracting twice the first identity from the second yields `(BALANCE)`. A
point has multiplicity one with owner `v` exactly when it is in `E_v`. ∎

This is a stronger simultaneous version of the stability inequality appearing
around THM-1028/1039: a hypothetical strict compact cover must make **all
thirteen** deletion cores near-tight relative to their deleted teeth, not just
the maximum-deletion core.

The balance law also identifies why pairwise/tournament telemetry repeatedly
stalls. Positive private mass is paid for by co-located multiplicity at least
three. The obstruction is intrinsically third-order; pair overlap totals do
not remember where the overlaps coincide.

There is no quantitative contradiction here from compactness alone. On `V0`,
the exact private mass is

```text
316140497/756575820 = 0.4178569928...,
```

whereas the sum of the thirteen Lipschitz component floors is only

```text
352684666109/30512702820600 = 0.0115586177....
```

The mass identity has a factor `36.1511...` of room over those guaranteed
components. Since `rho(V0)=104/19<13`, the crown inequalities plus the mass
law plus compactness cannot force a tight deletion. Cover14 (or an equally
strong cross-modulus input) must do actual work.

## 2. Farey-toothpick regeneration blocker

The dilation ray has more structure than merely `M(d{1,...,12})=1/13`.
Its explicit Farey-maximizer set contains

```text
T_d = { n/(13d) : 0 <= n < 13d and 13 does not divide n }.
```

There are `12d` points: `d` affine copies of the twelve base teeth. At each
one, multiplication by `n mod 13` permutes the twelve nonzero residues, so
the core minimum is exactly `1/13`.

> **Theorem B (regenerated AP cores cannot occur in a primitive compact
> `INVcov` counterexample).** If
>
> ```text
> M(d{1,...,12} union {v}) < 1/13,
> ```
>
> then `13d` divides `v`. Consequently, if this thirteen-set is primitive,
> covers every modulus `2,...,14`, and has `rho=max/second<13`, no such strict
> inequality is possible.

**Proof.** Put `N=13d`, `g=gcd(v,N)`, and `Q=N/g`. Strict inequality would
force `v` to be dangerous at every point of `T_d`. As `n` ranges modulo `N`,
the values `vn/N` run through each point of the `Q`-grid exactly `g` times.
The number of strict `1/13`-dangerous points on that grid is

```text
b(Q)=1+2 floor((Q-1)/13).
```

If `Q>=2`, then

```text
13 b(Q) <= 2Q+11 < 12Q,
```

so fewer than `12N/13=12d` numerator classes can be dangerous. They cannot
contain all of `T_d`. Therefore `Q=1`, or `13d|v`.

Now primitivity forces `d=1`. The core `{1,...,12}` carries no multiple of
`14`, so Cover14 forces `14|v`; together with `13|v`, this gives `v>=182`.
Thus `rho>=182/12=91/6>13`, contradicting compactness. ∎

This proof is the promised use of the ladder's toothpick self-similarity: one
must kill all `d` copies, not only the twelve representatives near zero.

## 3. What `n=12` rigidity would and would not finish

Suppose the full `n=12` tight classification were known:

```text
|W|=12 and M(W)=1/13  implies  W=d{1,...,12}.              (R12)
```

For a primitive Cover14 compact strict cover, Theorem B rules out every tight
deletion supplied by `(R12)`. Theorem A therefore forces the family into the
all-loose crown:

```text
M(V\{v})>1/13 and
M(V\{v})-1/13 < max(V\{v})/(13v)   for every v in V.       (ALL-LOOSE)
```

So proving the sporadic branch empty at `n=12` is important but does **not by
itself** close the compact LRC(14) inverse theorem. It eliminates the
tight-core branch and exposes the next exact obligation.

> **Crown-collapse obligation.** Prove that every primitive Cover14 family
> with `rho<13` and `M<1/13` has at least one tight deletion core.

`(R12)` plus crown collapse plus Theorem B would give a contradiction. If one
wants to use only the new shallow `A_12` carrier, the extraction must be
stronger: the tight deletion must contain no multiple of `13`. A generic
`u->0+` sheet then gives twelve edges all incident to sheet zero; coverage of
the other twelve sheets forces their inverse-residue endpoints to be all of
`F_13^*`. Thus the core is exactly a shallow full-residue packet and the
`A_12` ballot/regeneration theorem applies.

The extraction is real new mathematics, not set algebra. It is false in the
purely topological/divisor-complete analogue. Exact examples are:

| threshold | family `V` | `M(V)` | minimum deletion value |
|---|---|---|---|
| `1/3` | `{1,3,4}` | `2/7` | `2/5` |
| `1/4` | `{1,3,4,5}` | `2/9` | `2/7` |
| `1/5` | `{1,4,5,6,7}` | `2/11` | `2/9` |
| `1/6` | `{1,2,5,6,7,8}` | `2/13` | `1/5` |
| `1/7` | `{1,4,5,6,7,11,16}` | `2/15` | `2/13` |

Each row is primitive, covers every modulus `2,...,n+1`, is compact
`max/second<n`, lies strictly below `1/n`, and has every deletion strictly
above `1/n`. Hence essential-region topology, minimality, divisor coverage,
and compactness do not imply crown collapse uniformly in `n`.

At `n=13`, `V0` is the exact all-loose model and misses precisely modulus
`14`. Therefore any successful `n=13` crown-collapse proof must consume the
specific interaction of Cover14 with the thirteen-grid, not merely quote a
general cover theorem.

The weakest convenient sufficient composition is smaller than full Cover14:

> **Proposition C (corrected extraction implication).** There is no thirteen-
> speed set `V` satisfying all of the following:
>
> 1. `M(V)<1/13`;
> 2. some deletion `W=V\{v}` has `M(W)=1/13` and the available `n=12`
>    rigidity theorem identifies `W=d{1,...,12}`;
> 3. `gcd(V)=1`;
> 4. `V` contains a multiple of `14`; and
> 5. `rho(V)<91/6`.

Indeed Theorem B gives `13d|v`; primitivity gives `d=1`; the core
`{1,...,12}` has no 14-carrier, so `14|v`; hence `182|v` and
`rho(V)>=182/12=91/6`. Full Cover14 is used only through item 4, and the
project's compact hypothesis `rho<13` is stronger than item 5. Thus the exact
missing input is solely item 2's **tight-deletion extraction** (plus whichever
tight-core classification one invokes). No Cover2--Cover13 clauses enter the
post-extraction argument.

## 4. The exact sheet carrier beyond the shallow chart

The `A_12` walk extends exactly, but the extension explains why the shallow
state space is not enough. Write `t=(j+u)/13`.

- If `w=13h+r` with `r!=0`, then in a generic chamber
  `k<wu<k+1`, runner `w` is dangerous on the sheet edge

  ```text
  {-k r^{-1}, -(k+1) r^{-1}} in F_13.
  ```

- If `r=0`, then `||wt||=||hu||` independently of `j`: the runner is a
  **curtain**, covering all thirteen sheets or none.
- At an exact wall, strict cover (`M<1/13`) and closed cover
  (`M<=1/13`) differ. A nonzero-residue speed contributes the central sheet
  only in the strict convention, but the union of the adjacent edges in the
  closed convention. Curtain boundary events likewise toggle strict versus
  closed membership.

Thus the shallow `A_12` word is the no-curtain, twelve-edge specialization.
An arbitrary compact `INVcov` row requires an **edge-curtain event word**, its
private/triple-overlap chamber complex, and the integer lift sidecar of
THM-1099. Forgetting any one of those loses a named part of the predicate.

The Farey-comb regeneration target remains exact on the tight shallow face;
the crown theorem says a strict compact counterexample need not lie on that
face at all.

## 5. Tournament and alternate-carrier audit

Runners are not the only plausible vertices. This session tested:

1. **private stalks** (components of the `E_v`), oriented by cyclic midpoint
   after cutting at zero;
2. **wall events** `(v,k)`, ordered by `k/v` with simultaneous walls
   contracted;
3. **sheet cuts** in `F_13`, carrying the root-current vector; and
4. **Cover14 obligations/lift sheets**, with hyperedges for shared carriers.

For `V0` there are 94 positive private stalk chambers. Their midpoint-order
tournament is tautologically transitive: score histogram `0,...,93`, zero
directed triangles, 94 singleton SCCs, and one Hamiltonian path. It records
cyclic order and nothing about the weighted triple-excess identity. The event
tournament retains chronology but loses simultaneous higher incidence if
reduced to pairwise order. The owner-cost tournament of THM-1099 similarly
forgets shared lifts.

The smallest faithful carrier found here is therefore not a tournament:

```text
edge-curtain mechanical word
  + private-stalk / >=3-overlap chamber complex
  + root labels and tooth addresses
  + Cover14 owner-lift congruence sidecar.
```

Tournament fingerprints remain useful telemetry for edge flips and ordering,
but the proof predicate lives in a labelled hypergraph. This also refines the
Fano/`chi_7` probe: if a seven-point quotient is attempted, its vertices should
be triple-overlap obligations or lift-compatible chambers, not bare runners.
The balance law is exactly the information a pair shadow must be required to
retain.

## 6. Scope correction to HYP-7665

The S113 prose says “compact covering `rho<13 => M>=1/13` is equivalent to
LRC(14).” The exact statement still available after THM-1099 is narrower:

```text
primitive + Cover(2..14) + M<1/13  implies  rho>=13.        (INVcov)
```

Its contrapositive is a sufficient compact-floor route to LRC(14); no reverse
implication from LRC(14) to this stronger `1/13` floor has been proved.
Primitivity and Cover14 are indispensable:

- without primitivity, `2{1,...,13}` is Cover14 and compact but has
  `M=1/14<1/13`;
- without Cover14, `V0` above is primitive and compact with `M=8/105<1/13`;
- the S113 script's `range(2,14)` tests only moduli `2,...,13`, despite the
  Cover14 wording.

So “equivalent” should be read as “the chosen sufficient residual” unless an
additional converse is supplied. This is consistent with THM-1099 section 6,
which explicitly records only the implication chain.

### Exact S113/S114 audit after the latest pull

The source-level scope issues are sharper than the prose suggests.

1. `lrc14_compact_descent_boxeph_S113.py` checks one displayed boundary row,
   not the “100+” or “15” families discussed in the reflection. Its
   `covering` function omits modulus 14. Its bounded `QMAX=250` evaluator is a
   search rather than the repository's complete pair-sum ruler, although it
   happens to be ample for the displayed speed-24 row. The variable
   `core=V[:-1]` removes speed 24 and is
   `{2,4,6,8,10,12,13,14,16,18,20,22}`, **not** the hidden dilated core
   `2{1,...,12}` obtained by removing 13. Hence the descent diagnostic and the
   “dilated-AP core” diagnostic refer to different deletions.
2. The newly pulled S114 lines 11--18 identify

   ```text
   M(V)<1/13 + covering -> V\{v_max} is a dilated AP
   ```

   with the twelve-speed equality theorem

   ```text
   M(C)=1/13 -> C is a dilated AP.
   ```

   The forward use of the equality theorem is missing exactly the assertion
   `M(V\{v_max})=1/13`; Theorem A allows it to be strictly larger. More
   flexibly choosing a deletion still requires crown collapse. No reverse
   embedding of an arbitrary tight `C` into an admissible strict Cover14 row
   is supplied either. Thus the claimed equivalence is not established in
   either direction by the cited text.
3. S114's “dilated-AP-core compact” bullet invokes THM-1013 too broadly.
   THM-1013's proved witness requires **every** speed to be at distance at
   least `d` from `13d Z`; an arbitrary extra speed need not satisfy this.
   For strict `M<1/13`, Theorem B gives the opposite conclusion
   `13d|v`. The primitive/14-carrier/ratio contradiction above is the correct
   scoped discharge.
4. THM-1143 itself has the honest scope: its final paragraph explicitly says
   shallow ballot rigidity does not close the compact residual without an
   additional bridge. That statement agrees with the crown analysis and
   should control the older S113/S114 equivalence wording.

### Counterexample-first hard tests

- The admissible Cover14 lifts of the exact `V0` chart are
  `7 -> 112+210k`. All six compact lifts (`k=0,...,5`, carriers through 1162)
  have exact value `3/20`, at the same witness `3/20`; none preserves the
  strict crown. This is positive evidence that the lift sidecar has the right
  sign, not a uniform proof.
- In the complete single-substitution neighbourhood of
  `{2,4,...,24} union {13}`, with candidate speed at most 500, there are 2,748
  legal primitive Cover14 compact rows. A sound `q<=200` witness rejects all
  but eleven from having `M<=1/13`; complete pair-sum evaluation shows those
  eleven are exactly the equality rows obtained by replacing `13` with
  `39,65,...,299`. Every survivor retains the regenerated core `2{1,...,12}`;
  there are zero below-threshold rows. This is finite evidence only.

## 7. Smallest remaining obligation

After integrating all three lenses, the compact frontier is:

```text
strict compact INVcov row
  -> tight deletion OR all-loose essential crown                 [proved]

tight deletion + full n=12 rigidity
  -> regenerated d[12] core
  -> 13d | extra speed
  -> primitive Cover14 forces extra >=182 and rho>13             [proved]

all-loose crown
  -> thirteen inequalities (CROWN)
     + private/triple mass balance
     + edge-curtain word
     + Cover14 lift sidecar                                      [exact residual]
```

Accordingly the smallest genuinely new target is **not another reformulation
of shallow ballot rigidity**. It is the `n=13`-specific crown-collapse (or a
direct contradiction to `(ALL-LOOSE)`) using the Cover14 lift congruence.
The `7->112` pair is the mandatory positive-control test: a proposed invariant
must change under that lift even though the old `8/105` Farey chart does not.

Exact audit:
`04-computation/lrc14_compact_essential_crown_bridge_codex_20260718.py` and
`05-knowledge/results/lrc14_compact_essential_crown_bridge_codex_20260718.out`.

Normal and optimized Python executions are byte-identical to the frozen
output. Hashes after the hardening audit:

```text
source  335c42ae4d08426df09935bec629e62724e4b9fe3f15e1f0cd13bc089d2fcc55
output  db54526f86997674183824f9c1eec5688662be7ae147a2e31ae687a11dacb0ec
```
