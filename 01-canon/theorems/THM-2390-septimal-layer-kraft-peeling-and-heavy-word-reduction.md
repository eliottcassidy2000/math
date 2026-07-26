---
id: THM-2390
title: "Septimal layer Kraft peeling and heavy-word reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a
  finite family of centred open danger combs of widths a_i/7,
  1<=a_i<=6, put W_e=sum_(nu_7(v_i)=e)a_i. If every W_e<7, the common
  closed safe set has Haar measure at least
  product_e(7-W_e)/7. The proof is a pointwise seven-root peel: at the
  least remaining valuation, higher layers are invariant and the danger
  words use at most W_e of the seven roots. In THM-2367's last
  k=2, (t,b)=(1,0) alternative, the two higher masks each form a
  singleton layer and the lower weight is eight. If no lower layer has
  weight at least seven, optimizing all partitions of eight gives the
  explicit uncovered floor 180/2401. Hence every survivor forces the
  guard and at least five of the six lower ordinary labels into one
  septimal layer. On a positive-measure family of generic fibres, weight
  seven gives an exact seven-root partition and weight eight gives exactly
  one doubled root. This is a reduction only: it does not empty the
  branch, decrement the 165-row ledger, land a target, or prove LRC(14).
source: codex-2026-07-26-septimal-layer-kraft-peeling
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
related:
  - THM-2095-p-adic-guard-ratio-scale-terminal
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2381-one-top-one-blocker-septimal-root-stalk-closure
  - THM-2382-saturated-septimal-seven-bin-root-fibre-closure
  - THM-2385-two-top-septimal-blocker-collision-reduction
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
script: 04-computation/lrc14_septimal_layer_kraft_peeling_thm2390.py
output: 05-knowledge/results/lrc14_septimal_layer_kraft_peeling_thm2390.out
script_sha256: 9b9765372f1a1b66d384b416afedbd658038f12ef30dc251f4676fbbcc66b2c2
output_sha256: a793f3d2cb55a6e9249edd1f0f4bb2d4c60833adcc485ebcc97a422fb1d3c847
hash_basis: working-tree bytes (LF)
---

# THM-2390 -- septimal layer Kraft peeling

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The closest proved mechanism is THM-2367's single top-layer count.
THM-2377 supplies the corrected spectral near miss: pairwise-distinct
lower valuations force a selected target tensor to be circulant, but a
single repeated layer can already carry every nonzero target colour.  The
present argument retains positivity instead of a chosen Fourier channel
and peels **every** septimal layer.

The result is an exact Kraft-type alternative:

```text
every valuation layer has danger weight at most six
  -> a quantitatively positive common-safe set;

an almost-everywhere cover
  -> some valuation layer has danger weight at least seven.          (1)
```

Aligned strict endpoints are a sharp local hostile to an equality count:
an open arc of length `a/7` can contain `a-1`, rather than `a`, points of
an aligned seven-grid.  This helps the safe lower bound and never harms it.

## 1. The general septimal peeling lemma

For a positive integer `v` and an integer `1<=a<=6`, put

```text
D_(v,a)={x in T:||vx||<a/14},

G_(v,a)=T minus D_(v,a)={x:||vx||>=a/14}.             (2)
```

Thus `D_(v,a)` is an open centred danger comb of Haar measure `a/7`.
Ordinary LRC danger masks have `a=1`; the scalar guard-danger mask has
`a=2`.

Let

```text
F={(v_i,a_i):1<=i<=s}
```

be any finite labelled family, with repetitions allowed.  For each
septimal valuation layer define its total danger width

```text
W_e=sum_(i:nu_7(v_i)=e)a_i.                           (3)
```

> **Septimal layer Kraft lemma.** If
>
> ```text
> W_e<=6                         for every e,          (4)
> ```
>
> then
>
> ```text
> mu(intersection_i G_(v_i,a_i))
>   >=product_(e:W_e>0)(7-W_e)/7
>   >0.                                                (5)
> ```

The factors in (5) are exact pointwise root-capacity factors.  Equality
is not asserted for an arbitrary family.

### One peel

Let `e` be the least valuation occurring in the current family and split

```text
P_e(x)=product_(i:nu_7(v_i)=e)1_(G_(v_i,a_i))(x),

P_>(x)=product_(i:nu_7(v_i)>e)1_(G_(v_i,a_i))(x).      (6)
```

For `k in F_7`, use the seven roots

```text
x_k=x+k/7^(e+1).                                      (7)
```

If `nu_7(v_i)>e`, then

```text
v_i(x_k-x) in Z,
```

so every higher factor is constant on the fibre:

```text
P_>(x_k)=P_>(x).                                      (8)
```

If `nu_7(v_i)=e`, write `v_i=7^e u_i`, with `u_i` a
seven-unit.  Its seven phases are

```text
v_i x+u_i k/7,                  k in F_7.             (9)
```

Multiplication by `u_i` permutes the seven-grid.  An open circle arc of
length `a_i/7` contains at most `a_i` points of that grid.  Hence the
union bound on the root labels gives, pointwise for every `x`,

```text
#{k:P_e(x_k)=1}>=7-W_e,                               (10)

(1/7)sum_k P_e(x_k)>=(7-W_e)/7.                       (11)
```

At an aligned endpoint a factor can lose one danger root; (10) only
improves.  No genericity or boundary reassignment is required for the
lower bound.

Average (6) over the seven translations.  Haar invariance and (8) give

```text
integral_T P_e P_>
 =integral_T P_>(x)(1/7 sum_k P_e(x_k))dx
 >=(7-W_e)/7 integral_T P_>.                         (12)
```

This is the conditioning order: remove the **least** remaining
valuation, while every unpeeled higher layer is literally invariant on
its root fibre.  Iterating (12) through the finitely many occupied
layers proves (5).  There are no lower layers left inside a peel, and no
Fourier series, limiting interchange, or independence assumption is
used. QED.

## 2. Specialization to the last `(t,b)=(1,0)` lane

Retain THM-2367's only-`c_3`-dominant notation:

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))>0,

nu_7(H)<M,

nu_7(c_3)>M,

k=#{j in {1,2}:nu_7(c_j)<=M}=2.                      (13)
```

In the remaining profile

```text
(t,b)=(1,0),                                         (14)
```

there is a unique top unit `q_*` with

```text
nu_7(q_*)=M,                                         (15)
```

while the four other `q_i` and both `c_1,c_2` have valuation below
`M`.  The nine danger factors of the scalar cover are

```text
(H,2),

(q_1,1),...,(q_5,1),

(c_1,1),(c_2,1),(c_3,1).                             (16)
```

The `q_*` layer has weight one.  The `c_3` layer is strictly higher and
also has weight one.  The lower factors

```text
(H,2), (q_i,1)_(q_i!=q_*), (c_1,1),(c_2,1)          (17)
```

have total weight

```text
2+4+2=8.                                             (18)
```

Suppose every lower layer has weight at most six.  If its nonzero layer
weights are the integer partition

```text
8=w_1+...+w_r,                  1<=w_j<=6,            (19)
```

then the Kraft lemma gives

```text
mu(common safe set)
 >=(6/7)^2 product_j(7-w_j)/7.                       (20)
```

The first factor is from the two distinct singleton layers (15) and
`nu_7(c_3)>M`.

### Exact optimization of the lower partition

If `a+b<=6`, merging two parts strictly decreases the product:

```text
(7-a)(7-b)/49-(7-a-b)/7=ab/49>0.                    (21)
```

Every partition of eight into at least three parts has two parts whose
sum is at most six, so repeated merging reduces the minimization to

```text
8=6+2,                 8=5+3,                 8=4+4. (22)
```

Their lower products are respectively

```text
5/49,                  8/49,                  9/49.  (23)
```

Thus the sharp partition bound used here is

```text
product_j(7-w_j)/7>=5/49,                            (24)
```

and (20) becomes

```text
mu(common safe set)>=
  (6/7)^2(5/49)
 =180/2401
 >0.                                                 (25)
```

But the scalar cover says that the common safe set of all nine factors is
null.  Therefore:

> **Heavy-layer reduction.** Every surviving `(t,b)=(1,0)` packet has a
> lower septimal layer `e<M` with
>
> ```text
> W_e>=7.                                             (26)
> ```

The six lower ordinary labels have total weight six.  Consequently the
guard must belong to this layer, and

```text
nu_7(H)=e,

#{v in {q_i:q_i!=q_*} union {c_1,c_2}:nu_7(v)=e}
  >=5.                                                (27)
```

Since the entire lower weight is eight, the heavy layer has weight
exactly seven or eight.  This is strictly stronger than the bare
pair-collision conclusion furnished by the contrapositive in THM-2377.

## 3. The forced terminal seven-root words

The reduction also identifies the finite word that remains.

Let `e` be the heavy layer.  There is at most one lower ordinary factor
outside it.  Peel that factor first if its valuation is below `e`; the
one-layer estimate costs `6/7>0`.  It follows from the null common-safe
set that the heavy and higher factors still have null common-safe set.

Now average the heavy product over

```text
x+k/7^(e+1),                    k in F_7.             (28)
```

All factors above `e` are constant on this fibre.  Their valuation layers
have weight one: they consist of the possible remaining lower ordinary
factor above `e`, the top unit `q_*`, and the higher blocker `c_3`.
Applying the Kraft lemma to those higher layers gives a common-safe base
set of measure at least

```text
(6/7)^3=216/343                                      (29)
```

when all three are present, and at least `(6/7)^2=36/49` when the
remaining lower factor was already peeled or does not exist.

On almost every base in this positive-measure set, the heavy danger words
must cover all seven roots; otherwise one root would be safe for every
factor.  Discard the finite set of aligned heavy endpoints.  Every
ordinary heavy word then has size one and the guard word has size two.
Therefore:

1. if `W_e=7`, the two-point guard word and five singleton ordinary
   words form an exact partition of `F_7`;
2. if `W_e=8`, the two-point guard word and six singleton ordinary
   words cover `F_7` with multiplicity one at six roots and multiplicity
   two at exactly one root.

The second case has one double, not a hidden triple: the total generic
incidence is eight and every one of the seven roots is covered.

These are labelled words.  The two guard roots are adjacent only in the
guard coordinate, and the ordinary singleton positions retain their
seven-unit slopes.  Replacing them by an unlabelled cycle or a cosmetic
tournament would discard the exact object needed by the next step.

## 4. Relation to the incoming blocker reflection

THM-2388 is a proved, exact candidate under independent audit, and is not a
dependency here.  For the six thirteen-unit guard/`q` masks, its
multiplicity reflection

```text
P(K-1)=-(K-1)                                         (30)
```

exchanges unit holes and unit overlaps across a thirteen-root fibre and
cages them by quotient blockers.  The present theorem does not use (30).
It says that any packet on which that reflection is applied has already
collapsed septimally to one of only two heavy words in Section 3.

The cheapest next test is therefore finite and cross-prime:

```text
labelled 7-bin partition / one-double word
  + thirteen-root hole/overlap reflection
  + THM-2388's exact 36/343 multiplicity-excess current. (31)
```

No implication from (31) is asserted here.  In particular, a repeated
septimal layer is not itself contradictory: THM-2377 gives a physical
same-layer carrier with all twelve nonzero target colours.  That carrier
is the canonical stopping boundary to any claim that collision alone
empties the branch.

## 5. Scope

This theorem proves a structural reduction of the sole remaining
septimal alternative.  It does not prove:

- that either heavy word in Section 3 is impossible;
- that the heavy word aligns with a canonical owner or target;
- that the thirteen-root reflection preserves a terminal phase;
- that any one of the `165` thirteen-adic profiles is empty; or
- LRC(14).

The ledger remains `165`.

## 6. Exact companion

The dependency-free exact companion:

- enumerates the boundary chambers of an open `a/7` arc on a
  seven-grid for `a=1,...,6`, obtaining generic count `a` and aligned
  count `a-1`;
- checks the pointwise root-capacity inequality `7-W` for every
  `W=0,...,6`;
- enumerates every integer partition of lower weight eight with parts
  at most six and finds the unique product minimum `5/49` at `(6,2)`;
- verifies the full floor `180/2401`;
- exhausts the labelled heavy-role choices, giving six weight-seven
  choices and one weight-eight choice; and
- exhausts all abstract size-two guard words and labelled singleton
  words on `F_7`.  Every full weight-seven word is a partition; every
  full weight-eight word has the unique one-double profile.

Run

```bash
python3 04-computation/lrc14_septimal_layer_kraft_peeling_thm2390.py
python3 -O 04-computation/lrc14_septimal_layer_kraft_peeling_thm2390.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_septimal_layer_kraft_peeling_thm2390.out
```

after LF normalization.  Every executable check raises explicitly under
optimized Python.

An independent hostile audit reconstructed the least-valuation peel,
checked the strict-endpoint direction, independently minimized every
capped partition of eight, verified the `W=7,8` fibre conclusions, and
replayed the normal, optimized, and stored transcripts with the declared
LF-normalized hashes.  It found no proof or scope defect.
