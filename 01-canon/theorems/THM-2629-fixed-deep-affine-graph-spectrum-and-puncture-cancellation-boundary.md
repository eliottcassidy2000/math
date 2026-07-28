---
id: THM-2629
title: "Fixed-deep affine graph spectrum and puncture-cancellation boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Retaining one nonzero translated deep-probe label before THM-2616's
  seven-clock unit test gives 4,156 nonzero labelled rows, each supported on
  all twelve nonzero deep labels.  Across the 84 base cells the fixed-deep
  unit relations have 72 common (q,r) pairs.  Exhausting all 169 affine
  graphs r=alpha q+beta with the lexicographic score (minimum cell coverage,
  common-section count) gives the unique optimum r=-q-1, with score (10,9),
  unit-set histogram 11^78,10^6, and common q-set
  {1,3,4,5,6,7,8,9,10}.  The graph maps the missing future sheet q=h=12 to
  the structurally missing deep sheet r=0.  It therefore exposes the exact
  completion boundary: filling the future C13 carrier without also filling
  deep r=0 cannot yield any full affine bijection.  This is a coefficient-
  level fixed-label atlas; rail choices vary by cell and no chronological,
  semantic, or ordered transition gluing is proved.
source: kind-pasteur-2026-07-28-fixed-deep-affine-graph
depends_on:
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
related:
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
script: 04-computation/lrc14_fixed_deep_affine_graph_spectrum_thm2629.py
output: 05-knowledge/results/lrc14_fixed_deep_affine_graph_spectrum_thm2629.out
script_sha256: 6dd02116b5847ee9c8c8a76bdff6f4231bde9f5c19e60f34952faed783c8c28e
output_sha256: 1366cbfe04c4e2a88bb68eb4ff100c08868bdb14144a93d5e4fb4c66487a34c2
hash_basis: LF-normalized bytes
---

# THM-2629 -- the best deep/future graph hides one puncture in the other

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2616 first keeps the present target section `q`, translated deep probe
`r`, and delayed physical digit `h` in one nonnegative common-`x` tensor,
then finds a large unit atlas on the numerical diagonal `q=h`.  Its unit test
sums all twelve nonzero `r` labels.  That leaves open a sharper question:

```text
does one affine deep label r=f(h) already carry the unit class?       (1)
```

The answer is unexpectedly structured.  Positivity itself is completely
blind to `r`: every nonzero labelled row contains all twelve deep labels.
After the seven-clock determinant, however, the opposite affine graph

```text
r=-h-1                                                     (2)
```

is uniquely best among all `169` affine rules.  Its constant-section set is
strictly larger than the aggregate diagonal's.  The improvement is not a
principal-action theorem.  It is a nonmonotone cancellation phenomenon in
the fixed-label Bockstein and it aligns the two missing sheets exactly.

Throughout, all residues lie in `F_13`, and `q=h` on the diagonal carrier.

## 1. Retain one deep label before taking the unit test

Use THM-2616's notation.  Its globally primitive diagonal tensor is

```text
a_(e,q;ell,r,q)=A_(e,q;ell,r,q)/26,       r!=0,            (3)
```

where `e` is one of the `162` positive middle rails and `ell in F_7` is the
future owner clock.  For one fixed nonzero deep label define

```text
Y^r_ell(e,q)=a_(e,q;ell,r,q) r^(-1) mod 13,

Y^r_(e,q)(z)=sum_(ell=0)^6 Y^r_ell(e,q) z^ell
                 in F_13[z]/(Phi_7).                       (4)
```

Call `(e,q,r)` a fixed-deep unit when (4) is a unit.  Division by the same
global content `26` is load-bearing: no rail, `q`, or `r` slice is
reprimitivized.  The scalar `r^(-1)` does not affect unitness, but it is kept
because (4) is the literal one-row summand of THM-2616's Bockstein.

For a base cell `c=(s,ell_4)`, unite only its two already proved rail choices:

```text
U_c={(q,r): some rail e over c makes Y^r_(e,q) a unit}.    (5)
```

This is still a labelled coefficient relation, not an identification of the
present probe with the delayed physical digit.

## 2. Support is maximally non-atomic

Before the unit test there are exactly

```text
4,156
```

nonzero labelled rows `(e,q,ell,h=q)`.  Every one has

```text
{r:A_(e,q;ell,r,q)>0}=F_13^*.                             (6)
```

Thus the deep-support histogram is exactly `12^4156`.  No support rule,
positive Markov left inverse, or singleton fibre can select `r` from this
carrier.  Any structure below comes from the clockwise cyclotomic class, not
from a hidden sparse physical row.

The exact fixed-deep unit counts per base cell are

| `#U_c` | 108 | 120 | 122 | 130 | 131 | 132 |
|:--|--:|--:|--:|--:|--:|--:|
| cells | 4 | 2 | 1 | 1 | 6 | 70 |

Their intersection has `72` pairs.  Written by target/future label, the
common allowed deep sets are

| `q=h` | common `r` |
|:--:|:--|
| 0,2,11,12 | empty |
| 1 | `1,...,11` |
| 3,6 | `1,...,12` |
| 4,5 | `7,...,12` |
| 7 | `1,2,3,5,6,7,8,9,10,11,12` |
| 8,9 | `1,...,6` |
| 10 | `1,2` |

This is a genuine strengthening of the aggregate description: it locates
which single deep labels survive every base cell.

## 3. The complete affine-graph spectrum

For every affine rule

```text
f_(alpha,beta)(q)=alpha q+beta,          alpha,beta in F_13, (7)
```

delete the forbidden value `f(q)=0` and put

```text
Q_c(alpha,beta)={q:(q,f_(alpha,beta)(q)) in U_c},

m(alpha,beta)=min_c #Q_c(alpha,beta),

n(alpha,beta)=# intersection_c Q_c(alpha,beta).           (8)
```

Exhaustion of all `169` maps gives

```text
max m=10.                                                  (9)
```

Exactly four maps attain (9):

| `(alpha,beta)` | cell-size histogram | common `q` | `n` |
|:--:|:--|:--|--:|
| `(-1,-1)` | `10^6,11^78` | `1,3,4,5,6,7,8,9,10` | 9 |
| `(5,0)` | `10^7,11^77` | `1,3,4,5,6,7,8,9` | 8 |
| `(2,0)` | `10^7,11^77` | `1,3,4,5,6,7,8,9` | 8 |
| `(-1,0)` | `10^9,11^75` | `3,4,5,6,7,8,9` | 7 |

Hence `(alpha,beta)=(-1,-1)` is the unique maximizer of the lexicographic
score `(m,n)`.  Its six exceptional cells are

```text
(2,2),(2,3),(7,2),(7,3):   q=1,...,10;

(6,5),(11,5):              q=1,3,4,5,6,7,8,9,10,11.      (10)
```

The other `78` cells retain every `q=1,...,11`.  At the weaker positivity
level the same graph retains all eleven labels on `81` cells and ten labels
on the remaining three.  Thus its unit losses are spectral rather than
support-theoretic.

Three controls show that the optimum is not forced by a trivial convention:

```text
r=q:       (m,n)=(8,4),

r=7q:      (m,n)=(9,6),

best constant r=2: (m,n)=(8,7).                           (11)
```

In particular, the natural identity and inverse-doubling guesses are both
strictly worse.

## 4. Exact puncture cancellation

On THM-2616's safe delayed-word carrier,

```text
H_future={1,...,11},                  R_deep=F_13^*.       (12)
```

The sharp graph (2) satisfies

```text
{-h-1:h=1,...,11}={1,...,11},

h=12 -> r=0,                         h=0 -> r=12.          (13)
```

Thus one absent future sheet, `h=12`, hides the absent deep sheet `r=0`,
while the other absent future sheet leaves `r=12` unused.  This is the exact
puncture-cancellation alignment exhibited by the affine optimum.

It also proves a sharp completion no-go independent of any particular guard
construction.  Every affine bijection `h -> alpha h+beta`, `alpha!=0`, maps
one and only one `h` to zero.  Therefore:

```text
full future C_13 carrier + permanently absent deep r=0
  => no full affine deep/future graph.                     (14)
```

Filling `h=0` can extend (2) through `r=12`; filling `h=12` cannot extend it
unless a second physical chart also restores `r=0`.  Complementing only the
future role repairs the old eleven-sheet puncture but cannot by itself create
a principal affine bitorsor.

Algebraically, (2) intertwines `h -> h+a` with the **opposite** action
`r -> r-a`.  This is a legitimate abstract opposite-action graph.  The
current theorem does not prove that the physical deep-probe action has that
orientation, nor that the constant `-1` is a chronological carry rather than
an extremal coefficient choice.

## 5. Exact evidence and stopping boundary

Run

```text
python 04-computation/lrc14_fixed_deep_affine_graph_spectrum_thm2629.py
python -O 04-computation/lrc14_fixed_deep_affine_graph_spectrum_thm2629.py
```

The companion rebuilds all four THM-2616 shards, the `162` rails, all `4,156`
nonzero labelled rows, every fixed-`r` multiplication determinant, all `84`
base relations, and all `169` affine graphs.  It checks the complete common
pair table, the sharp four-map bank, the six exceptional cells, the three
controls (11), and the elementary full-affine puncture obstruction (14).
Every branch is exact integer or finite-field arithmetic; no floating point
or random fixture is used.  Normal and optimized runs are required to
byte-match the stored transcript.

An independent hostile audit rebuilt the fixed-`r` classes on the same
global `/26` primitive lattice, replayed normal and optimized modes against
the stored output, re-exhausted all `169` affine maps, and checked the six
exceptional cells and both puncture identities.  It also verified the
load-bearing scope: `r` is a translated-probe coefficient label and the
realizing rail may vary with the base cell, so the opposite affine action is
not being asserted as physical chronology.

The result selects a labelled deep-probe index inside a coefficient tensor.
The realizing rail may vary from cell to cell, and the seven-clock unit is an
aggregate over `ell`.  No common chronological stalk, outgoing-root equals
incoming-root identity, positive ordered seven-edge path, semantic owner,
relation residue, or endpoint current follows.  It removes no scalar row and
does not prove LRC(14).  QED.
