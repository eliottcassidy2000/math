---
id: THM-1274
title: CENTERED PROTRUSION ORIENTS THE SLOWEST TWO-CYCLE AND FORCES A FIVE-OWNER RETURN OR TERMINAL TAIL
status: PROVED (the THM-1267 protrusion sign fixes the binary descent digit; the protrusion-facing blocker wall continues through only five owner labels; distinct tooth-instance reuse gives a literal closest return of at most five edges and strict factor greater than 3/2; otherwise the word reaches the carrier endpoint within five tooth occurrences; the exceptional six-owner path either has a closest 6/5 return or is the whole six-tooth word, forcing d6<=7c; exact mirrored guardrails; optimization-safe exact referee; sorry-free Lean arithmetic/combinatorial core). This is a return-or-endpoint landing theorem, not an endpoint tax, six-comb noncoverage, or LRC(14)
source: codex-2026-07-19-S82 protrusion/return composition
depends_on: [THM-1248, THM-1252, THM-1253, THM-1256, THM-1262, THM-1264, THM-1266, THM-1267]
related: [THM-1244, THM-1283, HYP-7870]
script: 04-computation/lrc14_centered_tail_return_or_terminal_thm1274.py
output: 05-knowledge/results/lrc14_centered_tail_return_or_terminal_thm1274.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCenteredTailReturnOrTerminal.lean
script_sha256: 22b4a8d460fe76fd26ce040f3713cda2dcf228b4fe86d1d888d2f20733ab8566
output_sha256: 8a01b41177225116c59591d91c7b1f2b6af984de0bdabdf12ff4acc3a7e2c29b
formalization_sha256: 16597fe2a2611d65edfedf729b74eb8bda9232f70b20e49b8f1460e9ce895bcb
---

# THM-1274 — centered-tail digit and five-owner return-or-terminal

## 1. Setup and the oriented protrusion

Let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)],
c<d_1<...<d_6,
```

and suppose the six strict fast danger combs cover `G`.  Put `a=d_1` and
let `t_a` be THM-1240's centered `a`-spoke.  THM-1267 supplies the complete
closed `a`-safe component `S_a` through `t_a`, a nonzero rounding error
`epsilon_a`, and a sign

```text
sigma=sign(epsilon_a) in {-1,+1}.                    (1)
```

The component protrudes through exactly the `sigma` endpoint of `G`.  Its
normalized exterior length and the slowest-fast ratio satisfy

```text
ell=|S_a\G|/|S_a|>11/270,
270a<=563c-1.                                        (2)
```

Choose THM-1253's deletion-minimal chronological cover by individual fast
teeth, and choose all centered blockers from this word.  Suppose the blocker
lasso contains the slowest-rooted two-cycle

```text
a -> b -> a,                 a<b.                    (3)
```

Write

```text
H = the selected b-tooth containing t_a,
L = the selected a-tooth containing t_b.             (4)
```

The ascent protection in THM-1252/1262 gives

```text
closure(H) subset int(G intersect S_a),              (5)
```

whereas `L` is strictly `a`-dangerous and therefore disjoint from `S_a`.
Both teeth lie in `G`.

## 2. The protrusion sign fixes the descent digit

Because `S_a` is shorter than `G`, meets it, and protrudes through only one
endpoint, `G\S_a` is one interval on the opposite side.  Equation (5) puts
`H` in `S_a`, while `L` lies in `G\S_a`.  Hence

```text
sigma=-1  =>  t_a<t_b,
sigma=+1  =>  t_b<t_a.                               (6)
```

For the descent edge `b->a`, THM-1248/1256 gives a binary relative digit
`delta_ba`, with

```text
delta_ba=0 iff t_a<t_b,
delta_ba=1 iff t_b<t_a.                              (7)
```

Combining (6)--(7) yields the exact **tail-digit law**

```text
delta_ba=(1+sigma)/2.                                (8)
```

Thus the left protrusion forces digit zero and the right protrusion forces
digit one.  This is not an extra binary assumption: it is forced by the
global position of the centered safe component.

Every coherent blocker edge is already phase/word aligned by THM-1256.  The
THM-1262 corridor from `H` toward the reverse target `L` therefore leaves the
wall of `H` opposite the protrusion.  The other wall of `H` is canonically
the **protrusion-facing wall**.

## 3. The protrusion-facing continuation uses only five owners

Let `e_sigma` be the endpoint of `G` through which `S_a` protrudes, and let
`w_sigma` be the protrusion-facing wall of `H`.  The whole closed segment

```text
C_sigma=[w_sigma,e_sigma]                             (9)
```

lies in `S_a intersect G`.  Follow the deletion-minimal tooth word from `H`
in that direction until it covers `e_sigma`: a chronological suffix for a
right protrusion and the reversed chronological prefix for a left protrusion.
Every selected tooth on this continuation meets the interior of `S_a`; an
`a`-danger tooth cannot do so.  Consequently this tail subword uses only the
five owner labels

```text
{d_2,...,d_6}.                                       (10)
```

This is the operation-level gain missing from the abstract return DAG: the
canonical blocker wall chooses a literal consecutive continuation and a
specific terminal endpoint.

Split by **tooth occurrence**, not merely by owner label.

### 3.1 Distinct-occurrence reuse gives a `3/2` return

If an owner repeats at two distinct word positions, choose a closest pair
`p<q`.  Adjacent selected teeth have different owners, so `m=q-p>=2`.
Closestness makes all internal owners distinct and different from the outer
owner.  Since (10) has only five labels,

```text
2<=m<=5.                                             (11)
```

Let `A` be the repeated outer speed, let `R>=1` be its tooth-address jump,
and let `s_r` be the internal owner speeds.  THM-1264's literal return
identity gives

```text
sum_(r=p)^(q-1) omega_r
  =(1/7)sum_(r=p)^(q-1)1/s_r-R/A>0.                  (12)
```

If `B=min_(p<r<q)s_r`, then

```text
A>((7R-1)/(m-1))B >=(6/4)B=(3/2)B.                  (13)
```

All seams in (12) were already paid by THM-1253.  The new information is a
literal address-holonomy cell attached to the protrusion-facing blocker wall,
not additional overlap mass.

### 3.2 No reuse reaches the endpoint in five occurrences

If no owner repeats, the five-label word (10) has at most five tooth
occurrences.  Therefore the same canonical continuation reaches `e_sigma`
without a return after at most five selected teeth, including `H`, and at
most four handoffs.  This is the terminal alternative:

```text
protrusion-facing continuation
  => closest return with factor >3/2
     or carrier endpoint in at most five occurrences.                (14)
```

The second branch is not discharged inside this theorem.  It is a sharply
localized endpoint obligation carrying `(sigma,epsilon_a,ell)`, the last
selected tooth, and the positive five-comb survivor outside `G` supplied by
THM-1267.  THM-1283 subsequently consumes this payload: the terminal tooth
crosses the adjacent carrier tooth in a new exterior gcd seam, and removing
that tooth segment from the survivor tail yields a strict residue-dependent
protrusion tax.
The six-cover hypothesis applies on `G`; it does not permit the selected word
to be continued through that exterior survivor.

## 4. The exceptional six-owner path is return or compact terminal

THM-1256's incidence-free/reuse-free two-cycle residual contains a literal
chronological path through all six owner labels:

```text
h_L -> b -> h_1 -> ... -> h_2 -> a -> h_R.           (15)
```

Its four wall seams plus the first internal crossing handoff form a located
six-label spanning tree, and its descent corridor carries the nested tariff.
Equation (15) now has a further exact split.

If the full minimal word contains any tooth beyond one occurrence of each of
the six labels, some owner occurs at two distinct positions.  A closest pair
has at most six edges and THM-1264/1266 gives

```text
A>(6/5)B.                                            (16)
```

If there is no such return, (15) is the entire word: exactly six teeth, one
per owner.  Write `n_i` for the number of selected teeth owned by `d_i`.
THM-1253's private-tooth count gives

```text
6=N=sum_i n_i,
n_i>=ceil(d_i/(7c))>=1.                              (17)
```

Every term in (17) is therefore one.  Hence

```text
d_i<=7c for every i,             d_6<=7c.            (18)
```

Thus the exceptional branch is no longer an untyped slowest two-cycle:

```text
distinct-occurrence closest return with factor >6/5,
or a six-tooth Hamiltonian terminal cell with d_6<=7c.                (19)
```

## 5. Corrected recursive trichotomy

The THM-1248 two-wall lasso split must now be typed by tooth instances.

1. A wall owner lies on the lasso.
2. An outside owner label repeats.
3. The lasso is the exceptional slowest-rooted two-cycle.

In (1)--(2), two distinct selected occurrences give a closest return and
the THM-1264 rank.  If one selected tooth supplies both incidences, there is
no return: it is a located multi-incidence/turn obligation and must retain
both wall germs.  Label reuse alone must never be promoted to address return.

In (3), equations (14) and (19) apply simultaneously.  The inner object is
the located double-star-plus-bridge tree and nested seam tariff from THM-1256;
the outer object is the oriented protrusion and its five-owner continuation.
These credits are algebraically distinct: (12) classifies already-paid inner
seams, while `ell>11/270` is functional survivor mass outside `G`.  They may
not be added as two copies of interval measure.

## 6. Exact guardrails

THM-1248's binary exceptional packet

```text
(c;a,b;outside)=(2;4,6;15,19,28,43)
```

has `epsilon_a=-1/2`, left protrusion, phases

```text
t_a=1/6<t_b=1/4,
```

and descent digit zero.  Circle reflection gives `epsilon_a=+1/2`, right
protrusion, `t_b=3/4<t_a=5/6`, and digit one.  Both packets are globally
lonely.  They prove that both signs in (8) occur arithmetically and that the
tail-digit law is a state refinement, not by itself a contradiction.

THM-1266's primitive five-rung row realizes aligned blocker two-cycles and
all current local return laws while leaving four complete fastest-safe tails.
It is not a six-cover.  Any purported proof of (14) which rejects that row
before using endpoint completion has silently assumed the missing global
input.  More precisely, for its cycle `254<->1805`, the slowest rounding
error is positive, the protrusion is right, and the descent digit is one.
The protrusion-facing suffix begins with owner word
`(1805,255,1805)`, realizing the return branch with the stronger immediate-
backtrack factor `1805>6*255`; reflection realizes the left/digit-zero side.

## 7. Tournament and alternate-vertex audit

The speed tournament on the five non-`a` labels is transitive, with score
histogram `(0,1,2,3,4)`, no directed cycles, five singleton SCCs, and one
Hamiltonian path.  It cannot distinguish a repeated tooth occurrence from
two wall incidences supplied by the same tooth.

Use instead the protrusion-facing selected tooth occurrences as vertices.
The pairwise observable is equality of owner label together with equality or
inequality of tooth address.  Chronological order is the switch/gauge and the
tail word is the tie Hamiltonian path.  A distinct-address equality closes a
return polygon; absence of equality terminates at `e_sigma`.  The faithful
state is

```text
(sigma,epsilon_a,ell;
 protrusion-facing wall and carrier endpoint;
 ordered tooth instances with owner and address;
 paid seams; first repeated owner or terminal tooth).                 (20)
```

This quotient preserves the return-or-endpoint predicate.  Projecting to
runners destroys occurrence identity; projecting to wall events destroys
the common return address; projecting to interval mass destroys the endpoint
orientation.  The challenged assumption is that owner reuse is automatically
a graph cycle.  It is a return only after the two tooth instances are proved
distinct.

## 8. Verification and formal boundary

The optimization-safe dependency-free exact referee contains no Python
`assert` nodes.  It checks `27,305` adjacent-distinct five-owner words,
including `26,980` closest-return words and `325` no-return terminal words;
the sharp five-edge/five-occurrence bounds; `400` exact return-factor rows
with minimum `3/2`; all `720` six-owner Hamiltonian terminal paths and
`21,600` legal one-extra-occurrence insertions; the unique six-positive-count
vector of total at most six; `50,500` exact private-count speed-cap rows; and
both mirrored tail-digit guardrails.  Normal and optimized outputs are
byte-identical.

The sorry-free Lean module checks the abstract interval-side/digit law, the
five-label repeat/terminal pigeonhole, the `3/2` return consumer, the
six-label extra-occurrence split, and the terminal private-count speed cap.
Identification of `S_a`, the protrusion-facing wall, extraction of the
literal selected subword to the carrier endpoint, and THM-1264's occurrence
matching remain the named paper topology providers.  There are no proof
placeholders or `native_decide` calls.

Frozen artifact hashes are

```text
source         22b4a8d460fe76fd26ce040f3713cda2dcf228b4fe86d1d888d2f20733ab8566
output         8a01b41177225116c59591d91c7b1f2b6af984de0bdabdf12ff4acc3a7e2c29b
formalization  16597fe2a2611d65edfedf729b74eb8bda9232f70b20e49b8f1460e9ce895bcb
```

THM-1274 itself does not prove that the terminal branch pays new mass.
THM-1283 now supplies that downstream endpoint consumer: it exports a proper
`c`--terminal-owner crossing and strengthens the protrusion cut by its exact
endpoint residual.  The remaining operation after the pair is narrower: make
a return rank feed its next blocker, or follow the terminal owner's safe side
with the mixed above/below owner set while retaining the residue and disjoint
internal invoice.  The same-tooth multi-incidence branch also needs a strict
torsion/turn consumer.  These are smaller typed obligations than the former
arbitrary blocker walk, but global LRC(14) remains open.  ∎
