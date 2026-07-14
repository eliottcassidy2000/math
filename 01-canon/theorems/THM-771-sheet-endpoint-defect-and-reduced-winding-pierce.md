---
id: THM-771
title: Sheet endpoint defect and reduced-winding pierce — exact capacity/slack/overlap identity at seven exceptions, strict event tears, and the correct scale-free winding cutoff
status: PROVED (elementary exact grid geometry; corrects the false raw-winding and KCL claims in THM-767)
source: codex-2026-07-14-S5 (independent audit of the seven-exception sheet wall)
depends_on:
  - THM-761   # sheet construction and fibre-exact core clearance
  - LRCUpTo13 # only for the six-speed core margin in Corollary C
corrects:
  - THM-767   # its balance/event-pierce core survives; raw-w mesh, integral two-value display, and KCL necessity do not
related: [THM-760, THM-765, HYP-6830, HYP-6835, MISTAKE-146]
verification: 04-computation/lrc14_r7_sheet_endpoint_defect_codex_S5.py
  (+ 05-knowledge/results/lrc14_r7_sheet_endpoint_defect_codex_S5.out)
---

# THM-771 — Sheet endpoint defect and reduced-winding pierce

## Frame and exact owner-incidence defect

Let `c>=2`, let `P` be a nonempty finite set of positive integer speeds, and
let `W={w_1,...,w_7}` be seven distinct positive integers, none divisible by
`c`.  At target `delta=1/14`, put

```text
t_k=(x+k)/c,
B_a(x)={k in Z_c : ||w_a t_k||<1/14},
d_k(x)=#{a:k in B_a(x)},
F(x)=#{k:d_k(x)=0}.
```

The strict inequality in `B_a` is part of the statement.  A sheet outside
all `B_a` has the closed clearance `>=1/14` required by LRC.

For every owner `a`, define

```text
g_a=gcd(w_a,c),       C_a=c/g_a,       u_a=w_a/g_a,
A_a=g_a*ceil(C_a/7).
```

Thus `C_a` is the effective sheet-grid order, `u_a` is the reduced winding,
and `A_a` is the maximum possible number of bad sheets for that owner.  Define

```text
Q(x)=sum_a (A_a-|B_a(x)|)                 (capacity slack),
Omega(x)=sum_k max(d_k(x)-1,0)            (overlap debt),
sigma=sum_a A_a-c                         (ramification surplus).
```

Then `sigma>=0` and the following identity is exact:

```text
                     F(x)=Q(x)+Omega(x)-sigma.             (1)
```

In particular, this is not a probabilistic union bound.  It is an exact
owner-labelled Euler defect of the seven-row incidence deck.

## A. Exact count and endpoint lattice in the unramified layer

Suppose `7|C_a`, and write `C_a=7m_a`.  Define the endpoint-transfer set

```text
E_a={x in R/Z : u_a*x = C_a/14 (mod 1)}.                 (2)
```

Then

```text
|B_a(x)| = c/7-g_a,   if x in E_a,
           c/7,       if x notin E_a.                    (3)
```

The set `E_a` is one coset of the `u_a`-grid: it has exactly `u_a` points and
cyclic mesh

```text
                         1/u_a=g_a/w_a.                  (4)
```

It is not, in general, a `1/w_a` grid.  The two nominal `+/-1/14`
endpoint families coincide after projection because their difference is
`C_a/7=m_a`, an integer.

If all seven owners satisfy `7|C_a` (equivalently `7g_a|c`), then `A_a=c/7`
for every owner and `sigma=0`.  Formula (1) becomes

```text
                         F(x)=Q(x)+Omega(x).              (5)
```

Consequently a fully burned non-event fibre is necessarily an exact
owner-labelled partition: every owner has `c/7` sheets and no sheet has two
owners.  At an event `x in E_a`, however,

```text
                         F(x)>=g_a>=1.                    (6)
```

Thus every boundary switch tears the strict cover.  An exiting owner and an
entering owner both sit at equality at the switching moment; equality is safe,
so it cannot be treated as a KCL handoff of a bad sheet.

## B. Component-length event pierce

Let

```text
K(P)={x in R/Z : min_{p in P}||p x||>=1/14}
```

and retain the all-unramified hypothesis `7|C_a` for every owner.  If `K(P)`
contains a closed circular interval `J` with

```text
                         |J|>=1/u_a                       (7)
```

for some owner `a`, then `J` meets `E_a`.  At that event, (6) supplies a free
sheet `k`.  At `t=(x+k)/c`, every `cp` has clearance `||p x||>=1/14` and every
exception has clearance `>=1/14`.  Therefore

```text
                         M(cP union W)>=1/14.             (8)
```

More quantitatively, let `b=max(P)` and suppose `M(P)>1/14`.  A maximizer of
the core lower envelope has a closed `1/14`-safe neighbourhood of length at
least

```text
                         2*(M(P)-1/14)/b.                 (9)
```

Hence (8) follows whenever

```text
 max_a u_a >= b/(2*(M(P)-1/14)).                         (10)
```

For the LRC(14) seven-exception split, `|P|=6`; settled lower-dimensional LRC
gives `M(P)>=1/7`.  Thus the scale-free sufficient condition is

```text
             max_a w_a/gcd(w_a,c) >= 7*max(P).           (11)
```

This improves the component-count bound involving `sum(P)` and, crucially,
uses reduced winding.  Multiplying both `w_a` and its gcd sheet multiplicity
by the same factor cannot make endpoint events denser.

## C. The ramification obstruction is exactly surplus paid by overlap

Condition `7|c` alone does not imply exact tiling.  If `7` divides `g_a` to
the full depth present in `c`, then `7` need not divide `C_a`; such an owner
can have `A_a>c/7`.  Formula (1) says precisely how this extra capacity can
hide all free sheets: a full cover has

```text
                         Q(x)+Omega(x)=sigma.             (12)
```

The exact primitive example

```text
c=21,  P={1,2,3,4,5,6},
W={1,2,3,4,7,49,56},  x=1/7
```

has bad-set sizes `(3,3,3,3,7,7,7)` and

```text
Q=0,  Omega=12,  sigma=12,  F=0.
```

All sheets are strictly bad on a nontrivial open `x`-chamber, but the cover
has many overlaps and is not an exact tiling.  This is the concrete
counterexample to replacing `7|C_a for every a` by the weaker `7|c`.

At the opposite extreme,

```text
c=7,  P={1,...,6},  W={1,4,5,6,8,9,10},  x=1/7
```

is a static zero-defect exact partition with `Q=Omega=sigma=F=0`, persisting
on the chamber `|x-1/7|<1/140`.  Its owner `10` has zero `14*gcd` mirror
capacity, so a local exact-tiling chamber does not imply the KCL inequality
claimed in THM-767.  At the event `x=1/8` of owner `4`, sheet `5` becomes free
and gives the exact witness `t=41/56` with clearance `1/14`, exactly as (6)
predicts.

## Proof

As `k` runs over `Z_c`, the phases of owner `a` are the translated `C_a`-grid

```text
{(u_a*x+j)/C_a mod 1 : j in Z_{C_a}},
```

each point repeated `g_a` times.  An open arc of length `1/7` contains at most
`ceil(C_a/7)` grid points, proving the capacity formula.  Also

```text
sum_a |B_a| = sum_k d_k,
|union_a B_a| = sum_k d_k - sum_k max(d_k-1,0).
```

Substitute `|B_a|=A_a-(A_a-|B_a|)` and `sum A_a=c+sigma`; subtract the union
size from `c`.  This gives (1).

If `C_a=7m_a`, the open arc has length exactly `m_a` grid spacings.  Away from
endpoint alignment it contains `m_a` grid points.  At alignment, both
endpoints lie on the grid and are excluded, leaving `m_a-1`.  Endpoint
alignment is

```text
u_a*x+j = +/- C_a/14 (mod C_a).
```

The two signs have the same fractional part because their difference is
`C_a/7=m_a`.  This proves (2)--(4).  Under the all-unramified hypothesis,
`sigma=0`, so (5) and (6) follow from nonnegativity of `Q` and `Omega`.

A coset of the `u_a`-grid meets every closed interval of length `1/u_a`, which
proves B.  Finally, if `x_*` maximizes the core envelope, its `b`-Lipschitz
property gives

```text
min_p ||p(x_*+h)|| >= M(P)-b|h|.
```

Taking `|h|<=(M(P)-1/14)/b` proves (9)--(11).  Part C is direct exact
evaluation, reproduced by the verification artifact. ∎

## Owner-splice and Tournament Analysis

When `Omega=0`, label every occupied sheet by its unique owner.  Assume in this
paragraph that the cyclic word has at least one occupied and at least one free
sheet.  If `R` is the number of cyclic monochromatic occupied runs and `S` the
number of adjacencies where one occupied owner-run splices directly into a
different owner-run, then

```text
                 R-S = number of cyclic free-sheet runs.                  (13)
```

This is the component analogue of the cardinal defect (1).  It is useful for
owner-word classification, but (1) is stronger because it also sees sheet
multiplicity and ramification.  The full and all-free circular words are handled
separately; no artificial Euler convention is imposed on those edge cases.

The verification tournament uses owners as vertices.  Its pairwise observable
is oriented private-sheet adjacency; the switch compares the two directions,
and ties use the owner-order Hamiltonian path.  Score histograms, directed
cycles, SCCs, edge flips, and Hamiltonian-path counts are reported for the two
counterexample rows.  The challenged assumption is explicit: owners, runners,
and sheets are all lossy vertex choices.  The quotient preserves cyclic private
splices but destroys `Q`, `Omega`, and `sigma`; therefore the full labelled
owner-by-sheet incidence deck, not the tournament, preserves the free-sheet
predicate.

## Correction to THM-767

THM-767 was promoted from a claim stub with three errors.  Its balance identity
and its fixed-event pierce under `7g_a|c` are valid.  However:

1. at integral `C_a/7`, the two count values are `c/7` and `c/7-g_a`, not the
   displayed floor and floor-plus-one pair;
2. the endpoint mesh is `g_a/w_a=1/u_a`, so its raw condition
   `max w_a>7*sum(P)` does not follow and is not scale invariant;
3. strict endpoint equality is safe, so coincident exit/entry boundaries do
   not furnish a bad-sheet KCL handoff.  A local static tiling chamber can also
   avoid every event and violate the proposed KCL capacity.

Use THM-771, not the original THM-767 claim, for the seven-exception endpoint
reduction.
