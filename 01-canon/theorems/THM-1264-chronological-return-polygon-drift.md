---
id: THM-1264
title: CHRONOLOGICAL RETURN POLYGONS PAY EXACT ADDRESS DRIFT -- a literal owner triangle forces a factor-three speed ascent
status: PROVED (all-length consecutive return-subword overlap identity; address-drift ratio sum; minimum-owner ascent; simple-return indexing and depth; six-fast factor-6/5 and triangle factor-3 compact heights; seam-digit triangle interpretation; optimization-safe exact interval census; sorry-free Lean arithmetic core). The theorem requires one literal consecutive return subword with shared intermediate tooth occurrences, not merely support edges appearing elsewhere
source: codex-2026-07-19 placed-Fano continuation
depends_on: [THM-1253, THM-1256, THM-1260, THM-1262]
related: [THM-841, THM-1156, THM-1252, THM-1254, THM-1266, THM-1274]
script: 04-computation/lrc14_chronological_return_polygon_thm1264.py
output: 05-knowledge/results/lrc14_chronological_return_polygon_thm1264.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCChronologicalReturnPolygon.lean
script_sha256: 46d6de8819556e0720c33fa27412d572bdc670c1743525d417b2bbfa97141d34
output_sha256: d5fc0f5e2e8e1e7fe9c0f6e96655812d452d0c5d359e3651c40a0a59d56757ba
formalization_sha256: 2a49505403d0d66ba6b911633b07015ce2dd1c36f0b16b5cd6688bd6b4b423b6
---

# THM-1264 -- chronological return-polygon drift

THM-1260 proves that one fully placed sharp fork and its next binary-side
blocker carry no speed-`chi_7` line law.  THM-1262 forces a coherent blocker
two-cycle to open through a third-owner corridor.  This theorem identifies
the first object on which a genuine metric Fano consumer exists: not three
owner-pair edges in an abstract support graph, but one **literal consecutive
owner-return polygon** in the chronological tooth word.

## 1. The literal return object

Work in THM-1253's deletion-minimal chronological cover by individual open
teeth.  Take one consecutive subword

```text
I_0,I_1,...,I_m,                 m>=2,                (1)
```

where

```text
I_q=((14n_q-1)/(14s_q),(14n_q+1)/(14s_q)).           (2)
```

Both endpoint sequences are strictly increasing.  Consecutive teeth have
the positive raw overlaps

```text
omega_q=R(I_q)-L(I_(q+1))>0,        0<=q<m.          (3)
```

Assume the subword returns to its initial owner:

```text
s_0=s_m=a,                     n_m=n_0+R,             (4)
R in Z_(>=1).                                          (5)
```

The return is positive because the two `a`-teeth are distinct and occur in
increasing chronological order.  The phrase **literal consecutive** is an
essential hypothesis.  The tooth `(s_q,n_q)` in the `q-1` handoff is the
same occurrence, with the same absolute address, as the tooth in the `q`
handoff.  Choosing the owner edges independently at separated word positions
does not provide (1)--(5).

## 2. Exact all-length polygon identity

For each edge, endpoint subtraction gives

```text
omega_q=n_q/s_q-n_(q+1)/s_(q+1)
        +1/(14s_q)+1/(14s_(q+1)).                    (6)
```

Sum (6) over `q=0,...,m-1`.  The centres telescope to

```text
n_0/a-n_m/a=-R/a.                                    (7)
```

Every internal half-width occurs twice.  The two endpoint half-widths also
combine because `s_m=s_0=a`.  Therefore

```text
sum_(q=0)^(m-1) omega_q
  =(1/7)sum_(q=0)^(m-1) 1/s_q-R/a.                  (8)
```

This is an equality of the actual raw overlaps already charged by
THM-1253, not a second mass invoice.  It identifies their missing functional
coordinate: total chronological overlap is the reciprocal-speed perimeter
minus the integral address holonomy.

The same proof works for any radius `1/(2N)`, with `1/7` in (8) replaced by
`1/N`; the radius-`1/14` form is kept here because its integer return cost is
exactly the LRC(14) constant used below.

## 3. Ratio ascent and the minimum internal owner

The left side of (8) is strictly positive.  Multiplying by `7a` gives

```text
sum_(q=1)^(m-1) a/s_q>7R-1.                          (9)
```

Put

```text
s_min=min_(1<=q<m) s_q.                              (10)
```

Each summand in (9) is at most `a/s_min`, so

```text
a>((7R-1)/(m-1))s_min.                               (11)
```

This is the exact return-polygon ascent.  Large address return makes the
factor stronger.  Long polygons spread the address invoice among more
internal owners, which is why word length cannot be discarded.

For `m=2`, (8) is THM-1252/1256's backtrack law in summed form:

```text
a>(7R-1)s_1>=6s_1.                                   (12)
```

Thus the old toothpick return is the two-edge member of one all-length
polygon hierarchy.

## 4. Simple returns: exact indexing, depth, and breadth

Call (1) a **simple owner return** when

```text
s_1,...,s_(m-1) are distinct and none equals a.       (13)
```

With at most `K` available owner labels, (13) has at most `K-1` internal
positions.  Hence

```text
m-1<=K-1,                         m<=K.               (14)
```

The indexing is by **edges**: a simple word with `K` distinct owners has at
most `K+1` tooth occurrences but exactly `K` edges, because the last
occurrence repeats the first owner.

For a seven-owner ambient alphabet, (11), `R>=1`, and `m<=7` give only

```text
a>s_min.                                              (15)
```

This is a strict ascent but there is no uniform multiplicative constant
larger than one, so the projective bound `d<2345c` alone gives no scale-only
height.  On a **fixed** seven-owner packet, however, orient each simple
return from a minimizing internal owner to its repeated outer owner.  The
speed strictly increases, so this return DAG has depth at most six.

In the actual slow-gap tooth word only the six fast labels occur.  Then

```text
K=6,                  m<=6,
a>(6/5)s_min.                                         (16)
```

If one deliberately ignores the stronger fixed-label DAG, the compact ratio
box gives the exact scale-only bound

```text
(6/5)^42<2345<(6/5)^43.                              (17)
```

Thus any **transported chain** of factor-`6/5` returns has height at most
`42`.  On the fixed six-speed packet, strict speed ascent already gives a
potential DAG depth at most five.

These are **depth** statements, not breadth statements.  Many return polygons
may point into the same larger owner, and a star can have arbitrarily many
tooth occurrences as scale grows.  No return-DAG argument removes that
breadth.  THM-1253's full occurrence invoice is what pays every branch and
every repeated seam.

If a recurrent owner word is not simple, choose a pair of equal owners at
minimum positive separation.  Its intervening owners are distinct and avoid
the endpoint owner, so the resulting subword is simple.  In the six-fast word
this **closest-return corollary** has `2<=m<=6`, forces the strict factor `6/5`,
and carries the scale rank `42` from (17).

This is a rank, not a recursion-closure theorem.  Nothing in (8) says that a
blocker bridge emitted at one return lands on a child return, so the potential
DAG edges need not concatenate along the actual proof operation.  THM-1266's
exact primitive five-ray star makes the obstruction sharp: one high owner can
return separately across each of five lower owners, giving height one and
breadth five.  All five rays are already paid occurrences in THM-1253's full
invoice.  Consequently (14) is the correct local return cell, but transport
to its child stalk is an additional theorem, not a consequence of owner-word
contraction.

## 5. The literal owner triangle: the first metric Fano line

Specialize to three edges:

```text
a -> b -> c -> a,                  R>=1.              (18)
```

Equations (8)--(11) become

```text
omega_ab+omega_bc+omega_ca
 =1/(7a)+1/(7b)+1/(7c)-R/a,                          (19)

a/b+a/c>7R-1>=6,                                     (20)

a>3 min(b,c).                                        (21)
```

Orient (18) from its smaller internal owner to `a`.  It is a strict ratio
ascent by more than three.  In the compact THM-1233 box,

```text
3^7=2187<2345<6561=3^8,                              (22)
```

so a chain of such triangle-return ascents has at most seven edges by scale;
on one fixed six-speed packet the owner-DAG depth at most five remains
stronger.

This is the first honest Fano-line consumer after THM-1260.  For a handoff
from `(u,M)` to `(v,N)`, define its positive numerator and normalized seam
digit by

```text
W=v(14M+1)-u(14N-1),
g=gcd(u,v),                  w=W/g.                   (23)
```

Then

```text
w==(u/g+v/g) mod 14.                                  (24)
```

The three edges in (18) are three **actual occurrence digits**, and (19)
couples their metric values to the common address return `R`.  THM-1260
shows why their endpoint `chi_7` colours alone cannot do this: its sharp
positive fibre `w=1` already realizes all four ordered colour pairs.

A triangle in the unordered owner-transition support is not enough.  If its
three edges occur at separated places, their intermediate addresses need not
match, there is no single `R`, and the telescope (7) is unavailable.  Such a
support triangle is telemetry, not an instance of (18).

## 6. Relation to the placed blocker frontier

THM-1262 says that an aligned blocker two-cycle exits its ascent target
through a protected third-owner bridge.  That produces a literal initial
path in the common word.  THM-1264 supplies the following exact consumer if
subsequent consecutive handoffs return to one of its earlier owners:

1. choose the first repeated owner along that literal continuation;
2. the minimum-distance return subword is simple;
3. apply (11), or (21) when it has exactly two internal owners;
4. orient the result as a strict speed ascent in the bounded return DAG.

What is still missing is a theorem that the exported bridge must return
before the word terminates, or else an endpoint/private-stalk invoice that
pays the nonreturning continuation.  THM-1264 proves the complete return
branch.  It does not manufacture the return from three abstract support
edges and does not prove that consecutive closest-return ranks compose.
THM-1266 records the exact five-ray obstruction to that composition.  No
remaining six-comb cover is closed here.

There is a second occurrence-level caution.  A repeated wall-owner **label**
does not automatically supply (1)--(5): both wall incidences may be carried
by the same selected tooth.  THM-1264 applies only after two distinct tooth
instances are located.  THM-1274 supplies exactly this split on the
protrusion-facing continuation of a slowest two-cycle.  Distinct instances
give a closest return with at most five edges and factor greater than `3/2`;
without a repeat the continuation reaches the carrier endpoint within five
tooth occurrences.  The latter endpoint obligation, and same-instance
multi-incidence, remain open consumers.

This also reroutes the Fano probe.  The right vertices are not seven runner
points but handoff occurrences with shared tooth addresses.  A useful line
is a closed chronological transport cell carrying `(W,g,w,R)`.  Its law is
the metric holonomy (19), not a point-colour parity.

## 7. Tournament and alternate-carrier audit

Chronological tooth positions form a transitive tournament, with the word
order as its unique Hamiltonian path.  The owner projection of (1) is a
directed return polygon.  Neither quotient alone proves (8): position order
forgets repeated ownership, while the owner polygon forgets occurrence
matching and absolute addresses.

Orient each literal return obligation from a minimizing internal owner to
the repeated outer owner.  The pairwise observable is the exact factor in
(11), with speed order as the tie Hamiltonian path.  On simple returns this
is an acyclic ascent relation.  Its SCCs are singletons and it has no directed
cycles; breadth is recorded by parallel occurrence obligations rather than
collapsed to one tournament edge.

We challenged runners, owner pairs, individual teeth, tooth boundaries,
raw overlap cells, seam digits, support triangles, Fano lines, return
polygons, and proof obligations as vertices.  The smallest faithful carrier
for this theorem is

```text
(one consecutive tooth subword;
 every shared occurrence (s_q,n_q);
 its raw overlaps and seam digits;
 repeated endpoint owner a and integral address return R).              (25)
```

It preserves exactly the telescoping predicate.  Projecting to any unordered
support graph destroys it.

## 8. Exact replay, formalization, and scope

The optimization-safe dependency-free referee constructs an exact bank of
`335` teeth and enumerates every literal minimal return chain of two through
five edges in that bank:

```text
edges       2      3       4       5
returns    42      8    1,310   8,601.                (26)
```

Across the `48,353` constituent seam occurrences it checks endpoint
monotonicity, THM-1253 two-step separation, positive overlap, the gcd sheet,
the seam-digit congruence (24), identity (8), inequalities (9)--(11), and all
eight literal triangle rows.  The bank contains returns through address jump
four.  It separately checks the edge indexing, factor-`6/5` height `42`,
fixed-owner depths, and factor-three height `7`.

Every proof-critical Python check raises an explicit `RuntimeError`; the
referee parses its own AST and rejects any Python `assert` node.  Normal and
optimized runs are byte-identical.

The sorry-free Lean module proves the abstract interval-chain telescope, the
all-length tooth formula (8), the positive ratio sum (9), the triangle
factor-three consequence (21), and both exact compact-height comparisons.
Extraction of the deletion-minimal word and the fact that a selected
consecutive pair has raw overlap (3) remain the paper topology provider from
THM-1253; there are no proof placeholders or `native_decide` calls.

The frozen artifact hashes are

```text
script         46d6de8819556e0720c33fa27412d572bdc670c1743525d417b2bbfa97141d34
output         d5fc0f5e2e8e1e7fe9c0f6e96655812d452d0c5d359e3651c40a0a59d56757ba
formalization  2a49505403d0d66ba6b911633b07015ce2dd1c36f0b16b5cd6688bd6b4b423b6
```

THM-1264 proves no new point-colour law and does not prove LRC(14).  It closes
the literal chronological return branch with an exact functional identity
and turns the first genuine three-owner line into a bounded metric ascent.
The remaining operation-level target is return-or-endpoint debt plus an
honest bridge-to-child transport for the third-owner bridges exported by
THM-1262.
