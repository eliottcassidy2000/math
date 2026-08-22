# The nine relation slots are six `K4` edges plus three XOR axes

**Status: research synthesis around a FINITE-EXACT unnumbered sidecar; not a
truth source.**  The exact census is in
`04-computation/lrc_relation_k4_xor_star_triangle_probe_20260815.py`.

## 1. The four-vertex object hiding in the relation

The THM-3479 relation has coefficient packet

```text
(c1,c2,c3,H,q1,q2,q3,q4,q5)
 -> (-27,-27,-27,20110798,-41,-27,-27,-27,38).       (1)
```

Exactly six slots have coefficient `-27`.  A complete graph on four vertices
has exactly six edges and exactly three perfect matchings.  Label its vertices
by the XOR group

```text
V4=F_2^2={00,01,10,11}.                              (2)
```

The three nonzero differences partition the six edges into three opposite
pairs.  One explicit dictionary is

```text
difference 01: (c1,q2),
difference 10: (c2,q3),
difference 11: (c3,q4),                              (3)
```

where the `c` edge is incident with `00` and the paired `q` edge is opposite
to it.  Thus `{c1,c2,c3}` is the star at `00`, while `{q2,q3,q4}` is its
complementary triangle.  The three remaining slots `{H,q1,q5}` label the
three matching/difference axes in some declared order.

The axis order in (3) is a chart, not a theorem.  The conclusions below that
matter—V4 symmetry and star/triangle parity—are unchanged by permuting the
three axes.

## 2. The relation stabilizer is exactly V4

`S4` acts on the six edges and on the three perfect matchings.  The edge
coefficients in (1) are constant, while the three axis coefficients

```text
20110798, -41, 38                                    (4)
```

are pairwise distinct.  Therefore an `S4` permutation preserves (1) exactly
when it fixes all three matchings.  Exact enumeration gives

```text
|S4|=24,
|Stab(relation)|=4,
Stab(relation)=V4                                    (5)
```

where `V4` acts by translating the four vertex labels in (2).  The six cosets
are the six permutations of the three axis colours.

This makes the XOR appearance intrinsic to the declared four-vertex model:
it is the kernel of `S4 -> S3` on the three perfect matchings, not a numerical
analogy with the number four.

## 3. Three exact word-current characters

For each opposite pair write

```text
s_delta=U(c_delta)+U(q_delta),
d_delta=U(c_delta)-U(q_delta).                       (6)
```

Translation preserves every `s_delta`.  It changes the three differences by
the nontrivial V4 character table

```text
             00  01  10  11
chi_01        +   +   -   -
chi_10        +   -   +   -
chi_11        +   -   -   +.                        (7)
```

Consequently the nine relation coordinates split, under this V4 action, as

```text
six trivial coordinates
  = three pair sums + H + q1 + q5,
plus one copy of each of the three nontrivial characters. (8)
```

The triple `(d_01,d_10,d_11)` is therefore an explicit XOR-character current
packet.  It is not graph-cycle `H^1`, an endpoint phase, or a physical LRC
current; it is a representation-theoretic sidecar on the relation
coordinates.

## 4. Star versus triangle is one Boolean parity

When no opposite pair ties, choose the `c` edge if `d_delta>0` and the `q`
edge if `d_delta<0`.  One edge is then selected from each perfect matching.
There are eight such transversals:

```text
four stars:     odd number of c edges,
four triangles: even number of c edges.              (9)
```

Every nonzero V4 translation flips exactly two choices, so the parity in (9)
is invariant.  Equivalently,

```text
d_01*d_10*d_11 > 0  <=> STAR,
d_01*d_10*d_11 < 0  <=> TRIANGLE.                   (10)
```

This is a Boolean/XOR orbit invariant, not a forced tournament: it selects a
three-edge subgraph of `K4` but does not orient all six edges.

## 5. The two transplants occupy opposite Boolean types

For the canonical pairing (3), the exact packets are

```text
U_full:  (-14, 2066, 742533)          -> TRIANGLE,
U_clock: (-655166, -656336, 81141)    -> STAR,
U_q27:   (25648, 7598382, 18279864576)-> STAR,
U_q51:   (65779, 7191885, 30105545898)-> STAR.       (11)
```

More strongly, all `3!=6` bijections between the c-triple and q-triple give
the same type for each row:

```text
U_full:  6/6 TRIANGLE,
U_clock: 6/6 STAR,
U_q27:   6/6 STAR,
U_q51:   6/6 STAR.                                  (12)
```

The robustness has an elementary order explanation.  For `U_full`, `c1` is
below every q-value while `c2,c3` are above every q-value.  For `U_clock`, two
c-values are below every q-value and one is above.  In both CRT lifts every
c-value is above every q-value.

The CRT lifts are descendants of the clock construction, not two independent
experiments.  Thus (12) is a sharp classifier of the known branches and a
possible explanation of the two-transplant split, but it does not prove that
triangle type causes complete target nonvanishing or star type causes a
common-centre realization.

An interlacing hostile shows the classification is not automatic.  With

```text
C=(0,4,8), Q=(2,6,10),                               (13)
```

the six pairings produce both stars and triangles.

## 6. The mod-13 current and the Archimedean bit are separate

In every canonical packet of (11), the three differences are diagonal and
nonzero modulo 13:

```text
U_full,U_q27,U_q51: (12,12,12),
U_clock:            (8,8,8).                         (14)
```

So the three nontrivial V4 characters carry one common 13-adic amplitude,
while their integer signs distinguish triangle from star.  The 13-adic
diagonal does not determine the Archimedean Boolean type: `U_full` and the CRT
lifts all have residue `(12,12,12)` but lie in opposite sign orbits.

This gives a concrete two-sidecar grammar:

```text
13-adic amplitude + V4 sign parity.                  (15)
```

It may be useful when the endpoint-factorization lane tries to transport the
relation coordinates to actual edge responses.

## 7. Relation to the Jacobian three-class discussion

The three nonzero V4 characters in (7) also index the conventional three
quadratic characters of a V4 resolvent.  That supplies an index dictionary,
not a cohomology map.  On the fixed Jacobian map the three cubic discriminants
all collapse to the same Kummer class `[-L]`; they occupy the diagonal in a
three-view class packet.  Here the values `d_delta` are rational integers
transforming in three distinct real character lines.

Identifying them would still confuse ordinary finite-group representation
data with étale `H^1`.  A lawful cross-frontier map would have to send the
three relation characters to explicit boundary valuations or Kummer cocycles
and prove compatibility with the V4 action.  None is constructed here.

## 8. Cheapest decisive next test

The useful hypothesis prompted by (12) is now finite and falsifiable:

```text
within lawful THM-3479 relation tuples,
triangle type predicts complete unrestricted target support,
star type predicts same-clock/common-centre realizability. (16)
```

The cheapest test is to enumerate additional exact relation tuples satisfying
the literal q=11 masks and positive delayed mass, cross-tabulate (10) against
the endpoint-bank and common-clock predicates, and include interlacing rows.
One counterexample kills (16); a positive census would identify which extra
phase or ancestry coordinate must be added before any physical theorem.

Grouped address noncancellation, all-unit projectors, ancestry, bispectrum,
scalar-row exclusion, and LRC(14) remain open.
