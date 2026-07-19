---
id: THM-1193
title: Every covered slow gap carries a canonical sum-beat killer cycle
status: PROVED pointwise necessary structure.  Each fast comb supplies an interior sum-beat point at which it and the slow carrier have distance greater than 1/4; a six-comb cover therefore induces a fixed-point-free killer map and hence a directed cycle.  The cycle is not sufficient for coverage: an exact guardrail packet realizes all six canonical obligations while leaving survivor length 879/10780.  No slow-gap noncoverage theorem or LRC(14) closure is claimed
source: codex-2026-07-18-S75/S77 canonical beat-cycle synthesis
depends_on: [THM-1170, THM-1176, THM-1192]
related: [THM-1156, THM-1175, THM-1179, HYP-7678]
script: 04-computation/lrc14_canonical_sum_beat_killer_cycle_codex_20260718.py
output: 05-knowledge/results/lrc14_canonical_sum_beat_killer_cycle_codex_20260718.out
---

# THM-1193 -- the canonical sum-beat killer cycle

Work at radius `1/14`.  For a positive integer `s`, put

```text
D_s={t in R : ||st||<1/14}.
```

Let `a<b` and let

```text
G_k=[(14k+1)/(14a),(14k+13)/(14a)]                  (1)
```

be the `k`th closed slow gap of `D_a`.  Its centre and half-length are

```text
c_k=(2k+1)/(2a),               3/(7a).              (2)
```

## 1. A canonical sum beat in every slow gap

Put `q=a+b`.  Choose the integer `p` nearest to `q c_k`, resolving a tie
toward the larger integer, and define

```text
t_b=p/q.                                               (3)
```

> **Lemma 1 (interior sum-beat witness).**  The point `t_b` lies in the
> interior of `G_k` and
>
> ```text
> ||a t_b||=||b t_b||>=b/(2(a+b))>1/4>1/14.           (4)
> ```

**Proof.**  Nearest-integer rounding gives

```text
|p-qc_k|<=1/2,
|t_b-c_k|<=1/(2q).                                    (5)
```

Because `b>a`, one has `q>2a`, and therefore

```text
1/(2q)<1/(4a)<3/(7a).                                 (6)
```

Equations (2), (5), and (6) put `t_b` strictly inside (1).  Moreover
`a c_k=(2k+1)/2` is a half-integer, while
`|a(t_b-c_k)|<=a/(2q)<1/4`.  Distance to the integer lattice is consequently

```text
||a t_b||=1/2-|a(t_b-c_k)|
          >=1/2-a/(2q)=b/(2q)>1/4.                   (7)
```

Finally `q=a+b` and `q t_b=p` give

```text
b t_b=p-a t_b,
```

so the two lattice distances agree.  This proves (4). ∎

The construction is canonical after the rounding gauge is fixed.  It is a
particularly strong member of THM-1192's sum-beat puncture set: neither member
of the defining pair is merely endpoint-safe; both have margin greater than
`1/4`.

There is also a useful noncanonical extension.  For every integer `r` with
`r/q in G_k`, the sum-beat identity gives

```text
||a r/q||=||b r/q||>=1/14.                            (8)
```

Thus every such lattice point must be punctured by a third comb under a cover.
Lemma 1 proves that this finite lattice set is never empty.

## 2. Coverage forces a functional cycle

Now let

```text
a<b_1<...<b_6                                           (9)
```

and suppose the six fast danger combs cover the slow gap:

```text
G_k subset union_(j=1)^6 D_(b_j).                     (10)
```

For each `i`, form the canonical point `t_i` from Lemma 1 using `b=b_i`.
By (4), neither `a` nor `b_i` is dangerous at `t_i`.  Hypothesis (10) therefore
forces some label `j!=i` with

```text
||b_j t_i||<1/14.                                     (11)
```

Choose the least such label and call it `f(i)`.  This gives a map

```text
f:{1,...,6}->{1,...,6},             f(i)!=i.          (12)
```

> **Theorem 2 (killer cycle).**  Every cover (10) induces a directed cycle
>
> ```text
> i_1 -> i_2 -> ... -> i_r -> i_1,       2<=r<=6,     (13)
> ```
>
> where the edge `i->j` is the exact beat-puncture inequality (11).

**Proof.**  Iterate `f` from any label.  Among seven successive labels in a
six-element set, one repeats.  The labels from its first occurrence to the
next form a directed cycle.  Since (12) forbids a fixed point, its length is
at least two. ∎

Writing `q_i=a+b_i` and `p_i=q_i t_i`, every edge has the integer form

```text
there exists n_i in Z such that
|b_j p_i-n_i q_i|<q_i/14.                             (14)
```

The cycle is therefore a finite simultaneous Diophantine obstruction attached
to one slow gap.  Unlike the transitive speed-order tournament, it is produced
by the cover predicate itself.

## 3. Exact guardrail: the cycle is necessary, not sufficient

Take

```text
a=3,   k=1,   (b_1,...,b_6)=(4,5,6,7,9,11).           (15)
```

The canonical `(p_i,q_i)` are

```text
(4,7),(4,8),(5,9),(5,10),(6,12),(7,14).               (16)
```

Every row has a third-comb killer.  With zero-based labels their exact killer
sets are

```text
{3}, {0,2}, {4}, {0,2}, {0,2}, {0,2}.                 (17)
```

The least-label map contains the two-cycle `0->3->0`.  Nevertheless exact
closed-tooth merging leaves survivor length

```text
879/10780>0.                                           (18)
```

Closing teeth only helps coverage, so (15) genuinely does not cover the slow
gap.  Thus Theorem 2 is a compulsory finite skeleton, not a replacement for
the full interval or all-pairs beat obligations of THM-1192.

## 4. Tournament and carrier audit

The proof-bearing directed relation is not a tournament: a beat obligation
can have several killers, and two labels can kill one another.  For comparison,
define the directed clearance

```text
kappa_ij=1/14-||b_j t_i||                              (19)
```

and the antisymmetric observable

```text
O(i,j)=kappa_ij-kappa_ji.                              (20)
```

Orient `i->j` when `O(i,j)>0`; break equality by increasing label.  On the
guardrail (15), this quotient tournament has score multiset

```text
1,2,2,2,3,5,
```

four directed triangles, SCC sizes `5,1`, and fifteen Hamiltonian paths.  The
tie Hamiltonian path used by the exact replay is `(1,0,3,2,4,5)`.  This is a
genuinely nontransitive fingerprint, unlike speed order, but it still does not
preserve coverage: (17) holds while (18) remains positive.

The challenged assumption is again that vertices should be runners.  The
faithful vertices here are the six **sum-beat proof obligations** `(p_i,q_i)`;
edges are labelled killer incidences (14), with slow-gap address `k` and comb
modulus retained as sidecars.  Passing to (20) destroys multiple killers,
absolute margins, the other integer points in (8), and the interval between
beats.  The right object is therefore a beat-obligation digraph/hypergraph,
not a naked tournament.

## Honest frontier

The theorem turns slow-gap coverage into a mandatory directed cycle of exact
modular inequalities and supplies an interior `>1/4` beat in every row.  It
does not rule out those cycles.  THM-1192 shows that retaining all pair beats,
instead of only the six canonical carrier beats, is computationally much more
selective.  The open step is an all-range contradiction from the combined
cycle, mixed-gcd, and harmonic-pressure data; LRC(14) remains open.

The exact replay verifies Lemma 1 on 104,005 phase-labelled rows, reproduces
the guardrail, and audits every one of the 720 vertex orders in its tournament
fingerprint.  Normal and optimized runs are byte-identical to the stored
output.  Frozen SHA-256 hashes are

```text
source  4d7940a6fddefbf7baf7fe144a3d9eba9ef4b49f0d3e2d10efc83a2e6b3d3cdf
output  1836b2d9e52141a213fa98e3c60c2346d200c49249cfc68b6615cf77cad1ba16
```
