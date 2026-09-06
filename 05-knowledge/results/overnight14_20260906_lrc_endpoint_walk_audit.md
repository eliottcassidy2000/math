# Independent audit of the endpoint-gcd walk bridge

**Verdict: PASS.** The exact ratio identity and conditional LRC endpoint
criterion in
[the producer report](overnight14_20260906_lrc_endpoint_walk.md)
are valid. No mathematical repair is requested. A qualifying walk is not
proved to exist in every remaining core; LRC(14) remains open.

## 1. Exact identity and divisibility

Let d=gcd(w₀,…,w_m), D=gcd(w₀,w_m), and let L be the lcm of the reduced
denominators of w_i/w₀. Directly,

    L=lcm_i(w₀/gcd(w₀,w_i))=w₀/d.

The first equality is the definition of reduced denominator; the second
holds prime by prime. This independently establishes the producer's clearing
factor without its accumulated-prefix argument.

For the reduced oriented edge ratios w_i/w_(i+1)=a_i/b_i, set A=∏a_i,
B=∏b_i, C=gcd(A,B). Then w_m/w₀=B/A has reduced denominator A/C.
Consequently D=w₀/(A/C), and

    D/d=L/(A/C)=J.

The endpoint denominator divides L, so J is an integer. Every prefix
denominator divides the corresponding product of a_i, which divides A;
therefore L|A. It follows that C/J=A/L is an integer, proving J|C. The
orientation is correct: the endpoint denominator is A/C, not B/C.

The formulas hold for a zero-edge walk, repeated vertices, and returns to
the starting vertex. J is exactly endpoint gcd divided by collective gcd.
For fixed endpoints and visited set it is independent of the particular
walk order, while C can grow under repeated excursions. For example,
6→2→3 has C=J=3, whereas6→2→6→2→3 has C=9 but J=3. Thus the weaker
bound using C is valid but can lose information solely through repetition.

## 2. Endpoint and physical-scale qualifications

For the LRC consumer, the walk must start at the physical global core
maximum tK and end at a **different** core coordinate. Its r is the number
of distinct visited labels, not the number of entries or edges. Repetitions
do not change the gcd of that subset. The inherited cap M_r is a cap on
the physical subset gcd d, while the twelfth endpoint test concerns

    gcd(K,w_m/t)=D/t=dJ/t≤M_rJ/t.

This exactly justifies J≤floor(tH/M_r), including equality, with
H=76,388,115. The ratios, C, L, and J are unchanged by common physical
scaling; the endpoint gcd must still be divided by t. The six thresholds
are correctly computed:

| Physical scale | r=5 | r=6 | r=7 |
|---|---:|---:|---:|
| t=1 | 6 | 2,390 | 848,756 |
| t=2 | 13 | 4,781 | 1,697,513 |

The physical caps and t≤2 are inherited restrictions on a hypothetical
strict failure. The conclusion is therefore obtained by contradiction:
a qualifying actual-entry walk forces an endpoint pair satisfying the
audited twelfth closure theorem. It does not assert that every safe row
must obey those caps. The actual graph, global maximum, physical box, and
W=V_dec assumptions of that supplier remain necessary inputs. No five-label
Bézout identity is asserted to have support three.

The path6→2→3 is an actual atlas hostile: its reduced sums are4 and5,
its collective gcd is1, and its endpoint gcd is3. At least three distinct
vertices are required for this failure **when the endpoints are distinct**:
with only two distinct vertices the collective gcd already equals the
endpoint gcd. Without that endpoint condition,6→2→6 has two distinct
vertices and J=3, but supplies no distinct endpoint pair for the closure
criterion. This is a scope clarification, not a defect in the stated
consumer, which already requires another core endpoint.

The monotone seven-vertex path729→243→81→27→9→3→1 has J=C=1. It is
a positive packet control, not a complete thirteen-speed entry or an
existence proof for all residual cores. Bounding J on some qualifying
maximum-started walk is the remaining structural obligation.

The final producer adds two distinct-vertex strict-gain controls. The path
6→4→1 has J=1<C=2. More substantially,18→12→3→9→6 has J=2≤6<C=12:
the exact five-subset test passes while replacing J by C would fail. The
referee independently reconstructs the actual11+2 graph on the declared
core(1,2,3,4,5,6,7,9,10,12,18) and pair(g,3g), g=36Q+1, and verifies the
physical box and full crossing dominance. This is a genuine safe equality
control with the path starting at the global core maximum. It demonstrates
a strict gain in the retained packet; it is not a new unsafe family or an
example excluded by all earlier unit-component closure criteria.

## 3. Independent reproduction and frozen inputs

[Independent source](../../04-computation/overnight14_20260906_lrc_endpoint_walk_audit.py)
and [output](overnight14_20260906_lrc_endpoint_walk_audit.out):

```
python -B 04-computation/overnight14_20260906_lrc_endpoint_walk_audit.py
python -B -O 04-computation/overnight14_20260906_lrc_endpoint_walk_audit.py
```

Both runs pass **117,721 always-active gates**, with byte-identical LF
output. No producer routine is imported. Direct rational vertex ratios
provide the denominators independently of the edge products. The complete
banks are5,832 three-entry words over1,…,18 and all13,542 actual atlas
walks on vertices1,…,14 with zero through four edges, with counts
14,70,382,2,040,11,036. Additional controls check common scales, distinct
visited-set cardinality, returned endpoints, cancellation inflation, the
seven-vertex positive example, and each threshold plus its first rejected J.

The producer's7,776-word bank consists of **five-entry** walks, since
6⁵=7,776. The final source/report explicitly state that convention and
qualify the minimal hostile by its distinct endpoints, as requested during
the audit. These are scope clarifications, not changes to the identity.

Pinned producer source:
`e38ce143de236a35d089f0e7d4110fae5200ea819868e7b925a705f50428664a`.
Pinned producer output:
`4178d510f86a4cdb4cca86afcb6405439e5179cddd924e48fffd1654a1729561`.
Referee source:
`7210275c77bd15af84fe6a6a5369ae0c4be3377563cff3bfaf2f69a71e622067`.
Referee normal/optimized output:
`89bcba554f4bdd58bd94f8bd4284aa17407f806e958274f423aa5d6650afeb4b`.
Input lookup supports the eventual separated source/results filing layout.
No Git, navigation, producer code, or prior frozen artifact was modified.

**Filing:** root read these proofs and audits and integrated the fourteenth
checkpoint. Reproduction commands are relative to the repository root.
