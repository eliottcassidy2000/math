---
id: THM-2642
title: "Cyclic difference-relation saturation and thick-holotopy no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a
  prime cyclic torsor V=F_p and a directed n-cycle whose ith
  translation-equivariant edge allows the nonempty increment set S_i, the
  number of c-clutched sections is exactly
  p(1_{S_1}*...*1_{S_n})(c).  Iterated Cauchy--Davenport makes every clutch
  class occur once sum_i(|S_i|-1)>=p-1.  For two edges the exact universal
  lower bound is p max(0,|S|+|T|-p), and it is sharp.  Thus two eleven-sheet
  C_13 relations have at least 117 sections for every carry; the complete
  exact census is 117:55,770, 130:22,308, 143:1,014 over all relation-pair
  and carry triples.  This is the relation-valued dual boundary of the
  fixed-branch holotopy principle: thick equivariant relations can service
  every affine clutch while their support alone cannot decode a private
  carry.  No LRC transition or row conclusion follows, because the proved
  eleven-sheet rows are static/coefficient fibres rather than a same-base
  positive cyclic relation system.
source: deep-energy-audit-2026-07-28-cyclic-relation-holotopy
depends_on: []
related:
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2637-derangement-character-fixed-branch-holotopy-principle
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
script: 04-computation/lrc14_cyclic_difference_relation_holotopy_thm2642.py
output: 05-knowledge/results/lrc14_cyclic_difference_relation_holotopy_thm2642.out
script_sha256: ea885d7d9b066f4b1dc6f5d3ad16d4337e75ac5a69984c1e55b618d4028a53db
output_sha256: 7dd3af62f5fd61c495f9cd6d08dea3d1b643de87ad2a148c282802e69c45bf8f
hash_basis: LF-normalized bytes
---

# THM-2642 -- thick relations saturate cyclic clutch classes

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2637 treats a thin local system: every edge is one permutation, so a
cycle has one holonomy element and a fixed branch constrains its character.
The natural dual object is a cycle of translation-equivariant **relations**.
It behaves in the opposite way.  Composition becomes convolution, and enough
local thickness forces every affine clutch class to occur.  The result below
isolates that saturation boundary and prevents a dense root-support atlas
from being mistaken for a private carry decoder.

## 1. Difference relations and clutched sections

Let

```text
V=F_p                                                     (1)
```

for a prime `p`.  Give a directed `n`-cycle vertices `0,...,n-1`.  On its
`i`th edge choose a nonempty allowed-increment set `S_i subset V`, and define
the translation-equivariant relation

```text
R_i={(x,x+s):x in V, s in S_i}.                           (2)
```

For `c in V`, a `c`-clutched section is a sequence

```text
x_0,x_1,...,x_n,
(x_(i-1),x_i) in R_i,             x_n=x_0+c.              (3)
```

Use unnormalized additive convolution

```text
(f*g)(c)=sum_(u in V) f(u)g(c-u).                         (4)
```

Then the exact section count is

```text
N(c)=p (1_(S_1)*...*1_(S_n))(c).                         (5)
```

Indeed put `s_i=x_i-x_(i-1)`.  Equation (3) is equivalent to

```text
s_i in S_i,                    sum_i s_i=c.               (6)
```

For every solution of (6), `x_0` is arbitrary and determines all other
vertices.  This is a bijection and proves (5), including its factor `p`.

The set of serviced clutch classes is therefore the Minkowski sum

```text
Supp N=S_1+...+S_n.                                      (7)
```

## 2. The exact saturation threshold

Cauchy--Davenport says that nonempty `A,B subset F_p` satisfy

```text
|A+B|>=min(p,|A|+|B|-1).                                 (8)
```

For completeness, in the nontrivial range `|A|+|B|-2<p`, if `A+B` had at
most `|A|+|B|-2` elements, enlarge it to a set `C` of that size.  The
polynomial

```text
prod_(c in C)(X+Y-c)                                     (9)
```

vanishes on `A x B`, while the coefficient of
`X^(|A|-1)Y^(|B|-1)` is
`binom(|A|+|B|-2,|A|-1)!=0 mod p`; the coefficient-grid lemma gives a
contradiction.  In the remaining range `|A|+|B|>p`, every translate of `-B`
meets `A`, so `A+B=F_p`.  This proves (8).

Iterating (8) gives

```text
|S_1+...+S_n|
 >=min(p,1+sum_i(|S_i|-1)).                              (10)
```

Consequently

```text
sum_i(|S_i|-1)>=p-1
       ==> N(c)>0 for every c in F_p.                    (11)
```

The threshold is sharp as a universal support statement.  At `p=13`, take

```text
S={0,1,2,3,4,5},             T={0,1,2,3,4,5,6}.          (12)
```

Their excess is `5+6=11=p-2`, and `S+T={0,...,11}` misses clutch `12`.
Increasing `T` to `{0,...,7}` gives excess `12=p-1` and fills `F_13`.

## 3. Two edges: a sharp multiplicity floor

For two edges and a fixed clutch `c`, the number of increment decompositions
per choice of base sheet is

```text
r_(S,T)(c)=(1_S*1_T)(c)=|S intersect (c-T)|.             (13)
```

Two subsets of a `p`-element universe have intersection at least the excess
of their total cardinality over `p`.  Hence

```text
r_(S,T)(c)>=max(0,|S|+|T|-p),
N(c)>=p max(0,|S|+|T|-p).                                (14)
```

Both bounds are sharp for every pair of nonempty sizes and every prescribed
`c`: choose `S` and `c-T` with the smallest possible intersection.  Thus
(14) is not merely an averaged or asymptotic estimate.

At `p=13`, if both relations have eleven allowed increments, then

```text
r_(S,T)(c)>=9,                 N(c)>=13*9=117             (15)
```

for every carry `c`.  Already modulo global translation there are at least
nine distinct two-edge continuations of every clutch.

There is also an exact universal census.  Write

```text
S=V\A,          T=V\B,                 |A|=|B|=2.        (16)
```

Then

```text
r_(S,T)(c)=9+|A intersect (c-B)|.                        (17)
```

The `78` two-subsets split into six undirected-difference classes of size
`13`.  For `1,014=6*13^2` ordered pairs `(A,B)` in the same class, the last
term of (17), as `c` varies, has census

```text
0:10,        1:2,        2:1.                            (18)
```

For the remaining `5,070` ordered pairs it has census

```text
0:9,         1:4.                                        (19)
```

Therefore over all `78^2*13=79,092` triples `(S,T,c)`, the exact
representation and section-count censuses are

```text
r:          9:55,770,       10:22,308,       11:1,014,
N=13r:    117:55,770,      130:22,308,      143:1,014.   (20)
```

## 4. Thin holonomy versus thick holotopy

If every `S_i` is a singleton, each edge is a translation permutation.
Exactly one clutch class `c=sum_i s_i` has sections, and it has precisely
`p` of them.  This is the ordinary local-system regime used by THM-2637:
the cycle carries a well-defined holonomy translation.

For thick relations, there is generally no single transport element.  The
cycle carries the relation-valued clutch spectrum (5).  Once (11) holds,
mere existence of a `c`-clutched section is true for **every** `c`; in the
eleven-by-eleven case it is true with at least nine continuations over each
base sheet.  Thus support-level fixed-branch service cannot identify a
private carry or prove character-trivial holonomy.  It has crossed from thin
holonomy to saturated relation-valued holotopy.

This is deliberately a support statement.  The multiplicities in (20) can
still vary with `c`, and their Fourier transform can retain information.
The theorem does not claim that weighted character data vanish; it says that
nonempty clutch support alone has become maximally nonselective.

## 5. LRC typing boundary

THM-2623 proves many coefficient rows with support size `11` or `12`, and the
reserved THM-2640 targets a carry refinement of those rows.  It is tempting
to place two such rows on consecutive edges and invoke (15).  That inference
is currently illegal.  Equations (2)--(5) require all of the following on one
physical object:

```text
one common base sheet x,
positive transition relations R_i,
translation equivariance in the same carry gauge,
lawful cyclic composition of target and source endpoints.                (21)
```

The proved rows are static/coefficient fibres after marginalization, not
same-base positive transitions satisfying (21).  Conversely, if a future
theorem does produce two eleven-sheet relations satisfying (21), THM-2642
will show that their support is too thick to decode chronology: every carry
will already have at least `117` sections.  One would then need weighted
phase, a private branch selector, or a thinner matched subrelation.

No principal transition, common endpoint pair, scalar row exclusion, or
LRC(14) conclusion follows.  The scalar ledger remains `165`.

## 6. Exact companion

Run

```bash
python3 04-computation/lrc14_cyclic_difference_relation_holotopy_thm2642.py
python3 -O 04-computation/lrc14_cyclic_difference_relation_holotopy_thm2642.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_cyclic_difference_relation_holotopy_thm2642.out.
```

The dependency-free referee uses explicit optimization-safe guards.  It

1. exhausts all `16,129` ordered nonempty subset pairs of `F_7`, checking
   Cauchy--Davenport, (5), and (14) for every clutch;
2. constructs sharp extremizers for all `169` nonempty size profiles over
   `F_13`;
3. checks the threshold and one-below-threshold controls (12);
4. exhausts all `6,084` eleven-sheet relation pairs and all `79,092`
   pair/carry triples, recovering (20); and
5. checks the singleton thin-holonomy boundary.

The LF-normalized SHA-256 hashes are declared in the frontmatter.

An independent hostile audit rederived the exact base-sheet factor in (5),
the iterated Cauchy--Davenport threshold and its one-below hostile, the sharp
intersection construction in (14), and the complete same/different
undirected-difference census (17)--(20).  It separately checked the
thin-versus-thick holotopy interpretation and the load-bearing LRC typing
boundary: the current eleven-sheet coefficient rows are not common-base
positive equivariant transitions.  Normal and optimized executions both
byte-match the stored transcript, and the declared LF-normalized hashes were
independently reproduced.

QED.
