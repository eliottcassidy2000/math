# Transposition rigidity closes all remaining cumulative signed-cycle gaps

**Status: PROVED with FINITE-EXACT base, INDEPENDENTLY AUDITED.** The new
analytic transposition argument proves the
Hamilton-layer minimum at every order at least nine. Together with the
independently audited K8 Hamilton base and the deletion lemma, it closes
every cumulative layer D>=7, with exactly single-edge equality. The K9
Walsh and independent path-attachment recurrences also agree exactly on
every one of their 268,435,456 characters; this is an additional control,
not an assumption of the analytic all-order proof. The canonical statement is
[THM-4427, all-cumulative signed-cycle gaps by transposition
rigidity](../../01-canon/theorems/THM-4427-all-cumulative-signed-cycle-gaps-by-transposition-rigidity.md)
after independent mathematical review and literal parity/compiler-mode audits.
The packet-genus lane and parent both checked the analytic argument.

## Inheritance, hostile, and board

The anchor is the cumulative signed-cycle frontier after
[THM-4416, even-graph cumulative D5/D6 spectral gap](../../01-canon/theorems/THM-4416-even-graph-cumulative-d5-d6-spectral-gap.md).
The closest proved mechanism is its balanced-deletion trichotomy, inherited
from [THM-4083, cumulative D3/D4 spectral gap](../../01-canon/theorems/THM-4083-even-graph-cumulative-d3-d4-spectral-gap.md).
The canonical hostile is the antibalanced signing: all even cycles have
positive sign, so a positive even-layer lower bound on every nonbalanced
class is false. The corrected near miss is MISTAKE-496's inclusion of the
balanced character in a spectral-gap minimum. Here both even-layer zero
classes must be removed before testing a single-layer gap.

The least-used sidecar is the **type** of a zero deletion: balanced or
antibalanced. Keeping this one bit makes two exceptional deletions rigid,
and the individual cycle layer becomes easier to induct than the cumulative
sum. A second hostile is the negative C5 plus a positive apex on K6, with
`c6=20<24`, so the finite Hamilton-layer premise cannot be silently assumed
at every length.

The live concept board was: cumulative deletion weights; full local cycle
profiles; even-layer zero types; Hamilton-cycle minimum distance; a deleted
vertex's attachment quadratic form. The local K8 census yielded 233 profile
vectors and the comparison `3c8>=4c6`, but scalar weighted induction remained
weak near n=9. Passing to the two-zero-class carrier supplies the missing
all-order mechanism. Vertex transpositions then expose a two-star difference
whose Hamilton count is explicit. Its rigid near-zero cases close every
larger Hamilton base, making the final result an all-order theorem.

## 1. General finite-base-to-all-orders transfer

Let H be an edge signing of K_n, modulo cut switching, and let c_k(H)
count unoriented simple negative k-cycles. Put

```text
e_(n,k)=(n-2)!/(n-k)!.
```

The balanced class B has all cycle products +1. Write A for the antibalanced
class, represented by every edge negative. Its triangle signs are all -1,
and its even cycle products are all +1. Global edge negation H -> -H
interchanges B and A, preserves every even c_k, and complements every odd
cycle sign. Define

```text
Z_k={B}       if k is odd,
Z_k={B,A}     if k is even.

E_k={single-negative-edge switching classes}          if k is odd,
E_k={single-edge classes and their global negatives}  if k is even.
```

**Transfer lemma.** Fix k>=3 and n0>=max(k,4). Suppose on K_(n0),

```text
c_k(H)>=e_(n0,k)  for every H outside Z_k,
```

with equality exactly on E_k. Then for every n>=n0,

```text
c_k(H)>=e_(n,k)   for every H outside Z_k,             (1)
```

again with equality exactly on E_k. In particular these are the full
nonzero single-layer minima, and Z_k is exactly the zero set in that range.

### Proof, including the two-type overlap

Induct on n>n0. A vertex v is exceptional if H-v belongs to Z_k.
For even k, two exceptional deletions cannot have opposite types. If H-u
were balanced and H-v antibalanced, their common induced signing on
K_(n-2) would have every triangle both positive and negative. Since n>=5,
that common graph has a triangle, a contradiction. This argument concerns
the actual overlapping triangle signs; the individual switching gauges
need not agree.

If two deletions are balanced, every negative triangle contains both
deleted vertices u,v. Thus all negative triangles are of the form uvx.
The product of the four triangle signs on a four-set is +1, so uvx and uvy
have the same sign. Since H is not balanced, they are all negative.
Triangle signs determine a switching class, hence H is the single-negative-
edge class uv. If two deletions are antibalanced, apply the same argument
to -H; this gives a globally negated single-edge class. Such classes have
exactly e_(n,k) negative k-cycles for the relevant parity.

Otherwise at most one deletion is exceptional, so at least n-1 deletions
obey the inductive bound. Each k-cycle survives exactly n-k vertex
deletions. Consequently

```text
(n-k)c_k(H) = sum_v c_k(H-v)
           >= (n-1)e_(n-1,k),

c_k(H) >= [(n-1)/(n-2)]e_(n,k) > e_(n,k).             (2)
```

Here n>n0>=k, so division by n-k is valid. This proves the inequality and
strictness outside E_k, completing the induction. QED.

The transfer requires the entire finite base, including its equality set.
It does not claim a Hamilton base for arbitrary k. The K6 hostile above is
the first needed check before any such generalization.

## 2. Exact Hamilton bases at eight and nine vertices

The finite bases are

| order n and layer k | switching classes | unoriented k-cycles | minimum outside Z_k | exact labelled ties |
|---|---:|---:|---:|---:|
| n=k=8 | 2^21=2,097,152 | 2,520 | 720=6! | 56=2 binom(8,2), exactly E_8 |
| n=k=9 | 2^28=268,435,456 | 20,160 | 5,040=7! | 36=binom(9,2), exactly E_9 |

No frustration, degree, connectivity, sign density, or relabelling filter
is imposed. All edges incident to vertex zero are switched positive; every
switching class has exactly one representative in this gauge. The remaining
`binom(n-1,2)` bits are enumerated completely. The balanced mask is zero;
the antibalanced root-gauged mask has every remaining bit one. The latter
is excluded only for the even layer, and is a positive hostile control for
the odd layer.

The transfer lemma therefore gives

```text
c8(H)>=(n-2)!/(n-8)!   for all n>=8 and H outside {B,A},
c9(H)>=(n-2)!/(n-9)!   for all n>=9 and H outside {B}.  (3)
```

Equality in the first line is exactly on globally signed single-edge
classes; equality in the second is exactly on single-edge classes.
Section 2A proves the corresponding Hamilton base analytically for every
length at least nine, so (3) extends to all subsequent fixed layers.

### Primary full Walsh path and integer width

A Hamilton cycle through the positive-gauged root is encoded by the edges
of its Hamilton path on the other n-1 vertices. Enumerate every path up to
reversal; its frequency vector has length `2^binom(n-1,2)` and entries
zero or one. An integer Walsh transform returns the signed Hamilton-cycle
sum for every switching class. Thus

```text
c_n(H)=(total_cycles-Walsh(H))/2.
```

Every intermediate Walsh entry is a signed sum of a subset of the cycle
frequency vector, so its absolute value is at most total_cycles. For n=9
this is 20,160, below 2^15. One in-place int16 array therefore suffices,
using exactly 536,870,912 bytes (512 MiB); no int16 overflow or wraparound
is used. Subtraction and count conversion occur in wider integers.

### Independent base audits

For n=8, the independent checker builds cycles from all full vertex
permutations, deduplicates their masks, and literally counts odd edge
parity for every one of the 2,097,152 classes and 2,520 cycles. It uses
no transform and imports no primary source. Its complete spectrum agrees
with the primary path: FNV-style rolling digest `e2dfba14125e7983`.
The explicit source seed is `1469598103934665603` and multiplier
`1099511628211`; this is a reproducibility checksum, not a cryptographic
certificate or the conventional FNV offset basis.

For n=9, delete the last vertex from each Hamilton cycle. The result is a
Hamilton path on K8 with distinct endpoints u,v. If t_u are the signs of
the new incident edges, with t_0=1 under the root gauge, then the exact
signed Hamilton sum is the quadratic polynomial

```text
Q_H(t)=sum_(0<=u<v<8) W_uv(H)t_u t_v,                 (4)
W_uv(H)=sum_(Hamilton paths u->v on K8) sign(path).
```

Each endpoint pair has 6!=720 paths; the 28 pairs partition all 20,160
Hamilton cycles. The independent checker constructs these path channels
directly, evaluates each on all 2^21 K8 classes, and then evaluates (4)
on all 2^7 attachments. This is a distinct factorization of the computation
and never builds the full K9 cycle-frequency vector. Its coefficient arrays
use 112 MiB and its attachment buffer has 128 entries. For the final full
per-character comparison, loading the independent primary reference adds
512 MiB, bringing the array total to 624 MiB. The map from parent
and attachment bits to full K9 edge bits is a bijection, so all 2^28 classes
are covered exactly once. All root/deletion incidences are retained.

The final K9 comparison checked every signed sum, not just minima or a
checksum: `full_reference_compared=1`. The independent keyed spectrum sum
is `b68ca4e0ac0ebf9b`, agreeing with the primary path; the primary rolling
digest is `402a35163eea3383`. K8 also passed all three paths, including
per-character equality of the independent attachment recurrence.

## 2A. Analytic Hamilton-layer rigidity at every order at least nine

Let n>=9, write `w(H)=c_n(H)` and `e=(n-2)!`. We prove

```text
min_(H outside Z_n) w(H)=e,
equality exactly on E_n.                              (H)
```

Here Z_n and E_n have the parity-dependent meanings of Section 1. The
proof has no finite enumeration or asymptotic exception.

### A transposition exposes a two-star signing

For distinct vertices u,v, let tau transpose their labels and let
`H+tau H` mean edgewise multiplication of signs, or xor of negative-edge
sets. For x outside {u,v}, the edges ux and vx are both negative in this
difference exactly when their signs in H disagree. All other edges are
positive. Thus H+tau H is a negative K_(2,r), where

```text
r=#{x outside {u,v}: sign(ux)!=sign(vx)},
s=n-2-r.
```

For every cycle, its negative indicator under edgewise multiplication is
the xor of the two original indicators. Relabelling preserves the number
of negative Hamilton cycles, so

```text
w(H+tau H)<=w(H)+w(tau H)=2w(H).                      (8)
```

### Exact negative Hamilton count of the two-star difference

For the signing K_(2,r) just described, put q=rs. Its negative Hamilton
count is

```text
T_n(r)=2rs[(n-2)(n-3)-2rs](n-5)!.                     (9)
```

To count it directly, root an oriented Hamilton cycle at u and divide by
two. The distinguished vertex v is adjacent to u in two possible positions,
at distance two in two positions, or has one of n-5 remaining positions.
The relevant edge product is the product of the R/S types of the neighbours
of u and v, with repeated neighbours cancelling.

- Adjacent u,v: after orientation division, the negative count is
  `2rs(n-4)!`.
- Distance two: the shared neighbour cancels; summing its choices gives
  the same count `2rs(n-4)!`.
- Both distances at least three: four distinct neighbour slots have an
  odd number of R entries. The count is
  `2rs[(r-1)(r-2)+(s-1)(s-2)](n-5)!`.

Adding these three contributions yields (9). The division by reversal is
valid for labelled simple cycles; no cyclic order or repeated-neighbour
configuration is dropped. The formula also holds at n=6, where its weaker
values preserve the known negative-C5-apex hostile.

An independent root derivation deletes u,v first. On a Hamilton cycle of
the other m=n-2 vertices, let X count R/S-crossing gaps. There are `2X`
negative insertions when u,v occupy the same gap and `2X(m-X)` when they
occupy distinct gaps. Using

```text
E X=2rs/(m-1),
E X^2=4rs(rs-1)/[(m-1)(m-2)]
```

and multiplying `2 E[X(m+1-X)]` by `(m-1)!/2` independently recovers (9).
Thus both the cycle-position count and the deleted-cycle insertion count
retain the contribution that naive independent insertion would lose.

### Every interior row-disagreement count costs more than twice the edge

Suppose `2<=r<=n-4`. Then

```text
2(n-4)<=q=rs<=(n-2)^2/4.
```

The function `F(q)=q[(n-2)(n-3)-2q]` is concave. To show
`T_n(r)>2e`, it therefore suffices to compare F with
`(n-2)(n-3)(n-4)` at the two real endpoints of this interval. Their excesses
are respectively

```text
(n-4)(n^2-13n+38),
(n-2)(n-4)(n^2-12n+28)/8.                            (10)
```

Writing n=9+t, the two quadratic factors become `t^2+5t+2` and
`t^2+6t+1`. Both are strictly positive for t>=0. Consequently (8) implies

```text
w(H)<=e => r(u,v) in {0,1,n-3,n-2} for every u!=v.    (11)
```

This is an exact strict threshold, not a large-n estimate. It fails as a
proof method at n=8 and below, which is why the separately audited K8 base
is needed.

### Classifying the near-constant disagreement matrices

Switch all edges at vertex zero positive. Let G be the simple graph of
negative edges on the other m=n-1 vertices. Applying (11) to (0,v) gives

```text
deg_G(v) in {0,1,m-2,m-1}.                            (12)
```

Let L be the vertices of degree at most one and R those of degree at least
m-2; they partition V(G). Write a=|L| and b=|R|. Every low-degree vertex
has at most one neighbour in R, so the number of L/R edges is at most a.
Every high-degree vertex misses at most one vertex in L, so it is at least
b(a-1). If a,b>=2 this forces

```text
b(a-1)<=a, equivalently (a-1)(b-1)<=1.
```

Hence a=b=2, impossible since m>=8. Thus one part has size at most one.

If all vertices are low, G is a matching plus isolated vertices. Two
disjoint matching edges would give a pair of endpoints with row-disagreement
r=2, forbidden by (11). Therefore G is empty or one edge. Complementing G
preserves the allowed disagreement condition, so the all-high cases are
the complete graph or the complement of one edge.

If exactly one vertex h is high, it meets all other vertices or all except
one. The first possibility is a full star. In the second, its missed vertex
is isolated because every other low vertex already uses its one edge at h.
For any neighbour v of h the row-disagreement is then m-3, which is excluded
by (11). Thus only the full star remains. Complementation gives the case
of exactly one low vertex.

In this root gauge, the six families are therefore

```text
empty, complete, one edge, complement of one edge,
full star, complement of a full star.                (13)
```

The empty and complete graphs represent B and A. A full star is the
root-gauged representative of a single negative edge incident to vertex
zero. Graph complementation realizes global edge negation followed by
restoring the root gauge. Thus (13) is exactly `{B,A}` together with all
globally signed single-edge classes.

### Finishing the Hamilton theorem

For even n, B and A have no negative Hamilton cycles, and the two signed
single-edge families both have e. For odd n, B has zero, a single edge has
e, A has `(n-1)e/2`, and a globally negated single edge has `(n-3)e/2`.
The last two numbers are strictly greater than e for n>=9. Therefore (11)
and (13) prove (H), including its full equality set.

Applying Section 1 with n0=k now proves, for every k>=9 and every n>=k,

```text
c_k(H)>=e_(n,k) outside Z_k, equality exactly on E_k.  (14)
```

Together with the K8 base this covers every fixed cycle layer k>=8.

## 3. Cumulative closure for every D>=7

Write

```text
S_D=sum_(k=3)^(D+1)c_k,
A_(n,D)=sum_(k=3)^(D+1)e_(n,k).
```

For any nonbalanced H outside A, combine THM-4416's all-order D6 bound
with the first line of (3):

```text
S7(H)=S6(H)+c8(H)>=A_(n,6)+e_(n,8)=A_(n,7).          (5)
```

The antibalanced class requires its own actual consequence, because c8=0.
There

```text
S7(A)=sum_(k=3,5,7) n!/[2k(n-k)!].
```

Its excess over A_(n,7), after writing n=8+t, is

```text
t^7/14 + 3t^6/2 + 73t^5/5 + 189t^4/2
 +1328t^3/3 + 2723t^2/2 + 244667t/105 + 1652 > 0     (6)
```

for every t>=0. Thus (5) holds strictly at A. Equality in (5) elsewhere
forces equality in the inherited D6 bound, hence H is a single-edge class.

For n>=9, adding the second line of (3) to the proved D7 bound gives

```text
S8(H)>=A_(n,7)+e_(n,9)=A_(n,8),                      (7)
```

again with equality exactly on single-edge classes. Both lower bounds are
attained because an individual edge belongs to e_(n,k) simple k-cycles.

More generally fix any D>=7 and n>=D+1. For H outside {B,A}, use the
inherited D6 bound and every single-layer bound (14), including k=8:

```text
S_D(H)>=A_(n,6)+sum_(k=8)^(D+1)e_(n,k)=A_(n,D).       (15)
```

The antibalanced class is treated by pairing each odd length k-1 with the
next even length k. Its odd-layer count is
`N_(n,k-1)=n!/[2(k-1)(n-k+1)!]`, and

```text
N_(n,k-1) >= e_(n,k-1)+e_(n,k)
iff n(n-1)-2(k-1)(n-k+2)>=0.
```

Writing n=k+t, the last expression is

```text
(k-1)(k-4)+t+t^2,                                  (16)
```

nonnegative for every even k>=4 and t>=0. The first pair k=4 is strictly
positive since n>=8. If D+1 is odd, the unpaired final odd layer satisfies
`N_(n,D+1)>e_(n,D+1)`, since their ratio is `n(n-1)/[2(D+1)]>1`.
Thus S_D(A)>A_(n,D) for every D>=7 and n>=D+1.

Equality in (15) forces equality in the D6 summand and hence exactly a
single-edge class. Such classes attain every summand. This proves the full
cumulative conjecture for every D>=7, with exactly binom(n,2) labelled
equality classes and one relabelling orbit. Combined with THM-4083 and
THM-4416, all previously open cumulative layers are closed.

In THM-4078's multiplicity-weighted cycle operator, the cumulative
unnormalized Laplacian gap is 2A_(n,D), with that labelled multiplicity.
No Booleanized quotient spectrum is substituted.

## 4. Connection contract and literature-inspired method

The source is a switching class with its triangle signs. The first target
is its list of vertex restrictions and their zero types; the map is induced
deletion. This preserves the exact k-cycle parity and loses the ambient
number of occurrences, restored by the n-k sidecar in (2). Remembering the
zero type restores the overlap test. This is why applying a cumulative
scalar bound alone is weaker here.

The independent computation maps a Hamilton cycle with a distinguished
vertex to its deleted Hamilton path plus an endpoint pair. The preserved
predicate is the complete signed cycle sum in (4). Dropping the endpoint
pair would lose the factors t_u t_v; the 28 labelled channels restore them.
The decisive test is full per-character agreement with the direct cycle
transform, including the balanced, antibalanced, and every single-edge class.

[Heule--Scheucher, arXiv:2403.00737](https://arxiv.org/pdf/2403.00737),
Sections 4--5, motivated seeking a compact local witness and independently
checking a different encoding. Here the finite local witness is a Hamilton
layer; the graph-type overlap and the deleted-path channels supply the
compression. The paper's planar orientation clauses and empty-hexagon
theorem are not imported as mathematical dependencies. No SAT solver or
package installation was needed: native integer Walsh arrays fit well below
the 2 GiB memory bound.

The all-D cumulative gap does not classify the complete cycle-profile cone,
near-minimizers, or every short individual Hamilton layer. It gives no
tournament H>=disc inequality, Boolean quotient gap, or LRC(14) theorem.

## 5. Reproduction

```bash
g++ -O3 -std=c++17 -mpopcnt 04-computation/overnight_hexagon_sep05_hamilton_walsh.cpp -o /tmp/overnight_hexagon_sep05_hamilton_walsh
/tmp/overnight_hexagon_sep05_hamilton_walsh 8 /tmp/overnight_hexagon_sep05_k8_spectrum.bin
/tmp/overnight_hexagon_sep05_hamilton_walsh 9 /tmp/overnight_hexagon_sep05_k9_spectrum.bin
g++ -O3 -std=c++17 -mpopcnt 04-computation/overnight_hexagon_sep05_hamilton_direct.cpp -o /tmp/overnight_hexagon_sep05_hamilton_direct
/tmp/overnight_hexagon_sep05_hamilton_direct
g++ -O3 -std=c++17 -mpopcnt 04-computation/overnight_hexagon_sep05_hamilton_extension.cpp -o /tmp/overnight_hexagon_sep05_hamilton_extension
/tmp/overnight_hexagon_sep05_hamilton_extension 8 /tmp/overnight_hexagon_sep05_k8_spectrum.bin
/tmp/overnight_hexagon_sep05_hamilton_extension 9 /tmp/overnight_hexagon_sep05_k9_spectrum.bin
g++ -O3 -std=c++17 -mpopcnt 04-computation/overnight_hexagon_sep05_transposition_audit.cpp -o /tmp/overnight_hexagon_sep05_transposition_audit
/tmp/overnight_hexagon_sep05_transposition_audit
python3 04-computation/overnight_hexagon_sep05_transposition_symbolic.py
python3 -O 04-computation/overnight_hexagon_sep05_transposition_symbolic.py
```

The exploratory full-profile source is
`04-computation/overnight_hexagon_sep05_d7_profile.py`; its local comparison
search is not needed by the final all-order argument. Its reported floating
ratio was only exploratory; the comparison `3c8-4c6>=0` is now checked
independently with exact integer arrays in that source.

The transposition checker enumerates every unoriented Hamilton cycle on
K5 through K11 and every split `0<=r<=n-2`. Literal edge parity, three
cycle-position contributions, and a separate enumeration of deleted-cycle
crossing-gap moments/insertion counts all agree with (9). It also enumerates
every switching class on K6, K7 and K8, obtaining exactly 32, 44 and 58
classes satisfying (11), equal to `{B,A}` and the two signed edge families.
The negative-C5/apex hostile gives 20 rather than 24, and two disjoint
negative edges fail the rigidity predicate. Thus the algebraic threshold,
equality classes and false overgeneralization all have explicit controls.

The symbolic checker uses rational polynomial identities for (9), (10),
(16), the deletion ratio and strict concavity. Positivity for all n>=9
comes from the exact shifted polynomials with positive coefficients; its
additional integer-endpoint checks through n=200 are finite corroboration,
not the all-order proof. Normal and optimized Python outputs are identical.
The independent C++ transposition audit also gives identical output under
`-O0` and `-O3`; all failure gates are explicit, not removable assertions.

Frozen output records are the same-stem `hamilton_walsh8.out`,
`hamilton_direct8.out`, `hamilton_extension8.out`, `hamilton_walsh9.out`,
`hamilton_extension9.out`, `transposition_audit.out` and
`transposition_symbolic.out` in this results directory. Large native-endian
temporary spectrum binaries are replay artifacts, not portable proof files;
regenerate them on the same host as the independent comparison.

### Frozen raw-LF SHA-256 manifest

Paths below are relative to the repository root. The source and output
hashes provide artifact integrity, not a substitute for the mathematical
proof or independent full-spectrum comparison.

```text
c186bf9ef61d6557798a30439a4c1d0c407c894edf13515032ea1819d033812e  04-computation/overnight_hexagon_sep05_hamilton_walsh.cpp
68a96abdac3f62538df92bdb0288e62e536e25f7819b0036a3965e583530a6be  04-computation/overnight_hexagon_sep05_hamilton_direct.cpp
5f4a13e5ad564069c420b79d17376485e614a992ff727ecf051e7e2f8e8b266f  04-computation/overnight_hexagon_sep05_hamilton_extension.cpp
7baa74c5779bdf6298d3892f76d26ff7c0988b111f2633333414b37dd369a5f3  04-computation/overnight_hexagon_sep05_transposition_audit.cpp
8ca1a13e0063922a5f1f3a1f51e0ac30b10f0389d7c4d05253810e889e8ff287  04-computation/overnight_hexagon_sep05_transposition_symbolic.py
5defbffc3cc6819a681b67b2ebfb147ae82c03ba5384e2d61fb95c705e1c0501  04-computation/overnight_hexagon_sep05_d7_profile.py
6da47f2104e3aa035592311b1c43e7794fee701633e1909cfd337810d8e70657  05-knowledge/results/overnight_hexagon_sep05_hamilton_walsh8.out
797c621ae68f29e0e5cd24419dcdc5fe268f34d7bdf201eb4b322b1a7617decd  05-knowledge/results/overnight_hexagon_sep05_hamilton_direct8.out
9268b5950eb0d17e3b67c47c350248c280984d3d1b2f85c224acac7e4d1d1ffb  05-knowledge/results/overnight_hexagon_sep05_hamilton_extension8.out
94d0d37617100bf1942013e49472fd112e1574bb587a997ff9f73075bc3f1c4d  05-knowledge/results/overnight_hexagon_sep05_hamilton_walsh9.out
dd426ab3a48042e80e7224afb533dd10cd744b804ed9c68214adbb53edfa83bc  05-knowledge/results/overnight_hexagon_sep05_hamilton_extension9.out
f99b3fea6ad7fad6393cd540fd4cc644bbb699dc77505cc795874eee5828680c  05-knowledge/results/overnight_hexagon_sep05_transposition_audit.out
e23cdb095dd3bef7f31bf7e3f8a2e7054c58c339ee9078c6cdd937110325011d  05-knowledge/results/overnight_hexagon_sep05_transposition_symbolic.out
```
