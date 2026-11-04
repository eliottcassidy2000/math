# Exact phase and named-port transport across a zero-ready component reset

**PROVED:** the polynomial identities and operation-specific scalar closures
below. **FINITE-EXACT:** the bounded controls in `check.json` and the named
21-vertex branch counts in `phase.json`. These results do not establish an
unbounded return induction or the full `{7,21}` spectrum conjecture.

The inherited mechanism is the canonical zero-ready decomposition in
[`triangle_return_170807_20261103/ZERO_READY.md`](../../triangle_return_170807_20261103/ZERO_READY.md):
a strong tournament with marked root r and `E_r(2)=H` has a unique in-neighbor
w, and deleting r gives `A -> B`, with A nonempty and B the last strong
component. The hostile inherited example is the marked 81/27 pair in
[`return6325_completion_20261103/halving_response.json`](../../return6325_completion_20261103/halving_response.json).
The corrected near miss is treating a scalar `(h_B,F_w)` target identity as
transport of its phase or named attachment data. The useful sidecar is the
ordered path-bin calculus in
[`polarized_path_calculus_20261103/CALCULUS.md`](../../polarized_path_calculus_20261103/CALCULUS.md).

Our live concepts are: the marked strong component; the zero-ready wrapper;
rooted polynomial phase; mixed spare-port responses; and the universal17
future-H observer. Their relation is operational: changing the component
changes the input to each of the latter four objects. Scalar equality alone
preserves none of their omitted coordinates.

## 1. Polynomial lift through the exact wrapper

For a tournament B marked at w, let `P(z)` count z ordered, possibly empty,
path bins spanning B. Let `E(z)` / `S(z)` require the first bin nonempty and
ending / starting at w. Let `Q(z)=P_(B-w)(z)`. These are integer-valued
polynomials, determined by their nonnegative-integer evaluations. For a
nonempty tournament P(0)=0; for the empty tournament P(0)=1.

Write `T(A,B,w)` for `A -> B` plus a new vertex r that loses only to w.
No strongness assumption is needed for the counting identities. Put
`a(z)=P_A(z)`. Then

```
P_T(z)   = z [a(z+1)E(z+1) + a(z)P(z) - a(z-1)S(z)],
E_r,T(z) = a(z)E(z) + a(z-1)P(z-1),
S_r,T(z) = a(z)P(z) - a(z-1)S(z).
```

Proof. A bin in `A -> B` has its A segment before its B segment. Therefore
its unrestricted polynomial is `a(z)P(z)`, its first-bin-ending-w polynomial
is `a(z)E(z)`, and its first-bin-starting-w polynomial is `a(z-1)S(z)`.
For the last expression, the first bin cannot contain A vertices. Insert r
at the beginning of a bin or immediately after w. The inherited native-birth
identity is `P_new(z)=z[E_old(z+1)+P_old(z)-S_old(z)]`. The first two cases for
a bin ending at r are its being appended after w, or being a singleton;
these give E_r. A bin starting at r must, after removal of r, start anywhere
except at w; this gives S_r. This proves all three identities for positive
integers z and hence polynomially.

For `A=TT_k`, k>=1, `a(z)=z^k`. Put

```
h=P(1), s=S(1), q=Q(1), F=E(2), e=E(1),
P2=P(2), S2=S(2), E3=E(3), Q2=Q(2).
```

The wrapper has

```
H = h + 2^k F,
(E_r(1),S_r(1),H(T-r),E_r(2),S_r(2))
    = (e,h,h,H,2^k P2-S2),
P_T(2) = 2[3^k E3 + 2^k P2-S2].
```

Consequently its rooted golden coordinates are exactly

```
(H, (2^k P2-S2-e)/2, (2^k P2-S2+e)/2).
```

The five component numbers `(h,F,e,P2,S2)` suffice to evaluate H and this
one rooted phase after one specified wrapper. E3 is additionally needed for
the wrapper's full P2. Sufficiency for these evaluations is not a closure
or a minimality claim.

## 2. Two genuine finite scalar closures, with distinct operation scopes

Here R adds x beaten only by w and keeps w marked; N performs the same birth
and marks x. For a component with at least two vertices, the nine-vector
`(h,s,q,F,e,P2,S2,E3,Q2)` updates by

```
R: (h+F-s, q, q, 2F-q, e,
    2(E3+P2-S2), S2+Q2, 3E3-Q2, 2Q2),
N: (h+F-s, h-s, h, h+F, e,
    2(E3+P2-S2), P2-S2, E3+P2, P2).
```

The supporting polynomial identities, in addition to the common
`P'=z[E(z+1)+P(z)-S(z)]`, are

```
R: E'=zE-Q(z-1),   S'=Q(z)+(z-1)S, Q'=zQ,
N: E'=E+P(z-1),    S'=P-S,         Q'=P.
```

For E_R, x starts any of the z bins. Prepending x to the singleton-w
first bin is forbidden; the offending covers are Q(z-1). For S_R, x either
follows the initial w, leaving Q(z) possibilities, or starts one of the
other z-1 bins, giving `(z-1)S`. Deleting w makes x a source, proving Q_R.
For N the first bin ending at x comes from a bin ending at w or a singleton
x. Starting at x excludes first-bin start w. Deleting x recovers B.
Evaluating at 1,2,3 proves the nine-vector closure. The restriction n>=2
ensures Q(0)=P(0)=0 in the stated scalar updates.

The parent's complementary marked-component reader uses a global sink D,
or O, a newborn that beats only w, while keeping w. Its polynomial laws are

```
D: P'=zP,                E'=(z-1)E,      S'=zS,          Q'=zQ,
O: P'=z[P-E+S(z+1)],     E'=(z-1)E+Q(z), S'=zS-Q(z-1), Q'=zQ.
```

These follow by reversing all arcs in the corresponding source/one-port
birth decompositions. They close the eight-vector
`(h,e,S2,F,Q2,q,P2,S3)`, where S3=S(3), as

```
D: (h,0,2S2,F,2Q2,q,2P2,3S3),
O: (h+S2-e,q,2S2-q,F+Q2,2Q2,q,2(P2-F+S3),3S3-Q2).
```

This is sufficient to transport the zero-wrapper golden coordinates for
arbitrarily many specified D/O reader steps. It does not automatically
transport a separate named spare vertex. Neither the nine-vector nor the
eight-vector is claimed minimal or closed under freely mixing all four
operations and arbitrary root changes. The complete polynomial quadruple
`(P,E,S,Q)` obeys all displayed operations; a finite combined scalar quotient
would need its own proof.

## 3. What a retained spare u actually adds

Assume A is nonempty, and let u differ from w in B. Use the following genuine ordered-bin counts on B:

- `A_uv(m)`: first bin starts at u and ends at v, with m-1 arbitrary bins;
- `W_uw(2)`: first bin ends at u, second starts at w;
- `D_wu(2)`: first and second bins end at w,u respectively;
- `G_wu(3)`: first and second bins end at w,u, with a third arbitrary bin;
- `SE_uw(3)`: first bin starts at u, second ends at w, third arbitrary.

With `a_j=P_A(j)`, the old spare's five fields after wrapping are

```
e_u(T)  = a1 e_u(B) + a2 D_wu(2),
s_u(T)  = a1 A_uw(2),
q_u(T)  = a1 H(B-u) + a2 E_w,B-u(2),
F_u(T)  = 2a2 F_u(B) - a1[A_wu(2)+W_uw(2)] + 2a3 G_wu(3),
S_u,T(2)= a1 S_u,B(2) + a2[A_uw(3)+SE_uw(3)].
```

Proof. For e_u, the new r either begins the sole path, or lies internally
immediately after w; removing it gives respectively one u-ending path or
two paths ending at w,u. For s_u, r must lie internally after w; its first
segment starts u and ends w. The deletion identity is the unrestricted
wrapper formula on B-u. For F_u, first put r at the beginning of either
bin, giving `2E_u(2)`, and subtract the two forbidden starts at w:
`A_wu(2)` and `W_uw(2)`. If r lies internally after w, splitting that bin
gives two choices with three bins, two of which end at w,u: `2G_wu(3)`.
For S_u(2), r may start only the second bin, giving S_u(2) minus covers
whose two starts are u,w. Those latter covers cannot accommodate the
nonempty A prefix. An internal r in the first or second bin gives
`A_uw(3)` or `SE_uw(3)` respectively. Finally distribute A across the bins
that can contain it. A bin starting in B cannot contain any A vertices;
all other bins may. This gives exactly the displayed factors a1,a2,a3.

In particular, `3^k` is an actual mixed-interface mode even where the scalar
H clock sees only `2^k`. The component clock and the golden root phase omit
these mixed data. This is a structural missing-coordinate statement, not a
numerical failure of the particular 170807 source.

For the inherited actual components, independent subset-DP fields followed
by these formulas reproduce all five fields at the zero root and spare6:

| Wrapped H | Component nine-vector | Spare6 fields `(e,s,q,F,S2)` |
|---:|---|---|
| 56935 | `(18977,2459,2461,18979,31,2769008,425477,1387531,428504)` | `(4449,186,30241,2185924,162311)` |
| 170807 | `(2195,1,53,42153,53,672920,2673,3358499,52480)` | `(13297,318,91025,7705208,489119)` |

## 4. A same-named phase obstruction and the 81/27 hostile control

The strong five-vertex tournament with rows `(20,1,2,7,14)` has at roots
0 and3 the same `(h,F)=(9,11)`, the same e=1 and P2=84, but S2 respectively
11 and23. Relabel the second graph by `[3,1,2,0,4]`, so both marked vertices
have physical name0 and both spare vertices keep physical name1. Apply the
same TT2 wrapper and give its new root the same name in both constructions.
Both outputs have H53, but the golden triples are

```
(53,162,163) and (53,156,157),
```

with distinct F49/charge coordinates `(2,2,2)` and `(6,3,1)`. Thus even the
same named ports, same unmarked component isomorphism class, same prefix,
and same `(h,F,e,P2)` do not determine the new rooted phase. The omitted
S2 is decisive. `check.py` verifies the full small graphs directly.

The inherited 81/27 graphs both have H462915. Marking physical vertex8,
the same native cut256 gives the known unequal responses1061775 and1395279.
Their corresponding F values are646002 and963549. The new TT1 zero-ready
wrapper therefore gives unequal H1754919 and2390013. Independent literal
subset endpoint DP verifies the wrapper H, P2, root and spare6 profiles;
`check.json` retains the exact current F49/F343 codecs. This control warns
against quotienting out marked incidence even in a highly symmetric graph.
The automorphism orders are inherited labels here, not newly recomputed.

## 5. The universal17 observer and a small-feedback exact engine

The universal two-port theorem is
[`adaptive_state_induction_20261103/COMMON_OBSERVER.md`](../../adaptive_state_induction_20261103/COMMON_OBSERVER.md),
not the older right-suffix17 observer without reversal. It covers twelve
commands encoded `3*incoming_mask+discard`, where discard2 keeps both old
ports. On an actual B with active pair(w,u), prepend k global sources using
command2, then add r beaten only by w and replace slot0 using command3.
This literally constructs the wrapper, with ordered pair(r,u), and

```
X_(T;r,u) = X_(B;w,u) A_2^k A_3.
```

The prefix names must be recorded in reverse birth chronology to recover
`TT_k -> B`. The 17 coordinates close future H for precisely these permitted
commands. The golden phase and arbitrary component reset are not asserted
to factor through X. Pairwise17 records also do not, by themselves, identify
a simultaneous three-port interface. Full actual named rows remain the
source of truth for operations beyond a stored quotient's scope.

`gap.py` supplies exact ordered-bin counts on an actual tournament when
deleting a supplied feedback set leaves a transitive tournament. It promotes
requested endpoint vertices to feedback. For each permutation of feedback
vertices and their allocation to ordered bins, the remaining chain vertices
are assigned to gaps. A gap records empty, nonempty with legal last edge,
or nonempty with illegal last edge. Its first inserted chain vertex must
receive an arc from the left feedback endpoint; its last must send an arc
to the right endpoint. An empty internal gap requires the direct feedback
arc. Chain vertices are processed in their unique transitive order.

Every ordered path-bin arrangement has exactly one feedback permutation,
bin allocation and gap assignment; conversely each accepted assignment
gives exactly one arrangement. This proves the algorithm, including empty
bins and constrained endpoint queries. The state set for a fixed order is
bounded by `3^(r+m)`, with r promoted feedback vertices and m bins. There is
no claim of finite-state universality when r grows under arbitrary resets.

For the parent's two actual21-vertex targets, feedback names are
`(39,31,55)` for +2 and `(39,31)` for +4. `phase.py` independently counts the
marked18-vertex component and full wrapper. Four pairs are retained:
`(39,6),(31,6),(50,51),(39,31)`. Each record includes endpoint/deletion
fields, full P2, ordered W2, and all seventeen universal static continuation
query values. Appending such static queries preserves a transitive
complement after the pair is added to feedback; no `2^21` subset table is
used. The exact graph names and source-file digest are retained.

The independent results are:

| Branch | Component `(h,F)` | Golden triple at root39 | F49/charge |
|---|---|---|---|
| +2 | `(49955,115617)` | `(512423,154166620,154166623)` | `(6,2,0)` |
| +4 | `(72701,109931)` | `(512425,226399221,226399222)` | `(1,0,4)` |

Both F343 current-observer codecs have degree3. The component eight-states,
ordered joint responses and continuation values are in `phase.json`.

F49/charge is the established linear recoding of the golden triple modulo7.
The separate F343 codec treats that triple in the basis `(1,z,z^2)` with
`z^3+z^2+1=0`. Its Frobenius orbit is checked exactly. These are current
observer codecs, not a tournament automorphism or inherited ancestral
spectral-clock transport.

## 6. Reproduction, scope and the remaining bridge

From the repository root:

```
python3 05-knowledge/results/zero_ready_component_returns_20261103/interface/check.py
python3 -O 05-knowledge/results/zero_ready_component_returns_20261103/interface/check.py
python3 05-knowledge/results/zero_ready_component_returns_20261103/interface/phase.py
python3 -O 05-knowledge/results/zero_ready_component_returns_20261103/interface/phase.py
```

`check.py` compiles two small C++ programs in a temporary directory and
cleans them. Its universe is all64 labeled order-four tournaments for1024
literal gap queries, 128 R/N updates, 128 O/sink updates, and64 independent
22-field subset-DP controls. Another68 static query graphs compare the gap
engine with an independent subset endpoint DP. It then checks the two inherited actual
components, the same-named phase collision, and the 81/27 wrapper. Source
pins are in the receipts. Checks use explicit exceptions and remain active
under `-O`. C++ arithmetic is confined to order<=17, where every three-bin
count is at most `binom(n+2,2)n! < 2^63`; the larger targets use Python
integers and the gap engine.

The constructive bridge is now precise. A component reader can carry an
honest closed golden interface while its full named graph supplies the
joint attachment record. A proof of repeated 3H+2 and3H+4 returns still
needs a theorem producing a legal next marked component from each output,
with both target equations, admissible named-port transport and a maintained
operation invariant. A finite collection of successful resets does not
supply that invariant. Assuming spectrum coverage may propose target
components; it cannot discharge this missing construction or its closure.
