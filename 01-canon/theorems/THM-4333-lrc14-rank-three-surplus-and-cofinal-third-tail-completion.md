---
id: THM-4333
title: "LRC(14) rank-three surplus and cofinal third-tail completion"
status: >
  PROVED RELATIVE TO THM-4231/4326/4150/4170/4331 + FINITE-EXACT + O2/O3
  INVARIANT + INDEPENDENT HOSTILE AUDIT; LRC(14) OPEN. On every pair in
  THM-4231's exact 181,194-pair remainder, every fixed-pool nine-body has
  retained rank-three safe mass strictly above 2/27. This gives
  1,060,514,892,450 eight-pool/two-outsider cores with the same strict
  floor, a uniform cofinal third-core-outside cutoff 3,370,132,808, and a
  new two-sheet reserve above 1/189 that admits any further speed s at the
  explicit bound s>=27(C-1), where C is the sum of the preceding twelve
  speeds. These are structured residual-pair families, not arbitrary-row
  entry or LRC(14).
source: root + hypergraph_scout + cofinal_transfer_scout / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4331-lrc14-safe-component-endpoint-denominator-odd-wall-escape
related:
  - THM-4329-lrc14-complete-thirty-label-fixed-outsider-and-thirty-two-label-pascal-chart
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
  - THM-4332-lrc14-fixed-pool-single-constraint-implication-rigidity
artifact_root: 05-knowledge/results/lrc14_rank3_surplus_third_tail_thm4333
artifact_manifest: 05-knowledge/results/lrc14_rank3_surplus_third_tail_thm4333/SHA256SUMS
artifact_manifest_sha256: 796c8fcd9a02875eb5dadfd206c599ea4265f651ac5da63eb9ed83296508f67e
audit: >
  PASS / ACCEPT WITH SCOPED INDEPENDENCE. Full O2/O3 runs give byte-identical
  181,194-row screen and 67,198-row exact-fallback ledgers. The closed
  verifier pins the inherited universe hash/order, every tick identity,
  positivity, body rank, search total, normalized controls, and the exact
  cofinal ceiling, normally and under optimized Python. A separate literal
  midpoint-wall program imports neither the event sweep nor branch-and-bound
  and unpruned-enumerates all C(30,9) bodies on the canonical (50,70)
  hostile, reproducing the minimum and least mask under O2/O3. Two symbolic
  audits independently checked the downstream discrepancy and two-sheet
  fibre arguments, component bounds, equality signs, and speed counts.
---

# THM-4333 -- rank-three surplus and cofinal third-tail completion

**PROVED RELATIVE TO THM-4231/4326/4150/4170/4331 + FINITE-EXACT + O2/O3
INVARIANT + INDEPENDENT HOSTILE AUDIT. LRC(14) REMAINS OPEN.**

## 1. Rank-m failure hypergraphs

Let `P={p_0,...,p_(n-1)}` be a labelled pool and let `Q` be a fixed outsider
set. Resolve the walls of `P union Q` on an integer grid of circumference
`D`. Retain open cells on which every label of `Q` is safe. For such a cell
`C`, let

```text
F(C)={i:||p_i x||<1/14 on C},
w_S=sum of the integer widths of retained cells with F(C)=S.  (1)
```

For a failure-rank cutoff `m`, define

```text
W_m=sum_(|S|<=m)w_S,
d_T^(m)=sum_(S superset T, |S|<=m)w_S,
L_m(B)=sum_(S intersect B=empty, |S|<=m)w_S.           (2)
```

The complete safe mass and its truncation satisfy

```text
0<=L_m(B)<=D mu(G_(B union Q)).                         (3)
```

Finite inclusion-exclusion gives the exact identity

```text
L_m(B)=W_m+
 sum_(empty!=T subset B, |T|<=m)(-1)^|T| d_T^(m).      (4)
```

Indeed, the coefficient of `w_S` is
`sum_(T subset B intersect S)(-1)^|T|`, which is one exactly when
`B intersect S` is empty. Equivalently,

```text
C_m(B)=W_m-L_m(B)
      =sum_(S:B intersect S!=empty, |S|<=m)w_S         (5)
```

is weighted maximum coverage and is monotone submodular. At a partial body
`A`, the exact marginal of a new vertex is

```text
Delta(v|A)=sum_(S:v in S, S intersect A=empty, |S|<=m)w_S. (6)
```

Future marginals can only decrease. Hence the sum of the largest required
current marginals is a sound completion upper bound for branch-and-bound;
pruning strictly below the incumbent is exact, while retaining equality
preserves the least-mask tie break. At the root,

```text
L_m(B)>=W_m-(sum of the |B| largest vertex degrees).   (7)
```

For `m=3`, write `d_i=d_{ {i} }^(3)`,
`d_ij=d_{ {i,j} }^(3)`, and `w_ijk=d_{ {i,j,k} }^(3)`. Then

```text
L_3(B)=W_3-sum_(i in B)d_i
           +sum_({i,j} subset B)d_ij
           -sum_({i,j,k} subset B)w_ijk.               (8)
```

The codegree `d_ij` includes the raw rank-two weight and every rank-three
weight containing `{i,j}`. The intrinsic object is a weighted hypergraph,
not a tournament: triple overlaps cannot be replaced by pair orientations.

## 2. Exact residual-pair surplus theorem

Retain the THM-4326 pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                       (9)
```

Let `E_rem` be exactly THM-4231's frozen `181,194` unordered outsider-pair
remainder. Its ordered file hash is

```text
9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1.
                                                               (10)
```

> **Rank-three surplus theorem.** For every `{q,r} in E_rem` and every
> `B in binom(P,9)`, the rank-three truncation on cells where `q,r` are safe
> satisfies
>
> ```text
> L_3(B)/D>2/27.                                       (11)
> ```
>
> In particular,
>
> ```text
> mu(G_(B union {q,r}))>2/27.                          (12)
> ```

The complete exact split is

```text
pairs                                      181,194
degree-bound positive                      113,996
exact maximum-coverage fallback             67,198
exact fallback positive                     67,198
fallback search nodes                   75,201,674
fallback prunes                         37,153,312.     (13)
```

The least normalized exact fallback is

```text
(q,r)=(50,70),
B mask=011cd402,
B={10,80,85,95,120,143,145,168,193},
D=91,205,797,082,400,
L_3(B)=7,799,459,654,598,
27L_3(B)-2D=28,173,816,509,346.                        (14)
```

This is the minimum only among the `67,198` exact fallbacks. The weakest
uniform certificate over the combined screen/fallback proof is instead the
degree bound on `(509,640)`:

```text
D=74,278,001,143,906,560,
coarse mass=5,502,084,265,759,472,
coarse ticks=272,887,692,624,
top-degree mask=311c6400.                              (15)
```

Thus every row in `(11)` has the explicit certified surplus

```text
mu(G_(B union {q,r}))-2/27
 >= delta_*=1895053421/13927125214482480>0.            (16)
```

This is a lower certificate, not an actual full-mass minimum.

By heredity, `(12)` also holds with `B` replaced by every `K subset P` of
size at most nine: extend `K` to a nine-body and then drop the added
constraints. The LRC-relevant eight-body layer contains exactly

```text
|E_rem| binom(30,8)
=181,194*5,852,925
=1,060,514,892,450                                     (17)
```

distinct ten-label cores `K union {q,r}`.

## 3. Cofinal third core outsider

For each `{q,r} in E_rem`, let `T_cert(q,r)>0` be the pairwise certificate
stored by the complete ledger: the degree-screen ticks on a coarse-positive
pair and the exact minimum ticks on a fallback pair. Thus every
`B in binom(P,9)` satisfies

```text
U=G_(B union {q,r}),
27L_3(B)-2D>=T_cert(q,r)>0.                            (18)
```

If `c(U)` is the number of positive-length components, THM-4170's interval
discrepancy estimate gives

```text
mu(U intersect G_s)>=(6/7)mu(U)-6c(U)/(49s).           (19)
```

Because the nine largest pool labels sum to `2061`,

```text
c(U)<=sum(B)+q+r<=2061+q+r.                            (20)
```

Combining `(18)`, `(19)`, and `(21)` with

```text
(6/7)(2/27)=4/63                                      (21)
```

shows that

```text
s>=ceil(27(2061+q+r)D/(7T_cert(q,r)))                  (22)
```

is sufficient for mass at least `4/63`. The complete pair ledger proves
that the maximum of these conservative pairwise bounds is attained by the
weakest certificate `(509,640)` and equals

```text
ceil(6386581705498394400/1895053421)
=3,370,132,808.                                        (23)
```

No floating-point ceiling is used.

Now take `K in binom(P,8)`, extend it by any `p in P\K`, apply `(22)` to
`B=K union {p}`, and finally drop `p`. This proves:

> For every `{q,r} in E_rem`, every `K in binom(P,8)`, and every integer
> `s>=3,370,132,808`,
>
> ```text
> mu(G_(K union {q,r,s}))>=4/63.                       (24)
> ```

The third outsider is automatically fresh. THM-4150 and Haar invariance
under common dilation therefore give a safe thirteen-speed row

```text
2d(K union {q,r,s}) union {a,b}                        (25)
```

for every positive integer `d` and every two distinct positive odd integers
`a,b`.

## 4. A two-sheet 1/189 reserve

There is a second consequence that appends the thirteenth speed after the
two physical tails rather than before them.

> **Two-sheet reserve lemma.** Let `H` be a nonempty finite positive set with
>
> ```text
> mu(G_H)>2/27.                                        (26)
> ```
>
> For every positive integer `d` and distinct positive odd integers `a,b`,
> put
>
> ```text
> U_phys=G_(2dH union {a,b}).                          (27)
> ```
>
> Then
>
> ```text
> mu(U_phys)>=1/2(mu(G_H)-4/63)>1/189.                 (28)
> ```

**Proof.** Under the doubling map `x -> y=2x`, both points in the fibre over
`y in G_(dH)` are body-safe. Let `F_(a,b)` be the quotient locus on which
the two odd tails spoil both lifts. The cross-comb calculation in THM-4150,
including the gcd pullback for nonprimitive tail pairs, gives

```text
mu(F_(a,b))<=4/63.                                     (29)
```

Outside `F_(a,b)`, at least one of the two fibre points is fully safe. If
`N(y)` counts safe fibre points, exact two-sheet disintegration gives

```text
mu(U_phys)=1/2 int N(y)dy
 >=1/2 mu(G_(dH)\F_(a,b))
 >=1/2(mu(G_H)-4/63).                                  (30)
```

Multiplication by `d` preserves Haar measure. Equation `(26)` and
`1/2(2/27-4/63)=1/189` prove `(28)`. **QED.**

Put

```text
C=2d sum_(h in H)h+a+b                                (31)
```

so `C` is the sum of the speeds in `(27)`. The global reserve corollary of
THM-4331 applies to `(28)` and shows that every positive integer

```text
s>=27(C-1)                                             (32)
```

satisfies

```text
G_(2dH union {a,b,s}) is nonempty.                     (33)
```

Here `C>=6`, so `(32)` gives `s>C`; the new speed is distinct from every
preceding speed. In the application below, `(27)` has twelve speeds.

Apply this with `H=K union {q,r}` from `(17)`. It gives a second explicit
thirteen-speed family. The appended `s` may have either parity. Restricting
it to odd integers produces exactly ten distinct even speeds and three
distinct odd tails.

The endpoint address can be much sharper than the global reserve. For the
hostile-derived core

```text
K={10,80,85,95,120,143,145,168},       (q,r)=(50,70),
```

and physical tails `(1,9)`, THM-4331 Section 6 independently reconstructs a
component that certifies every appended integer `s>=120`. The global reserve
bound is `52,407` for this row, and the earlier discrepancy bound was
`52,434`. This is one VERIFIED-EXACT addressed control, not a uniform
replacement for `(23)` or `(32)`.

## 5. Audit and scope

The primary program reuses THM-4326's wide event-state merger, but its
rank-three aggregation, zeta codegrees, `2/27` target, and optimizer are new.
The two full optimized builds agree byte-for-byte on both ledgers. The closed
verifier checks the frozen remainder order and hash, split equality, tick
identities, positive flags, body ranks, search totals, normalized controls,
and `(23)`. The separate `(50,70)` program reconstructs literal midpoint
walls and flat-enumerates all `14,307,150` bodies without importing the
event sweep or optimizer.

```text
source:       pair-safe wall cells over the fixed pool
target:       every nine-body on E_rem and its eight-body descendants
map:          retain pool-failure ranks <=3
preserved:    exact retained mass and zeta overlaps through rank three
destroyed:    rank>=4 mass, cyclic address, owners, arrival, parity entry
sidecar:      exact grid, full pair ledger, normalized ticks/D
hostiles:     exact fallback (50,70); weakest certificate (509,640)
decisive test: 27L_3(B)-2D>0.                          (34)
```

This theorem covers only `E_rem`. THM-4231's proof outside that remainder
establishes `4/63`, not the stronger `2/27`, so no all-pair rank-three claim
is inherited. The cutoff in `(23)` is sufficient, not minimal. The theorem
does not say rank-three cells exhaust the safe set, identify the actual
safe-mass minimizer, control every third outsider or every odd triple, give
arbitrary-row entry, or prove LRC(14).

## 6. Reproduction

From the repository root:

```text
python 05-knowledge/results/lrc14_rank3_surplus_third_tail_thm4333/verify_pair_rank3_packet.py
python -O 05-knowledge/results/lrc14_rank3_surplus_third_tail_thm4333/verify_pair_rank3_packet.py
```

The packet README gives replay instructions and an independent hostile audit.
