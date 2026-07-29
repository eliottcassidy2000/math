---
id: THM-2881
title: Rank-impossible pair-residual closure of thirty-eight THM-741 roots
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.
  Thirty-eight non-flood THM-741 roots
  on which the scalar top-two/two-tail gate fails are closed uniformly.
  A globally sealed exact two-comb union cap gives at most one tail above
  a root-specific cutoff W<=4628; exact pair-aware one-tail and all-head
  recursions then have zero failures.  Normal and optimized replays agree.
  Global THM-741 and LRC(14) remain open.
source: root/lrc-rank-impossible-overlap-2026-07-29
depends_on:
  - THM-731-covering-middle-order-x-integral-autocorrelation-discrepancy
  - THM-732-disc-v-bernoulli-edge-pair-dedekind-form-exact-certificates-far-element-tail
  - THM-738-near-AP-three-slot-closure-all-1001-bodies-in-1-14
related:
  - THM-741-near-AP-four-slot-closure-all-2002-bodies-in-1-14
  - THM-2209-sharp-quadratic-reversed-peel-and-joint-fourier-ledger
verification:
  - 04-computation/lrc14_thm741_rank_impossible_pair_residual_closure_codex_20260729.py
  - 05-knowledge/results/lrc14_thm741_rank_impossible_pair_residual_closure_codex_20260729.out
---

# THM-2881 — rank-impossible pair-residual closure of 38 THM-741 roots

## 1. Statement and exact universe

Let `B_38` be the following lexicographically ordered collection of
nine-subsets of `{1,...,14}`:

```text
1234578(11)13; 1234579(11)13; 123478(10)(11)13;
123479(10)(11)13; 12348(10)(11)1314; 1235789(12)13;
123578(10)(11)13; 123578(11)(12)13; 123579(10)(11)13;
1236789(10)13; 123679(10)(11)13; 123789(10)(12)13;
12378(10)(11)(12)13; 12378(10)(11)1314; 1245678913;
1246789(10)13; 124689(10)1314; 12469(10)(12)1314;
12489(10)(12)1314; 1256789(10)13; 126789(10)(12)13;
12689(10)(12)1314; 134578(10)(11)13; 13458(10)(11)1314;
1346789(10)13; 134789(10)(12)13; 13478(10)(11)1314;
13489(10)(12)1314; 1456789(10)13; 1456789(11)13;
14589(10)(12)1314; 146789(10)(11)13; 146789(10)1314;
14689(10)(12)1314; 156789(10)(11)13; 234578(10)(11)13;
23478(10)(11)1314; 24689(10)(12)1314.
```

Parentheses only separate two-digit labels.  The verifier contains the same
literal tuples; the SHA-256 of their newline-joined CSV representation is

```text
81e0807f389f0b10ac0ff4bfd9e3c84f86b003edc9b5469d1c581ca6e25b0c4d. (1)
```

For every `E in B_38` and every four distinct integers

```text
15<=a<b<c<d,
```

the thirteen-speed family `E union {a,b,c,d}` has positive lonely-time
measure.  By the proved
[`THM-738-near-AP-three-slot-closure-all-1001-bodies-in-1-14.md`](THM-738-near-AP-three-slot-closure-all-1001-bodies-in-1-14.md),
the same holds when one or more of the four added speeds is at most `14`.
Consequently all `38` whole THM-741 roots in `B_38` are closed uniformly.

The set `B_38` is the independently reproduced non-flood rank-impossible
subatlas on which the scalar two-head/two-tail margin is nonpositive.  That
classification explains why these rows were selected; the theorem itself
reconstructs and checks every type condition it uses.

## 2. Root coverages and the missing coordinate

For a root `E`, write

```text
G_E = its good set,             m_E=|G_E|,
r_E = number of components of G_E,
c_E(w)=|G_E intersect D_w|.                            (2)
```

Put `s=99/70`.  The THM-731/732 covariance-discrepancy chain gives the
strict estimate

```text
c_E(w)<u_E(w):=m_E/7+s r_E/(7w).                       (3)
```

The exact `15..600` scan orders the coverages
`q_1>=q_2>=q_3>=q_4`.  On every row,

```text
q_4>m_E/7,
s r_E/[7(q_4-m_E/7)]<601,                              (4)
```

so `(3)` proves that these are the global top four values over every
`w>=15`.  The verifier also checks

```text
m_E-q_1-q_2-q_3<=m_E/7,
q_1+q_2>=5m_E/7.                                      (5)
```

Thus the scalar fourth-rank and scalar two-head/two-tail gates genuinely fail
on every row.  The missing coordinate is not another single coverage.  It is
the exact two-comb union

```text
U_E(a,b)=|G_E intersect (D_a union D_b)|.              (6)
```

## 3. The global pair cap

First, every row satisfies

```text
q_1<4m_E/7.                                            (7)
```

The smallest margin in `(7)` is

```text
4m_E/7-q_1=5633/490490
```

at `E=123578(10)(11)13`.  Therefore

```text
w>s r_E/[7(4m_E/7-q_1)]
  ==> q_1+c_E(w)<5m_E/7.                              (8)
```

The largest crossing threshold in `(8)` is

```text
57081024/189619<2501
```

at `E=12689(10)(12)1314`.

It remains to inspect pairs with both speeds at most `2500`.  Order these
speeds by decreasing `c_E`.  Maintain the largest exact pair union already
seen.  A pair whose single-coverage sum is no larger than that incumbent
cannot improve it, because `U_E(a,b)<=c_E(a)+c_E(b)`.  This gives an exact
branch-and-bound maximum after only `313` paid pair carriers across all
`38` roots.  Every paid carrier is rebuilt directly as the eleven-speed good
set `G(E union {a,b})`.

Let `H_E` be that exact finite-head maximum and define

```text
Ubar_2(E)=max(H_E, q_1+u_E(2501)).                     (9)
```

Equations `(3)`, `(4)`, and `(8)` show that `(9)` is a global upper bound for
`U_E(a,b)` over all distinct `a,b>=15`.  Exact arithmetic gives

```text
delta_E:=5m_E/7-Ubar_2(E)>0                            (10)
```

on every row.  The smallest margin is

```text
delta_E=56731/40612572
```

again at `E=123578(10)(11)13`.

## 4. At most one tail

Choose

```text
W_E=floor(2s r_E/(7delta_E)).                          (11)
```

Then

```text
Ubar_2(E)+2u_E(W_E+1)<m_E.                             (12)
```

The maximum cutoff is

```text
W_E=4628
```

at `E=123578(10)(11)13`; the smallest strict margin in `(12)` is

```text
58481/313325992980.
```

Suppose a four-speed extension contains at least two speeds greater than
`W_E`.  Treat any other two speeds as a pair and the chosen two tails as
single combs.  By `(3)`, `(9)`, and `(12)`, their union inside `G_E` has
measure strictly less than

```text
Ubar_2(E)+2u_E(W_E+1)<m_E.
```

This argument covers the `2+2`, `1+3`, and `0+4` head/tail splits at once:
the globally capped pair may itself contain zero, one, or two tail speeds.
Thus every possible obstruction has at most one speed above `W_E`.

## 5. Exactly one tail

Assume exactly one speed `t` exceeds `W_E`, and order the other three by
decreasing root coverage:

```text
c_E(a)>=c_E(b)>=c_E(c).                               (13)
```

The pair union `U_E(a,b)` is computed exactly.  Unless

```text
U_E(a,b)+c_E(c)+u_E(W_E+1)>=m_E,                      (14)
```

the branch is already closed.  The old scalar prefilter admits `82,346`
triples.  Pair-aware branch-and-bound reduces them to `328` exact pair
carriers and only `686` literal three-comb residuals

```text
G_{E;a,b,c}=G_E minus (D_a union D_b union D_c),
h=|G_{E;a,b,c}|,       r=#components.                  (15)
```

On `(15)`, THM-732 gives, uniformly for `t>W_E`,

```text
|G_{E;a,b,c} minus D_t|
 >6h/7-sr/[7(W_E+1)].                                 (16)
```

Every one of the `686` values in `(16)` is positive.  The minimum is

```text
110133/8472100
```

for the twelve-speed carrier

```text
{1,4,5,6,7,8,9,11,13,19,20,48},
h=30721/1452360,       r=16.                           (17)
```

Hence the exactly-one-tail chamber is empty.

## 6. All four speeds in the head

Order a head quadruple by decreasing individual root coverage.  Its top pair
is again evaluated exactly.  The top-three residual is constructed only when
the exact pair union plus the two remaining single coverages can reach
`m_E`.  A fourth comb is tested only when its root coverage is at least the
literal top-three residual measure.

Across all `38` rows this output-sensitive recursion needs

```text
381 exact pair carriers,
827 exact triple carriers,
343 literal quadruple obligations.                    (18)
```

Every one of the `343` quadruples is rebuilt both by nested exact subtraction
and by direct thirteen-comb union.  All survivors are positive and none is
tight.  The minimum is

```text
67811/8621550
```

at

```text
{1,3,4,5,7,8,10,11,13,17,18,23,25},
```

whose survivor has eight components.  This closes the all-head chamber.

Sections 4--6 exhaust the pure tail.  If one added speed is at most `14`, the
completed family has at least ten speeds in `{1,...,14}` and THM-738 applies.
This proves the whole-root statement.

## 7. Verification, controls, and scope

The canonical ledger digests are

```text
pair       2929656c496ab00fe415790fcedda8e5c627186891eb32721cd0c48a155d9ac2
one-tail   298f679473c4faadbd4c7171c858e82c8ed90bbfa4b9023a77ab95abfde52ffc
all-head   42cd555dd945fe0d1e4487fb5820c23fb2d3a077ec1c675675a946a285db1232
consequence f17bb1d243cf7657948a4da530f130618df34881ee5a318d3ec90f4481e04529.
```

The script and stored-output hashes are

```text
script 05aab2283513ded0747c5708a0bfd9d9ba8d163e4622b5e3638bddd6a162a6d1
output 4ee45dbe8e087f5c1bef9f00fca9916fa709d5a857aaf79ca85d1a0538db2923.
```

Reproduction:

```bash
python3 04-computation/lrc14_thm741_rank_impossible_pair_residual_closure_codex_20260729.py --workers 6
python3 -O 04-computation/lrc14_thm741_rank_impossible_pair_residual_closure_codex_20260729.py --workers 6
```

Both commands reproduce the stored report byte for byte.  Load-bearing
controls include the fixed core hash, literal body digest, root-rank and
scalar-failure type checks, strict tail crossings, independent direct
reconstruction of every paid pair/triple/quadruple carrier, and positive
hostiles in both the one-tail and all-head chambers.

An independent proof audit read the full verifier and separately checked:

1. the global `q_1` seal and completeness of the pair branch-and-bound;
2. the exhaustive `0/1/>=2` tail partition and strict floor choice in `(11)`;
3. the unique coverage-rank orientation, monotone early breaks, and retained
   equality cases in both finite recursions; and
4. the literal `38`-body membership and digest `(1)`.

It found no mathematical defect.  Its only wording correction was adopted:
the `343` count is the number of exact literal obligations after safe screens,
not the larger scalar-danger census.

The theorem closes exactly the stated `38` roots.  It does not promote the
separate ranked-head computation, close the other THM-741 roots, prove global
THM-741, or prove LRC(14).  Its reusable mechanism is the two-comb residual
carrier: a scalar rank can be impossible while the exact pair union still
lies uniformly below its projective limit threshold. ∎
