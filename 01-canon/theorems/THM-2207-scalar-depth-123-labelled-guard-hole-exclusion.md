---
id: THM-2207
title: "Scalar depth-(1,2,3) labelled guard-hole exclusion"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch, the
  actual blocker valuation profile (1,2,3) is empty. On the primitive
  13^4 guard-safe annulus, an exact labelled branch certificate determines
  the five largest conditional capacities for all 1,014*78 typed shallow
  pairs and all 13,182 unit sign classes. The unique minimum deficit is
  1,608. A second one-owner Ky-Fan certificate closes 79,091 of the 79,092
  pairs and reduces its sole exception to a direct margin of 1,676. Direct
  full-residue controls reproduce both named rows. Together with THM-2198
  and THM-2204/2205, this excludes every unique-deepest profile through
  depth three. THM-2203 bounds every remaining scalar deepest valuation by
  4..19, leaving 1,136 finite profiles. None is settled here, and this is
  not a proof of LRC(14).
source: codex-klein-2026-07-24-scalar-depth123-labelled-capacity
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2204-scalar-depth-223-thirteen-lift-capacity-law
  - THM-2205-scalar-depth-113-exact-lift-capacity-exclusion
script: 04-computation/lrc14_scalar_depth123_branch_certificate_thm2207.py
output: 05-knowledge/results/lrc14_scalar_depth123_branch_certificate_thm2207.out
script_sha256: 8e5bd4bb0f6e9072fd9f1992a0d2e454d63e92f476395ca021f773fcf811f2df
output_sha256: 6710488fdb3b9b1f0a31f8b8ed771f079e75b4e7dadfb79d8863a7d929df3890
hash_basis: working-tree bytes (LF)
---

# THM-2207 -- scalar depth-`(1,2,3)` exclusion

Put

```text
D_a={t in R/Z:||at||<1/14},
C_H={t in R/Z:||Ht||>1/7}.                           (1)
```

In the scalar `5+3` branch of THM-2192 and THM-2198, suppose

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere. The coefficients `H,q_1,...,q_5` are positive
thirteen-units, the three actual blockers are positive multiples of
thirteen, and after relabelling

```text
1<=lambda_1<=lambda_2<lambda_3,
lambda_j=nu_13(c_j).                                 (3)
```

This theorem excludes the last profile of deepest valuation three:

```text
(lambda_1,lambda_2,lambda_3)=(1,2,3).                (4)
```

The finite carrier is the mixed-depth labelled root-capacity vector. It
retains the correlation between each unit label, each guard hole, and the
two different blocker scales. No thirteen-lift family is replaced by its
average.

## 1. The mixed primitive layer

Assume (4) and put

```text
N=13^4=28561,                    Q=13^3=2197.         (5)
```

Multiplication of primitive numerators by `H modulo N` normalizes the guard
to one and replaces each terminal coefficient by its product with
`H^(-1) modulo N`. This bijection preserves all three thirteen-valuations.

Define

```text
U_N={z mod N:13 does not divide z and 7||z||_N>N},
||z||_N=min(z mod N,-z mod N).                       (6)
```

The exact cardinality is

```text
|U_N|=18830.                                         (7)
```

The depth-three blocker is safe throughout `U_N`. Its normalized
coefficient is `13^3w` with `w` a unit modulo thirteen, and at a primitive
numerator

```text
13^3w*z/13^4=wz/13
```

is a nonzero thirteenth root. Its norm is at least `1/13>1/14`. The exact
companions check this directly for all six unit sign classes modulo
thirteen.

Every primitive phase `r modulo Q` has thirteen primitive roots

```text
z=r+kQ modulo N,                    k in F_13.        (8)
```

There are `phi(Q)=2028` such phases. Let

```text
B(r)={k in F_13:r+kQ belongs to C_1},
h(r)=|B(r)|.                                         (9)
```

Exact root counting gives

```text
h(r) in {9,10},
#{r:h(r)=9}=1450,               #{r:h(r)=10}=578,
sum_r h(r)=18830.                                    (10)
```

A depth-one coefficient has the form `13u`. Its danger bit is constant on
the roots (8):

```text
A_u(r)=1_[14||ur||_Q<Q].                             (11)
```

The unit part `u modulo Q`, up to sign, ranges through

```text
phi(Q)/2=1014                                        (12)
```

classes. A depth-two coefficient has the form `13^2v`; its phase bit is

```text
B_v(r)=1_[14||v(r mod 13^2)||_(13^2)<13^2],          (13)
```

where `v modulo 13^2`, up to sign, ranges through

```text
phi(13^2)/2=78.                                      (14)
```

The two blocker types cannot be interchanged. Their complete typed pair
universe therefore has

```text
1014*78=79092                                        (15)
```

rows.

For a unit sign class `q modulo N`, define the guard-surviving root count

```text
w_q(r)=#{k in B(r):14||q(r+kQ)||_N<N}.               (16)
```

The root-window law of THM-2198 gives

```text
w_q(r) in {0,1,2}.                                   (17)
```

There are `phi(N)/2=13182` unit sign classes. The two threshold sets

```text
E_q^1={r:w_q(r)>=1},             E_q^2={r:w_q(r)=2}  (18)
```

therefore determine every conditional capacity.

For a mixed blocker pair put

```text
P_(u,v)={r:A_u(r)=B_v(r)=0}.                         (19)
```

Its full residual size and a unit class's conditional capacity are

```text
R(u,v)=sum_(r in P_(u,v))h(r),

C_q(u,v)=sum_(r in P_(u,v))w_q(r)
        =|P_(u,v) intersection E_q^1|
         +|P_(u,v) intersection E_q^2|.              (20)
```

This is an exact disintegration of the full torsion masks.

## 2. The exact branch certificate

For a fixed unit class put

```text
F_q=sum_r w_q(r),             X_q(S)=sum_(r in S)w_q(r).
```

With `A={r:A_u(r)=1}` and `B={r:B_v(r)=1}`,
inclusion--exclusion gives

```text
C_q(u,v)
 =F_q-X_q(A)-X_q(B)+X_q(A intersection B).           (21)
```

Since `w_q(r)<=2`, and the residual lies in each one-owner complement,

```text
C_q(u,v)<=U_q(u,v),

U_q(u,v)=min(
 F_q-X_q(A),
 F_q-X_q(B),
 F_q-X_q(A)-X_q(B)+2|A intersection B|
).                                                   (22)
```

For each typed pair, the certificate visits unit classes in decreasing
order of the exact integer bound (22), evaluates capacities from the two
bitsets (18), and retains the five largest `(capacity,-label)` pairs. It
stops only when the next bound is strictly smaller than the fifth retained
capacity. All cutoff ties are therefore evaluated.

Across all `79092` pairs this requires

```text
940857 exact candidate evaluations,
average 11.895729024427 per pair,
maximum 27 on any pair.                              (23)
```

No pruning decision uses floating point, an unchecked integer width, or an
optimization-sensitive assertion.

### General owner-overlap lemma

The bound is not special to two owners. For arbitrary active phase masks
`A_1,...,A_k`, put

```text
m(r)=#{i:r in A_i},             R=(union_i A_i)^c.
```

The pointwise identity

```text
1_(m=0)=1-m+(m-1)_+
```

gives the exact correction

```text
sum_(r in R)w_q(r)
 =F_q-sum_i X_q(A_i)+sum_r(m(r)-1)_+ w_q(r).         (24)
```

Whenever `w_q<=2`,

```text
sum_(r in R)w_q(r)
 <=F_q-sum_i X_q(A_i)+2sum_r(m(r)-1)_+.              (25)
```

It is also bounded by every one-owner complement capacity. Thus owner
overlap multiplicity plus the two thresholds (18) gives a rigorous
candidate cutoff for later rooted problems. It does not identify which
labels attain that cutoff.

## 3. Exact deficit and hostile row

Order the capacities decreasingly:

```text
C_(1)(u,v)>=...>=C_(13182)(u,v).
```

The exhaustive result is

```text
R(u,v)-sum_(i=1)^5 C_(i)(u,v)>=1608                 (26)
```

on the full primitive annulus for every typed pair. The minimum is unique
in the least-positive sign convention:

```text
(u,v)=(799,46),                  R(u,v)=13526,

((C_(i),q_i))_(i=1)^5
 =((2604,5193),(2472,10386),(2292,7773),
   (2288,10388),(2262,7775)).                        (27)
```

The five capacities sum to `11918`, leaving

```text
13526-11918=1608.                                    (28)
```

Two digests freeze the complete search:

```text
pair table:
79b9b75f3732e47b43c2bba726906250ce0c1069b90d8952c33caa8a8364570f

branch trace:
2699c79c62d4c0e5805529b26dc4ce7ae5ed1f0146be25714b898ea42176802a.
                                                               (29)
```

The primary companion independently enumerates all `18830` full torsion
points and all `13182` unit classes on the hostile row, reproducing
(27)--(28). It also checks direct/image parity for every shallow class and
THM-2204's lift-family sum on all `1014` hostile-residual families.

## 4. Independent one-owner certificate

A second exact companion gives a cheaper proof with a different
intermediate object. Define

```text
C_q^1(u)=sum_(r:A_u(r)=0)w_q(r),
C_q^2(v)=sum_(r:B_v(r)=0)w_q(r),                    (30)
```

and let `K_1(u),K_2(v)` be the sums of their five largest entries. Residual
containment gives Ky-Fan monotonicity:

```text
sum_(i=1)^5 C_(i)(u,v)<=min(K_1(u),K_2(v)).           (31)
```

The one-owner scan proves

```text
R(u,v)-min(K_1(u),K_2(v))>0
```

for `79091` of the `79092` pairs. The least positive coarse margin is `36`
at

```text
(u,v)=(1,1),             R=13310,
min(K_1,K_2)=13274.                                  (32)
```

There is exactly one nonpositive coarse row:

```text
(u,v)=(1098,84),         R=13310,
min(K_1,K_2)=13328,          coarse margin=-18.       (33)
```

Its exact top five are

```text
((2396,6),(2330,14275),(2324,14278),
 (2296,5),(2288,12)),                                 (34)
```

with sum `11634` and exact margin

```text
13310-11634=1676.                                    (35)
```

Direct strict inequalities on the full residue universe independently
reproduce both (27) and (33)--(35). This companion is

```text
04-computation/lrc14_scalar_depth123_labelled_guard_hole_thm2207.py
05-knowledge/results/lrc14_scalar_depth123_labelled_guard_hole_thm2207.out

LF SHA256:
3c32b7c55197b5ff96ef892df1b29b61fc79df71494d1980decb34bc870032d8
41dd14abaed66e37719eaefc216646082937aa70a76781e154dbed078604e1b3.
                                                               (36)
```

An independent referee reconstructed the one-owner proof from definitions,
audited the primary pruning rule and tie handling, recomputed the hostile
row, and replayed both normal and optimized companions. The verdict was
`ACCEPT / PROVED + VERIFIED-EXACT`.

## 5. Exclusion of `(1,2,3)`

If (2) held, the depth-three blocker would contribute nothing on `U_N`.
For the actual shallow pair, the five unit masks would have to cover every
point counted by `R(u,v)`. Duplicate unit sign classes are redundant in a
union, and the union of at most five masks has size at most the sum of the
five largest individual capacities. Equation (26) says that this sum is
strictly smaller than the residual.

Hence one primitive residue is safe from all five unit masks and both
shallow blockers. It is safe from the depth-three blocker by Section 1 and
lies in the guard-danger set by construction. Powers of thirteen are
coprime to seven and fourteen, so no relevant torsion equality occurs. The
missed point satisfies every inequality strictly and thickens to an
uncovered open interval, contradicting the almost-everywhere cover (2).
This proves that (4) is empty.

## 6. Consequence and structural boundary

The unique-deepest ordering (3) leaves only

```text
depth 2: (1,1,2);
depth 3: (1,1,3), (1,2,3), (2,2,3).                 (37)
```

THM-2198 excludes the first, THM-2205 excludes `(1,1,3)`, this theorem
excludes `(1,2,3)`, and THM-2204 excludes `(2,2,3)`. Therefore every
surviving scalar `5+3` branch has

```text
4<=lambda_3<=19,                                     (38)
```

where the upper bound is THM-2203. That theorem counts `1140` profiles
before the four exclusions, so `1136` remain.

The exact carrier ledger is

```text
source:
  one depth-one owner, one depth-two owner, and a unit label;
map:
  disintegrate the depth-four torsion layer into thirteen-root phases;
target:
  owner phase masks and threshold sets E_q^1,E_q^2;
preserved:
  exact residual size and every labelled conditional unit capacity;
destroyed by scalar averaging:
  phase-aligned guard-hole correlation and the top-five order statistic;
cheapest decisive test:
  the sole one-owner exception (1098,84).             (39)
```

THM-2204's lift-family sum is only the zeroth correlation moment. The
hostile top five lie in five different lift families. THM-2201's full
triangular Hasse jets retain more information than (18), but the present
profile needs only the two thresholds because every root count is at most
two.

The remaining depth-`4..19` ledger is finite but not feasibly flat: its
largest torsion layer and owner-current sidecars remain enormous. This
theorem gives no uniform higher-depth top-five recurrence, does not empty
the primitive rank-twelve locus, and does not prove LRC(14). QED.
