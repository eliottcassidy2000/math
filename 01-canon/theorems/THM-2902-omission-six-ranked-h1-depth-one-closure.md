---
id: THM-2902
title: "Omission-six ranked H1 depth-one closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  The two one-hard omission-six bodies
  (1,2,3,4,5,10,12) and (1,2,3,4,5,12,13) close by the ranked r=4
  complement cap, globally sealed finite H1 cores of sizes 49 and 6, and
  exact earliest-label depth-one certificates.
source: codex-lrc-rank-impossible-overlap-2026-07-29
depends_on:
  - THM-2892-eight-body-five-slot-heavy-triangle-closure
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2899-all-root-ranked-suffix-scalar-census
related:
  - THM-2901-all-hard-exact-global-pair-cap-and-route-partition
verification:
  - 04-computation/lrc14_j6_omission6_ranked_h1_depth1_closure_thm2902.py
  - 05-knowledge/results/lrc14_j6_omission6_ranked_h1_depth1_closure_thm2902.out
---

# THM-2902 -- omission-six ranked H1 depth-one closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Statement

Let

```text
E_1=(1,2,3,4,5,10,12),
E_2=(1,2,3,4,5,12,13).                                  (1)
```

For `i=1,2`, no six distinct external speeds `w>=15` have danger combs
whose union covers `G_(E_i)`.  Together with THM-2892's eight-body chamber,
each body in `(1)` is therefore a whole seven-body/six-slot root closure.

Both roots contain the five-label low prefix `{1,2,3,4,5}` and omit `6`.
They are the two ranked-`r=4` one-hard targets on the omission-six boundary
isolated after THM-2899.

## 2. Unique scalar-hard suffixes

THM-2896 gives `E_1` and `E_2` ordered gates of sizes `10` and `6`.
THM-2899 closes every marked suffix except the rank-one apex `22` on each
root.  Put

```text
C_i=G_(E_i) minus D_22,             h_i=|C_i|.            (2)
```

The excluded prefix is exactly `(22)`.  The exact hard-carrier data are

```text
                 h             components       q1+q2+q3+q4
E_1          689/2310              28              18371/72765
E_2         4127/15015             22           3380627/15135120. (3)
```

For each carrier, write

```text
R_4=q_1+q_2+q_3+q_4,       lambda=h-R_4.                  (4)
```

The strict ranked-complement margins are

```text
6h_1/7-R_4 = 232/72765,
6h_2/7-R_4 = 26443/2162160.                              (5)
```

Thus THM-2893/2897 with four complementary singleton slots forces all five
labels of every hypothetical cover into

```text
H_1(C_i)={w allowed:c_(C_i)(w)>=lambda_i}.                (6)
```

## 3. Global H1 cores

If `C_i` has `r_i` components, THM-735(ii), through the hash-pinned
THM-2899 engine, gives

```text
c_(C_i)(w)<h_i/7+(99/70)r_i/(7w).                        (7)
```

Since `lambda_i-h_i/7` equals the positive margin in `(5)`, every member
of `(6)` satisfies

```text
w<=ceil((99/70)r_i/(7(6h_i/7-R_4)))-1.                   (8)
```

The two exact cutoffs and globally sealed cores are

```text
E_1: cutoff 1774, |H_1|=49,
  (26,39,27,28,16,78,18,65,32,54,21,52,56,40,69,63,130,
   91,46,120,108,42,45,41,68,25,286,220,81,92,156,97,107,
   82,125,182,264,275,121,260,250,195,104,242,162,175,50,
   23,237);

E_2: cutoff 363, |H_1|=6,
  (27,28,18,16,21,40).                                  (9)
```

The order in `(9)` is decreasing literal coverage on `C_i`, with speed as
the tie-breaker.  The verifier scans every allowed speed through the
cutoff and independently recomputes all `2107` vector entries by the
scalar coverage primitive.  Strictness in `(7)` excludes every omitted
speed from `(6)`.

## 4. Earliest-label depth-one closure

Give a hypothetical five-set in `H_1(C_i)` to its unique earliest label
`x` in the order `(9)`.  The other four labels lie in the strict suffix.
On

```text
C_(i,x)=C_i minus D_x
```

the union of those four labels is at most the sum of their four individual
coverages.  The parent-carrier coverages certify `43` earliest labels for
`E_1` and one for `E_2`; the last four labels on each core are impossible
as earliest choices because their suffix has size below four.

Only three earliest labels require literal residual recomputation.  Their
exact local certificates are

```text
E_1, x=26:
  margin 1957/116424,
  top four 28:11/196, 27:7538/135135,
           16:169/3080, 18:50/1001;

E_1, x=39:
  margin 64553/2522520,
  top four 27:2696/45045, 18:311/6006,
           28:379/7644, 40:5773/120120;

E_2, x=27:
  margin 228223/15135120,
  top four 28:71/1323, 21:9668/189189,
           40:4091/83160, 16:1871/39312.                 (10)
```

All three margins are strict.  Across parent and local certificates, the
smallest positive depth-one margin is

```text
E_1: 7949/3027024,             E_2: 22991/5045040.       (11)
```

Hence no five labels cover either carrier in `(2)`.  The unique hard suffix
on each root closes; all later suffixes were already scalar-closed by
THM-2899, proving the external statement in Section 1.

## 5. Root recomposition and scope

The two bodies in `(1)` are disjoint from the four THM-2895 roots, the
THM-2898 root, the five THM-2899 roots, and the five THM-2901 roots.
Consequently the proved union now contains `17` roots, leaving

```text
3432-17=3415                                             (12)
```

in the official seven-body residual.

This theorem proves only the two named roots.  It does not promote the
partial `9/32` all-root H1 scout, claim that every ranked-`r=4` branch has
a depth-one closure, close the seven-body rung, or prove LRC(14).

## 6. Verification

The verifier hash-pins the promoted THM-2899 engine, reconstructs both
complete gates and their unique hard suffixes, independently reconstructs
both literal carriers, globally seals `(9)`, and locks every parent/local
earliest-label certificate.  It contains no Python `assert`.  Ordinary and
optimized outputs are byte-identical.  The semantic digest is

```text
7b173f9c5817ab6e73ddabba48bddf0f7d4714f7384ce034b23911e2056b5e96.
```

Canonical artifacts:

```text
04-computation/lrc14_j6_omission6_ranked_h1_depth1_closure_thm2902.py
SHA-256 91bb54803bddd72b516b8c21dc51d3bd3f38eb6be6e2e0beade98b6676012f3f

05-knowledge/results/lrc14_j6_omission6_ranked_h1_depth1_closure_thm2902.out
SHA-256 cff46490d4a904947ec0fbe0cedfa59484c6b7974923656ba2459a55781192d7
```
