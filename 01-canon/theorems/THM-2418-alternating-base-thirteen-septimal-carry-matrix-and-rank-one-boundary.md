---
id: THM-2418
title: "Alternating base-thirteen septimal carry matrix and rank-one boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. For
  R=13^k, the physical-to-terminal seven-root map is the affine
  permutation l=floor(Ry)+(-1)^k r modulo seven, and the carry is the
  alternating sum of the first k base-thirteen digits. Haar averaging
  gives an explicit seven-state kernel K_k=P^k of rank seven; its six
  charged singular values are exactly 1/R. A terminal-only rational
  cylinder profile is nonflat exactly when every six nonzero septimal
  colours survives, with attenuation 1/R. This does not close the real
  word problem: a flat centred-comb terminal word is exact, a
  one-cylinder word makes the terminal matrix rank one, and for every k
  an even centred BV-two source set of mass at least 7/13 makes the
  complete source-weighted carry matrix rank one. No canonical
  THM-2305 source correlation, all-91-unit address, row exclusion, or
  LRC(14) conclusion is proved.
source: codex-2026-07-26-septimal-digit-cocycle
depends_on: []
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2409-unfiltered-septimal-source-completion-and-word-phase-boundary
  - THM-2414-thirteen-skew-septimal-word-transport-and-local-stopping-atlas
script: 04-computation/lrc14_septimal_carry_matrix_thm2418.py
output: 05-knowledge/results/lrc14_septimal_carry_matrix_thm2418.out
script_sha256: 697e6f6e5de8536efa563f1ff123bcaedea0da9914108efe57ff6bb6f2a59096
output_sha256: d9daa1a359426a3ddd697c05d5d2d0d408af707aaee93f02f61c30d928281193
hash_basis: working-tree bytes (LF)
---

# THM-2418 -- the septimal word phase is an alternating digit carry

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2409 isolates the real-word obstruction: modulo thirteen the
transported word is neutral, while modulo seven the same coordinate
translation moves it. Candidate THM-2414 finds the corresponding
one-step affine digit in a last-lane local atlas. The complete
base-thirteen iteration is elementary and exact:

```text
physical septimal root
  -> alternating base-thirteen carry
  -> full-rank seven-state Haar kernel
  -> charged singular value 13^(-k),

but

source-prefix weighting
  -> exact rank-one collapse is possible.                       (1)
```

The abstract theorem below is independent of the candidate status of
THM-2414. Its application boundary is the main result.

## 1. The exact affine carry cocycle

For a seven-unit integer `A`, width `a in {1,2}`, and real base `y`,
put

```text
W_(A,a)(y)
 ={r in F_7: ||A(y+r)/7||<a/14}.                               (2)
```

Fix `k>=1` and write

```text
R=13^k,              epsilon=(-1)^k,

z={Ry},              M_k(y)=floor(Ry).                         (3)
```

For an unreduced lift, multiplication by `R` sends the physical root
`r` to

```text
l=Rr=epsilon r                         modulo seven.           (4)
```

Reducing `Ry=z+M_k` translates the terminal root by `M_k`. Therefore

```text
W_(RA,a)(y)
 =epsilon (W_(A,a)(z)-M_k(y))                 in F_7.           (5)
```

Here affine operations act elementwise on the root set. Equivalently,
a physical root `r` is transported to

```text
l=M_k(y)+epsilon r                         modulo seven.        (6)
```

Write the first `k` base-thirteen digits of `y` as

```text
d_j=floor(13{13^(j-1)y}),              1<=j<=k.                (7)
```

The ordinary prefix integer and its septimal reduction are

```text
M_k
 =sum_(j=1)^k d_j 13^(k-j),

M_k
 =sum_(j=1)^k d_j (-1)^(k-j)              modulo seven.        (8)
```

Thus the word phase is not an independent colour. It is the alternating
digit carry of the real base transition.

For `c in F_7`, define the pointwise affine permutation matrix

```text
A_(c,epsilon)(r,l)
 =1_(l=c+epsilon r).                                          (9)
```

Equation (6) says that every individual prefix cylinder acts by exactly
one such permutation.

## 2. The cheapest exact seven-state transfer

Disintegrate

```text
y=(z+n)/R,                 0<=n<R,              0<=z<1.       (10)
```

Then `M_k(y)=n`. Haar averaging over the `R` prefix cylinders gives

```text
K_k(r,l)
 =1/R #{0<=n<R: n=l-epsilon r mod 7}.                          (11)
```

The digit recurrence

```text
c_(j+1)=-c_j+d_(j+1)                         modulo seven     (12)
```

shows that

```text
K_k=P^k,                                                       (13)
```

where the one-digit matrix is

```text
             1
P=          --- *
             13

  [2 2 2 2 2 2 1
   2 2 2 2 2 1 2
   2 2 2 2 1 2 2
   2 2 2 1 2 2 2
   2 2 1 2 2 2 2
   2 1 2 2 2 2 2
   1 2 2 2 2 2 2].                                           (14)
```

Let

```text
Pi=J/7,

A(r,l)=1_(r+l=6 mod 7),                                      (15)
```

where `J` is the all-one matrix. Since `A^2=I` and `A Pi=Pi`,
equations (11)--(14) give

```text
K_k
 =Pi+R^(-1)(I-Pi),                  k even,

K_k
 =Pi-R^(-1)(A-Pi),                  k odd.                    (16)
```

Consequently,

```text
K_k K_k^T
 =Pi+R^(-2)(I-Pi).                                          (17)
```

Thus `K_k` has rank seven. Its constant singular value is one, and all
six charged singular values are exactly

```text
1/R=13^(-k).                                                  (18)
```

The same fact appears in the raw digit Fourier sum. For
`e in F_7^*`,

```text
1/R sum_(n=0)^(R-1) zeta_7^(e n)
 =
  1/R,                         k even,

 -zeta_7^(-e)/R,               k odd.                         (19)
```

In particular, every raw carry colour is nonzero, but its size is
already only `1/R`.

## 3. Terminal-cylinder conditioning

Let `Q:T->{0,1}` be a positive-mass rational finite step word. Define
its seven terminal-cylinder masses by

```text
q_l
 =integral_0^1 Q((z+l)/7)dz
 =7 mu(Q intersection [l/7,(l+1)/7)).                         (20)
```

Then `q in Q_+^7` and

```text
sum_l q_l=7mu(Q).                                             (21)
```

In the **terminal-only** quotient--Haar prefix weight, with no
source-owner factor--the exact filtered transfer and its source-row
survival vector are

```text
T_Q=K_k diag(q),

M=K_k q.                                                       (22)
```

Put `h=floor(R/7)`. Every entry of `K_k` is at least `h/R`, so

```text
M_r
 >=7h mu(Q)/R
 >0                              for every r.                  (23)
```

Since `K_k` is invertible,

```text
rank(T_Q)=#support(q),                                        (24)

M is nonconstant iff q is nonconstant.                        (25)
```

For the normalized septimal transform

```text
qhat(e)=1/7 sum_l q_l zeta_7^(e l),                            (26)
```

equation (16) gives, up to reflection and a unit phase,

```text
|Mhat(e)|=|qhat(epsilon e)|/R,            e!=0.               (27)
```

Rationality makes (25) an all-colour statement. If one
`qhat(e)=0` with `e!=0`, the degree-at-most-six polynomial

```text
sum_(l=0)^6 q_l X^l
```

is divisible by `Phi_7`, so all `q_l` are equal. Hence

```text
q nonconstant
  iff M nonconstant
  iff Mhat(e)!=0 for every e in F_7^*.                         (28)
```

This is the exact positive survivor. For example, any positive terminal
word missing one complete seventh-cylinder fires all six terminal-only
carry colours.

## 4. Two sharp terminal boundaries

Write

```text
D_v={x in T:||vx||<1/14}.                                    (29)
```

First take

```text
Q=D_7.
```

Each seventh-cylinder contains mass `1/49`, so

```text
q=(1/7,...,1/7),

M=(1/7,...,1/7)                                               (30)
```

for every `k`. Thus a nonconstant centred comb word can have a perfectly
flat terminal-cylinder profile and no charged terminal-only carry mode.

At the opposite boundary, take

```text
Q=D_2 intersection D_1^c
  =(13/28,15/28).                                             (31)
```

This interval lies entirely in terminal cylinder `l=3`, and

```text
q=(1/2)e_3,

rank(T_Q)=1.                                                   (32)
```

Every normalized surviving row is the same point mass `e_3`. At `k=1`,
the unnormalized row masses are

```text
1/13,1/13,1/13,1/26,1/13,1/13,1/13.                          (33)
```

Thus full raw carry rank does not prevent a terminal word from reducing
the filtered transfer to rank one.

## 5. A real, even, large-mass source rank-one hostile

Terminal conditioning is not the main obstruction. The actual THM-2305
current also contains a source-owner weight before the terminal word.
Carry state alone cannot control that weight.

For every `k`, define a rational prefix-cylinder set `S_k` as follows.

If `k` is even, then `R=1 mod 14`. Remove the reflection-fixed central
prefix cylinder:

```text
S_k
 =T minus [(R-1)/(2R),(R+1)/(2R)).                            (34)
```

If `k` is odd, then `R=13 mod 14`. Remove the first and last three
prefix cylinders:

```text
S_k=[3/R,1-3/R).                                              (35)
```

In both cases `1_(S_k)` is real and even, has exactly two circular
jumps, and

```text
mu(S_k)
 =
  1-1/R,                         k even,

  1-6/R,                         k odd,

mu(S_k)>=7/13.                                               (36)
```

Among the retained prefix indices, every carry residue occurs exactly
`h=floor(R/7)` times. Therefore the source-weighted transfer is

```text
K_(S_k)=(h/R)J.                                               (37)
```

It has rank one and annihilates the whole six-dimensional charged
subspace. With an arbitrary terminal profile,

```text
T_(S_k,Q)=(h/R)J diag(q),                                    (38)
```

still has rank one; every normalized source row is identical.

This hostile is stronger than a small exceptional cylinder. Its mass is
uniformly positive, its circular variation is exactly two, and it is
centred/even. It also explains why a BV or Perron approximation does not
repair the carry colour. The unweighted signal in (19) is `1/R`, while
deleting only one or six prefix cylinders--mass `O(1/R)`--kills that
signal exactly. Signal and admissible BV error live at the same scale.

The sets `S_k` depend on the selected clock. They are formal rational
source-prefix hostiles, not asserted to be canonical exclusive-owner
sets.

## 6. Canonical scope and the missing sidecar

For a THM-2305 terminal stratum `Q_(j,sigma)`, positivity proves only

```text
sum_l q_l>0.                                                   (39)
```

It does not prove that its seven cylinder masses are nonconstant.
Candidate THM-2414 supplies one strict local `W=8` stopping atlas, but
explicitly does not claim that atlas is a scalar cover or a
positive-mass global terminal-cylinder profile.

More importantly, the genuine THM-2305 current contains `G_j` or `E_j`.
That source weight depends on the full prefix `n`, tail `z`, and root
label `r`; it need not factor through (22). The rank-one construction
(34)--(38) proves that positivity, evenness, bounded variation, and
large mass do not replace the missing source--terminal correlation.

THM-2409 proves an all-colour theorem for `Q=1` because its source
profile has both a nonzero total and an anchored zero. With a real
transported word, the lawful co-shift can make every inserted-source
term vanish. The carry matrix supplies no nonzero total for that
co-shifted coefficient profile.

The missing sidecar is therefore one of:

1. a canonical source--terminal correlation excluding rank one;
2. a fixed-cylinder imbalance for the actual THM-2305 word together
   with a source-factor intertwiner; or
3. a common phase/endpoint mechanism that works at the natural `1/R`
   charged scale.

No canonical source current is completed, no all-`91`-unit relation
address is produced, no row is excluded, the ledger remains `165`, and
LRC(14) remains open.

## 7. Exact companion

The dependency-free exact companion:

- verifies (5)--(8) on `18,040` rational word instances;
- exhausts every base-thirteen digit word through depth five;
- checks `P^k`, the direct prefix counts, (16)--(19), and the Gram
  identity through depth eight;
- exhausts all `127` nonempty Boolean terminal-cylinder supports;
- checks the flat and one-cylinder terminal hostiles; and
- verifies the even centred rank-one source constructions, exact carry
  counts, variation, and mass floor through depth eight.

Run:

```bash
python3 04-computation/lrc14_septimal_carry_matrix_thm2418.py
python3 -O 04-computation/lrc14_septimal_carry_matrix_thm2418.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_septimal_carry_matrix_thm2418.out
```

Every truth-bearing finite check raises explicitly, so optimized mode
executes the same audit.
