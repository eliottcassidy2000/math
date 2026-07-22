---
id: THM-2086
title: "Relative-Hunter Fourier channel decomposition and a lacunary rank-seven cone"
status: >
  PROVED. Every guard-restricted THM-2081 pair edge is the sum of five
  sevenths of its global pair overlap, two explicit mixed-fold corrections,
  and one absolutely convergent genuine three-frequency relation channel.
  This gives an exact spanning-tree margin identity and isolates the
  all-height residual as a short-relation-lattice sum. Independently, a BV
  mixing estimate and the sharp odd-guard mixed-overlap spectrum prove that
  rank-seven guard containment is impossible when the largest speed B obeys
  sum(Q minus {B})+6h<(17/1078)B. The same channel split closes every
  divisor-complete hereditary rank-seven terminal with 7 dividing h, by an
  exact bipartite-tree margin of at least 5/294. When 7 does not divide h it
  also closes the maximal profile having five 7-divisible speeds, with margin
  at least 23/1911. These are genuine all-height branches; they do not close
  the remaining relation templates or LRC(14).
source: codex-2026-07-22-LRC-relative-Hunter-Fourier
depends_on:
  - THM-1221
  - THM-1234
  - THM-2080
  - THM-2081
  - THM-2083
related:
  - THM-2082
  - THM-2087
  - THM-2088
script: 04-computation/lrc14_relative_hunter_fourier_lacunary_codex_20260722.py
output: 05-knowledge/results/lrc14_relative_hunter_fourier_lacunary_codex_20260722.out
script_sha256: 67fbb82f4aa07a4835e8c9673c8bbd0698e8f9330d04aae2cc4dd5a78a8c1426
output_sha256: dbdc5e5742c819150e41e4b85931fb4849878095f6625a671745174ccd2f7767
hash_basis: working-tree bytes (LF)
---

# THM-2086 -- the relative-Hunter Fourier channel and lacunary cone

Put

```text
D_q={t in R/Z:||qt||<1/14},
E_h={t in R/Z:||ht||<1/7},       C_h=(R/Z) minus E_h.  (1)
```

For distinct positive speeds `p,q`, retain THM-2081's quantities

```text
I_p=measure(D_p intersect E_h)=2/49+epsilon_p,
rho_pq=measure(D_p intersect D_q),
w_pq=measure(D_p intersect D_q intersect C_h).         (2)
```

The first part of the theorem identifies exactly what the rank-one code wheel
of THM-2082 forgot: the genuine three-frequency relation channel. The second
part uses a different old tool, bounded-variation mixing, to close an
unbounded cone without classifying that channel term by term.

## 1. Absolutely convergent three-frequency formula

Use the Fourier convention

```text
hat f(n)=integral_(R/Z) f(t) exp(-2 pi i n t) dt.       (3)
```

For the base danger and guard-complement indicators, set

```text
s_0=1/7,       s_n=sin(pi n/7)/(pi n),                 n!=0,
u_0=5/7,       u_n=-sin(2 pi n/7)/(pi n),              n!=0. (4)
```

Then

```text
w_pq=sum_(a p+b q+c h=0) s_a s_b u_c.                 (5)
```

The series in (5) is absolutely convergent. Its part with all three indices
nonzero is

```text
R_h(p,q)=sum_(a p+b q+c h=0; abc!=0) s_a s_b u_c.      (6)
```

Consequently the exact edge decomposition is

```text
w_pq=(5/7)rho_pq-(epsilon_p+epsilon_q)/7+R_h(p,q).     (7)
```

### Proof of convergence and (5)

Take Fejer means of the three base indicators. They remain between zero and
one, converge in `L2`, and their triple products converge in `L1`. Integrating
the finite Fourier products leaves exactly the characters satisfying
`ap+bq+ch=0`.

It remains to justify passage to the infinite coefficient sum. The origin is
finite. On any one-zero axis the relation lattice is primitive one-dimensional
and the terms are `O(k^-2)`. On the genuine locus, (4) gives

```text
|s_a s_b u_c|<=1/(pi^3 |abc|),
c=-(pa+qb)/h.                                          (8)
```

Thus it suffices to prove

```text
sum_(a,b!=0;pa+qb!=0) 1/(|ab(pa+qb)|)<infinity.        (9)
```

When `a,b` have the same sign,
`pa+qb>=2 sqrt(pqab)`, and the sum is bounded by a constant times
`zeta(3/2)^2`. For opposite signs write `b=-v`, with `a,v>0`. In the ranges
`pa>=2qv` and `qv>=2pa`, the larger term controls the difference and the sum
is bounded by `sum log(a)/a^2` or its symmetric copy. In the remaining range
`v` is comparable to `a`; summing `1/|pa-qv|` along the relevant arithmetic
progression costs `O(log a)`, again leaving `sum log(a)/a^2`. This proves
(9), hence absolute convergence. Dominated convergence on the Fejer weights
proves (5).

### Proof of the split

The origin together with the `c=0` channel contributes

```text
(5/7)rho_pq.                                           (10)
```

On the `b=0` axis, excluding the origin, the contribution is

```text
(1/7)[measure(D_p intersect C_h)-5/49]
=(1/7)[1/7-I_p-5/49]
=-epsilon_p/7.                                         (11)
```

The `a=0` axis similarly contributes `-epsilon_q/7`. The remaining terms are
(6), proving (7). QED.

## 2. Exact spanning-tree margin identity

For a labelled seven-set `Q={q_1,...,q_7}` and any spanning tree `T`, let
`d_T(i)` be the degree of vertex `i` and put

```text
tau_T=sum_(ij in T)w_ij,
Delta=2/7-sum_i I_i=-sum_i epsilon_i.                  (12)
```

Summing (7) over the six tree edges gives

```text
tau_T-Delta
 =(5/7)sum_(ij in T)rho_ij
  +sum_i(1-d_T(i)/7)epsilon_i
  +sum_(ij in T)R_h(q_i,q_j).                          (13)
```

Since THM-2081's `tau_h(Q)` is the maximum of `tau_T`, any tree for which the
right side of (13) is positive contradicts `G_Q subset E_h` and supplies that
much safe mass outside the guard.

Equation (13) is the exact all-height ledger. Choosing THM-1221's global
maximum pair tree supplies the uniform positive bulk

```text
(5/7)sum_(ij in T)rho_ij >=(5/7)(15/154)=75/1078.      (14)
```

Only the explicit degree-weighted fold correction and the genuine rank-two
channel sum can spend it. Bounding `R_h` termwise is generally wasteful: a
small resonant triple can make one `R_h(p,q)` strongly negative while its
`epsilon` terms compensate. The combined identity (13), not `R_h` alone, is
the correct classifier.

### Exact hostile boundary

On THM-2081's minimum-margin packet

```text
Q=(1,9,10,11,13,14,24),               h=23,           (15)
```

one maximum restricted tree has speed edges

```text
(1,14),(14,24),(24,10),(24,9),(24,11),(24,13).         (16)
```

Its degree vector in the displayed order of `Q` is

```text
(1,1,1,1,1,2,5).                                      (17)
```

The three terms in (13) are exactly

```text
global-pair bulk        100421/1177176,
degree-fold correction -16117/4512508,
genuine channel sum    -2833331/203062860.             (18)
```

They sum to

```text
561797/8288280>0,                                      (19)
```

the independently atomized THM-2081 margin. The six genuine edge values are

```text
R(1,14)  =-11/2254,       R(14,24)=-11/27048,
R(24,10)=-71/118335,      R(24,9) =-281/284004,
R(24,11)=-7663/2082696,   R(24,13)=-4181/1230684.      (20)
```

The full height-24 replay verifies (13) on all `1,322` scalar survivors and
all `2,982` distinct guard/two-speed edges.

## 3. A BV high-frequency mixing lemma

For arbitrary positive `B,q,h`,

```text
|measure(D_B intersect D_q intersect C_h)
 -(1/7)measure(D_q intersect C_h)|
 <=(q+h)/(3B).                                         (21)
```

### Proof

Let `g=1_(D_q intersect C_h)`. Its jumps are among the `2q+2h` boundary
points of the two combs, so

```text
Var(g)<=2(q+h).                                        (22)
```

The standard BV Fourier estimate and (4) give, for `n!=0`,

```text
|hat g(nB)|<=(q+h)/(pi |n|B),
|s_n|<=1/(pi |n|).                                     (23)
```

Expand `1_(D_B)` in Fejer means and pair with `g`. The zero mode is the term
subtracted in (21), while the nonzero modes have total absolute value at most

```text
(q+h)/(pi^2 B) sum_(n!=0)1/n^2=(q+h)/(3B).             (24)
```

Fejer convergence and the null endpoints complete the proof. QED.

## 4. The odd-guard overlap spectrum

For odd `h` and any positive integer `q`, THM-2080 sharpens as follows:

```text
I(q,h)>=1/35,                                          (25)
```

except for the two possible strict-low ratios

```text
q=6h:       I=1/42,
q=h/11:     I=2/77,                                    (26)
```

where the second row is present only when integral. Among a set of distinct
speeds, each exception occurs at most once.

### Proof

Write the reduced THM-2080 variables as

```text
a=h/gcd(h,q),                    b=q/gcd(h,q).          (27)
```

Here `a` is odd. Its fold bound gives

```text
I(q,h)>=2/49-1/(4ab)>=1/35             when ab>=21.    (28)
```

For coprime `a,b` with odd `a` and `ab<=20`, exact reduction of the fold has
only four rows at or below `1/35`:

```text
(a,b)=(1,5):1/35,       (1,6):1/42,
      (3,5):1/35,       (11,1):2/77.                   (29)
```

The two strict rows translate exactly to (26); the other two satisfy (25)
with equality. QED.

## 5. Two apex-modular branches

### Every apex-divisible guard

Let `Q` be seven distinct positive speeds. Assume it is divisor-complete
through `7` and hereditarily primitive. If `h` is odd and

```text
7 divides h,                                             (29a)
```

then

```text
G_Q is not contained in E_h,
measure(G_Q minus E_h)>=5/294.                           (29b)
```

Thus every rank-seven dyadic-tower terminal obstruction has

```text
7 does not divide h.                                    (29c)
```

### Proof

Partition

```text
L={p in Q:7 does not divide p},
H={q in Q:7 divides q},                 m=|H|.           (29d)
```

Divisor completeness gives `m>=1`. Hereditary primitivity gives `|L|>=2`:
if there were only one nonmultiple of seven, deleting it would leave gcd
divisible by seven. Hence

```text
1<=m<=5.                                                (29e)
```

For `p in L`, reduce the THM-2080 variables for `(p,h)`. The guard numerator
is an odd multiple of seven, so its residue modulo fourteen is seven. In the
fold function,

```text
F(1/2,y)=0,                                             (29f)
```

and therefore

```text
epsilon_p=0.                                            (29g)
```

For a cross pair `p in L`, `q in H`, every nonzero pair relation
`ap+bq=0` forces `7|a`, so `s_a=0`. Only the zero mode survives and

```text
rho_pq=1/49.                                            (29h)
```

Likewise every genuine relation `ap+bq+ch=0` forces `7|a`; hence every term
in (6) vanishes and

```text
R_h(p,q)=0.                                             (29i)
```

The edge split (7) is now exact on every cross edge:

```text
w_pq=5/343-epsilon_q/7.                                 (29j)
```

Choose any spanning tree of the complete bipartite graph `K_(L,H)`. It has
six edges. If `d_q` is the degree of a high vertex `q in H`, then

```text
sum_(q in H)d_q=6,
Delta=-sum_(q in H)epsilon_q.                           (29k)
```

Therefore

```text
tau_T-Delta
 =30/343+sum_(q in H)(1-d_q/7)epsilon_q.                (29l)
```

Every coefficient lies in `[1/7,6/7]` and their sum is `m-6/7`. THM-2080's
sharp floor `I_q>=1/42` says

```text
epsilon_q>=1/42-2/49=-5/294.                            (29m)
```

Using (29e),

```text
tau_T-Delta
 >=30/343-(5/294)(m-6/7)
 >=30/343-(5/294)(5-6/7)
 =5/294.                                                (29n)
```

THM-2081 converts this into (29b). QED.

This closure is genuinely arithmetic, not merely a special case of the
lacunary estimate below. It uses the apex prime seven twice: the mixed fold
axis becomes independent and the sine coefficient kills the entire genuine
cross channel.

### Maximal apex multiplicity when `7` does not divide the guard

Continue with a divisor-complete hereditary seven-set, but now assume

```text
7 does not divide h,
#{q in Q:7 divides q}=5.                                (29o)
```

Then

```text
measure(G_Q minus E_h)>=23/1911>0.                      (29p)
```

For a high speed `q` divisible by seven, the reduced speed residue in
THM-2080 is zero. Thus

```text
epsilon_q=0.                                            (29q)
```

For two high speeds `p,q`, every genuine relation
`ap+bq+ch=0` forces `7|c`. The complement coefficient `u_c` vanishes, so

```text
R_h(p,q)=0,
w_pq=(5/7)rho_pq.                                       (29r)
```

THM-1234 gives total global pair mass at least `44/273` on these five high
vertices. In a uniformly random labelled spanning tree of `K_5`, each edge
occurs with probability `2/5`. Hence some high-vertex tree has

```text
sum rho_pq >=(2/5)(44/273)=88/1365,
sum w_pq   >=88/1911.                                   (29s)
```

Attach the two low vertices by arbitrary nonnegative restricted edges. This
produces a spanning tree on all seven vertices with the same lower bound.
Only the two low vertices contribute to `Delta`; THM-2080 gives each
`epsilon>=-5/294`, and therefore

```text
Delta<=5/147,
tau_h(Q)-Delta>=88/1911-5/147=23/1911.                  (29t)
```

This proves (29p). Combined with the preceding branch and divisor
completeness, the live rank-seven modular profile is now exactly

```text
7 does not divide h,
1<=#{q in Q:7 divides q}<=4.                            (29u)
```

QED.

## 6. An explicit unbounded lacunary terminal cone

Let `Q` be seven distinct positive speeds, let `h` be odd, and put

```text
B=max(Q).                                               (30)
```

If

```text
sum_(q in Q minus {B})q+6h < (17/1078)B,               (31)
```

then

```text
G_Q is not contained in E_h.                           (32)
```

More quantitatively,

```text
measure(G_Q minus E_h)
 >=17/3234-[sum_(q!=B)q+6h]/(3B)>0.                    (33)
```

### Proof

Use the six-edge star centered at `B`. Applying (21) to every edge gives

```text
tau_star-Delta
 >=I_B+(6/7)sum_(q!=B)I_q-8/49
   -[sum_(q!=B)q+6h]/(3B).                             (34)
```

To minimize its first line under (25)--(26), place the smallest exceptional
value `1/42` on the coefficient-one center, the next value `2/77` on a
coefficient-`6/7` leaf, and use `1/35` on the other five leaves. This only
relaxes the actual distinct-speed constraints, so

```text
I_B+(6/7)sum_(q!=B)I_q-8/49
 >=1/42+(6/7)(2/77+5/35)-8/49
 =17/3234.                                             (35)
```

Equations (34)--(35) prove (33), and THM-2081 proves (32). QED.

No divisor-completeness, hereditary-primitivity, or terminal relative-height
assumption is needed. In particular, THM-2082's frozen family

```text
{1,...,s-1,360360*32^j},       7<=s<=10, h=13          (36)
```

has its rank-seven member `s=7` detected immediately: its far speed overwhelms
the fixed side mass in (31). The theorem is not being applied to the higher
ranks `s=8,9,10`; their common `3/31` exit remains the relevant certificate.
Thus the frozen family proves the scalar/code filters cannot create height,
while restored relative incidence does create an all-height exit on the live
rank-seven lane. The two results fit rather than conflict.

## 7. Transfer to the other active problems

Equation (7) is the LRC analogue of the unique-versus-coincident channel split
in the Gaussian Moment work. The axes are determined by positive marginal
data; only `R_h` contains genuine three-way cancellation. This suggests that
GMC tied-face arguments should likewise separate marginal/axis channels from
the genuinely multi-charge radial sum before attempting asymptotics. It does
not prove arbitrary-coefficient DvdK1: the GMC channel weights are complex and
have no positive Hunter inequality.

For code and tournament problems, (7) identifies the lost sidecar precisely.
Hamming weight or a score sequence can preserve marginal support while
destroying labelled incidence; `R_h` is what first sees the relation lattice
among three labels. In the `[72,36,16]` gate that role is played by labelled
cocircuit-design incidence, not merely `A_16`.

THM-2083 says every remaining rank-seven obstruction has a bounded relation.
The present theorem makes the complement constructive in a second direction:
large fastest-speed separation closes by BV mixing, while the residual is
both nonlacunary and short-relation-rich. That intersection is substantially
smaller than either condition alone.

## 8. Assumption challenge and Tournament Analysis

Candidate vertices considered here include runners, gaps, fixed circle
sections, section boundaries, wall events, residues, Fourier modes, matroid
cocircuits, and proof obligations. The useful vertices remain THM-2081's
restricted events `D_q intersect C_h`. The pairwise observable is their exact
intersection weight `w_pq`; the switch `w_pq>w_rs` yields a priority relation,
with residue label as tie gauge. But its tournament fingerprints do not
determine any term in (13): orientations discard rational weights, fold signs,
and the ternary relation lattice supporting `R_h`.

The fastest-speed star in Section 6 is an honest Hunter tree, not a tournament
claim. Its faithful sidecar is `(I_q; rho_pq; R_h(p,q))`, or equivalently the
labelled three-frequency relation hypergraph. The challenged assumption was
that either global pair overlaps or scalar mixed overlaps should suffice.
Their exact missing joint coordinate is (6).

## Exact referee

The companion replays the complete height-24 rank-seven bank, checks (13) on
all scalar survivors and the exact hostile decomposition (18)--(20), exhausts
the odd reduced spectrum (29), audits both modular identities
(29g)--(29t), and tests (21) by exact boundary atomization on a broad finite
triple bank. It uses explicit runtime checks; normal and `python -O`
transcripts byte-match and end in `PASS`. QED.
