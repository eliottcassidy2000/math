---
id: THM-2418
title: "Alternating base-thirteen septimal carry matrix and rank-one boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For
  R=13^k, the physical-to-terminal seven-root map is the affine
  permutation l=floor(Ry)+(-1)^k r modulo seven, and the carry is the
  alternating sum of the first k base-thirteen digits. Haar averaging
  gives an explicit seven-state kernel K_k=P^k of rank seven; its six
  charged singular values are exactly 1/R. A terminal-only rational
  cylinder profile is nonflat exactly when all six nonzero septimal
  colours survive, with attenuation 1/R. This does not close the real
  word problem: a flat centred-comb terminal word is exact, a
  one-cylinder word makes the terminal matrix rank one, and the single
  fixed real/even BV-two source interval [3/13,10/13) of mass 7/13
  makes the complete source-weighted carry matrix rank one for every
  clock k>=1. More generally, one rational depth-K carry histogram
  classifies all deeper clocks as all-six-colour or flat, and an
  arbitrary rational finite-step source has an exact finite
  scale-periodic charged classifier. No canonical
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
script_sha256: 0cd8532b58aa28a4ab4ed7abbc454c6b1c07064205f3a29d2b49ba732f679b47
output_sha256: c72ea8192cecffc45b182d94e7e7e29981c3646b198edbb68dc7af76724a293f
hash_basis: working-tree bytes (LF)
---

# THM-2418 -- the septimal word phase is an alternating digit carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

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

For a seven-unit integer `A`, width `a in {1,2}`, and a circle point
represented by `0<=y<1`, put

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
digit carry of the canonical base transition. For an arbitrary real
lift one must first replace `y` by `{y}` (or retain the additional
integer-gauge term); equation (8) is not a digit expansion of
`floor(Ry)` outside the canonical interval.

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

## 5. One fixed real, even, large-mass source kills every clock

Terminal conditioning is not the main obstruction. The actual THM-2305
current also contains a source-owner weight before the terminal word.
Carry state alone cannot control that weight.

The obstruction is universal.  Let `B>p` be coprime positive integers
and put `R=B^k`.  For an integer `0<=A<A+p<=B`, the fixed first-level interval

```text
G_(B,p,A)=[A/B,(A+p)/B)                                    (34)
```

is, at every depth `k>=1`, exactly the union of the `pB^(k-1)`
consecutive prefix cells

```text
n=A B^(k-1),...,(A+p)B^(k-1)-1.                            (35)
```

Every residue modulo `p` occurs `B^(k-1)` times, so every entry of the
source-weighted `p`-state carry kernel is `1/B`:

```text
K_(B,p,A)=(1/B)J_p.                                        (36)
```

When `B,p` are odd, choosing `A=(B-p)/2` makes the interval inversion
invariant up to endpoints, BV-two, and of mass `p/B`.

Specialize now to `(B,p,A)=(13,7,3)`.  The single fixed circular
interval

```text
G=[3/13,10/13)                                              (37)
```

Up to endpoints, `G=-G`, so `1_G` is real and even.  It has exactly two
circular jumps and

```text
mu(G)=7/13.                                                 (38)
```

At depth `R=13^k`, the disintegration (10) sees `G` as the union of the
prefix cells

```text
n=3*13^(k-1),...,10*13^(k-1)-1.                            (39)
```

This is one block of `7*13^(k-1)` consecutive integers.  Every class
modulo seven therefore occurs exactly `13^(k-1)` times, independently
of `z`.  The source-weighted carry kernel is

```text
K_G=(1/13)J                                                (40)
```

for **every** `k>=1`.  It has rank one and annihilates the entire
six-dimensional charged subspace.  With an arbitrary terminal profile,

```text
T_(G,Q)=(1/13)J diag(q),                                   (41)
```

has rank at most one, and has rank one exactly when `q` is not the zero
profile.  Whenever a row has positive mass, every normalized source row
is identical.

The mechanism has a complete rational family.  If `A,s,t` are integers
with `s>=0` and `t>=1`, and

```text
G_(A,s,t)=[A/13^s,(A+7t)/13^s),
             0<=A<A+7t<=13^s,                               (42)
```

then for every `k>=s` its prefix block has length
`7t*13^(k-s)`.  Hence every carry class occurs
`t*13^(k-s)` times and

```text
K_(G_(A,s,t))=(t/13^s)J.                                   (43)
```

The special choice `(A,s,t)=(3,1,1)` is simultaneously fixed across
all clocks, even, BV-two, and of mass `7/13`.  Thus positivity,
evenness, bounded variation, large mass, and consistency of the source
set across all clocks still do not repair the charged carry colour.
Only a source--terminal correlation excluding equal carry counts can do
so.  The interval (34) is a formal rational source hostile, not asserted
to be a canonical exclusive-owner set.

There is also a hostile satisfying the basic ordinary-danger and
one-sheet source geometry.  For the ordinary speed `u=2`, put

```text
J=[81/169,88/169).                                         (43a)
```

This is the seven-cell depth-two block `81,...,87`, is
inversion-invariant almost everywhere, has variation two and mass
`7/169`, and lies strictly inside the central danger component

```text
(13/28,15/28) subset D_2,                                  (43b)
```

with clearance `71/4732` at both endpoints.  It lies in the single
thirteen-branch `[6/13,7/13)`, and `T_13` maps it bijectively to
`[3/13,10/13)`.  Thus it occupies one predecessor sheet.  At every
depth `k>=2`, each carry residue occurs `13^(k-2)` times, so

```text
K_J=(1/169)J_7.                                             (43c)
```

More generally, a sufficiently fine nonwrapping seven-cell block can
be chosen inside an open component of any ordinary danger comb.  If
that block and the canonical nonwrapping representative of its
reflection are disjoint, their union is even, BV-at-most-four, occupies
at most two predecessor sheets, and has kernel `(2/13^K)J_7`.
Disjointness and the canonical prefix lift are load-bearing: wrapped or
partially overlapping reflected blocks need not have uniform carry
counts.
Consequently danger support and one-/two-sheet sparsity alone do not
repair the charged colour.  The example is still not an exclusive
owner packet and does not encode the remaining cover factors.

## 6. One cylinder histogram classifies every deeper clock

Let a fixed rational source `G` be constant with rational values `g_m`
on the depth-`K` base-thirteen cells

```text
[m/13^K,(m+1)/13^K),             0<=m<13^K.                 (44)
```

Define its seven carry-bin sums

```text
b_r=sum_(m=r mod7) g_m,                    r in F_7.          (45)
```

For `k=K+d` and `e in F_7^*`, put

```text
C_(k,e)
 =integral_T G(y) zeta_7^(e floor(13^k y))dy.                (46)
```

Disintegrating each cell in (44) into `13^d` descendants gives

```text
C_(K+d,e)
 =13^(-(K+d)) theta_(d,e)
   sum_m g_m zeta_7^(e(-1)^d m),                             (47)

theta_(d,e)
 =1,                              d even,

 =-zeta_7^(-e),                   d odd.                     (48)
```

Indeed the descendant digit sum is

```text
sum_(n=0)^(13^d-1) zeta_7^(e n)
 =1                         for d even,

 =-zeta_7^(-e)              for d odd.                       (49)
```

Because the `b_r` are rational, irreducibility of `Phi_7` gives the
exact tail dichotomy:

```text
b_0=...=b_6
 iff C_(k,e)=0 for every e!=0 and every k>=K;

b nonconstant
 iff C_(k,e)!=0 for every e!=0 and every k>=K.               (50)
```

Thus one finite depth-`K` histogram completely classifies the infinite
clock tail.  The fixed hostile (37) has `K=1` and `b_r=1` for all `r`.
Conversely, any nonuniform rational histogram is a positive all-six
source-colour certificate at every later clock.  This is a source-only
statement; it still must be aligned with the terminal word/current.

## 7. Every fixed rational source--terminal pair has a finite scale classifier

The cylinder hypothesis has an exact weighted finite extension.  Let
`G` be a rational finite-step function whose endpoints have common
denominator

```text
D=13^K D_0,                 gcd(D_0,13)=1.                   (51)
```

First take terminal weight `Q=1`.  For `k=K+d`, refine to the grid of
denominator `13^k D_0`.  The sequence

```text
j -> zeta_7^(e floor(j/D_0))                                (52)
```

has period `7D_0` and zero sum over one period.  Therefore the
scaled charged coefficient

```text
13^k C_(k,e)                                                 (53)
```

depends only on

```text
13^d mod 7D_0.                                               (54)
```

It is consequently periodic in `d` with period dividing
`ord_(7D_0)(13)`.  Equivalently, after removing the uniform mean from
the seven carry masses, the scaled deviation vector is periodic.

At each fixed depth the seven masses are rational, so one vanished
charged colour again forces the mass vector to be uniform and hence
all six colours to vanish.  A finite exact scan of one period therefore
classifies every deeper clock of any fixed rational finite-step source.
Unlike (50), different period classes may alternate between flat and
all-six survival.

That last possibility really occurs.  At denominator `39`, take the
Boolean source supported on cells

```text
{0,5,8,9,13,16,18}.                                        (55)
```

Here `ord_21(13)=2`.  At extra depths `d=0,1,2,3`, the unscaled
carry-count vectors are respectively

```text
(1,1,1,1,1,1,1),
(15,19,18,17,12,4,6),
(169,169,169,169,169,169,169),
(2199,2203,2202,2201,2196,2188,2190).                       (56)
```

The last vector is the second plus the uniform vector `(2184,...,2184)`.
Thus the even class is flat while the odd class has all six charged
colours, exactly realizing the scale-period alternative.

Now let `Q` be any fixed integrable terminal profile and define

```text
C_(G,Q)(k,e)
 =integral_T G(y) Q({13^k y})
    zeta_7^(e floor(13^k y))dy.                              (57)
```

Write the step function in jump form

```text
G=g_0+sum_j gamma_j 1_[a_j,1),
a_j=A_j/(13^K D_0),                 0<a_j<1.                 (58)
```

For `R=13^k`, put `R a_j=m_j+tau_j`, with `m_j` integral and
`0<=tau_j<1`.  Disintegrating `y=(n+z)/R` gives, almost everywhere,

```text
R C_(G,Q)(k,e)
 =g_0 mu(Q) sum_(n=0)^(R-1) zeta_7^(en)
  +sum_j gamma_j [
      mu(Q) sum_(n=m_j+1)^(R-1) zeta_7^(en)
      +zeta_7^(e m_j) integral_[tau_j,1) Q(z)dz
    ].                                                       (59)
```

This includes `tau_j=0`: the last integral then supplies the `n=m_j`
term.  If `N=13^(k-K)`, the residue of `N` modulo `7D_0` determines
`R mod7`, `m_j mod7`, and `tau_j`.  Hence

```text
13^k C_(G,Q)(k,e)
```

depends only on `13^(k-K) mod 7D_0` and is periodic with period dividing
`ord_(7D_0)(13)`.  The denominator of `Q` does not enter this period:
`Q` is never rescaled after the `z={Ry}` disintegration.

When `G,Q` are rational finite-step functions, the seven weighted carry
masses

```text
a_(k,r)
 =integral_T G(y)Q({13^k y})
    1_(floor(13^k y)=r mod7)dy                               (60)
```

are rational.  Therefore one vanished nonzero colour is equivalent to
all seven `a_(k,r)` being equal, and hence to all six charged colours
vanishing.  The periodic histogram object is

```text
13^k (a_(k,r)-sum_s a_(k,s)/7),                             (61)
```

not the raw scaled masses.  In particular `G=Q=1` is the immediate
hostile to raw-mass periodicity.

This weighted theorem directly covers a fixed transported terminal word
`Q({13^k y})` only when the remaining source factor is one fixed scalar
rational step function `G(y)`, independent of the clock, role, and root
gauge.  A fixed finite root-indexed family can be treated entrywise, but
that produces only a periodic matrix and does not force a common phase
or a nonzero class.  Clock-selected Perron words, lawful co-shifts, and
prefix/tail/root-dependent owner packets need a separate alignment
lemma.  Even a surviving class here is a carry coefficient, not yet an
all-`91`-unit terminal current.

## 8. Canonical scope and the missing sidecar

For a THM-2305 terminal stratum `Q_(j,sigma)`, positivity proves only

```text
sum_l q_l>0.                                                   (62)
```

It does not prove that its seven cylinder masses are nonconstant.
Candidate THM-2414 supplies one strict local `W=8` stopping atlas, but
explicitly does not claim that atlas is a scalar cover or a
positive-mass global terminal-cylinder profile.

More importantly, the genuine THM-2305 current contains `G_j` or `E_j`.
That source weight depends on the full prefix `n`, tail `z`, and root
label `r`; it need not factor through (22). The rank-one construction
(34)--(43) proves that positivity, evenness, bounded variation, and
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

## 9. Exact companion

The dependency-free exact companion:

- verifies (5)--(8) on `18,040` rational word instances;
- exhausts every base-thirteen digit word through depth five;
- checks `P^k`, the direct prefix counts, (16)--(19), and the Gram
  identity through depth eight;
- exhausts all `127` nonempty Boolean terminal-cylinder supports;
- checks the flat and one-cylinder terminal hostiles; and
- verifies the fixed even rank-one source through depth eight and
  spot-checks five centred odd `(B,p)` instances of the universal
  analytic block law through four depths;
- verifies the fixed one-sheet `D_2` source through depths two to seven;
- exhausts all `8,192` Boolean depth-one source profiles and checks the
  fixed-cylinder tail formula; and
- checks the eventual finite scale-periodicity on rational endpoint
  denominators coprime to thirteen, including the period-two
  flat/all-six Boolean control (55)--(56); and
- checks the weighted fixed-source/fixed-terminal period, including a
  nonconstant denominator-five terminal profile and the distinction
  between raw scaled masses and centred scaled defects.

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

## 10. Independent hostile audit

The promoted base theorem was independently reconstructed from residue
counts.  With `A(r,l)=1_(r+l=6 mod7)` and `Pi=J/7`, the audit obtained

```text
P=(2J-A)/13,
K_(2m)=Pi+13^(-2m)(I-Pi),
K_(2m+1)=Pi-13^(-(2m+1))(A-Pi),
K_k K_k^T=Pi+13^(-2k)(I-Pi).
```

It also reproduced the flat terminal and one-cylinder hostiles, the
rational all-or-flat gate, and the corrected domain `y in [0,1)`.
Normal and optimized transcripts and the then-current hashes passed.

A second hostile audit independently verified the fixed source
`G=[3/13,10/13)` at every depth: its prefix block is inversion-invariant
and contains each carry residue `13^(k-1)` times, so
`K_G=(1/13)J`, with mass `7/13` and variation two.  It also checked the
nonwrapping aligned interval family (42), the one-sheet `D_2` hostile,
and the partial-overlap/wrapping boundaries.

A separate audit reconstructed the fixed-cylinder and weighted
source--terminal formulas in Sections 6--7, including endpoint typing,
the fact that `Q` does not enlarge the scale period, the rational
all-or-flat implication, and the centred-defect rather than raw-mass
periodicity.  Final normal and optimized executions byte-match the
stored twenty-eight-line transcript; the LF hashes in the frontmatter
were independently reproduced.
