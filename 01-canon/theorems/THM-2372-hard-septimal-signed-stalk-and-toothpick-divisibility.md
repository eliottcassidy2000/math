---
id: THM-2372
title: "Hard septimal signed stalk and toothpick divisibility"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the
  hard W=k=1 septimal lane of THM-2367, the six lower masks have a
  mean-zero integer defect F supported on the top unit and two high
  blocker combs, while its 7^M-root average vanishes. Exact moment
  identities give integral F^3=6P_3>0 and hence a fully nonzero additive
  colour triangle; at M=1 its 30 terms split into 12 same-chi_7 and 18
  mixed-chi_7 ordered triangles, but positivity need not select the
  Fano sector. More decisively, all possible low-blocker/lower-unit wall
  handoffs lie on one 1/13-density stalk. Exact absorber incidence then
  forces 14c_* to divide one of the two high blockers, and the hard
  valuation gap upgrades this to 98c_* divisibility. The associated
  minimal lawful high-target carrier renormalizes exactly to the
  isolated graft of THM-2367 and is necessarily circulant. This is a
  nested-toothpick constraint and a sharp one-factor no-go, not a
  scalar-row exclusion or target landing; the ledger remains 165 and
  LRC(14) remains open.
source: codex-2026-07-25-hard-septimal-signed-stalk
depends_on:
  - THM-1156-tooth-seam-chi7-bipartition
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
related:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2366-retained-probe-target-covariance-and-subthirteen-budget-bridge
  - THM-2369-complete-line-target-dirichlet-and-balanced-observable-no-go
script: 04-computation/lrc14_hard_septimal_signed_stalk_thm2372.py
output: 05-knowledge/results/lrc14_hard_septimal_signed_stalk_thm2372.out
script_sha256: cf5e26f938dc352be0fce9c5c2b620b0a78d21983e4854de92981559076f7b75
output_sha256: cc06096ff60c0e4f8adc93845dca59f181b785831bf404a40cf5ddd477d482cf
hash_basis: working-tree bytes (LF)
---

# THM-2372 -- hard septimal signed stalk and toothpick divisibility

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2367 turns the hardest surviving septimal role chamber into an exact
six-mask tiling away from three absorber combs. Two different structures
then appear:

```text
physical side:
  a mean-zero signed defect is supported on three combs;

event side:
  all low-blocker handoffs collapse onto one repeated endpoint stalk;

consequence:
  one high blocker is a fourteen-fold nested copy of the low blocker. (1)
```

The first side gives a finite additive-colour cubic probe. The second is
stronger arithmetically and supplies the first exact ladder
self-similarity law in this branch.

## 1. The hard `W=k=1` lane

Retain a canonical first-depth-one scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3)          (2)
```

almost everywhere, where

```text
D_v={x:||vx||<1/14},

C_H={x:||Hx||>1/7},

E_H={x:||Hx||<1/7}.                                 (3)
```

The inherited row typing gives:

- `H` is odd;
- `H,q_1,...,q_5` are units modulo thirteen; and
- every blocker is divisible by thirteen.

Put

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))>0,

nu_7(H)<M.                                          (4)
```

Assume the `k=W=1` alternative of THM-2367. Thus there are unique labels

```text
q_* among the q_i,              c_* among the c_j
```

such that

```text
nu_7(q_*)=M,

nu_7(q_i)<M                         for q_i!=q_*,

nu_7(c_*)<M,

nu_7(a),nu_7(b)>M                  for {a,b}={c_j}\{c_*}. (5)
```

No ordering by numerical size is intended. The word “high” below refers
only to the septimal valuation in (5).

## 2. The signed lower defect and its missing root colour

Define

```text
A=D_(q_*) union D_a union D_b,

L
 =1_(E_H)
  +sum_(q_i!=q_*)1_(D_(q_i))
  +1_(D_(c_*)),

F=L-1.                                               (6)
```

There is one guard mask of mass `2/7` and five ordinary danger masks of
mass `1/7`, so

```text
integral_T F=0.                                      (7)
```

THM-2367's seven-bin equality says

```text
F=0                         almost everywhere on A^c. (8)
```

This is the exact signed-defect form of the nonlinear scalar cover.

Let

```text
K=7^M,

(P_K F)(x)=1/K sum_(j=0)^(K-1)F(x+j/K).             (9)
```

Every absorber speed is divisible by `K`, so `A` is `1/K`-invariant.
Every lower speed has septimal valuation below `M`. Its `K`-root word
therefore contains exactly its expected number of occupied roots, away
from finitely many endpoints. Hence

```text
P_K F=0                         almost everywhere.  (10)
```

In the Fourier convention

```text
f^hat(n)=integral_T f(x)exp(-2 pi i n x)dx,
```

equations (7), (9), and the interval transform give

```text
F^hat(0)=0,

F^hat(Kn)=0                         for every n,

F^hat(n)
 =1_(H|n) sin(2 pi(n/H)/7)/(pi(n/H))
  +sum_(v in {q_i:q_i!=q_*} union {c_*})
     1_(v|n) sin(pi(n/v)/7)/(pi(n/v))              (11)
```

for `n!=0`. Thus `F` has six explicit lower frequency lattices, physical
support in three higher combs, and no `K`-multiple colour.

## 3. Exact moments and the additive seven-colour cubic

Let

```text
P_r
 =sum_(S subset {the six lower masks}, |S|=r)
    measure(intersection_(B in S)B).                (12)
```

Let `Z={L=0}`. Since `L` is integer-valued, `F=L-1`, and
`integral L=1`, the pointwise binomial identities give

```text
integral F_+=measure(Z),

||F||_1=2 measure(Z),

||F||_2^2=2P_2,

integral F^3=6P_3.                                  (13)
```

All six centred masks share an interval of length

```text
lambda
 =2 min(
      1/(7H),
      min_(v in {q_i:q_i!=q_*} union {c_*})1/(14v)
    )>0.                                            (14)
```

There are `binom(6,3)=20` triples. Therefore

```text
P_3>=20 lambda,

integral F^3>=120 lambda>0.                         (15)
```

On the same core `F=5`. Equation (7) then forces

```text
measure(Z)>=5 lambda.
```

Both the core and `Z` lie in the essential support of `F`, hence in `A`
up to null sets. In particular

```text
measure(A)>=6 lambda.                               (16)
```

The cubic has a finite root-stalk form. For `y in T`, put

```text
f_y(j)=F(y+j/K),                       j in Z/KZ,

f_y~(r)=1/K sum_j f_y(j)exp(-2 pi i rj/K).          (17)
```

Equation (10) is `f_y~(0)=0`. Finite Fourier inversion, followed by
integration in `y`, gives the absolutely finite identity

```text
integral F^3
 =sum_(
    r,s in Z/KZ;
    r!=0, s!=0, r+s!=0
  )
   integral f_y~(r)f_y~(s)f_y~(-r-s)dy.             (18)
```

There are `(K-1)(K-2)` ordered terms, so some fully nonzero additive
colour triangle has real part, and hence magnitude, at least

```text
120 lambda/((K-1)(K-2)).                            (19)
```

When `M=1`, the `30` ordered triples split under the Legendre character
`chi_7` as follows:

```text
6 permutations of (1,2,4),             all residues;

6 permutations of (3,5,6),             all nonresidues;

18 mixed-colour triples.                             (20)
```

Thus (18) is a genuine fully nontrivial additive `F_7` probe, but it
does **not** force the current into the same-colour sector. The two
six-term sectors in (20) are the precise still-open `chi_7` selector.
They are Fourier-mode triangles, not yet the runner-labelled Fano lines
of THM-1156.

This distinction is necessary: THM-2366 and THM-2369 show that an
uncoloured signed mixed current can remain inverse-character flat and
carry zero target energy.

## 4. Every lower-blocker handoff lies on one stalk

Write

```text
c_*=13u.                                             (21)
```

Consider a wall of `D_(c_*)` with sign `epsilon in {+1,-1}`:

```text
x=(14A+epsilon)/(14c_*).                            (22)
```

An opposite wall of a lower ordinary danger `D_q` has phase
`-epsilon/14`. Write

```text
g=gcd(u,q),                 u=ga,              q=gd.
```

The wall equality is

```text
d(14A+epsilon)=13a(14B-epsilon)                    (23)
```

for some integer `B`. Reduction modulo fourteen gives

```text
d congruent a                   (mod 14),
```

which is the exact seam criterion of THM-1156 in this gauge. Reduction
modulo thirteen, using `13` not dividing `d`, gives

```text
A congruent -epsilon            (mod 13).           (24)
```

Substituting (24) into (22) yields the universal stalk law

```text
u x congruent -epsilon/14       (mod 1).            (25)
```

For each sign, (25) contains only `u=c_*/13` of the `c_*` walls. The
union of **all** possible handoffs with all four lower `q` labels
therefore occupies at most

```text
1/13
```

of the signed `c_*` event layer. It is not four independent `1/13`
budgets: all four labels reuse the same stalk.

The guard cannot provide an additional handoff. A guard wall has
numerator `2` in fourteenth units, so the opposite-wall compatibility is

```text
14 gcd(c_*,H) divides H+2c_*.
```

The right side is odd because `H` is odd, while the left side is even.

This is the event-level refinement of the tooth-seam bipartition. The
exact-seam graph is still `chi_7`-bipartite by THM-1156, but every edge
incident to `c_*` is supported on the same repeated endpoint stalk.

## 5. Exact absorber incidence

Let `w` be any of the three absorbers `q_*,a,b`, and put

```text
g=gcd(c_*,w),              C=c_*/g,              m=w/g. (26)
```

By (5),

```text
7|m,                       7 does not divide C.     (27)
```

On either fixed sign branch of the `c_*` walls, exact counting on the
`C`-grid gives the fraction lying in the **closed** absorber comb
`{||wx||<=1/14}`:

```text
h_0(C)
 =[2 floor((C-1)/14)+1]/C             if m is even,

h_7(C)
 =2 floor((C+6)/14)/C                 if m is odd.  (28)
```

Only the parity of `m`, equivalently `m mod 14 in {0,7}`, remains after
the invertible permutation of the `C`-grid.

For the top unit `q_*`, thirteen does not divide `q_*` while thirteen
divides `c_*`. Hence `13|C`, and (28) gives

```text
h(C)<=2/13.                                         (29)
```

For a high blocker `w`, suppose

```text
14c_* does not divide w.                            (30)
```

If `C=1`, condition (30) forces `m` odd, and `h_7(1)=0`. If `C=2`,
coprimality forces `m` odd, and `h_7(2)=0`. For every `C>=3`, the two
elementary floor bounds in (28) give

```text
h(C)<=1/3.                                         (31)
```

## 6. The nested-toothpick theorem

The signed current in THM-2367 vanishes away from the absorber closure.
Consequently every `c_*` event atom not in an absorber must be cancelled
by an opposite lower event. This is a statement about the
distributional derivative of the almost-everywhere identity (8), not a
pointwise endpoint convention. Apply the union bounds (25), (29), and
(31) to one fixed sign branch. If neither high blocker were divisible
by `14c_*`, the fraction which could be serviced would be at most

```text
1/13+2/13+1/3+1/3
 =35/39
 <1.                                                (32)
```

This is impossible. Therefore:

> **Nested-toothpick conclusion.** In every hard `W=k=1` lane,
>
> ```text
> 14c_* divides a             or             14c_* divides b. (33)
> ```

The hard role inequalities sharpen (33). Since

```text
nu_7(c_*)<=M-1,

nu_7(a),nu_7(b)>=M+1,
```

the selected quotient has at least two factors of seven. Together with
its evenness in (33), this gives

```text
98c_* divides a             or             98c_* divides b. (34)
```

This conclusion is insensitive to the thirteen-adic depth of `c_*`.
Writing

```text
c_*=13^alpha u_*,

w=13^beta u_w,                  13 does not divide u_*u_w
```

for the selected high blocker, (34) forces `beta>=alpha` and

```text
u_w=98h u_*,

w/c_*=98h 13^(beta-alpha)                   for some h>=1. (35)
```

After quotienting by the common core `u_*`, the two teeth are a fixed
thirteen-adic ladder pair with a ninety-eight-fold hard-role scale jump
(the event-current argument itself first forces the factor fourteen).

Two role consequences are immediate:

- if `c_*=c_1`, then one of `c_2,c_3` has the form (35);
- if `c_*=c_2` on a strict row with `nu_13(c_2)>1`, the depth-one
  blocker `c_1` cannot be a multiple of `c_2`, so necessarily

  ```text
  98c_2 divides c_3.                                 (36)
  ```

  On a repeated-first row, the honest alternative remains
  `98c_2|c_1` or `98c_2|c_3`.

The local `N=49` chamber of THM-2367 fails this new law:

```text
98*195
```

divides neither `16562` nor `215306`. This independently explains why
that positive open chamber cannot globalize.

## 7. The nested minimal carrier renormalizes and loses drift

The forced toothpick is not itself the missing H-drift mechanism. In
fact it exactly preserves the bad functional form.

Put

```text
d_0(x)=1_(||x||<1/14),              g=1-d_0,

u_L(x)=1_(||x||>=L/14),             L in {1,2}.      (37)
```

Let `c,q,R` be positive integers with

```text
13|c,                  13 does not divide q,

nu_7(c)<nu_7(q),       7|R,                 d=Rc.    (38)
```

Consider the minimal nested lawful carrier

```text
K_L(r,t)
 =integral_T
   d_0(cx)
   d_0(dx-r/13)g(dx-t/13)
   u_L(qx+t/13) dx.                                  (39)
```

Write

```text
h=gcd(c,q),                  C=c/h,             Q=q/h,

g_0=gcd(Q,R),                n_0=Q/g_0,

ell_0=R/g_0,                 D=C ell_0.              (40)
```

Then

```text
13|D,             13 does not divide n_0,

gcd(D,n_0)=1.                                         (41)
```

The exact renormalization is

```text
K_L(r,t)
 =1/7 integral_T
    d_0(y-(r-t)/13)g(y)
    A_(D,L)(n_0 y+n_0 t/13)dy,                       (42)

A_(D,L)(Y)
 =1/D sum_(j=0)^(D-1)u_L((Y+j)/D).
```

To prove it, first average the `c` roots after Haar reduction:

```text
P(z)=d_0(z)A_(C,L)(Qz).
```

The `R`-root average

```text
B_RP(Y)=1/R sum_(j=0)^(R-1)P((Y+j)/R)
```

has the exact Fourier identity

```text
B_RP(Y)=1/7 A_(D,L)(n_0Y)                 almost everywhere. (43)
```

Indeed

```text
P^hat(N)
 =sum_m d_0^hat(N-Qm) u_L^hat(Cm).
```

At `N=Rn`, both `Rn` and `Qm` are seven-multiples. The factor
`d_0^hat(Rn-Qm)` vanishes unless `Rn=Qm`, in which case it equals
`1/7`. Writing `n=n_0k` gives exactly the Fourier coefficients on the
right of (43). This is an `L^2` convolution justified by
Cauchy--Schwarz; no conditionally rearranged double series is used.
Changing variables `y=Rz-t/13` gives (42).

Equation (42) is precisely `1/7` times THM-2367's isolated lawful graft
for the reduced pair `(D,n_0)`. Therefore the toothpick criterion applies
with no new detection argument:

```text
K_L(r,t)=Phi(r-t)

 iff 7|D

 iff nu_7(R)>nu_7(Q).                                (44)
```

The last equivalence uses `nu_7(c)<nu_7(q)`, so `7` does not divide `C`.

Now take `c=c_*`, `q=q_*`, and let `d` be the high blocker selected by
(34). Then

```text
nu_7(Q)=M-nu_7(c_*),

nu_7(R)=nu_7(d)-nu_7(c_*)>nu_7(Q).                  (45)
```

Thus the forced nested minimal high-target carrier is **always
circulant**. The wall gate creates an intrinsic low-to-high incidence,
but promoting that incidence directly to the high target kills the
one-factor H-drift by exact root self-similarity.

This locates the first possible escape. In `P=d_0 A_(C,L)(Q·)`, the
root average lives in septimal colour zero, and the nonconstant
zero-colour coefficients of `d_0` vanish. To evade this calculation, a
successful carrier must retain additional lower-factor information
which was lost in `P`, with nonzero mod-seven colours able to form a
zero-sum pair before the thirteen-target colour is read. This is the
concrete bridge between the toothpick renormalization and the cubic
probe (18).

## 8. Exact boundary and next use

The theorem does not exclude the nested branch. Indeed arbitrary
thirteen-adic profiles admit formal pairs of the shape (35). What has
changed is the object:

```text
before:
  an unstructured unique low blocker plus two high absorbers;

after:
  a common-core nested toothpick pair with quotient
  98h*13^(beta-alpha), whose minimal high-target carrier is
  exactly circulant.                                (46)
```

The highest-leverage next test is therefore the first **two-low-factor**
term in the lawful owner tensor, rather than another isolated graft.
One must decide whether its nonzero septimal colour pair survives the
cover-anchored co-shift and the thirteen-target projection. A downward
simplified co-shift would kill the zero-drift branch by the positive
deletion cone; equality would expose the exact residual character phase.

The additive cubic (18) supplies a second route. A theorem selecting one
of the same-`chi_7` sectors in (20), or converting its signed current
into THM-2369 charged target energy, would finally connect the Fano probe
to target landing. Neither selector is proved here.

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 9. Exact companion

The dependency-free exact companion:

- checks all polynomial moment identities for lower multiplicities
  `0,...,6`;
- enumerates all `30` nonzero additive `F_7` triangles and their
  `12+18` `chi_7` split;
- hostile-tests the universal handoff stalk and exact seam count over a
  finite bank;
- verifies that odd guards never provide an opposite handoff;
- checks the two exact absorber formulas and their `2/13` and `1/3`
  ceilings over a finite bank;
- verifies the `98`-fold nested normal form on all `165` valuation
  profiles;
- checks the exact root-average renormalization (42)--(44) for both
  `L=1,2`, including equal-valuation drift and strict-valuation
  circulancy controls; and
- checks that THM-2367's local hard chamber violates (34).

Run

```bash
python3 04-computation/lrc14_hard_septimal_signed_stalk_thm2372.py
python3 -O 04-computation/lrc14_hard_septimal_signed_stalk_thm2372.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_hard_septimal_signed_stalk_thm2372.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independently audited. QED.
