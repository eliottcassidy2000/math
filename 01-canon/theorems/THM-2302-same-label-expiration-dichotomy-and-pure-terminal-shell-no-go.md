---
id: THM-2302
title: "Same-label all-clock dichotomy and pure terminal-shell no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED MASS/CLOCK/SELECTOR CORE.
  On every strict
  first-depth-one scalar row, the shallow label selected by THM-2293's
  Perron-weighted quadratic shell has exclusive-owner mass at least
  227189785662847/26519324819222476800. At every fixed time k, either that
  same label makes a genuine blocker-to-other-blocker handoff of mass at
  least 84157587746251753/341303710423393276416000, or its ancestry
  multiplicity and the same current blocker-only residual share an exact
  nonzero Fourier atom n with 1<=n<=4SS_j-1. This includes the selected
  owner's prescribed expiration and the common shell clock b+1. At the
  latter clock every marked source atom has thirteen-adic valuation at
  least b+1, whereas every THM-2293 shell-graph vertex has exact valuation
  b, so the scalar marked atom cannot give incidence. On the positive-return
  arm, retaining the unsquared root phase instead marks an actual grade-b
  vertex in the shell-selected character and pairs it at one bounded gauge
  index with a named current blocker. Every such marked vertex has a bounded
  incident c_3-multiple edge after the coprimality color is dropped; a
  seven-periodic hostile control proves that endpoint recurrence cannot
  restore the missing gcd(m,91)=1 color. An exact one-sheet rooted
  middle-owner carrier has every nonzero root energy equal to 1_(D_s^c),
  making the entire quadratic covariance functionally pure in the terminal
  c_3 channel. This hostile carrier satisfies only the rooted middle-owner/
  deep-exclusion interface, not the full scalar LRC cover. No profile is
  excluded and LRC(14) remains open.
source: codex-2026-07-25-same-label-connected-covariance
depends_on:
  - THM-2273-shallow-owner-flow-and-deep-successor-gap-spread
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
related:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2289-finite-spectral-certificates-bound-null-exclusive-owner-handoff-times
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2300-small-owner-multipliers-force-same-character-relation-multiples
  - THM-2301-essential-affine-arrangement-and-visible-rank-six-address-bank
  - THM-2303-terminal-component-phase-current-and-defect-rank
script: 04-computation/lrc14_same_label_shell_no_go_thm2302.py
output: 05-knowledge/results/lrc14_same_label_shell_no_go_thm2302.out
script_sha256: 1a3425a691303f873ca6bea31ec1ab248647b19fdd59e3eb96673d9a07a30f5e
output_sha256: b7734398a3b2436d4b057b5f829827f109a1c5116bb3a9ec298d673e5b252aaf
hash_basis: working-tree bytes (LF)
---

# THM-2302 -- same-label all-clock dichotomy and the pure terminal shell

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED MASS/CLOCK/SELECTOR
CORE.**

THM-2293 produces an edge between two Fourier atoms of one shallow
exclusive-owner label. THM-2296 produces a marked source atom paired with
current blocker service, but its owner was selected by a separate mass
pigeonhole. This leaves two apparent losses:

```text
the shell edge and the current-service atom may have different labels;

even after the labels agree, the marked atom need not be a vertex of the
shell graph.                                                        (1)
```

The first loss can be removed quantitatively. The second is not merely a
weakness of the constants: at the common shell clock the two objects occupy
disjoint thirteen-adic grades.

There is a further functional obstruction. A one-sheet rooted owner can have
quadratic energy exactly equal to the deepest safe-comb indicator. Squaring
then erases every trace of the owner selector, and a connected projection
which removes the pure deepest-coordinate channel removes the complete
certificate.

## 1. The shell-selected shallow label

Use the strict scalar notation

```text
T(x)=13x mod 1,

D_a={x in R/Z:||ax||<1/14},
C_H={x in R/Z:||Hx||>1/7},

c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

2<=b<c,                         5<=c<=19.             (2)
```

The guard and five unit speeds, as well as `u_1,u_2,u_3`, are
thirteen-units, and the live scalar cover is

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),

A_0 subset D_(c_1) union D_(c_2) union D_(c_3)       (3)
```

almost everywhere.

THM-2273's two shallow pieces are exactly the full exclusive-owner sets

```text
E_j
 =A_0 intersection D_(c_j)
       minus union_(h!=j)D_(c_h),          j in {1,2}. (4)
```

Put

```text
f_j=1_(E_j),          e_j=measure(E_j),

h_j=P^b f_j,

P f(y)=(1/13)sum_(r=0)^12 f((y+r)/13).               (5)
```

For the thirteen roots `x_r(y)=(y+r)/13`, define

```text
N_(j,kappa)(y)
 =sum_(r=0)^12 h_j(x_r(y))zeta^(-kappa r),

V_j(y)=sum_(kappa=1)^12 |N_(j,kappa)(y)|^2,

nu_j=integral V_j,                 zeta=exp(2*pi*i/13). (6)
```

THM-2293 proves

```text
0<=V_1<=C_1=42,          0<=V_2<=C_2=22,

nu_1+nu_2>=nu_0
 =227189785662847/58436012221844400.                 (7)
```

It selects a label `j in {1,2}` satisfying

```text
nu_j/C_j>=nu_0/(C_1+C_2)=nu_0/64.                   (8)
```

Fix this same label throughout the theorem. THM-2293 then supplies
`kappa in F_13^*`, a multiplier `m` with

```text
0<|m|<=578982,               gcd(m,91)=1,            (9)
```

and two actual Fourier atoms `A,A'` of `f_j` such that

```text
nu_13(A)=nu_13(A')=b,

A/13^b congruent A'/13^b congruent kappa mod 13,

A-A'=m c_3.                                          (10)
```

Thus the selected shell graph

```text
Gamma_(j,kappa)={
 A:
 f_j_hat(A)!=0,
 nu_13(A)=b,
 A/13^b congruent kappa mod 13
}                                                    (11)
```

has an edge.

## 2. A uniform mass floor for the same label

The weighted root energy itself recovers a source-mass floor. Put

```text
a_r=h_j(x_r(y)),                Sigma=sum_r a_r.     (12)
```

Finite root Parseval gives the exact identity

```text
V_j
 =13sum_r a_r^2-Sigma^2
 =sum_(r<t)(a_r-a_t)^2.                              (13)
```

Every `a_r` lies in `[0,1]`. Hence

```text
(a_r-a_t)^2<=a_r+a_t,

V_j<=12 Sigma.                                       (14)
```

The constant twelve is sharp at one root of weight one and twelve roots of
weight zero. Integrating the root sum and using that `P` preserves mass,

```text
integral Sigma
 =13 integral h_j
 =13e_j.                                             (15)
```

Equations (14)--(15) give the exact universal comparison

```text
nu_j<=156e_j.                                        (16)
```

Since `C_j>=22`, equations (8) and (16) imply

```text
e_j
 >=22nu_0/(64*156)
 =:e_*
 =227189785662847/26519324819222476800.              (17)
```

This is stronger by a factor `13/12` than the coarser estimate obtained
from `V_j<=13 Sigma`.

The target floor is label-independent. Put

```text
R_j=A_0 minus D_(c_j).                               (18)
```

Every positive-depth guard/blocker intersection has mass at most `10/91`.
Indeed THM-2273's parity-sharp cap has its global maximum at odd depth one.
Since

```text
measure(A_0)>=961/6930,
```

one has, for either shallow label,

```text
measure(R_j)
 >=961/6930-10/91
 =:beta
 =2593/90090.                                        (19)
```

The resulting same-label handoff floor is

```text
delta_*
 =e_* beta
 =84157587746251753
  /341303710423393276416000.                         (20)
```

## 3. One dichotomy at every fixed clock

Let

```text
S=H+sum_(i=1)^5q_i+sum_(h=1)^3c_h,

S_j=H+sum_(i=1)^5q_i+c_j.                            (21)
```

For every integer `k>=0`, define

```text
g_(j,k)=P^k f_j,

rho_j(k)
 =measure(E_j intersection T^(-k)R_j)
 =integral g_(j,k)1_(R_j).                           (22)
```

Then the following alternative holds at **every fixed clock**:

```text
rho_j(k)>=delta_*,

or there exists an integer n with

  1<=n<=4SS_j-1,

  f_j_hat(13^k n)!=0,
  1_(R_j)_hat(n)!=0.                                 (23)
```

### Proof

Suppose the first arm fails. Equations (17), (19), and (20) give

```text
rho_j(k)
 <delta_*
 <=measure(E_j)measure(R_j).                         (24)
```

Parseval therefore gives a nonzero covariance

```text
sum_(n!=0)
 g_(j,k)_hat(n) overline(1_(R_j)_hat(n))
 !=0.                                                (25)
```

The sum is absolutely convergent by Cauchy--Schwarz.

Up to null endpoints, `E_j` is a nine-factor circular rectangle pullback.
Its boundary has at most `2S` points. Every jump of `P^k f_j` is an image
of a source jump, so

```text
#Jump(g_(j,k))<=2S.                                  (26)
```

The target `R_j` is a seven-factor pullback and has at most `2S_j` jumps.
Apply THM-2296's bilinear endpoint-Prony lemma to the two step functions.
The endpoint-difference exponential sum has at most

```text
(2S)(2S_j)=4SS_j
```

nodes. Since (25) is nonzero, one common positive atom occurs by index
`4SS_j-1`. Finally,

```text
g_(j,k)_hat(n)=f_j_hat(13^k n),                      (27)
```

which proves (23). QED.

Whenever the first arm holds, it is a genuine blocker handoff. At a source
point of `E_j`, only `c_j` is available among the three blockers. At a
target point of `R_j`, the guard is active, all five unit masks and `c_j`
are unavailable, and the global cover forces one of the other two blockers.

The all-clock statement includes two useful specializations:

```text
k=lambda_j+1:
  the selected owner's prescribed expiration;

k=b+1:
  THM-2293's common quadratic shell clock.            (28)
```

For `j=2` these clocks agree. For `j=1`, the prescribed expiration is
`2`, whereas `b+1` is a later common shell clock. Conflating these clocks
would be an error.

## 4. Exact grade nonincidence at the common shell clock

Assume the spectral arm of (23) at

```text
k=b+1
```

and put

```text
B=13^(b+1)n.                                         (29)
```

Then

```text
f_j_hat(B)!=0,

nu_13(B)=b+1+nu_13(n)>=b+1.                          (30)
```

Every vertex of every THM-2293 shell graph on this label has exact
valuation `b`, by (10)--(11). Therefore

```text
B notin Gamma_(j,kappa')
            for every kappa' in F_13^*.              (31)
```

In particular the common atom from the same-label Prony theorem cannot be
an endpoint of the shell edge (10). The obstruction is stronger than
failure to select the right multiplier: after normalization to depth `b`,
the marked atom is in the zero root-character channel, whereas every shell
vertex lies in a nonzero character.

At the actual prescribed expiration, the valuation ledger is:

```text
j=2:
  k=b+1, so the same absolute nonincidence (31) holds;

j=1:
  k=2, and a marked atom can have shell grade b only if
  nu_13(n)=b-2.                                       (32)
```

Even in the second line, incidence with (10) additionally requires the
correct normalized residue and the exact absolute atom. Neither follows
from endpoint-Prony.

The clock `k=b` exposes the minimal lawful repair. If its common atom `n`
is a thirteen-unit, then `13^b n` is at least a marked vertex in the
grade-`b` graph

```text
Gamma_(j,n mod 13).
```

It still need not lie in the particular character graph selected by
THM-2293. Thus the remaining positive target is a rooted current-service
pairing which forces both an unramified atom and the same character, not a
scalar push to `b+1`.

## 5. The return arm marks the right grade and character

The last sentence of Section 4 describes the scalar projection. On the
positive-return arm at the common clock, the unsquared rooted word supplies
the missing grade and character automatically.

Assume

```text
rho_j(b+1)>=delta_*>0.                               (33)
```

Retain

```text
G_j=P^b f_j.                                         (34)
```

THM-2293 gives the exact middle-owner support split

```text
support(G_1) subset D_(u_2)^c,

support(G_2) subset D_(u_2).                         (35)
```

Thus every active thirteen-root vector of `G_1` has at most twelve
nonzero entries, and every active vector of `G_2` has at most two. The
entries are nonnegative rational multiples of `13^(-b)`.

For `a in {1,...,12}`, define the raw rooted word and its periodic gauge:

```text
M_a(y)=sum_(r=0)^12 G_j(x_r(y))zeta^(-ar),

N_a(y)=exp(-2*pi*i*a*y/13)M_a(y).                    (36)
```

Every nonzero proper rational root vector has `M_a(y)!=0`. Indeed, if

```text
sum_r v_r zeta^(-ar)=0
```

with rational `v_r`, the coefficient polynomial has degree at most twelve
and vanishes at a primitive thirteenth root. It must be a rational multiple
of

```text
Phi_13(X)=1+X+...+X^12.
```

All thirteen coefficients would then be equal. Proper support and
nonzeroness rule this out.

Partition `R_j` measurably by the first available one of the two other
blockers, exactly as in THM-2299. One named target `R_(j,t)` satisfies

```text
rho_(j,t)
 :=integral_(R_(j,t))P G_j
 >=rho_j(b+1)/2
 >0.                                                  (37)
```

Wherever the integrand in (37) is positive, the rooted vector is nonzero
and proper. Hence, for every nonzero character `a`,

```text
integral_(R_(j,t))|N_a(y)|^2dy>0.                   (38)
```

Put

```text
W_a=1_(R_(j,t))N_a.                                 (39)
```

Then

```text
<N_a,W_a>
 =integral_(R_(j,t))|N_a|^2
 >0.                                                  (40)
```

Parseval already proves that some integer gauge index `h` satisfies

```text
(N_a)_hat(h) overline((W_a)_hat(h))!=0.              (41)
```

There is also a uniform endpoint bound. The step amplitude `M_a` has at
most `2S` jumps: every jump is the image under `T` of one jump of `G_j`.
The named target uses no endpoint outside the same nine scalar banks, so
the step amplitude `1_(R_(j,t))M_a` has at most `4S` jumps.

For a piecewise exponential function

```text
U(y)=exp(-2*pi*i*a*y/13)A(y)
```

with step amplitude `A`, distributional differentiation gives

```text
2*pi*i(h+a/13) U_hat(h)
 =sum_(x in Jump(A))
    Delta_A(x)exp(-2*pi*i(h+a/13)x).                 (42)
```

Here the jump at the identified endpoint is counted in the covariant
section `A`; multiplication by the gauge makes `U` genuinely periodic, so
there is no unrecorded boundary term. Also

```text
h+a/13!=0                         for every h in Z,
```

including `h=0`.

Apply (42) to `N_a` and `W_a`, multiply one identity by the conjugate of
the other, and combine equal endpoint-difference nodes. The resulting
sequence in `h` is a nonzero exponential sum with at most

```text
(2S)(4S)=8S^2                                       (43)
```

nodes; it is nonzero by (41). A nonzero `L`-node exponential sum cannot
vanish at `L` consecutive integers. Therefore (41) holds for some

```text
0<=h<=8S^2-1.                                       (44)
```

Unlike THM-2296's ordinary step-function lemma, no zero at index zero is
used here. Index `h=0` is allowed and already gives the valid nonzero root
frequency `a`.

The gauge was chosen so that

```text
(N_a)_hat(h)
 =13 (G_j)_hat(a+13h)
 =13 f_j_hat(13^b(a+13h)).                          (45)
```

Consequently

```text
A_(a,h):=13^b(a+13h)
```

is an actual source atom with

```text
nu_13(A_(a,h))=b,

A_(a,h)/13^b congruent a mod 13.                    (46)
```

Take `a=kappa`, where `kappa` is the character of THM-2293's shell edge
in (10). Equations (39), (44)--(46) produce:

```text
an actual marked vertex A_(kappa,h) in Gamma_(j,kappa),

and a nonzero signed named-current-service coefficient
  (W_kappa)_hat(h)

at the same bounded gauge index h.                   (47)
```

Thus, on the return arm, the same-label and same-character vertex problem
is solved. What remains is incidence of this marked vertex with a
**unit-colored** terminal edge.

## 6. Marked degree is automatic only after erasing the unit color

There is a general endpoint-recurrence lemma.

> **Marked-degree lemma.** Let `f` be a real step function with at most `J`
> jumps, let `f_hat(A)!=0`, and let `c` be a positive integer. Choose
> `epsilon in {+1,-1}` so that `A` and `epsilon c` have the same sign.
> Then some
>
> ```text
> 1<=m<=J
> ```
>
> satisfies
>
> ```text
> f_hat(A+epsilon m c)!=0.                           (48)
> ```

Indeed, `A+epsilon mc` never vanishes. Distributional differentiation
turns

```text
C_m=(A+epsilon mc)f_hat(A+epsilon mc)
```

into an exponential sum over at most `J` endpoint nodes. Since

```text
C_0=A f_hat(A)!=0,
```

the sequence is not identically zero. If all `C_m` vanished for
`m=1,...,J`, its first `L<=J` nodes would give `L` consecutive zeros and
the Vandermonde lemma would force the sequence to vanish identically, a
contradiction. A constant step function cannot have the assumed nonzero
atom at `A!=0`, so `J>=1`; collisions among endpoint nodes only replace
`J` by a smaller positive `L`. Choosing the sign away from zero also
prevents the derivative factor from introducing an exceptional index.

Apply the lemma to `f=f_j`, the marked atom in (47), and `c=c_3`.
The source has at most `2S` jumps. Hence some

```text
1<=m<=2S
```

gives a second source atom

```text
A_(kappa,h)+epsilon m c_3.                           (49)
```

Since `c>b`, adding a multiple of `c_3` preserves exact valuation `b` and
the normalized residue `kappa`. Thus every marked vertex has bounded
positive degree in the **enlarged** shell graph in which all nonzero
terminal multipliers are allowed.

This does not recover THM-2293's load-bearing color

```text
gcd(m,91)=1.                                         (50)
```

The failure is sharp for endpoint recurrence. Let `F` be any nonconstant
`1/7`-periodic finite union of intervals with `F_hat(7)!=0`, for example

```text
F=union_(r=0)^6 [r/7,(r+1/4)/7).                    (51)
```

Then

```text
F_hat(n)=0                  unless 7 divides n.
```

If `c` is prime to seven, the marked atom `A=7` obeys

```text
F_hat(7+mc)!=0       only if 7 divides m.             (52)
```

Thus every incident edge visible to this recurrence can have a forbidden
seven-divisible multiplier. This is a universal step-function hostile
control, not an LRC row.

The return arm has therefore reached the exact boundary:

```text
same label
  + same nonzero root character
  + grade-b marked source vertex
  + named current blocker service at the same gauge index
  + bounded incident uncolored c_3 edge,

but no forced gcd(m,91)=1 incident edge.              (53)
```

The missing object is unit-colored degree or an equivalent phase-sensitive
third-order coefficient, not graph connectivity without colors.

THM-2301's proved rank-six pivot atlas gives one sharply typed candidate
sidecar. Of its `18` possible pivots with exactly one blocker column,
exactly `6` contain the selected source blocker `j`; only those six leave
both other blocker columns in the three-column complement. A theorem
forcing one of these owner-`j` type-one pivots, and pairing its complement
phase with (47), would retain the whole target-blocker word. The current
atlas guarantees only some unit pivot among all `83` types, so this is a
precise next target, not a consequence of THM-2301.

## 7. A one-sheet rooted-owner carrier with purely terminal energy

The grade obstruction concerns the actual LRC source and target. We now
give a separate exact hostile carrier for the proposed connected/cumulant
repair. Its scope is deliberately smaller.

Fix

```text
s=13^d u,          d>=0,          13 does not divide u,

Y=D_s^c.                                             (H1)
```

Let

```text
B_0={x in R/Z:||x||<1/26},

G=B_0 intersection T^(-1)Y.                         (H2)
```

The interval `B_0` has length `1/13`, so

```text
T:B_0 -> R/Z
```

is one-to-one away from its identified endpoints and preserves every target
point exactly once. Moreover

```text
B_0 subset D_1,

T:G -> Y is one-to-one and onto almost everywhere,

measure(G)=measure(Y)/13=6/91.                       (H3)
```

For each `y`, label the thirteen roots by `x_r(y)=(y+r)/13` and define

```text
m_y(r)=1_G(x_r(y)),

M_kappa(y)=sum_r m_y(r)zeta^(-kappa r).              (H4)
```

For `y in Y`, exactly one root lies in `G`; for `y notin Y`, none does.
If `sigma(y)` is the selected root label, then for every nonzero character

```text
M_kappa(y)
 =zeta^(-kappa sigma(y))1_Y(y),

|M_kappa(y)|^2=1_Y(y).                               (H5)
```

Thus the phase remembers the owner selector, but squaring erases it
completely. Summing the twelve nonzero characters gives

```text
V(y)=12 1_(D_s^c)(y).                                (H6)
```

The exact nonzero Fourier coefficients are

```text
Fourier(|M_kappa|^2,n)=0,              s does not divide n,

Fourier(|M_kappa|^2,ms)
 =-sin(pi*m/7)/(pi*m),                m!=0.          (H7)
```

Consequently, with `Q=13^(d+1)`,

```text
sum_(Q does not divide n)
 Fourier(|M_kappa|^2,n)
 overline(Fourier(1_(D_s),n))

 =-sum_(m!=0,13 does not divide m)
    (sin(pi*m/7)/(pi*m))^2
 <0.                                                   (H8)
```

This is precisely the signed whole-class terminal covariance used by
THM-2293. Every term in it is already a multiple of `s`.

The construction can be lifted to the same depth pattern. For any `b>=2`,
put

```text
E=T^(-b)G,

c_2=13^b,               c_3=13^(b+1)s.              (H9)
```

Then

```text
E subset D_(c_2),
E intersection D_(c_3)=empty,

P^b1_E=1_G.                                          (H10)
```

Thus (H5)--(H8) are the actual Perron-weighted root energies of a positive
indicator source with middle-owner containment and deepest exclusion.
Pulling `ms` back through `T^(b+1)` gives exactly

```text
13^(b+1)ms=m c_3.                                   (H11)
```

Let `Pi_s` be conditional expectation onto the closed subspace of functions
of `sx`. Equation (H6) gives

```text
(I-Pi_s)V=0.                                         (H12)
```

Therefore any connected or cumulant-style projection whose defining
property is to remove the pure deepest-coordinate channel can remove the
**entire** quadratic certificate under the presently retained rooted
axioms. Positivity of a connected remainder cannot be proved from
proper-root activation plus deepest-successor exclusion alone.

This is a functional no-go, not a word-by-word claim about THM-2293's
Jackson polynomial. It does not prove that the actual bounded relation
remainder in that theorem vanishes.

### Exact scope of the hostile carrier

The carrier (H1)--(H10) satisfies:

```text
positive rooted middle-owner support;
one-sheet properness on every active fibre;
all twelve nonzero root characters;
deepest-successor exclusion;
activation on every deepest safe gap;
the exact quadratic shell covariance and terminal pullback.            (H13)
```

It does **not** supply:

```text
the guard and five unit conditions defining A_0;
c_1 exclusivity;
the three-blocker scalar cover;
pairwise-distinct integer speed data for the whole row;
the actual current-service target R_j;
or a hypothetical LRC counterexample.                                  (H14)
```

In particular, (H8) cannot be cited as a hostile full-cover row. Its role
is sharper: it proves that the rooted-owner/deep-exclusion interface
exported to the quadratic shell does not itself contain a nonzero connected
channel.

## 8. The correct composition object

The same-label theorem removes one quotient loss:

```text
THM-2293 shell edge
  and
THM-2296-style current-service atom

can be taken on the same shallow exclusive-owner label.                 (47)
```

At the common shell clock, however, the result is a two-layer graded
fork:

```text
small-return / scalar spectral arm:
  the marked ancestry atom has grade at least b+1 and is disjoint from
  every shell graph;

positive-return / rooted arm:
  every nonzero character has an actual grade-b marked vertex paired
  with named current service at the same bounded gauge index.           (54)
```

On the rooted arm, the marked-degree lemma supplies an incident edge after
forgetting the multiplier color, while THM-2293 supplies a
`gcd(m,91)=1` edge somewhere in the same character graph. Neither result
forces the unit-colored edge to touch the marked vertex.

This is not naturally a tournament. Shell edges reverse by replacing `m`
with `-m`; current service is time-oriented and target-colored; and the
multiplier color is essential. The faithful carrier is a rooted
edge-colored graph together with a directed service incidence, or
equivalently a small graded hypergraph.

The first potentially sufficient object is third-order:

```text
a gcd(m,91)=1 shell edge
  whose one endpoint is paired with named current service,
```

equivalently a rooted bispectrum/marked-degree coefficient retaining the
unsquared root phase. The pure carrier (H5)--(H12) shows why merely centering
the squared energy need not create it.

THM-2303 supplies the faithful phase variable: the complex current on each
terminal handoff component and relative phase transport between parallel
components. The present theorem identifies where that current must land:
on a `gcd(m,91)=1` edge incident to the marked grade-`b` vertex. Neither
the unsigned energy nor an uncolored incidence edge is enough.

## 9. Frontier effect and loss ledger

The proved connection is

```text
source:
  THM-2293's Perron-selected shallow label and actual shell edge;

target:
  a quantitative all-clock blocker handoff or a bounded exact
  ancestry/current-service common atom on that same label;

map:
  convert quadratic root energy to an owner-mass floor, then apply the
  bilinear endpoint-Prony lemma at a chosen clock; on the return arm,
  retain the unsquared root gauge and apply endpoint-difference Prony;

preserved:
  the shallow owner label, actual exclusive source, current blocker-only
  target, exact clock, finite endpoint bank, Fourier frequency, and on the
  return arm the shell character plus one named target blocker;

destroyed at k=b+1:
  by the scalar spectral projection, the nonzero root character and exact
  grade-b vertex;

restored on the return arm:
  a grade-b marked vertex in the shell-selected character, a signed named
  current-service coefficient at the same bounded gauge index, and a
  bounded incident edge after forgetting the multiplier color;

sharp hostile mechanism:
  a one-sheet selector whose phase is erased by modulus-square, leaving a
  function only of the deepest successor coordinate, together with a
  seven-periodic endpoint sequence whose incident multipliers all have
  forbidden color;

needed sidecar:
  a unit-colored rooted current-service degree, or an equivalent nonzero
  third-order phase/bispectral coefficient.                             (55)
```

The theorem excludes none of the `150` strict profiles. The fifteen
repeated-first profiles lie outside THM-2293's two-shallow shell input.
The scalar ledger remains `165`, and LRC(14) remains open.

## 10. Exact companion and audit

The companion freezes:

1. the exact `150=120+15+15` strict-profile census;
2. `nu_0`, the sharp energy-to-mass coefficient `156`, `beta`, `e_*`, and
   `delta_*`;
3. the pointwise pair inequality behind (14);
4. the common-clock valuation separation; and
5. endpoint-free one-sheet selector cells at
   `s=1,2,13,26,169`; and
6. the `8S^2-1` rooted common-index ledger, the `2S` uncolored degree
   ledger, and the first visible multiplier `7` in the periodic hostile
   color test.

Reproduce with

```bash
python3 04-computation/lrc14_same_label_shell_no_go_thm2302.py
python3 -O 04-computation/lrc14_same_label_shell_no_go_thm2302.py
```

Normal and optimized outputs are byte-identical to the stored transcript.
An independent audit derived (17)--(23), caught the distinction between the
selected owner's prescribed expiration and the common clock `b+1`, improved
the coarse coefficient `169` to the sharp pointwise coefficient `156`, and
verified the exact rooted-carrier scope. QED.
