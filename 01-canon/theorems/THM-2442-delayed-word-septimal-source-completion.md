---
id: THM-2442
title: "Parity-BV restoration of delayed-word septimal source completion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  source-deletion branch of THM-2407 run first with terminal word one,
  the nonflat source phases occur in opposite pairs. For any fixed
  positive canonical terminal word, deleting its source-safe factor and
  lawfully translating that factor through C_7 gives seven word means of
  total mass six times the deleted word, so at most one phase has zero
  mass. Along either parity of R=13^k, BV mixing separates the fast
  source-coshifted word from the slow source-inserted packet. Hence the
  limiting source-phase target vector has phase zero equal to zero and
  another nonzero entry. For every sufficiently large clock the actual
  rational vector is nonconstant; irreducibility of Phi_7 over
  Q(zeta_13) forces all six lawful source/word skew-diagonal colours.
  Deep diagonal cancellation then gives a fixed-frequency coefficient
  with nonzero source residue mod 7, nonzero first-target colour mod 13,
  and a deep multiplier coprime to 91. This restores a positive literal
  terminal word in the deletion branch, but it does not make every
  coordinate a 91-unit, prove owner-conditioned graft drift, remove a
  scalar row, or prove LRC(14).
source: codex-2026-07-26-delayed-word-source-completion
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2407-owner-or-source-deletion-target-current-dichotomy
  - THM-2409-unfiltered-septimal-source-completion-and-word-phase-boundary
  - THM-2445-twenty-four-cell-graft-owner-conditioning
related:
  - THM-2418-alternating-base-thirteen-septimal-carry-matrix-and-rank-one-boundary
  - THM-2421-all-clock-septimal-ancestry-endpoint-event-detector
  - THM-2441-septimal-ancestry-event-period-collapse
  - THM-2448-right-endpoint-cospan-transition-atlas
script: 04-computation/lrc14_delayed_word_source_completion_thm2442.py
output: 05-knowledge/results/lrc14_delayed_word_source_completion_thm2442.out
script_sha256: 042efba1790288fee4affa61ccf70aa4969591e464e144f711d55b9a7fe75b4b
output_sha256: e69de343254fe35f8562867ed4e19af073989ad0daec7e5990686322102c92fb
hash_basis: working-tree bytes (LF)
---

# THM-2442 -- delayed words do not kill source completion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2409 completes the missing shallow source coordinate when the
terminal word is one, but stops at a real delayed word: an exact-address
source twist must move the transported source factor by the skew phase
`R ell`, and the seven danger translates no longer partition the
fixed-word current.

The missing operation is to take the clock large *after* retaining this
lawful skew.  On either parity subsequence, the seven shifted terminal
words are fixed finite-step functions.  BV mixing then separates them
from the slow endpoint packet.  Two elementary sidecars prevent the
limiting coefficient from vanishing:

```text
unfiltered nonflat source phases occur in a reflection pair;

a positive source-deleted terminal word can have zero mass in
at most one of its seven source-safe translates.                 (1)
```

The pair cannot fit into the one exceptional phase.  This repairs the
precise word-phase boundary of THM-2409 while preserving the actual
source/word skew diagonal.

## 1. The unfiltered deletion branch

Use the owner pivot and notation of THM-2407.  Thus `j` is the
stationary positive shallow source, `U_(s,t)` is the present packet with
only the source-danger factor deleted, and

```text
Delta_r(x)=d(c_3 x-r/13).
```

Run the THM-2407 dichotomy first with terminal word `1`.  This theorem
concerns its source-deletion branch.  Insert the source danger at its
seven septimal phases:

```text
d_(j,ell)(x)=d(c_j x-ell/7),             ell in F_7,

H^unf_ell(r,s,t)
 =integral_T U_(s,t)(x)d_(j,ell)(x)Delta_r(x) dx.   (2)
```

Put

```text
o_ell(s)=sum_r H^unf_ell(r,s,0),

z_ell(b)
 =(1/13)sum_s o_ell(s)zeta_13^(bs).                 (3)
```

Every `o_ell(s)` is rational.  THM-2409 gives, for every `b!=0`,

```text
z_0(b)=0,

sum_ell z_ell(b)=uhat(b)!=0.                        (4)
```

The first identity says that the actual owner phase is target-flat; the
second says that deleting the source leaves every nonzero first-target
colour.

## 2. Reflection supplies two nonflat phases

Every danger/safe base factor is centred and even.  The repaired target
dipoles are linear in `(s,t)`.  Therefore the change of variable
`x -> -x` gives the exact table symmetry

```text
H^unf_(-ell)(r,s,t)
 =H^unf_ell(-r,-s,-t).                              (5)
```

In particular,

```text
o_(-ell)(s)=o_ell(-s),

z_(-ell)(b)=z_ell(-b)=conjugate(z_ell(b)).           (6)
```

Some `o_lambda` with `lambda!=0` is nonflat.  Otherwise every
`z_ell(b)` in (4) would vanish.  By (5), `o_(-lambda)` is also
nonflat, and the two phases are distinct because seven is odd.

Prime-cyclic rational rigidity now removes a quantifier mismatch.  For
a rational thirteen-vector, being nonflat is equivalent to having
every nonzero `C_13` Fourier colour nonzero.  Hence

```text
z_lambda(b)!=0,
z_(-lambda)(b)!=0             for every b!=0.       (7)
```

The same opposite pair works for all twelve target colours.

## 3. A positive terminal word misses at most one phase

Fix any positive canonical word `Q_(j,sigma)` from THM-2305/2349.  It
contains the source-safe factor exactly once.  Delete only that factor
and write, almost everywhere,

```text
Q_(j,sigma)(y)
 =T_sigma(y) g(c_j y),

tau=measure(T_sigma)>0.                             (8)
```

For `epsilon in {+1,-1}` and `ell in F_7`, define the lawfully shifted
terminal word

```text
Q^(epsilon)_ell(y)
 =T_sigma(y)g(c_j y-epsilon ell/7),

q^(epsilon)_ell=measure(Q^(epsilon)_ell).           (9)
```

The seven danger translates partition one, so the seven safe
translates sum to six.  Consequently

```text
sum_ell q^(epsilon)_ell=6 tau,

0<=q^(epsilon)_ell<=tau.                           (10)
```

At most one of the seven numbers in (10) can vanish.  Indeed, two zero
entries would leave at most five entries of size `tau`, contradicting
their total `6 tau`.  Also

```text
q^(epsilon)_0=measure(Q_(j,sigma))>0.               (11)
```

There is a useful strengthening.  Every canonical scalar comb is
centred, so `T_sigma` is even.  The change of variable `y -> -y`
therefore gives

```text
q^(epsilon)_ell=q^(epsilon)_(-ell).                 (11a)
```

A nonzero zero-mass phase would occur with its distinct opposite,
contradicting the at-most-one-zero conclusion, while (11) excludes zero
at the fixed phase.  Thus in fact

```text
q^(epsilon)_ell>0             for every ell.        (11b)
```

In particular, apply this to the pair `{lambda,-lambda}` in (7).  For
each parity `epsilon`, one may choose

```text
lambda_epsilon in {lambda,-lambda}
```

such that

```text
q^(epsilon)_(lambda_epsilon)>0.                    (12)
```

Equations (7) and (12) retain one phase which is simultaneously
target-nonflat and visible to the lawfully shifted positive word.

## 4. BV mixing restores the full skew word

Let

```text
R=13^k,
epsilon=R mod 7=(-1)^k.
```

The full source-coordinate twist must translate the left present
source, bare source, and transported-word source occurrences.  At the
indicator boundary the two slow source occurrences collapse to the
single factor `d_(j,ell)`.  Thus the lawful delayed table is

```text
H^k_ell(r,s,t)
 =integral_T
    U_(s,t)(x)d_(j,ell)(x)Delta_r(x)
    Q^(epsilon)_ell(Rx) dx.                        (13)
```

This is not the fixed-word shortcut rejected by THM-2409.  Its word
phase is exactly

```text
R ell=epsilon ell             mod 7,                (14)
```

the necessary and sufficient skew-diagonal descent law.

For fixed parity, the functions on both sides of (13) range over finite
families of rational finite-step functions.  The two-BV estimate used
in THM-2365 gives, uniformly in all `ell,r,s,t`,

```text
H^k_ell(r,s,t)
 =q^(epsilon)_ell H^unf_ell(r,s,t)+O(1/R).          (15)
```

More explicitly, if

```text
A_(ell,r,s,t)
 =U_(s,t)d_(j,ell)Delta_r,

V_A=max Var(A_(ell,r,s,t)),
V_Q=max Var(Q^(epsilon)_ell),
```

then the absolute error in (15) is at most

```text
V_A V_Q/(12R).                                     (16)
```

Both parity classes contain arbitrarily large clocks.  Define the
actual target marginal

```text
z^k_ell(b)
 =(1/13)sum_(r,s)
   H^k_ell(r,s,0)zeta_13^(bs).                     (17)
```

Equations (4), (7), and (15) imply, for each fixed parity,

```text
z^k_0(b) -> 0,

z^k_(lambda_epsilon)(b)
 ->q^(epsilon)_(lambda_epsilon)
     z_(lambda_epsilon)(b)!=0                     (18)
```

for every `b!=0`.  There are only twelve target colours and two
parities.  Hence, for every sufficiently large `k`, the seven-vector

```text
(z^k_ell(b))_(ell in F_7)                          (19)
```

is nonconstant for every `b!=0`.

This is the step which the exact finite partition alone could not
supply: the word need not partition at finite clock, but its correctly
shifted mean cannot erase both members of the reflection pair.

### 4.1 The error has an exact periodic `1/R` form

The qualitative BV limit has a sharper arithmetic functional form.
Fix one phase and table entry, abbreviate its slow factor by `A`, and
its parity-fixed terminal word by `Q`.  Put all endpoints of `A` over a
common denominator

```text
D=13^K D_0,                    gcd(D_0,13)=1,
```

and write

```text
A=disjoint_union_i [A_i/D,B_i/D),
S=sum_i(B_i-A_i)=D measure(A).
```

For `R=13^K N`, inverse-branch disintegration gives

```text
R integral A(x)Q(Rx)dx
 =integral_Q N^A_(R,total)(y)dy,                   (18a)
```

where `N^A_(R,total)` counts all prefixes, without septimal colouring.
For arbitrary positive reduced clocks `N,N'`, if

```text
N'=N+D_0 u,
```

then the two integer prefix endpoints of component `i` shift by
`uA_i,uB_i`.  Therefore, pointwise in `y`,

```text
N^A_(R',total)(y)
 =N^A_(R,total)(y)+uS.                             (18b)
```

Since

```text
(R'-R)measure(A)=uS,
```

equations (18a)--(18b) give the exact covariance

```text
R'[
  integral A(x)Q(R'x)dx-measure(A)measure(Q)
 ]

=R[
  integral A(x)Q(Rx)dx-measure(A)measure(Q)
 ].                                               (18c)
```

Consequently, along each fixed word parity and residue class,

```text
integral A(x)Q(Rx)dx
 =measure(A)measure(Q)+C_class/R                  (18d)
```

with an exact rational `C_class`.  The class is periodic in the delayed
exponent with period dividing `ord_(D_0)(13)`; adjoining the word parity
still leaves only finitely many classes.  This is the uncoloured
counterpart of THM-2441's centred seven-ancestry period.  The denominator
of `Q` does not enter the clock modulus because `Q` lives in the fixed
terminal coordinate after (18a).

The identity is pointwise for the displayed half-open representative.
For the original strict-open comb product it holds almost everywhere,
which is all (18a) requires.  If `D_0=1`, interpret the residue period
as one.

Applied to every entry of (13), (18d) upgrades the `O(1/R)` in (15) to
an exact finite-class `1/R` tail.  Thus the threshold in (18) is not
merely existential: for a fixed scalar instance and word it can be
found by checking finitely many rational constants.  It is still not a
uniform numerical clock over all speed tuples.

## 5. Every lawful source colour survives

All entries of (13) are rational.  Therefore

```text
z^k_ell(b) in K=Q(zeta_13).
```

For `kappa in F_7`, put

```text
Z^k(kappa,b)
 =(1/7)sum_ell z^k_ell(b)zeta_7^(kappa ell).        (20)
```

The fields `Q(zeta_7)` and `Q(zeta_13)` have intersection `Q`, so
`Phi_7` is irreducible over `K`.  If (20) vanished for one
`kappa!=0`, the polynomial

```text
P_(k,b)(X)=sum_(ell=0)^6 z^k_ell(b)X^ell
```

would be a scalar multiple of `Phi_7`; all seven coefficients would be
equal.  This contradicts (19).  Thus, at every sufficiently large
clock,

```text
Z^k(kappa,b)!=0
       for every kappa!=0 and every b!=0.           (21)
```

Because (13) shifts the word by (14), `kappa` is the full
exact-address source residue, not THM-2409's forbidden off-diagonal
probe.

## 6. Deep landing and exact-frequency extraction

Retain all three thirteen-colours:

```text
B^k(kappa,alpha,b,h)
 =1/(7*13^3) sum_(ell,r,s,t)
    H^k_ell(r,s,t)
    zeta_7^(kappa ell)
    zeta_13^(alpha r+b s+h t),

J^k(kappa,alpha,b)=sum_h B^k(kappa,alpha,b,h).
                                                               (22)
```

At `t=0`, equations (20)--(22) give

```text
J^k(kappa,0,b)=Z^k(kappa,b)/13!=0.                (23)
```

Every table in (13) retains the moving deepest-safe factor.  Hence

```text
H^k_ell(t,s,t)=0
```

and source-colour-retained diagonal cancellation gives

```text
sum_alpha J^k(kappa,alpha,b)=0.                    (24)
```

For each `kappa,b!=0`, equations (23)--(24) force some
`alpha!=0` and some `h` with

```text
B^k(kappa,alpha,b,h)!=0.                           (25)
```

The absolutely convergent `m`-then-frequency expansion of THM-2365
applies to the finite source transform.  Some exact fixed-frequency
term in (25) has

```text
source relation residue kappa!=0             mod 7,
first target colour b!=0                     mod 13,
deep multiplier m=alpha                     mod 13. (26)
```

The centred deepest-danger coefficient vanishes at every nonzero
multiple of seven.  Thus every live term has

```text
gcd(m,91)=1.                                      (27)
```

The actual positive terminal word `Q_(j,sigma)`, its complete lawful
source translate, and the deep/target colours are all retained.

### 6.1 Same-coefficient form and the THM-2445 ghost

The word-restoration argument is coefficientwise.  Let

```text
z_ell in Q(zeta_13),              ell in F_7,
```

be any fixed deep/target coefficient of seven unfiltered lawful source
translates, and suppose

```text
z_0=0,                    sum_ell z_ell!=0.          (27a)
```

For the fixed positive canonical word in Section 3, (11b) gives
`q^(epsilon)_ell>0` for every phase.  Hence the limiting coefficient
vector

```text
(q^(epsilon)_ell z_ell)_ell
```

has phase zero equal to zero and at least one nonzero entry.  The same
mixing and `Phi_7` argument proves that every nonzero lawful source
colour survives at that **same fixed deep/target coefficient** for all
sufficiently large clocks.

This applies directly to THM-2445 Section 5.  In its circulant-companion
branch, equations (22)--(23) there give (27a) at one fixed eligible
coefficient

```text
(alpha,b,h),             alpha!=0,
```

while retaining the labelled `(c_3,q_*)` graft and common partial bare
endpoint.  Therefore THM-2445's same-coefficient all-six-source
completion admits any fixed positive canonical delayed word, with the
same `(alpha,b,h)`, at every sufficiently large row-dependent clock.
In THM-2445's alternative `D(O)>0`, the ordinary target-neutral BV
limit of THM-2365 Section 7 already restores the delayed word.

Thus the unique ownerless ghost in THM-2445 has no remaining
**delayed-word** branch.  This corollary still does not restore the
missing right-endpoint factors: its endpoint remains THM-2445's common
partial bare endpoint.

## 7. Scope and sharp stopping boundary

The theorem removes one specific debt:

```text
THM-2407 source-deletion branch
 + THM-2409 unfiltered source completion
 + any fixed positive canonical terminal word
 -> full delayed-word source-coordinate completion.            (28)
```

The reflection pair and the positive-word assumption are
load-bearing.  Without reflection, a unique nonflat source phase could
coincide with the unique zero word phase.  Without `tau>0`, (10) has no
content.  Without rational prime-cyclic rigidity, nonconstancy need not
make every charged source colour survive.  Holding the word fixed
instead of applying (14) remains gauge-invalid.

This is not yet the owner-conditioned no-cancellation theorem for the
surviving `nu_7(c_3)<=M` graft.  It does not handle the owner branch of
THM-2407, make every relation coordinate nonzero modulo both primes,
identify the extracted coefficient with a preselected triangle, remove
one of the `165` scalar rows, or prove LRC(14).

After composition with THM-2445, the residual graft debt is narrower:
transport one of its twenty-three repair/blocker-labelled cells through
the semantic terminal-word/repair alignment, and restore the fully
masked right endpoint without losing the typed coefficient.  The ghost
cell's delayed-word phase is discharged by Section 6.1. THM-2448
subsequently makes the missing-factor expansion finite and exact; it
does not supply semantic alignment or positive same-root service.

Its preserved object is stronger than a scalar drift statistic:

```text
fixed positive terminal word
 + exact source owner
 + lawful present/word skew phase
 + first-target colour
 + deep 91-unit multiplier.                         (29)
```

This is a direct transplantable sidecar for the remaining noncirculant
graft branch.

## 8. Exact companion

The dependency-free companion checks the finite algebra behind the
proof:

- every bounded seven-word mean vector of total `6 tau` has at most
  one zero;
- every nonempty reflection-closed subset of `F_7^*` meets every
  six-phase visible set;
- the unique-zero/unique-nonflat hostile becomes possible if reflection
  is deleted;
- the degree-six `Phi_7` kernel consists exactly of constant
  coefficient vectors; and
- positive rational controls retain a nonconstant word-weighted source
  vector for every possible missing phase.

Run

```text
python3 04-computation/lrc14_delayed_word_source_completion_thm2442.py
python3 -O 04-computation/lrc14_delayed_word_source_completion_thm2442.py
```

and compare both transcripts byte-for-byte with

```text
05-knowledge/results/lrc14_delayed_word_source_completion_thm2442.out.
```

An independent hostile audit reconstructed the source/word shift from
THM-2334's gauge and verified that it is exactly `omega=R ell`, checked
the reflection signs and the stronger even-word conclusion
(11a)--(11b), and audited the finite-clock use of cyclotomic rigidity.
In particular, phase zero is used only in the limit (18), never asserted
to vanish at finite clock.  A separate audit of (18a)--(18d) replayed
the half-open prefix count, endpoint qualification, denominator typing,
parity dependence, and exact `1/R` covariance.  Normal and optimized
executions are byte-identical to the stored transcript.  QED.
