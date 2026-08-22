---
id: THM-3476
title: "Rule 30 source P-adic reconstruction, depth-four transverse-jet barrier, and slack Pascal atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The full transverse
  Hasse-jet tower is a faithful P-adic
  reconstruction of the radial Rule 30 source, but no fixed jet depth is
  uniformly sufficient: one physical depth-four strip has target
  coefficients whose first live jet grows linearly with target time.  On
  every target coefficient, slack order is exactly twice the repeated-root
  order of its Green-selected distance packet.  No
  Rule 30 prize or unrestricted complexity lower bound is claimed.
source: root-rule30-next-targets-20260815
audit: >
  An independent hostile audit rederived the complete evaluation kernel,
  exact P-adic/transverse order equality, Pascal/Lucas tensor atlas, and the
  physical depth-four family; it also extended the Green-kernel checks beyond
  the companion universe.  A separate cross-carrier derivation and focused
  re-audit check the ballistic ramification-two atlas.  Ordinary and optimized
  companion runs match the stored output byte-for-byte: ACCEPT.
depends_on:
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
related:
  - THM-2043-period14-parity-hasse-jet-completeness
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-3466-factorial-face-stokes-and-keller-boundary-current
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum
script: 04-computation/rule30_transverse_jet_barrier_thm3476.py
output: 05-knowledge/results/rule30_transverse_jet_barrier_thm3476.out
script_sha256: f2a0c2f759910e88258fa8c8cfe69b3b987a7ff081fb3d1155fbd3af941f5fe8
output_sha256: 3fbd07442c10784a926b3bc9015db62aa2de726c8389449ff6ceb4cb4fb8778b
hash_basis: raw bytes
---

# THM-3476 -- Rule 30 source P-adic reconstruction, depth-four transverse-jet barrier, and slack Pascal atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The transverse slack marker introduced in THM-3471 is not merely a way to
repair one three-strip cancellation.  Its complete Hasse-jet tower is an
exact local coordinate on the whole radial source.  The coordinate is
faithful only as an unbounded tower: a genuine fixed source-depth strip has
an infinite family of target coefficients invisible to every jet below an
order proportional to the target time.

This is a representation theorem and a sharp representation-specific
barrier.  It is not nonperiodicity, balance, or a general lower bound for
computing the fixed-seed Rule 30 center.

## 1. Inheritance, board, and conventions

The closest proved mechanism is THM-3471's marked whole-strip compiler

```text
R_u(z,q)=z^(u+1) A_u(W_q)/(1+qz),                    (1)
q^2 W_q^2+(1+qz)W_q+z^2=0,                          (2)
```

where `u` is radial source depth, `d` is source distance, `v` is Green
transport slack, and

```text
alpha_u(d)=F_(u+d)(d),
A_u(X)=sum_(d>=0) alpha_u(d)X^d,
R_u(z,q)=sum_(d,v>=0) alpha_u(d)K(d+v,v)
                z^(u+1+2d+v)q^v,                   (3)
K(n,r)=[x^r](1+x+x^2)^n.                            (4)
```

The canonical hostile is the exact boundary circuit
`R_0(z,1)+R_1(z,1)+R_2(z,1)=0`: scalarizing `q` destroys a nonzero labelled
current.  The corrected near miss is to retain the marker and audit how many
of its jets are actually needed.  The least-used sidecars are the complete
integer slack and the source polynomial before Green evaluation.

The live concept board is:

1. the two-variable radial source series;
2. the quadratic Motzkin evaluation kernel;
3. transverse Hasse jets at `q=1`;
4. Lucas/Pascal residue coordinates;
5. a physical depth-four strip; and
6. coefficient extraction as a moving observer; and
7. the ballistic distance-to-slack ramification.

All algebra below is over `F_2`.  Every target coefficient contains only
finitely many source events, so all displayed coefficientwise slack sums are
well defined.

## 2. The full source and its exact evaluation kernel

Package every radial source depth before transport as

```text
S(z,X)=sum_(u,d>=0) alpha_u(d)z^u X^d
       in A:=F_2[[z,X]].                              (5)
```

Summing (1) over `u` gives the full marked radial response

```text
R(z,q)=z S(z,W_q)/(1+qz).                            (6)
```

At `q=1`, write `W=W_1` and put

```text
P(z,X)=X^2+(1+z)X+z^2.                               (7)
```

The small root `W=z^2+O(z^3)` satisfies `P(z,W)=0`.  In `A`,

```text
P(z,X)=(X+W)(X+W+1+z).                               (8)
```

The second factor has constant term one and is a unit.  Evaluation
`ev_W:A -> F_2[[z]]`, `X |-> W`, therefore has the exact kernel

```text
boxed: ker(ev_W)=(P).                                (9)
```

This identifies THM-3471's unmarked cancellation ideal, rather than merely
exhibiting one element of it.

## 2.1 Transverse order equals `P`-adic order

Put `q=1+epsilon`, write `W_epsilon=W_(1+epsilon)`, and use (2).  Direct
subtraction from (7) gives

```text
P(z,W_epsilon)
 =epsilon z W_epsilon+epsilon^2 W_epsilon^2
 =epsilon z W_epsilon(1+epsilon W_epsilon/z).         (10)
```

The last factor is a unit and `zW_epsilon` is nonzero with transverse order
zero.  Hence `P(z,W_epsilon)` has exact `epsilon`-order one.

Changing variables from `X` to `Y=X+W` is an automorphism of `A`.  Thus every
nonzero `T in A` has a unique finite `P`-adic order.  If
`T=P^r T_0` with `P` not dividing `T_0`, then (9)--(10) give

```text
ord_epsilon T(z,W_epsilon)=r=ord_P T.                (11)
```

The prefactor `z/(1+(1+epsilon)z)` in (6) is injective and has transverse
order zero, so it does not alter (11).

For `m>=0`, let `Jet_m` retain Hasse orders `0,...,m`, equivalently reduce
modulo `epsilon^(m+1)`.  Equation (11) proves the exact finite kernel

```text
boxed: ker(Jet_m o R)=(P^(m+1)).                     (12)
```

Since `A` is `P`-adically separated,

```text
intersection_(m>=0) (P^(m+1))={0}.                   (13)
```

Consequently the **entire** transverse jet tower reconstructs the complete
source series faithfully.  Formal packets `P^(2^r)` show sharply that every
prescribed finite initial segment of that tower can nevertheless vanish.
These formal packets are algebraic controls; they are not asserted to be
physical Rule 30 source packets.

## 2.2 Every fixed-strip finite jet stays in one quadratic field

Expand

```text
W_epsilon=sum_(j>=0) w_j epsilon^j,
w_0=W.                                                (14)
```

Coefficients of (2) give, for `j>=1`, with `w_(-1)=0`,

```text
(1+z)w_j
 =z w_(j-1)
  +[j even](w_(j/2)^2+w_(j/2-1)^2).                 (15)
```

Induction places every `w_j` in the same quadratic field
`F_2(z,W)`.  THM-3468 makes each fixed `A_u` rational, so every finite Hasse
jet of every fixed strip (and finite strip sum) also lies in that field.
This is a closure theorem for fixed `u`; it says nothing about the unbounded
sum over all source depths.

## 3. The old boundary circuit is exactly one `P`-factor

For the forced depths `u=0,1,2`, THM-3471 proves

```text
A_0=X^2/(1+X),
A_1=X^2/(1+X^2),
A_2=X/(1+X^2).                                       (16)
```

Their source polynomial is

```text
S_[0,2](z,X)=A_0+zA_1+z^2A_2
            =X P(z,X)/(1+X^2).                       (17)
```

The quotient is not divisible by `P`.  Equations (11)--(12) therefore say
that the three-strip circuit has transverse order exactly one.  Its scalar
specialization vanishes, its first jet is nonzero, and no higher-order
mystery is hidden in that particular cancellation.  This recovers
THM-3471's first-jet identity as the first case of the complete kernel
theorem.

## 4. Hasse jets are exactly a Pascal atlas of slack residues

Let a fixed target coefficient be

```text
C(q)=sum_(v>=0) a_v q^v in F_2[q],                   (18)
```

and fix a power of two `M=2^m`.  Define its first `M` Hasse values and its
slack-residue parity histogram by

```text
J_j=[epsilon^j]C(1+epsilon)
   =sum_v binom(v,j)a_v                 (0<=j<M),
E_r=sum_(v congruent r mod M) a_v       (0<=r<M).     (19)
```

Lucas' theorem gives `binom(v,j)=binom(v mod M,j) mod 2`, hence

```text
J_j=sum_(r=0)^(M-1) binom(r,j)E_r.                   (20)
```

The finite Pascal matrix is self-inverse over `F_2`.  Explicitly,

```text
boxed: E_r=sum_(j=0)^(M-1) binom(j,r)J_j.            (21)
```

Indeed the product coefficient is
`binom(s,r) sum_h binom(s-r,h)`, which is zero in `F_2` unless `s=r`.

Thus the first `M` jets recover **exactly** the odd/even occupancy of every
slack class modulo `M`.  They do not recover the quotient `floor(v/M)`, an
integer multiplicity before reduction mod two, or the ordering of events
inside a residue class.  This is the precise preservation/loss boundary of
the finite jet atlas.

### 4.1 The higher blocks recursively recover quotient carries

The apparent loss of the quotient has an exact multiscale organization.  Put

```text
v=r+Mh,
j=s+Mk,                0<=r,s<M.                     (21a)
```

Lucas factorization across the `m` low binary digits gives

```text
binom(r+Mh,s+Mk)=binom(r,s)binom(h,k) mod 2.          (21b)
```

Consequently the jet block `j=kM,...,kM+M-1`, after the same inverse Pascal
transform in its low index `s`, is

```text
E_(r,k)=sum_(h>=0) binom(h,k)a_(r+Mh).                (21c)
```

For each fixed residue `r`, the sequence `(E_(r,k))_(k>=0)` is exactly the
Hasse/Pascal transform of the quotient-slack sequence
`(a_(r+Mh))_(h>=0)`.  Iterating (21a)--(21c) yields a dyadic tree of slack
digits: the first block sees residues, the next blocks see their quotient
carries, and the full tower is the tensor-product Pascal transform on all
binary digits.  This is another proof of faithfulness, now coefficientwise
rather than through the global ideal (12).

For the physical pair in (23), the two exponents agree in the low `m` bits
and have quotient indices differing by one.  The residue block cancels;
the first quotient-carry coordinate, at jet `M`, detects the pair.

## 4.2 Ballistic transport is a ramification of exact index two

There is a smaller faithful carrier hidden behind the slack tower.  Fix a
source depth `u` and target `t`, put `R=t-u-1`, and retain only the source
distances selected by the Green kernel:

```text
beta_d=alpha_u(d)K(R-d,R),
B_(u,t)(X)=sum_(0<=d<=R/2) beta_d X^d.               (21d)
```

The corresponding target coefficient of the marked strip is

```text
C_(u,t)(q)=[z^t]R_u(z,q)
           =sum_d beta_d q^(R-2d)
           =q^R B_(u,t)(q^(-2)).                    (21e)
```

The last expression is polynomial because every selected `d` is at most
`R/2`.  Since

```text
q^(-2)+1=q^(-2)(q+1)^2,                              (21f)
```

reversal and the monomial factor are units at one.  For every nonzero packet,

```text
boxed:
 ord_(q+1) C_(u,t)=2 ord_(X+1) B_(u,t).              (21g)
```

This is an exact ramification-index-two bridge between the repeated-root
phase/Hasse filtration and the transverse slack filtration.

The complete coefficient atlas is equally explicit.  Let

```text
delta=R mod 2,
L=(R-delta)/2,
H_(u,t)(X)=X^L B_(u,t)(X^(-1)).                      (21h)
```

Then `C(q)=q^delta H(q^2)`.  If

```text
h_j=[eta^j]H(1+eta),
J_i=[epsilon^i]C(1+epsilon),                         (21i)
```

Frobenius gives, for every `j>=0`,

```text
boxed: J_(2j)=h_j,       J_(2j+1)=delta h_j.         (21j)
```

Thus half of the slack tower is a forced duplicate (when `R` is odd) or zero
(when `R` is even).  The distance-packet Hasse tower is the compressed
carrier; the slack tower is its ballistic twofold dilation.  In particular,
every first-live slack order is even.

For the depth-four family below, the selected packet is

```text
B_(4,t_M)(X)=X^6(1+X)^(M/2).                         (21k)
```

Equation (21g) predicts its live slack order `M`, and because `R` is odd,
(21j) sharpens (24) to

```text
J_M=J_(M+1)=1.                                       (21l)
```

This is coefficientwise ramification, not global `P^M` divisibility of the
whole physical source series.  It is also stripwise: summing source depths of
different parity mixes the two values of `delta` in (21h).  Source-depth
parity (or the full `u` label) is therefore a necessary sidecar before using
the compressed atlas.

## 5. A physical depth-four family defeats every fixed jet bound

The formal `P`-power controls above have a physical analogue after target
coefficient extraction.

### Theorem 5.1 (depth-four transverse-jet barrier)

For `M=2^m`, `m>=5`, set

```text
t_M=5M/4+10.                                         (22)
```

Then the genuine Rule 30 depth-four strip satisfies

```text
boxed:
 [z^(t_M)]R_4(z,q)
   =q^(M/4-7)(1+q^M)
   =q^(M/4-7)(1+q)^M.                                (23)
```

Consequently its Hasse jets at `q=1` obey

```text
D_q^[j] [z^(t_M)]R_4(z,q)|_(q=1)=0  for 0<=j<M,
D_q^[M] [z^(t_M)]R_4(z,q)|_(q=1)=1.                  (24)
```

The first detecting order grows linearly with target time:

```text
M=4(t_M-10)/5.                                       (25)
```

In the residue atlas (19)--(21), the two live slack exponents differ by
exactly `M`.  They occupy the same residue class and cancel in every jet below
order `M`.

### 5.1 The exact source word at depth four

Let `a_s(j)` be the centered Rule 30 row, and use the edge offsets

```text
b_r(s)=a_s(s-r),
ell_r(s)=a_s(-s+r),                                  (26)
```

with negative offsets zero.  The triangular edge recurrences are

```text
b_r(s+1)=b_r(s)+(b_(r-1)(s) OR b_(r-2)(s)),
ell_r(s+1)=ell_(r-2)(s)+ell_(r-1)(s)
             +(1+ell_(r-1)(s))ell_r(s).              (27)
```

At source time `s=d+4`, the prefixes close on the following exact cycle.
The right word lists `b_0,...,b_4`, the left word lists
`ell_0,...,ell_5`, and the last column is `alpha_4(d)`:

```text
d mod 8   right   left     alpha_4
0         10001   110010   0
1         11101   110111   1
2         10010   110010   0
3         11111   110111   0
4         10000   110010   0
5         11100   110111   1
6         10011   110010   1
7         11110   110111   1.                        (28)
```

Applying (27) to the last row returns the first, so this is an all-time
induction rather than a fitted prefix.  Hence

```text
alpha_4(d)=1 iff d mod 8 is in {1,5,6,7}.            (29)
```

### 5.2 The Motzkin kernel leaves exactly two events

Put `R=t_M-5=5M/4+5`.  At source distance `d`, the forced slack is
`v=R-2d`; Green palindromy makes its coefficient

```text
K(R-d,R).                                            (30)
```

The number `R` is odd.  If `d` is odd, then `R-d` is even, and the even
Frobenius power has no odd-degree coefficient.  By (29), only
`d=8j+6` remains.

Write `M=32*2^r` and `A=5*2^r-1`.  Three binary digit reductions give

```text
K(8N+7,8Q+5)=K(N,Q-1),                              (31)
```

so (30) becomes `K(A-j,A)`.  For
`0<=j<=floor(A/2)`,

```text
K(A-j,A)=1 iff j is in {0,2^(r+1)}.                  (32)
```

For `r=0`, this is the direct row `A=4`.  In the induction step, odd `j`
again gives an even power and odd degree, while `j=2h` reduces by the odd/odd
digit rule to the preceding value of `r`.  Thus (32) is universal.

The two surviving distances and slacks are

```text
d=6,       v=5M/4-7,
d=M/2+6,   v=M/4-7.                                  (33)
```

Their XOR is exactly (23), completing the proof.

## 6. What the linear jet barrier does and does not prove

The family (23) proves that no fixed transverse jet depth can recover all
target coefficients, even on the single physical strip `u=4`.  It is a
genuine moving-observer obstruction: the whole fixed strip has a rational
source word and all of its finite jets lie in one quadratic field, but the
jet order needed after extracting a target coefficient is unbounded.

It is **not** an information-theoretic lower bound.  The sparse polynomial in
(23) is described by its two exponents using only `O(log t_M)` bits.  A
representation retaining exact integer slack sees it immediately, while the
residue/Hasse representation needs order `M`.  Therefore the theorem rules
out bounded-jet carriers, not compressed sparse carriers, nonrectangular
macroblocks, uniform algorithms, or advice models.

Nor does infinite jet faithfulness make a finite Rule 30 prize argument.
Prize 1 would require control of the unbounded physical source completion;
Prize 2 requires the marked temporal address and correlation, not only slack
residues; Prize 3 requires a fixed computation model and a lower bound outside
this declared representation.

## 7. Connection ledger

### 7.1 Period-14 Hasse completeness

THM-2043 proves that a complete finite Hasse coordinate can still be globally
blind without magnitude, owner height, and first-exit data.  Here the same
warning appears internally: the first `M` coordinates are complete for slack
**modulo `M`**, while the quotient of slack by `M` is precisely what separates
the two physical events in (23).

### 7.2 Factorial boundary Stokes

THM-3466 and THM-3471 show that an unlabelled boundary current can vanish
while a transverse face/jet remains nonzero.  Equation (12) identifies the
entire vanishing ideal; (23) adds the hostile that the first live transverse
layer itself can move outward linearly with the observer.

### 7.3 Berggren three-branch circuits

The ternary Green kernel and the three forced strips form a labelled circuit,
not a ternary ancestry tree.  Summing slack passes from ternary words to an
endpoint parity and destroys path order.  The exact source, target, and
sidecar are:

```text
source:     radial source series S(z,X),
map:        X |-> W_q followed by target extraction,
preserved:  center response and slack-residue parity,
destroyed:  exact slack quotient and ternary path ancestry,
sidecar:    full q-polynomial, or a proved sufficient moving jet order. (34)
```

No Pythagorean, tournament, LRC, factorial-conjecture, or Jacobian-conjecture
conclusion transfers through this quotient.

## 8. Exact companion and proof boundary

The companion

```text
python3 04-computation/rule30_transverse_jet_barrier_thm3476.py
python3 -O 04-computation/rule30_transverse_jet_barrier_thm3476.py
```

checks with explicit optimization-stable gates:

1. direct centered Rule 30 rows against independent left/right edge
   recurrences through the full physical source universe;
2. the exact period-eight depth-four word;
3. polynomial ternary powers against the digit recursion for all
   `0<=n<=192` and every coefficient;
4. (31) on a `64 x 64` low-digit universe and (32) for `0<=r<=7`;
5. the physical family (23)--(24) and paired jet (21l) for `5<=m<=11`, by
   both Green engines;
6. the Pascal/Lucas tensor atlas for moduli `2,...,256`;
7. all `645` nonzero ballistic packets with `0<=u<=8` and `u+1<=t<=96`,
   including their valuations and full even/odd jet atlas; and
8. the symbolic characteristic and boundary-source identities.

The finite universes audit the implementation.  Equations (8)--(13),
(20)--(21), and the inductions in Sections 5.1--5.2 are the universal proofs.

## 9. Next exact targets

The theorem leaves three sharply typed targets:

1. classify the physical source's `P`-adic profile across unbounded `u`, not
   merely formal packets or one extracted coefficient;
2. combine the slack atlas with THM-3471's calibrated terminal phase arc, so
   transport quotient and physical owner are retained simultaneously; and
3. test sparse/nonrectangular cocycles that retain exact slack without paying
   for all Hasse layers.

These are possible routes toward the Rule 30 problems, not consequences of
the present theorem.  All three prizes remain open in this repository.
