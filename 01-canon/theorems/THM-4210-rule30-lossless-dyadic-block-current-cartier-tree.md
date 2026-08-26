---
id: THM-4210
title: "Rule 30 lossless dyadic block-current Cartier tree"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT. The isolated-seed
  collision current has an exact all-scale dyadic block representation. Its
  two physical children are the even Cartier channel at the two temporal
  parities; the odd Cartier channel is the transverse spatial-coset current
  needed for lossless reconstruction. The odd channel is exactly null for
  the fixed marked center after decimation, but it re-enters physical current
  generation through the Rule-30 quadratic coupling. The resulting
  infinite-state two-child automaton gives exact dyadic-defect and Prize-1/
  Prize-2 interfaces. Ambient endpoint, character-bank, and arbitrary-tail
  hostiles show why bounded truncations do not decide either prize without
  an all-scale physical-admissibility theorem. No Rule-30 prize is solved.
source: codex/session-rule30-debruijn-reset-20260826
depends_on: []
related:
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification
  - THM-3500-rule30-dyadic-section-cut-defect-and-cross-depth-valuation-carrier
  - THM-4048-rule30-periodicity-balance-and-model-firewalls
  - THM-4064-rule30-cyclotomic-kernel-character-and-c60-alias-obstruction
  - THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation
  - THM-4206-rule30-characteristic-address-contrast-deck-entropy-decomposition
script: 04-computation/rule30_dyadic_block_current_cartier_tree_thm4210.py
output: 05-knowledge/results/rule30_dyadic_block_current_cartier_tree_thm4210.out
script_sha256: e682fcbc2f33e78e590e484ac7a9a6229824ca606c1c9e0190735d6d93c19b97
output_sha256: 06bf5f020e1b8389fe443643aa6797f117500db64353899767febcdb863b446f
independent_audit_script: 04-computation/rule30_dyadic_block_current_cartier_tree_thm4210_independent_audit.py
independent_audit_output: 05-knowledge/results/rule30_dyadic_block_current_cartier_tree_thm4210_independent_audit.out
independent_audit_script_sha256: f1bde45597487c6f0b1e072bb3dcf8ac12a513086e2f432790d4890989b77bc1
independent_audit_output_sha256: e28a40ff220210f20a13aa8eb1b1234137d0e586ef01497b7757a449f9126119
hash_basis: raw LF bytes
---

# THM-4210 -- the lossless dyadic block-current Cartier tree

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.** The universal proofs
below are algebraic. The exact computations audit their physical Rule-30
instances; they do not promote a finite bank to an asymptotic statement.

The main correction is load-bearing. The odd Cartier channel is necessary to
reconstruct the full current, but it does not drive the fixed even-lattice
center directly. It is a transverse physical-admissibility sidecar, not a
third dyadic child and not by itself a Prize-1 or Prize-2 observable.

## 1. Physical block currents

Work over

```text
L=F_2[w,w^(-1)].
```

Write the isolated-seed Rule-30 row and its oriented adjacent-`11` current as

```text
A_t(z)=sum_j a_t(j)z^j,
N_t(z)=sum_j a_t(j)a_t(j+1)z^j,
q(z)=z^(-1)+1+z.                                      (1)
```

Since

```text
f(l,c,r)=l+c+r+cr       in F_2,
```

the exact row law is

```text
A_(t+1)=qA_t+N_t.                                      (2)
```

For a dyadic scale `h=2^e` and `0<=r<h`, define

```text
S_n^(h,r)
 =sum_(u=0)^(h-1) q^(h-1-u)N_(r+nh+u),                (3)

pi_h(sum_k f_k z^k)=sum_j f_(jh)w^j,

Y_n^(h,r)=pi_h A_(r+nh),
J_n^(h,r)=pi_h S_n^(h,r).                              (4)
```

The variable in `q` is understood from context. Frobenius gives

```text
q(z)^h=z^(-h)+1+z^h.                                  (5)
```

Iterating (2) for one block and applying `pi_h` proves

```text
boxed: Y_(n+1)=q(w)Y_n+J_n.                            (6)
```

The physical row at time `r` is supported in `[-r,r]`. Since `r<h`, its only
site divisible by `h` is the origin, so

```text
boxed: Y_0^(h,r)=c_r delta_0.                          (7)
```

Let

```text
K(d,k)=[x^k](1+x+x^2)^d.                              (8)
```

The central coefficient of `q(w)^d` is one for every `d`: for even `d`,
Frobenius reduces it to depth `d/2`; for odd `d`, the two noncentral terms of
the last `q` cannot meet the even support of `q^(d-1)`. Solving (6) therefore
gives the full fixed-node profile

```text
boxed:
c_(r+nh)
 =c_r+sum_(s<n)sum_j
      J_s^(h,r)(j)K(n-1-s,n-1-s-j).                  (9)
```

This is THM-3463's one-block collision-current formula at every node of the
2-kernel tree. The nonlinear difficulty is the physical sequence `J`, not
the ternary Green transport in (9).

## 2. Exact scale doubling

Split a `2h` block into its two consecutive `h` blocks. Equation (3) gives

```text
S_n^(2h,r)
 =q^h S_(2n)^(h,r)+S_(2n+1)^(h,r),                   (10)

S_n^(2h,r+h)
 =q^h S_(2n+1)^(h,r)+S_(2n+2)^(h,r).                 (11)
```

For `i in {0,1}`, define Laurent Cartier extraction by

```text
C_i(sum_k f_k w^k)=sum_j f_(2j+i)w^j.                (12)
```

At a fixed parent node, put

```text
F_n=qJ_n+J_(n+1),
E_n=C_0F_n,
O_n=C_1F_n.                                           (13)
```

Applying `pi_(2h)=C_0 pi_h` to (10)--(11) proves the two-child law

```text
boxed:
J_n^(2h,r)=E_(2n),
J_n^(2h,r+h)=E_(2n+1).                                (14)
```

The corresponding row decimations are

```text
Y_n^(2h,r)=C_0Y_(2n)^(h,r),
Y_n^(2h,r+h)=C_0Y_(2n+1)^(h,r).                       (15)
```

Thus the two temporal parities of the **even spatial channel** are exactly
the two physical dyadic children. Temporal parity and spatial parity are
different coordinates.

## 3. Lossless current lift

Every Laurent polynomial has the unique Cartier decomposition

```text
F_n(w)=E_n(w^2)+wO_n(w^2).                            (16)
```

Hence

```text
J_(n+1)=F_n+qJ_n.                                     (17)
```

Equations (14), (16), and (17) prove the bijection

```text
(J_n)_(n>=0)
 <-> (J_0,(J_n^left)_(n>=0),(J_n^right)_(n>=0),
          (O_n)_(n>=0)),                              (18)

J_n^left=E_(2n),             J_n^right=E_(2n+1).      (19)
```

The inverse interleaves the two child sequences to recover `E`, uses (16)
to recover each `F_n`, and then applies (17) successively. This is a linear
bijection on the ambient space `L^N_0`. On the physical current locus it is
an injective coordinate change onto its physical image; arbitrary target
tuples need not obey Rule 30's nonlinear current constraint.

Iterating (18) to finite depth is lossless when one retains the leaf current
sequences, `J_0` and `O` at every internal node. Since

```text
J_n^(1,0)=N_n,                                        (20)
```

this is an exact all-scale encoding of the original collision-current
sequence. It has three outputs at an internal node, but only `left` and
`right` are physical descendants. Calling `O` a third child is a type error.

### 3.1 Carrier covariance

For `tau_kF=w^kF`, even translations preserve the two channels:

```text
C_i tau_(2a)F=w^a C_iF.                               (21)
```

Odd translations exchange them, with the address shift

```text
C_0 tau_(2a+1)F=w^(a+1)C_1F,
C_1 tau_(2a+1)F=w^a C_0F.                             (22)
```

For reflection `iota F(w)=F(w^(-1))`,

```text
C_0 iota F=iota C_0F,
C_1 iota F=w^(-1)iota C_1F.                           (23)
```

These are symmetries of the Laurent carrier. The marked origin breaks odd
translation covariance, while reflection sends physical Rule 30 to its
mirror rule rather than preserving the named orbit. Shifting the block
sequence once sends `(J^left_n,J^right_n)` to
`(J^right_n,J^left_(n+1))`; a bare left/right swap omits the odometer carry.

## 4. The odd channel is transverse and center-null

Combining two copies of (6) gives

```text
Y_(n+2)=q^2Y_n+F_n.                                   (24)
```

Since

```text
q(w)^2=w^(-2)+1+w^2                                  (25)
```

has only even powers, put

```text
Z_n=C_0Y_n,                 W_n=C_1Y_n.               (26)
```

Cartier extraction in (24) yields two decoupled parity-coset laws:

```text
boxed:
Z_(n+2)=qZ_n+E_n,
W_(n+2)=qW_n+O_n.                                     (27)
```

For a physical node, the support bound also gives

```text
Z_0=c_r delta_0,             Z_1=c_(r+h)delta_0,      (28)

c_(r+nh)=Z_n(0).                                      (29)
```

Therefore `E`, not `O`, is the complete two-step forcing seen by the marked
even spatial coset. Once the two initial even rows are fixed, changing `O`
cannot change the center trace. The odd channel is nevertheless necessary
for the full row and current.

The minimal ambient kernel witness is

```text
J_0=0,       E_n=0 for all n,
O_0=1,       O_n=0 for n>0.                           (30)
```

Its lossless inverse is

```text
J_n=q^(n-1)w                    (n>=1),                (31)
```

so already `J_2=1+w+w^2` has a nonzero marked-current coefficient. But the
row difference generated from `Y_0=0` is

```text
Y_(2k)=q^(2k-2)w                (k>=1),
Y_(2k+1)=0.                                             (32)
```

Every nonzero row in (32) has odd spatial support, so its center is always
zero. This proves both sides of the correction: omitting `O` destroys the
current, while claiming that `O` directly changes the fixed center is false.

There is also a small physical control. The nodes `(h,r)=(2,0)` and `(8,3)`
both have

```text
c_r=1,              J_0=1,              E_0=0,       (33)
```

but their `O_0` values are respectively zero and one. Their center profiles
both begin `101010`. This is **FINITE-EXACT**; it does not claim equality of
their full even-channel sequences.

### 4.1 Why the transverse channel remains physically load-bearing

The linear decimation (27) does not make the physical nonlinear current
autonomous on `Z`. If

```text
A_t(w)=Z_t(w^2)+wW_t(w^2),                            (34)
```

and `odot` denotes coefficientwise multiplication, then the physical current
has parity channels

```text
C_0N_t=Z_t odot W_t,
(C_1N_t)(j)=W_t(j)Z_t(j+1).                           (35)
```

Thus future current generation couples the even and odd row cosets at every
site. An autonomous recursion for the physical `E` tree must retain the
transverse row/current data or prove a genuine elimination theorem. The
statement that `O` is center-null in (27) is not such an elimination.

## 5. The exact infinite-state 2-kernel automaton

For an abstract node state

```text
Xi=(b,(J_n)_(n>=0)),                                  (36)
```

form `E` by (13) and define

```text
T_0(Xi)=(b,(E_(2n))_(n>=0)),
T_1(Xi)=(b+J_0(0),(E_(2n+1))_(n>=0)).                 (37)
```

The physical root is

```text
Xi_(0,0)=(1,(N_n)_(n>=0)).                            (38)
```

Equations (7), (14), and (29) prove inductively that reading the low binary
digits of `r` reaches

```text
Xi_(e,r)=(c_r,(J_n^(2^e,r))_(n>=0)).                  (39)
```

This is an exact LSD-first two-child automaton for the entire center 2-kernel,
but its state contains an infinite Laurent-polynomial sequence. It is not a
finite-state theorem.

Define the dyadic collision-flip word

```text
delta_e(r)=J_0^(2^e,r)(0),       0<=r<2^e.             (40)
```

The first block of (6) and (7) gives

```text
boxed: delta_e(r)=c_(2^e+r)+c_r.                       (41)
```

The two child labels are read directly from the even current:

```text
boxed:
delta_(e+1)(r)=E_0^(2^e,r)(0),
delta_(e+1)(r+2^e)=E_1^(2^e,r)(0).                    (42)
```

The anchored dyadic-difference transform is lossless. For every `D`,

```text
(c_0,...,c_(2^D-1))
 <-> (c_0,{delta_e(r):0<=e<D, 0<=r<2^e}),             (43)
```

with inverse

```text
boxed:
c_n=c_0+
    xor_(e: bit_e(n)=1) delta_e(n mod 2^e).            (44)
```

Indeed, adding the nonzero binary digits from low to high makes each term in
(44) telescope one edge `c_r -> c_(r+2^e)` from (41).

## 6. Exact Prize-1 and Prize-2 interfaces

Put

```text
Gamma_delta(n)=xor_(e: bit_e(n)=1)delta_e(n mod 2^e). (45)
```

By (44), eventual `T`-periodicity from `N` is exactly

```text
Gamma_delta(n+T)=Gamma_delta(n)       for every n>=N. (46)
```

There is a simpler finite-witness criterion. Write

```text
T=2^a u,        u odd,
L=ord_u(2),     with L=1 when u=1.                    (47)
```

If the center is `T`-periodic from `N`, then for every `e>=a` and every
`0<=r<2^e` with `2^e+r>=N`,

```text
boxed: delta_(e+L)(r)=delta_e(r).                      (48)
```

Indeed, the two sides differ by

```text
c_(2^(e+L)+r)+c_(2^e+r),                              (49)
```

and the time gap is divisible by `T`. Hence it would suffice for Prize 1 to
prove

```text
for every N,T, some e>=a with 2^e>=N and some r<2^e
has delta_(e+L)(r)!=delta_e(r).                        (50)
```

Such a mismatch is a finite certificate. Along the arithmetic progression
between the two unequal endpoints, at least one adjacent `T`-step relation
fails; THM-4048 then supplies a nonzero inverse-boundary bit no later than
depth `2^(e+L)+r`. Criterion (50) is sufficient, not known necessary for an
arbitrary nonperiodic sequence.

Hypothetical eventual periodicity also forces the within-shell relation

```text
delta_e(r+T)=delta_e(r)                               (51)
```

whenever `r>=N` and `r+T<2^e`. Thus a periodic center would make the dyadic
cocycle finite-patterned both within shells and across scales.

For Prize 2, put `s_n=2c_n-1`. Equation (41) gives the exact signed shell
identity

```text
boxed:
sum_(r<t)s_(2^e+r)
 =sum_(r<t)s_r(-1)^(delta_e(r)).                       (52)
```

Combining (52) with THM-4048's dyadic-shell equivalence shows that Prize 2 is
exactly the obligation

```text
boxed:
max_(0<=t<=2^e)
 |sum_(r<t)s_r(-1)^(delta_e(r))|=o(2^e).               (53)
```

Marginal balance of `delta_e` is insufficient. The load-bearing statistic is
its signed correlation with the inherited prefix `s_r`.

Neither (50) nor (53) is proved uniformly here. Both named prizes remain
**OPEN**.

## 7. What finite truncations cannot prove by themselves

The following stopping results concern the universal driven-Rule-150
envelope of (6), not the one physical Rule-30 current. Their role is to type
the missing physical-admissibility theorem.

### 7.1 Growing spatial cone

Equation (3) makes the marked dyadic defect

```text
delta_e(r)
 =sum_(u=0)^(h-1)[z^0]q^(h-1-u)N_(r+u).              (54)
```

At time offset `u`, it sees exactly the sites `|j|<=h-1-u`. The bound is
sharp in the ambient current category: an impulse `z^(h-1)` at `u=0`
reaches the center through the unique all-left endpoint path. Hence no fixed
spatial coefficient window determines (54) at every scale.

More generally, for `P_m=w^m`, the response trace

```text
R(P_m)_n=[w^0]q^nP_m                                  (55)
```

has its first nonzero bit exactly at `n=m`. The traces in (55) are pairwise
distinct. Therefore no finite set-valued static encoding of arbitrary
Laurent impulses can determine their complete future marked-center response.
This is stronger than a no-go for linear summaries, but it is still an
ambient statement.

### 7.2 Fixed character and Hasse banks

Fix finitely many quotient jets

```text
L/(1+w^(M_i))^(d_i)                                   (56)
```

and a coefficient window `|j|<=R`. Choose `2^k>=max_i d_i`, put

```text
H=lcm_i(M_i)2^k,
P=w^b(1+w^H),                 b>R.                    (57)
```

In characteristic two, `1+w^H` is divisible by every ideal in (56). Thus
`P` vanishes in all retained jets and in the coefficient window. These
quotients remain blind after multiplication by `q`. Nevertheless,

```text
[w^0]q^bP=1,                                           (58)
```

because the nearer monomial arrives by its unique endpoint path before the
farther one can arrive. No fixed finite cyclic-character/Hasse bank plus a
fixed coefficient window is therefore target-faithful on the ambient
carrier.

The primitive-cubic failure occurs on the physical orbit already:

```text
J_1^(2,1)=w^(-1),
J_3^(2,1)=w^(-3)+1+w^2.                               (59)
```

The two polynomials are congruent modulo `w^3-1`, but their marked constant
coefficients are respectively zero and one. This is **FINITE-EXACT**. It
shows that the full cubic character packet is a current sensor, not a marked
origin observer.

### 7.3 Every finite formal prefix has opposite ambient continuations

Let `A_0,...,A_(H+1)` be any finite physical Rule-30 row prefix. Continue it
by any chosen finite-support rows `B_t`, and define the unobserved forcing by

```text
D_t=B_(t+1)+qB_t.                                     (60)
```

Then (6), every block identity, and every Cartier identity above hold for
the resulting ambient driven system. The currents through time `H` remain
the physical adjacent-`11` currents.

One continuation may take every later row to be `delta_0`, giving the
eventually constant-one, unbalanced center. Another may take the later rows
to be

```text
B_(H+2+k)=t_k delta_0,
t_k=popcount(k) mod 2,                                 (61)
```

giving a balanced, non-eventually-periodic Thue--Morse tail. Balance follows
from `t_(2k)=t_k`, `t_(2k+1)=1-t_k`, which cancels every adjacent sign pair.
For nonperiodicity, if a period `T` held, compare `2^m-T` with `2^m` for
arbitrarily large `m`; the first word is the `m`-bit complement of `T-1`, so
its digit-sum parity changes when the parity of `m` changes.

Thus every finite physical prefix of the formal carrier has ambient
continuations with opposite answers to both asymptotic questions. This does
**not** prove that no finite certificate or finite-state theorem can decide
the named physical orbit. It proves that the unobserved forcing must be tied
to the physical quadratic law (35); the Frobenius/Cartier identities alone
cannot supply the missing uniform argument.

The same arbitrary-trace realization shows that the ambient algebra imposes
neither automaticity nor rational-coefficient P-recursiveness: it realizes
every binary center sequence, whereas each of those finite-description
classes is countable. This is a formalism firewall, not a classification of
the physical center.

## 8. The finite-state conditional and exact state-count scope

Suppose the physical reachable node states (39) admit a finite quotient
`Sigma` such that

1. the marked bit `c_r` is constant on each quotient state; and
2. both transitions `T_0,T_1` in (37) descend to `Sigma`.

Reading the binary digits of `n` least significant first then computes `c_n`
with `O(log n)` transitions. Equivalently, the center has a finite 2-kernel
and is 2-automatic. This would be a genuine query-complexity consequence,
but no such quotient is proved.

Conversely, infinitely many pairwise output-distinguishable physical node
states would prove an infinite 2-kernel. Since every eventually periodic
sequence is 2-automatic, that would prove Prize 1. It remains **OPEN**.

THM-3463's `25,830` distinct length-16 node-profile prefixes are a finite
lower bound, not this infinite quantifier. They force at least `25,830`
states in an LSD-first residual automaton. For a conventional `q`-state
MSD-first binary DFAO with leading-zero normalization, a fixed low-digit
suffix changes only the output labeling on its `q` states, so it can induce
at most `2^q` 2-kernel profiles. The same bank therefore gives only

```text
q>=ceil(log_2 25830)=15.                               (62)
```

in that convention. Do not quote the LSD bound as an MSD state bound.

## 9. Canon comparison and non-duplication

The source/target ledger is

```text
source:     isolated-seed rows and actual adjacent-11 collision currents
target:     dyadic coarse rows; two physical child currents; transverse O
map:        block Duhamel sum, lattice sampling, and Laurent Cartier split
preserves:  every coarse row/current coefficient; the marked 2-kernel;
            the full parent current when J_0 and O are retained
destroys:   E-only descent drops the odd spatial row/current coset;
            scalar delta drops every unmarked current coefficient
sidecar:    O plus boundary data for lossless current reconstruction;
            physical quadratic admissibility for autonomous future descent
```

This theorem does not duplicate THM-4204 or THM-4206.

- THM-4204 maps a binary word to a saturated 41-state **inverse cyclic
  spatial** monoid. Its finite-ring Haar entropy bridge uses complete-row
  inverse fibre sizes. The present map follows the **forward deterministic
  named-seed** collision current. The monoid erases `E`, `O`, current
  location, and every post-reset continuation; there is no inverse map from
  its saturated state to this tree.
- THM-4206's source is a finite spacetime family under Haar input and its
  target is a characteristic-address contrast deck. Here the source is one
  deterministic orbit and the target is a Laurent current carrier. Neither
  THM-4204's row entropy nor THM-4206's pivot entropy transfers to (53).
- THM-3488's two Cartier channels encode target-fixed slack/source-depth
  parity and Hasse jets. Here Cartier acts on coarse **spatial** parity after
  two block steps. The syntax is shared; the carriers are not.
- THM-3500's Mersenne operator is lossless within one depth but does not
  transport its actual phase current across depths. Equations (10)--(18)
  supply such a law for the different block current `J`, not for
  THM-3500's current.
- THM-3463 is the closest proved mechanism: (9) extends its single dyadic
  jump to every fixed 2-kernel profile, while (18) identifies the exact
  cross-scale coordinate its scalar defect did not retain.

At a primitive cube root `omega in F_4`, `q(omega)=0`, so (6) gives

```text
Y_(n+1)(omega)=J_n(omega).                             (63)
```

THM-4204 has a primitive-cubic factor on its `N=4h` inverse-ring geometry.
Equation (59) proves that the shared annihilator is not a center observer;
without a carrier map on the same physical geometry, it is not an
intertwiner.

## 10. Exact verification

The primary companion uses centered Laurent supports and direct Rule-30 rows.
It independently checks the center with the packed recurrence, then verifies
the coarse law, both child laws, lossless inversion, parity-coset dynamics,
dyadic words, shell prefixes, ambient hostiles, and the physical cubic alias.
The no-import audit uses scalar coefficient dictionaries and a separate
ordinary-row evolution.

Run from the repository root:

```bash
python3 04-computation/rule30_dyadic_block_current_cartier_tree_thm4210.py
python3 -O 04-computation/rule30_dyadic_block_current_cartier_tree_thm4210.py
python3 04-computation/rule30_dyadic_block_current_cartier_tree_thm4210_independent_audit.py
python3 -O 04-computation/rule30_dyadic_block_current_cartier_tree_thm4210_independent_audit.py
```

The finite universes and comparison counts are printed in the stored
transcripts. Universal claims are the proofs above. No finite sign pattern,
state saturation, automaticity conclusion, or named prize is inferred.
