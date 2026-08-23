---
id: THM-3734
title: "Automorphic Cohn diagonal binomial divided-power towers"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT AUDIT.  On every diagonal
  constant right-SL2 slice of the automorphic Cohn matrix, all polynomial
  lower-row and upper-row closures are classified.  They occur in two paired
  integer-depth binomial towers (including their constant depth-one
  boundaries).  Every closed row is the gradient of an explicit unimodular,
  noncoordinate divided-power potential, but its complementary Hamiltonian
  debt has no polynomial solution.  Thus no final opposite left shear
  completes any member to a Jacobian matrix.  The proof is uniform in the
  depth; exact checks through depth 12 are positive controls, not the source
  of the quantifier.
source: root + jc_sparse_direct_search / 2026-08-22
audit: >
  PENDING.  Two independent derivations found the lower tower (with shifted
  depth indices), and the second derivation also rederived the reflected
  upper tower and the homogeneous charge equations.  A final hostile audit
  of the written theorem, signs, exact transcript, and boundary scope is
  still required before audit promotion.
depends_on:
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
related:
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
  - THM-3721-automorphic-cohn-one-right-shear-nonentry
  - THM-3723-automorphic-cohn-c3-two-right-resonance
  - THM-3725-automorphic-cohn-opposite-two-right-hyperbolic-resonance-nonentry
  - THM-3726-automorphic-cohn-constant-sl2-orbit-broughton-classification
script: 04-computation/jc2_automorphic_cohn_diagonal_binomial_towers_thm3734.py
output: 05-knowledge/results/jc2_automorphic_cohn_diagonal_binomial_towers_thm3734.out
script_sha256: 363fc28b1077ce209fe747b8db3c16d7c699aa9ca2280b9f886e89f0acae11b1
output_sha256: 7d3e372b8eecf143e881a40aa9f939614f04b43d3a8851288d91972593f8616e
semantic_sha256: bd170c83f89be6a626f759285b19a220da42c4eb9f2132876327309de8226ee2
hash_basis: raw LF bytes
---

# THM-3734 -- the diagonal slice is exactly two finite binomial towers

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT AUDIT.**  The constant-
parameter theorem THM-3726 leaves a deliberate boundary: a constant right
matrix followed by a genuinely polynomial exposed left shear.  On the
diagonal right orbit, that boundary is now complete.  Surprisingly, it is not
empty at the first closure equation.  It contains two infinite towers of
smooth noncoordinate polynomials, governed by finite binomial powers and
converging formally to exponential kernels.  None survives the second
closure equation.

Work over a characteristic-zero field containing the displayed square roots
(in particular over `C`).  Set

```text
M_0=[4T^2,2XT-1;1+2XT,X^2],
R_a=diag(a,a^{-1}),                 A=a^2 !=0,
N_a=M_0R_a=[alpha;beta],            det N_a=1,          (1)
z=XT.
```

Here a row `(U,V)` is *closed* when `partial_T U=partial_X V`.  Use
`J(f,g)=f_Xg_T-f_Tg_X`.

## 1. Complete lower-row classification

There is a polynomial `h` for which `beta+h alpha` is closed if and only if
there is an integer `r>=1` such that

```text
A=(r+1)/(2r),
K_r^L(z)=(1+2z/r)^r,
H_r^L(z)=[K_r^L(z)-1-2z]/(4z^2),
h=X^2 H_r^L(z).                                      (2)
```

The quotient in `(2)` is polynomial.  At `r=1`, it is zero, so `(2)` includes
the constant depth-one boundary `A=1,h=0`.  For `r>=2`, `H_r^L` has degree
`r-2`.

The closed row is the gradient of

```text
Q_r^L=a q_r^L,
q_r^L= r[(1+2XT/r)^(r+1)-1]/[2(r+1)T]
     = X Phi_r^L(z),                                  (3)
Phi_r^L(z)=r[(1+2z/r)^(r+1)-1]/[2(r+1)z].
```

Although `(3)` is written as a divided difference, it is a polynomial.

## 2. Complete upper-row classification

There is a polynomial `h` for which `alpha+h beta` is closed if and only if
there is an integer `r>=1` such that

```text
A=r/[2(r+1)],
K_r^U(z)=-(1-2z/r)^r,
H_r^U(z)=[K_r^U(z)+1-2z]/z^2,
h=T^2 H_r^U(z).                                      (4)
```

Again `H_1^U=0`, and `deg H_r^U=r-2` for `r>=2`.  This time the potential is

```text
Q_r^U=a^{-1} q_r^U,
q_r^U= -rT[1-(1-2z/r)^(r+1)]/[2(r+1)z]
     = T Phi_r^U(z),                                  (5)
Phi_r^U(z)=-r[1-(1-2z/r)^(r+1)]/[2(r+1)z].
```

The two diagonal slopes satisfy `A_L(r)A_U(r)=1/4`; the towers are reflected
by the `X/T` charge duality.

## 3. Why these are all the first closures

Direct differentiation of `(1)` shows that the lower closure equation,
after multiplication by `a`, is

```text
4AT^2 h_T -(2z-1)h_X +2(4A-1)Th +2(A-1)X=0.          (6)
```

Give `X^iT^j` charge `i-j`.  The linear operator on `h` in `(6)` lowers
charge by one, while the forcing term has charge `+1`.  Hence the forced
piece of a polynomial solution has charge `+2`, so write it as
`h=X^2H(z)`.  Equation `(6)` becomes

```text
[z+(4A-2)z^2]H' +[2+(8A-6)z]H +2(A-1)=0.             (7)
```

If `H=0`, then `(7)` gives `A=1`, the `r=1` member.  Otherwise let
`n=deg H`.  The constant equation gives `H(0)=1-A`, while the unique top
coefficient is

```text
[4A(n+2)-2n-6] lead(H).                               (8)
```

It must vanish.  Therefore

```text
A=(n+3)/[2(n+2)].                                     (9)
```

Putting `r=n+2` yields exactly `(2)`.  The remaining triangular recurrence
has a unique solution, and substitution of `(2)` verifies it.

Other charge components would have to solve the homogeneous part of `(6)`.
For nonnegative charge `q`, write `h=X^qG(z)`; then

```text
[z+(4A-2)z^2]G' +[q+(8A-2-2q)z]G=0.                 (10)
```

The least `z`-order forces `q=0` and that least order to be zero.  On the
resonance `(2)`, the resulting solution is a scalar multiple of

```text
(1+2z/r)^(-(r+2)),                                    (11)
```

which is not a nonzero polynomial.  For negative charge `-p`, write
`h=T^pG(z)`; its homogeneous equation is

```text
[1+(4A-2)z]G' +(4Ap+8A-2)G=0,                        (12)
```

whose solution on `(2)` has exponent
`-[(r+1)p+r+2]`.  Thus it too has no nonzero polynomial solution.  There are
no homogeneous additions, proving the word *complete* in the lower claim.

For the upper row, the full equation is

```text
A(1+2z)h_T-X^2h_X+2(A-1)Xh+2(4A-1)T=0.              (13)
```

Its forced charge is `-2`.  Substitution `h=T^2H(z)` gives

```text
[Az+(2A-1)z^2]H' +[2A+(6A-2)z]H+2(4A-1)=0.          (14)
```

For nonzero degree `n`, the top coefficient is

```text
[A(2n+6)-(n+2)] lead(H),                              (15)
```

so `A=(n+2)/[2(n+3)]`, which is `(4)` with `r=n+2`.
The reflected homogeneous equations have the same two negative exponents
`-(r+2)` and `-[(r+1)p+r+2]`; hence no polynomial additions occur.  This
proves the upper classification.

## 4. Every first closure is unimodular but has no mate

Equations `(2)--(5)` obey

```text
Phi_r^L+z(Phi_r^L)'=K_r^L,
Phi_r^U+z(Phi_r^U)'=K_r^U.                            (16)
```

These identities prove directly that `dQ_r^L=beta+h alpha` and
`dQ_r^U=alpha+h beta`.  Each gradient is unimodular without a gcd
calculation: its determinant with the untouched complementary row is still
`det N_a=1`.

For a final upper left shear after the lower closure, one would need

```text
J(f,q_r^L)=4(r+2)/(r+1) T.                            (17)
```

The Hamiltonian map against `q_r^L`, which has charge `+1`, raises charge by
one.  Projecting `(17)` to the only relevant charge forces
`f=T^2F(z)`.  Then

```text
J(T^2F(z),q_r^L)
 =-T{z Phi_r^L F' +2K_r^L F}.                         (18)
```

If `F` has degree `m`, the highest coefficient inside braces is

```text
[m+2(r+1)] lead(Phi_r^L) lead(F),                     (19)
```

because `lead(K_r^L)=(r+1)lead(Phi_r^L)`.  It cannot
vanish in characteristic zero.  A nonzero constant right side is therefore
impossible.

For a final lower shear after the upper closure, the required equation is

```text
J(f,q_r^U)=-(r+2)/(r+1) X.                            (20)
```

Now charge forces `f=X^2F(z)`, and

```text
J(X^2F(z),q_r^U)
 =X{z Phi_r^U F' +2K_r^U F}.                          (21)
```

The same coefficient `(19)`, with `Phi_r^U`, excludes `(20)`.  Consequently
there are no polynomials `f,g` for which either

```text
E_+(f)E_-(h)M_0R_a,             E_-(g)E_+(h)M_0R_a   (22)
```

is a Jacobian matrix in the corresponding tower.

## 5. Cyclotomic fibres and the exponential limit

The first-closure potentials have a concrete geometric anatomy.  Over an
algebraic closure,

```text
q_r^L=0: X=0, or XT=(r/2)(zeta-1),
q_r^U=0: T=0, or XT=(r/2)(1-zeta),                    (23)
```

where `zeta` ranges over the nontrivial `(r+1)`-st roots of unity.  Thus the
zero fibre is one affine axis plus `r` disjoint hyperbolas.  The flanks are
simple because the characteristic is zero.  The depth `r=2` member is a
`C_3` fibre, while `r=8` is the corresponding `C_9`/ternary-depth-two fibre.
The recurring cyclic-three patterns are therefore literal finite fibres of
the diagonal Cohn closure equation, not merely analogies.

Coefficientwise as `r` tends formally to infinity,

```text
K_r^L -> exp(2z),             K_r^U -> -exp(-2z),
Phi_r^L -> [exp(2z)-1]/(2z),
Phi_r^U -> [exp(-2z)-1]/(2z).                         (24)
```

So the binomial towers are finite divided-power shadows of two exponential
solutions.  Polynomiality occurs exactly when the closure recurrence
terminates; mate failure occurs at the opposite, highest-degree boundary.

This theorem does **not** prove `JC(2)` and does not exclude non-diagonal
constant right matrices with polynomial exposed parameters, nonconstant
right matrices, or longer alternating words.  It instead identifies the
counterexample-directed escape condition: a successful continuation must
break the single-charge unique-top-edge mechanism, for example by coupling
charges or by producing a second top edge capable of genuine cancellation.

Reproduce the exact positive, hostile, and optimized checks with

```bash
python3 -B 04-computation/jc2_automorphic_cohn_diagonal_binomial_towers_thm3734.py
python3 -B -O 04-computation/jc2_automorphic_cohn_diagonal_binomial_towers_thm3734.py
```

The companion derives both unspecialized closure PDEs; checks the potentials,
debts, resonance perturbations, cyclotomic flank identities, and unique top
edges for `r=1,...,12`; rejects polynomial mate profiles through degree eight
at every tested depth; and verifies the first eight exponential coefficient
limits.  Equations `(6)--(21)` carry the all-depth proof.  **QED.**
