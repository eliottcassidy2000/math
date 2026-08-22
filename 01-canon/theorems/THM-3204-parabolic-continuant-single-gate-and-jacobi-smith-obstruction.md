---
id: THM-3204
title: "Parabolic continuant single gate and Jacobi Smith obstruction"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.
  A tridiagonal matrix with unit sub-diagonal has cyclic cokernel, so every
  arithmetic gate it produces is ONE divisibility, never a higher p-rank.
  Its constant form is the Chebyshev minor `U_n(alpha/2)`, and over `F_p` the
  vanishing set of the continuant is an arithmetic progression whose common
  difference is exactly `p` if and only if the step is parabolic
  (`alpha^2=4c`); otherwise the difference divides `p^2-1` and is never `p`.
  Consequently a non-cyclic Smith profile forbids continuant form: THM-3182's
  weighted Gauss--Manin reset `(1,p,p)` is NOT equivalent to any unital Jacobi
  matrix, while the scalar reset `(1,1,p)` is not obstructed.  For THM-3183's
  factorial exterior transfer every single step can be made parabolic, but no
  two CONSECUTIVE steps can, the obstruction being the exact nonzero integer
  `8i^3+20i^2+12i-1`.
audit: >
  The exact symbolic/integer companion verifies the cyclic-cokernel theorem on
  225 unital Jacobi matrices with a non-unital non-cyclic hostile, the deleted
  minor identity, the Chebyshev identity and its parabolic value through size
  11, the parabolic closed form on 78 instances, the complete vanishing-period
  dichotomy over ALL `(alpha,c)` with `c` a unit for `p=3,5,7,11,13`, the two
  reset Smith profiles on a 27-row parameter bank, the reduced factorial
  parabolic condition, the elimination identity that proves the consecutive
  obstruction, and single-step parabolic witnesses whose neighbours stay
  non-parabolic.  Normal and `-O` replay are byte-identical.  Independent
  immutable audit is pending.
source: death-star-imw-pi-n-s2-2026-08-02
depends_on: []
related:
  - THM-3182-factorial-gauss-manin-rank-one-reset-and-two-transverse-smith-bands
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
  - THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall
  - THM-3205-odd-primary-lambda-algebra-engine-and-toda-gate-family
  - THM-1880-the-a-b-functional-frame-chebyshev-pell-companions
external:
  - "S. O. Ivanov, R. Mikhailov, J. Wu, On nontriviality of homotopy groups of
    spheres, arXiv:1506.00952v1; Homology Homotopy Appl. 18 (2016) 337-344.
    Full source entry and convention warnings:
    05-knowledge/reference/CORE-PAPERS-HOMOTOPY.md.
    Its Lemma 3 is the parabolic instance of section 3 below."
script: 04-computation/parabolic_continuant_single_gate_thm3204.py
output: 05-knowledge/results/parabolic_continuant_single_gate_thm3204.out
script_sha256: 55830b22729c56015ecb0e9815ff1ff0798d57ab0f28cfadc6f4fc3216624ece
output_sha256: d3c7b70a8f62fb933184bd0546a7d521ff20cc53edac6f07ae19cb0b93ea0f88
hash_basis: LF-normalized bytes
---

# THM-3204 -- parabolic continuant single gate and Jacobi Smith obstruction

**PROVED + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.**

THM-3183 and THM-3186 put a Jacobi continuant at the centre of the factorial
exterior transfer and left two questions open in an invariant form: why the
two integral frames have different Smith profiles, and what kind of arithmetic
wall a continuant can produce at all.  This theorem answers the second
question completely and the first one structurally.  The trigger was an
external proof (Ivanov--Mikhailov--Wu) whose entire arithmetic content is one
`tridiag(1,-2,1)` determinant.

## 1. Unital Jacobi matrices have cyclic cokernel

Let `R` be a principal ideal domain.  Call an `n x n` matrix `J` a **unital
Jacobi matrix** when it is tridiagonal,

```text
J_(i,i)=alpha_i,  J_(i,i+1)=gamma_i,  J_(i+1,i)=beta_i,               (1)
```

and every sub-diagonal entry `beta_1,...,beta_(n-1)` is a unit.

**Theorem 1.** For `n>=2` the minor of `J` on rows `2,...,n` and columns
`1,...,n-1` equals `beta_1 beta_2 ... beta_(n-1)`.  Hence the first `n-1`
invariant factors of `J` are units, and if `det J != 0` then

```text
coker(J) = R/(det J),   invariant factors (1,...,1,det J).             (2)
```

The cokernel of a unital Jacobi matrix is **cyclic**.

*Proof.* Index the submatrix by `i,j in {1,...,n-1}` with entry `J_(i+1,j)`.
Tridiagonality makes it zero unless `|i+1-j|<=1`, that is unless
`j in {i,i+1,i+2}`; so the submatrix is upper triangular, and its diagonal
entry at `i` is `J_(i+1,i)=beta_i`.  Its determinant is `prod beta_i`, a unit,
so the gcd `d_(n-1)` of all `(n-1)`-minors is `1`, forcing
`d_1=...=d_(n-1)=1`.  The product of all invariant factors is `det J`. QED

The same argument with rows `1,...,n-1` and columns `2,...,n` gives the
statement for unit **super**-diagonal.

**Consequence (single gate).** Any obstruction computed as the cokernel of a
unital Jacobi matrix is governed by one number, `det J`.  It can never be a
`(Z/p)^2` obstruction, and localizing at `p` gives a kernel and cokernel of
dimension exactly `min(1, v_p(det J))` over the residue field, with the whole
`p`-torsion cyclic of order `p^(v_p(det J))`.

## 2. The converse is an obstruction

Smith invariant factors are invariants of the equivalence class `UJV`,
`U,V in GL_n(R)`.  So:

**Theorem 2.** If `M in M_n(R)` has `d_(n-1)(M)` a non-unit, then `M` is not
equivalent to any unital Jacobi matrix -- in particular to no continuant
transfer, companion form, or three-term-recurrence matrix with unit
sub-diagonal.

This settles the invariant status of THM-3182's two frames.  At `n=p-1`,
`d=p+s`, with `svDelta` a unit, that theorem's reset profiles are

```text
scalar   S_(p-1): (1,1,p)   so d_2 = 1  -> NOT obstructed,
weighted G_(p-1): (1,p,p)   so d_2 = p  -> OBSTRUCTED.                 (3)
```

Hence THM-3183's observation that the *scalar* frame carries a Jacobi block
`K_n` while the weighted Gauss--Manin frame does not is **not** an artefact of
the chosen coordinates: no coordinate change over the DVR can give `G_(p-1)`
a unit-sub-diagonal tridiagonal form.  The "one elementary modification of the
output lattice" recorded in THM-3183 section 2 is exactly the failure of
`d_2=1`.  The companion verifies `(3)` on a 27-row bank of `(p,s,v)`.

## 3. Chebyshev normal form and the parabolic dichotomy

For constant band data the leading principal minors of `tridiag(1,alpha,1)`
obey `D_n=alpha D_(n-1)-D_(n-2)`, `D_0=1`, `D_1=alpha`, so

```text
det tridiag(1,alpha,1)_n = U_n(alpha/2),                              (4)
```

`U_n` the Chebyshev polynomial of the second kind.  More generally, for
sub-diagonal/super-diagonal product `c` and diagonal `alpha`,

```text
D_n=alpha D_(n-1)-c D_(n-2),   discriminant  Delta=alpha^2-4c.        (5)
```

Call the step **parabolic** when `Delta=0`.

**Theorem 3.** Let `p` be an odd prime, `alpha,c in F_p`, `c != 0`.  The set
`{n>=0 : D_n=0}` is a nonempty arithmetic progression `{q-1,2q-1,3q-1,...}`,
and

```text
q = p                      if Delta = 0   (parabolic),
q = ord(x_+/x_-) | p^2-1   if Delta != 0,                             (6)
```

where `x_+,x_-` are the roots of `x^2-alpha x+c`.  Since
`gcd(p,p^2-1)=1`, the two cases are exclusive: **the vanishing progression of
a unital constant continuant has common difference `p` if and only if the step
is parabolic.**

*Proof.* If `Delta=0` the double root is `r=alpha/2`, which is nonzero because
`c=r^2!=0`, and `D_n=(n+1)r^n`; its zeros are `n = -1 mod p`.  If `Delta!=0`
then `x_+ != x_-` are nonzero and
`D_n=(x_+^(n+1)-x_-^(n+1))/(x_+-x_-)`, so `D_n=0` iff `rho^(n+1)=1` for
`rho=x_+/x_-`.  If `Delta` is a nonzero square then `rho in F_p^*` and
`ord(rho) | p-1`; if `Delta` is a nonsquare then `x_-=x_+^p`, so
`rho=x_+^(1-p)` and `rho^(p+1)=x_+^(1-p^2)=1`, giving `ord(rho) | p+1`.  In
both cases `ord(rho)` divides `p^2-1` and is at least `2`. QED

The parabolic locus of `(4)` is `alpha = +-2`, the two endpoints of the
Chebyshev interval, where `U_n(+-1)=(+-1)^n(n+1)` attains its extreme value
and the trigonometric solution degenerates to a linear one.  The companion
verifies `(6)` **exhaustively** over every `(alpha,c)` with `c` a unit for
`p=3,5,7,11,13`.

**Reading.** A "`p | n + const`" gate is the *signature* of a parabolic
continuant.  If a lane's arithmetic walls are not congruences of that shape,
its transfer is not uniformly parabolic, and no reparametrization of a
Jacobi/continuant model will make them congruences.

For comparison, THM-1880's transitive-tournament Chebyshev--Pell frame has
step matrix `[[1,x],[x,1]]` with discriminant `4x^2`, so its only parabolic
member is `x=0`, where the transfer is the identity.  That frame therefore
contains **no** nontrivial parabolic member, and its cotangent spectra are the
generic (non-degenerate) branch of `(6)`.

## 4. The factorial exterior transfer is nowhere consecutively parabolic

Take THM-3183 section 3 and THM-3186 section 2 notation:

```text
a_i=2(i+1)(2i+1)v,  b_i=i(i+1)Delta,  Delta=1-4dv,
alpha_i=a_i d,      beta_i=b_i d,
C_r=alpha_(n+r)C_(r-1)+d beta_(n+r)C_(r-2).                           (7)
```

The step discriminant is `alpha_i^2+4 d beta_i = 4d^2(i+1)P_i` with

```text
P_i(v,d)=(i+1)(2i+1)^2 v^2-4 i d v+i.                                 (8)
```

So with `d != 0` step `i` is parabolic exactly on the conic `P_i=0`.  Each
single step is achievable: for `i>=1` and any `v != 0`,

```text
d=[(i+1)(2i+1)^2 v^2+i]/(4 i v)                                       (9)
```

makes step `i` parabolic, and the companion checks that step `i+1` then is
not.  That is no accident.

**Theorem 4.** For every `i>=1`,

```text
(i+1)P_i - i P_(i+1) = -(8i^3+20i^2+12i-1) v^2.                      (10)
```

Since `8i^3+20i^2+12i-1 >= 39 > 0` for `i>=1`, the system `P_i=P_(i+1)=0`
forces `v=0`.  THM-3183's standing hypothesis makes `2vDelta` invertible, so
**no two consecutive steps of the factorial exterior transfer are
simultaneously parabolic.**

*Proof.* `P_i` is affine in `d` with coefficient `-4iv`, and its constant term
is `i`; in `(i+1)P_i-iP_(i+1)` both the `d`-terms and the constant terms cancel
identically, leaving
`v^2[(i+1)^2(2i+1)^2-i(i+2)(2i+3)^2]`, and expanding gives `(10)`. QED

**What this closes and what it does not.** It closes the hope of importing the
single-congruence mechanism of section 3 into the factorial lane by choosing
parameters: the closed form `C_r ~ (r+1)` needs a *uniform* parabolic step, and
two consecutive parabolic steps are already impossible.  It does **not** prove
that the factorial walls have period dividing `p^2-1`; `(6)` is a
constant-coefficient statement, and `(7)` has `i`-dependent coefficients.  It
does explain why THM-3176/3183's wall polynomials `H=24p+109`,
`J=256p^4-27648p^3-365600p^2-1528800p-2096649` and `K` are genuine arithmetic
polynomials rather than congruence conditions, and it says that looking for a
`floor(s/2)` staircase law of the form "`p` divides depth plus a constant" is
looking for a parabolic signature the transfer does not have.

## 5. Scope and exact evidence

Proved: Theorems 1--4 and the identities `(2),(4),(6),(8),(10)`.  Not proved:
any statement about `NC(2)`, `GMC(2)`, `JC(2)`, `LRC(14)`, the empirical
`floor(s/2)` Euclidean-depth staircase, or the non-constant analogue of `(6)`.

Run

```text
python 04-computation/parabolic_continuant_single_gate_thm3204.py
python -O 04-computation/parabolic_continuant_single_gate_thm3204.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
exact integer and symbolic arithmetic only.  It carries one non-unital
non-cyclic hostile (`invariant factors (3,3)`), an exhaustive period sweep
rather than a sample, and positive controls for the single-step parabolic
witnesses.  There is no floating point, random sampling, imported executable,
or assertion-sensitive test.

**QED.**
