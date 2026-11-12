# Correct the representative before correcting its numerical value

**Session278 cross lane. PROVED identities / FINITE-EXACT examples.**
No claim that JC(2), Chan–Pak, or a general exact-correction problem is solved.
The new contribution is a defect compiler and its quantitative application, not
novelty of integration by parts, Laguerre polynomials, or Taylor's theorem.

## Inheritance and the surprising connection

Source pin: `e07f05b9fee2`. S25's
`kps_s25_euler_field_bridge_law_20261111/README.md`, Lemma A1, supplies the
exact-Euler identity `[W^m h]=(D)_m[h]` **when** `xi(W)=W` and
`div(xi)=delta`, with `D=xi+delta`. We retain that named input, not re-audit
its JC applications. Its hostile Broughton examples show why quotient powers
must not be manipulated as powers in a quotient algebra.

Two older underused connections are:

* `07-reflections/single-edge-zero-mode-resolvent-update-codex-20260803.md`
  and THM-3305: a low-rank update needs the active observer and a finite head,
  not just the tail denominator. Here the exact correction needs the actual
  defect multiplication operator, not just its scalar class.
* `07-reflections/watson-estimates-for-gmc2-and-the-laguerre-polya-boundary-boxeph-S202.md`:
  factorial corrections can be asymptotic rather than convergent; cancellation
  invalidates top-term dominance. The Laguerre polynomials below are ordinary
  finite Laguerre polynomials, not that reflection's Laguerre–Pólya class.
* THM-4757, `THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`,
  retains all admissible signed insertions in a Hamiltonian path. It is an
  inherited exact insertion expansion, not a new result here. Its analogue
  here is a **typed** defect insertion: composition with D cannot be replaced
  by commuting scalar multiplication. This is a method bridge, not an
  identification of Gaussian moments with tournament path counts.

Session277 recovers exact integers by interval/residue intersection. This lane
improves the interval *before* applying any reconstruction: change an observable
by an exact annihilator, then apply its numerical correction in the better
representative. The example shows over 5,000-fold improvement at fixed Taylor
order. Neither an integral nor its finite decimal display is called an integer.

## 1. The defective-Euler compiler

Work over a characteristic-zero field in a polynomial algebra with a fixed
volume form. Let

    M_W = {div(v) - v(W): polynomial vector fields v}.

This is a vector subspace, **not assumed an ideal**. Brackets denote its
vector-space quotient. Fix a vector field xi of constant divergence delta,
and put `D=xi+delta`, `R=xi(W)-W`. Multiplication by R is an operator on
polynomials; we never assert it descends to the quotient.

**Theorem 1.** Set P_0=identity, P_1=D-R, and, for m>=2,

    P_m = P_(m-1) o (D+(m-1)I-R) + (m-1) P_(m-2) o R.        (1)

Then `[W^m h]=[P_m h]` for every polynomial h and m>=0.
For R=0 this recovers the inherited rising factorial. For R nonzero, (1)
is an exact finite correction, without a smallness hypothesis.

**Proof.** Apply the defining annihilator to the vector field
`W^(m-1)h xi`. The product rule gives

    div(W^(m-1)h xi)-W^(m-1)h xi(W)
      = W^(m-1)(D+m-1-R)h
        +(m-1)W^(m-2)Rh-W^m h.                              (2)

Taking classes and induction proves (1). The m=1 case has no second term.
Crucially, induction is applied to an actual polynomial argument, not to
multiplication of quotient classes. QED.

This proof supplies a literal witness, not just a membership assertion.
Define V_0(h)=0 and

    V_m(h) = -W^(m-1)h xi
             +V_(m-1)((D+m-1-R)h)+(m-1)V_(m-2)(Rh).          (3)

Omit the last term for m=1. Then
`W^m h-P_m h = div(V_m(h))-V_m(h)(W)` identically. The supplied program
checks these full polynomial identities, with rational coefficients.
Complexity may grow with m and the support of R; no uniform bounded-state
compiler or efficient general JC decision algorithm follows.

## 2. Homogeneous defects: finite signed Laguerre correction

Suppose W=W_0+z, xi(W_0)=W_0, xi(z)=2z, and xi(h)=e h. Put a=e+delta.
Then R=z and:

**Theorem 2.** For every m>=0,

    [W^m h] = [ h sum_(j=0)^m (-1)^j binom(m,j)
                            (a+j)_(m-j) z^j ].              (4)

Equivalently the polynomial on the right is `m! L_m^(a-1)(z)h`.
The special-function name is classical; (4) uses it as an exact correction
inside the full deformed Stein quotient. It is not a Gaussian assertion
about the undeformed potential W_0.

**Proof.** In (1), D(z^j h)=(a+2j)z^j h. Substituting the coefficients
in (4) verifies (1), using Pascal's identity; m=0,1 initialize it. Alternatively
its exponential generating series is

    (1-t)^(-a) exp(-zt/(1-t)),

whose coefficient recursion is (1). This is a formal power series identity;
no convergence or integration interchange is used. QED.

For a=1, the first polynomials are

    1;
    1-z;
    2-4z+z^2;
    6-18z+9z^2-z^3;
    24-96z+72z^2-16z^3+z^4.

These signed finite corrections are exact at every parameter value. Their
terms need not decrease in magnitude. Alternating coefficients alone do not
license an alternating-series remainder bound on (4).

### Several homogeneous defects at once

The same compiler handles interactions rather than adding independent error bars.
Suppose `W=W_0+sum_i z_i`, `xi(W_0)=W_0`, `xi(z_i)=lambda_i z_i`,
and `xi(h)=e h`. Then the following formal series is exact coefficientwise:

    sum_(m>=0) [W^m h] t^m/m!
      = [h (1-t)^(-a)
           exp(sum_i z_i (1-(1-t)^(1-lambda_i)))].            (8)

Only finitely many monomials in the z_i contribute at each t degree. Mixed
coefficients retain all cross-defect interactions; independent scalar error
corrections would discard them. Formula (8) is valid for arbitrary
characteristic-zero scalar weights lambda_i, using formal binomial series.

For a short proof, apply the annihilator to `h exp(tW) xi`. This gives the
formal differential system

    d_t F_h(t) = a/(1-t) F_h(t) - sum_i (lambda_i-1) F_(z_i h)(t).

Here `F_h=sum_m [W^m h]t^m/m!`, and the eigenweight of z_i h is
`e+lambda_i`. Differentiating the right side of (8) verifies this entire
system, including every multiplied argument z^J h. Its initial conditions
are `[z^J h]`; coefficient induction uniquely determines the formal solution.
This proof never multiplies an equality of quotient classes by z_i.

For lambda_i=1 that perturbation contributes zero Euler defect, exactly as
required. For integer lambda_i>1, each factor
`(1-t)^(1-lambda_i)-1` has nonnegative coefficients; thus the coefficient
of a fixed defect monomial z^J in (8) has sign `(-1)^|J|` when a>0.
This is an exact grading of the corrections, not a sign for the sum of the
omitted grades.

## 3. A second, independent correction axis with rigorous signs

Let X,Y be independent standard real Gaussians, g>=0, and

    W_g=(X^2+Y^2)/2+g X^4,
    L_g(f)=E[f(X,Y) exp(-g X^4)].

The full density is proportional to exp(-W_g); polynomial divergences
integrate to zero because the boundary terms vanish. Set xi=(x dx+y dy)/2,
so delta=1, z=g x^4. Theorem 2 gives

    L_g(W_g^m)=L_g(m! L_m(z)).                              (5)

For any nonnegative polynomial f (or nonnegative function with all displayed
Gaussian moments through order N+1 finite), define

    S_N(f)=sum_(j=0)^N (-g)^j E[f X^(4j)]/j!.

Pointwise Taylor's theorem for exp(-u), u>=0, gives the exact sign and bound

    0 <= (-1)^(N+1)(L_g(f)-S_N(f))
       <= g^(N+1) E[f X^(4N+4)]/(N+1)!.                    (6)

For fixed nonzero nonnegative polynomial f, the upper constant 1 in (6)
is sharp as g tends down to zero: the ratio of the signed remainder to
the displayed next-term bound tends to 1 by dominated convergence.
Thus S_N and S_(N+1) bracket L_g(f). For a signed even polynomial, split
its monomials into positive and negative coefficients and combine their
bounds with reversed endpoints for negative coefficients. This is rigorous
but does not pretend that the individual moment errors are independent or
that the resulting bound is globally optimal.

All Gaussian moments used are exact integers:
`E[X^(2k)]=(2k-1)!!`, odd moments zero. Therefore for rational g every bound
in (6) is rational and reproducible without floating integration.

At g=1/1000 and N=14, compare two certificates for the SAME integral:
(A) apply (6) directly to the nonnegative W_g^m;
(B) first use (5), then bound the signed polynomial monomial by monomial.
The exact rational intervals and decimal displays are in `output.txt`.

| m | direct interval width | Euler-reduced width | width improvement |
|---:|---:|---:|---:|
|1|7.7845020e-16|1.0820062e-16|7.1945075|
|2|2.8203527e-14|7.6203501e-16|37.010803|
|4|4.1461464e-11|7.5038797e-14|552.53370|
|6|7.0473872e-8|1.3922049e-11|5062.0331|

In these examples the reduced interval is contained in the direct interval;
this is a finite observation, not a theorem for all parameters or orders.
The reduction removes expensive radial powers in favor of the actual small
quartic defect. It improves the *absolute* error, not merely a relative display.

**Hostile boundary.** For f=1 the ratio of successive perturbation magnitudes is

    g(4j+1)(4j+3)/(j+1),

which tends to infinity for every g>0. Increasing N indefinitely does NOT
converge at fixed g. The compiler (1) is finite and exact; the independent
perturbation series is generally divergent. One may choose the best certified
finite bracket or switch to a different enclosure, but cannot claim arbitrary
precision from (6) alone.

## 4. A general certificate-intersection principle

Let L be a linear functional, K a known subspace of ker(L), and
f_1,...,f_r representatives with f_i-f in K. Suppose rigorous enclosures
[L_i,U_i] are available for each L(f_i). Then

    max_i L_i <= L(f) <= min_i U_i.                          (7)

This elementary principle is **not new**. The work lies in compiling useful
representatives and certifying K. Equations (1)-(5) do so nontrivially here.
For an integer-valued target, combining (7) with Session277's exact residue
and a uniqueness check is valid. Here L_g is real-valued, so no integer
rounding conclusion is drawn.

The corresponding tournament task is to find a physically realizable family
of interfaces giving exact annihilator/telescoping relations while retaining
the prescribed target. THM-4757 gives a signed forest expansion and the other
Session278 lane treats actual block states. A Gaussian correction polynomial
is NOT a construction of a tournament and its negative coefficients cannot
be interpreted as negative numbers of physical witnesses.

## 5. What remains open, and what can be used now

* Proven usable instrument: exact defective-Euler reduction plus explicit
  divergence witnesses, including nonhomogeneous defects through (1).
* Proven quantitative architecture: rewrite by exact annihilators, then
  intersect signed certified residual bounds in multiple representations.
* Demonstrated benefit: >5,000-fold reduction of an absolute interval width
  at the same finite expansion order for a specified tilted Gaussian moment.
* For Keller potentials, the native Euler field already has R=0. A nonzero
  defect arises after choosing an approximation, projection, or a different
  field. One must prove that choice improves the target while retaining the
  necessary source algebra. No implication from (1) to JC(2) is claimed.
* It remains open whether these corrections give a uniform tame representative
  for the JC obstruction or physically reusable Chan–Pak target interfaces.
  Finite identities and asymptotic corrections do not supply those quantifiers.

Reproduce: `python3 .../cross/probe.py`. The script checks exact polynomial
certificates through m=6 for quartic and m=4 for sextic perturbations, the
Laguerre formula on three eigenpolynomials, rational signed bounds, and
intersection consistency. No external numerical quadrature is used.
