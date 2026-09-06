# Independent audit of the coefficient-prism response theorem

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS; FINITE-EXACT independent
reconstruction.** The complete [proof](long_frontier2_sep06_rectangle.md),
[producer](../../04-computation/long_frontier2_sep06_rectangle.py),
[output](long_frontier2_sep06_rectangle.out), and
[coefficient certificate](long_frontier2_sep06_rectangle_certificate.json)
pass without a mathematical correction. The no-interlacer scope, repeated
beta-root boundary, original-zero condition, and nonvacuity are valid.

## 1. Types and complete original carriers

The actual object proved here is the stated real coefficient model. Its
parameters a,b,f are not silently identified with arbitrary factorial
parameters. The conclusion is Q(-s)<0 at original zeros of P. It does not
assert negativity at arbitrary phases, and it uses neither the real-rooted
beta predicate nor interlacing in the initial coefficient-prism theorem.

An independent source formed O,E and the three ordinary carrier numerators
directly. It expanded the products `t*O^2+E^2` and
`t^2*(beta^2+2*t*c*d)`, then took their matching coefficients with the
appropriate distinct shifts. This avoids the producer's array convolutions
and binomial-28 shortcut. The degree ranges give exactly Q indices -1
through 8. The coefficient at -1 is 28 and contributes -28 to s Q(-s).
No negative-support term is discarded.

The ordinary first carrier gives P(-s)=2002 p(s) with the printed quartic.
In particular P has degree four when f>0 and degree three when f=0;
its constant coefficient is nonzero. Elimination at s>0 divides only by
a nonzero phase and exactly substitutes the original equation for f.
The independent response reconstruction recovers all 19 monomials of H
and its complete multidegree (2,2,8).

## 2. Exhaustion of every original zero

At each fixed endpoint, p is affine in all three parameters. Its extrema
on the closed prism therefore occur at the eight corners. All 56 strict
endpoint signs were independently recomputed. For f>0 the last interval
has its positive sign at infinity from the positive quartic leading
coefficient. Four disjoint sign changes give four distinct positive roots
of p and exhaust its degree. Every root is therefore simple, and no
complex or additional real root remains. Under f=0, the first three sign
changes exhaust the actual cubic degree because b>0. No finite fourth
phase is inferred at that boundary.

The response charts are also rigorous on the closed coefficient rectangle.
Each finite transform is a multihomogeneous Bernstein expansion after
clearing positive denominators. Its endpoint limits retain the corresponding
highest-power faces, and every such face has strictly positive coefficients.
Thus a coordinate endpoint at infinity in a chart cannot create a hidden
zero. The tail transform has the same positive rectangle faces and a
strictly positive polynomial in w>=0. These charts cover all the original
phases just established. Because s>0, the identity
`s Q(-s)=-(14/11)H` has the required strict sign.

The independent coefficient path used the binomial tensor formula

```text
[u^i](lo+hi*u)^p(1+u)^(D-p)
  =sum_j binom(p,j) lo^(p-j) hi^j binom(D-p,i-j),
```

and the ordinary translation formula for the tail. Exact Fraction arithmetic
recovered all 324 saved chart coefficients individually and checked each
strictly positive. This is a complete polynomial identity check, not a
sampling check on chart values.

## 3. Full beta-root scope, including repeated and zero roots

Assume B has five nonnegative real roots, counted with multiplicity.
Rolle's theorem with multiplicities puts all four derivative roots in the
interlacing intervals. The derivative at zero is b>0, so its least root
gamma is strictly positive even if B has a zero root. Since the derivative
at 2/5 is uniformly negative, some derivative root lies in (0,2/5), hence
the least one does too. It lies in [beta1,beta2]. The quintic product is
nonnegative on that interval; if the first two roots coincide it is zero
there. Consequently B(gamma)>=0 in every repeated-root case as well.

The parameter inequalities correctly give F_(a,b)<=F_upper on [0,2/5]
and F'_(a,b)(2/5)<=F'_upper(2/5)=-81/10. The independent exact Bernstein
calculation recovers all six positive coefficients of 5-F_upper. Positivity
also holds at 2/5 by its positive leading coefficient. Therefore

```text
0<=f<=F_(a,b)(gamma)<=F_upper(gamma)<5.
```

The first inequality is the product of the nonnegative beta roots. This
places the entire admissible beta domain in the already-proved prism,
without C or D hypotheses. The proof does not borrow a tail estimate
whose hypotheses require interlacers.

The separate smaller prism is genuinely three dimensional: the affine
corner values at 0,1/10,1,3,5,7 have six strict alternating signs for
every a,b in the closed rectangle and f in [1/2,3/2]. The resulting five
positive roots exhaust the quintic degree and are distinct. Thus the
region contains a nonempty open set of positive-root B polynomials, not
just one identified centre or a repeated-root surface.

## 4. Hostiles and frozen reproduction

The f=6 control retains an original zero, and the independent interval
transform recovers all nine negative H coefficients there. Its response
is positive and it is outside the nonnegative-root beta domain by the
proved product barrier. At the genuine centre f=1, s=4 is independently
checked to be off the original zero, and the full ordinary carrier gives
`4Q(-4)=350398552675052>0`. These controls test different necessary
conditions and are not conflated with admissible counterexamples.

The [independent audit source](../../04-computation/long_frontier2_sep06_rectangle_audit.py)
and [output](long_frontier2_sep06_rectangle_audit.out) have **500
always-active gates**. Normal and optimized independent runs agree. Each
also reruns the 457-gate producer normally and under optimization, checking
both raw output and regenerated certificate against the frozen bytes.

```text
producer source 4fb5fe880fbc030b3461889d021f7fa24b6376b23cdc1251afd99b0327591ea5
producer output 572d5f8192b0a567efd95b6db25d23f8dcd7277477dacfd27d75a98134668bad
certificate cd1af1f580b51fc1b556a7cdf44f28203e464a3ba41e3d03e55363278816d412
audit source 85f79d0371b952415803d6a1691aec9984db5e5c923113972e5f06f0c2ff150c
audit output 214608be5b58dcf48837c3c505e87084289390ad04008296b77d41c75fa560da
```

No numerical root decision, broad census, or producer import supplies the
independent reconstruction. The candidate is safe to promote within its
declared model scope; maximality of the prism and general actual Laurent
transport remain outside the result.
