# Symmetry cannot refute a counting conjecture

**Session:** death-star-fc-archimedes-port, 2026-08-03.
**Trigger:** the owner asked for the Archimedes mechanism of THM-3290 to be
ported to the Factorial simplex lane, and supplied three papers as inspiration:
Ivanov--Mikhailov--Wu on `pi_n(S^2)` (arXiv:1506.00952), the Wegner rectangle
refutation (arXiv:2606.17854), and *The Maxwell Conjecture is False*
(arXiv:2607.27197).
**Outcome:** THM-3300, THM-3301.  This reflection records provenance and
reasoning; the theorem files are the truth surface.

## 1. The port succeeded as an identification and failed as a mechanism

The first thing that mattered was reading THM-3018 properly rather than
trusting a remembered summary.  My prior note said the homogeneous Factorial
case was "proved for all `n` by maximum-modulus Laplace".  THM-3018 says the
opposite: that closure is **AUDIT-REQUIRED**, because maximum-modulus
localization does not control oscillatory cancellation at an interior saddle.
So HFC was open on both sides and the port was live territory, not a settled
question.  Correcting that was worth more than an hour of clever work.

The port then produced something better than a mechanism transfer: an
identification.  With `z_i` independent standard complex Gaussians normalized
by `E|z_i|^2=1`,

```text
L_fac(x^alpha) = alpha! = E[ prod_i |z_i|^(2 alpha_i) ],
```

because `|z_i|^2` is exactly `Exp(1)`.  The moment map then sends the uniform
measure of `S^(2n-1)` to the uniform measure of `Delta_(n-1)`, so

```text
HFC(n) = the U(1)^n-invariant part of the Gaussian moment problem in 2n
         real variables.
```

The Factorial lane and the GMC/NC2 lane are not neighbouring problems that
might be bridged.  The first is a **subalgebra** of the second.  That is the
kind of statement worth having even when it kills the plan that produced it.

And it does kill it.  THM-3290's engine is a coefficient extraction in an
angular variable -- available only because the restricted polynomial had a pole
`.../z^2`.  Torus-invariant polynomials have no angular variable, so the
extraction degenerates to the identity and the spherical average is just the
simplex average again.  There is nothing left to be flat.  (A smaller check
makes the same point immediately: THM-3290 lives in real dimension 3, which is
odd, so it is not of the form `2n` at all.)

## 2. Running the mechanism in reverse, and the surprise

The obvious repair is to stop fighting the torus and *use* it: pick `P` with a
nonzero weight, and every moment `L(P^m)` vanishes for free, with no flatness
argument at all.  That looked like a shortcut to a whole new counterexample
family.

It is not, and the reason is worth stating carefully.  If `P` is a
`chi`-eigenvector then `P^m` is a `chi^m`-eigenvector, so `L(P^m)=0` whenever
`chi^m != 1` -- which for an infinite-order character is every `m>=1`.  But for
any fixed `Q`, decompose `Q` into its finitely many isotypic pieces; then
`Q_psi P^m` is a `psi chi^m`-eigenvector and dies unless `psi = chi^(-m)`.
Only finitely many `psi` occur, so only finitely many `m` survive.

```text
L(Q P^m) = 0 for all m >> 0.
```

That is exactly the Mathieu/vanishing conclusion.  **The symmetry mechanism
produces `P` with all moments zero and simultaneously certifies that `P` is not
a counterexample.**  It cannot refute GMC, GVC, FC, or the Image Conjecture --
ever.

The companion makes this concrete rather than abstract: for `P=x^w` and six
test `Q`, the nonvanishing exponents are *exactly* the weight-cancelling ones,
at most one per pair.  And THM-3290's `P`, checked directly, carries six
distinct torus weights `{-2,0,2,4,6,8}` -- it is as far from an eigenvector as
its degree permits.  By the theorem it had to be.

The mechanical reason is now visible.  A character acts on `Q P^m` by the same
multiplication as on `P^m`, so it cannot distinguish them.  THM-3290's collision
distinguishes them because `Q` shifts the coefficient index `delta` while the
flatness order tracks `m`: the *gap* is what is uniform, and `Q` moves it.  Any
counterexample must have that asymmetry between `P` and `QP`.

## 3. The torsion dichotomy, and a link back to my own earlier work

When `chi` has finite order `e`, part (a) only reaches `m` outside `eZ`.  The
simplex's affine automorphism group is exactly `S_n` -- finite -- so on
`Delta_(n-1)` the cyclic route never even achieves total vanishing.  For `n=3`
the surviving obligation is `<g^(3k)>=0`, with the exact evaluation

```text
< (s_1 + omega s_2 + omega^2 s_3)^m > = (n-1)! m!/(m+n-1)! * [n divides m].
```

I then closed the small cases: cyclic eigenvector families of degree 1, 2 and 3
on the triangle are excluded outright (explicitly; by an exact `Q(omega)` gcd;
and by a mod-`p` resultant certificate).

Writing this down, the pattern matched something I proved two sessions ago in a
different lane.  THM-3204's continuant dichotomy says a semisimple transfer has
an eigenvalue ratio of finite multiplicative order and gives *periodic*
vanishing, while the parabolic (unipotent) transfer has no ratio and gives a
single congruence of difference `p`.  That is the same statement in
multiplicative dress: **torsion in the acting character turns total
annihilation into periodic gaps.**  Two lanes, one phenomenon, and I did not
see it until the FC port forced the group-theoretic framing.

## 4. What the three inspiration papers actually contributed

Honest accounting, because it is easy to over-read a reading list.

The Maxwell and Wegner refutations are **illustrations, not evidence** for
anything proved here.  Neither is cited by THM-3300 or THM-3301 and neither is
needed.  What they did was supply a frame: both refute a conjectured bound that
is a *naive count* -- Maxwell's `(n-1)^2` is a Bezout/Morse-type count of
critical points of an electrostatic potential, beaten by 24 with five charges;
Wegner's `tau <= 2 nu - 1` is a piercing-versus-packing count, beaten toward
ratio 2 by 64 rectangles.  Reading them next to my own lane suggested the
question "what does a refuting construction have to *not* be", and Theorem 1 is
the answer in my setting: it cannot be an eigenvector.  The heuristic
generalization -- that extremal counterexamples to counting conjectures tend to
be deliberately asymmetric, because symmetry forces the naive count -- is a
**tangent**, not a theorem, and is recorded as such.  I have proof of it in
exactly one setting.

IMW is the complementary case and was already mined (THM-3204/3205/3210): a
*non*-vanishing result where a residue pattern suggested gaps, whose engine is
the parabolic continuant.  Its role this session was to supply the
torsion-versus-torsion-free vocabulary of section 3.

## 5. Candidate META-PATTERNS card (not promoted)

`META-PATTERNS.md` is at the shared startup byte budget -- there were 49 bytes
free when I routed this work, and the frontier is at exactly its 450-line cap.
Parked here.

> **Compute the isotypic decomposition before hunting a counterexample**
>
> **Trigger:** searching for a counterexample to a statement of the form
> "all moments vanish implies the perturbed moments eventually vanish"
> (GMC, GVC, FC, Image Conjecture, any Mathieu-subspace question).
> **Action:** decompose the candidate under whatever compact group preserves
> the functional.  If it is concentrated in one isotypic component, stop --
> THM-3301 proves it satisfies the conjecture.  A counterexample must mix
> components, and its mechanism must be an order or degree collision, which
> treats `P` and `QP` differently, rather than a character, which cannot.
> **Counterindication:** the argument needs the character to have infinite
> order.  With a finite group of exponent `e` you do not even get total
> vanishing -- only `m` outside `eZ` -- so a finite symmetry is neither a
> counterexample source nor a certificate.
> **Evidence:** THM-3301 (proved, plus the exact exceptional-set prediction on
> 18 pairs) and THM-3300 (the simplex has only `S_n`, so the cyclic route
> stalls at `m in 3Z` for `n=3`); THM-3204 is the same dichotomy in
> multiplicative dress.

## 6. Honest remaining frontier

`FC(n)`, `HFC(n)`, `SFC` and `JC` are untouched -- THM-3300 and THM-3301
constrain *mechanisms*, not truth values.  THM-3018's Laplace closure is still
AUDIT-REQUIRED, so HFC remains open on both sides.  The cyclic exclusions cover
degree at most three on `Delta_2` only; degree 4 has a 5-dimensional eigenspace
and my resultant code handles two parameters, not four, so that case is open
and is the obvious next computation.  The search criterion the bridge yields --
look for torus-invariant members inside known GMC counterexample families -- has
not been run against anything except THM-3290's family, which it excludes by
weight inspection alone.
