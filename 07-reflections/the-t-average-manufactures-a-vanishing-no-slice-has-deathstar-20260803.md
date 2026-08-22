# The t-average manufactures a vanishing that no slice has

**Session:** death-star-gvc3-counterexample, 2026-08-03.
**Trigger:** the owner supplied a three-variable "GVC(3) counterexample" as a
bare list of polynomials plus a closed form, together with a two-block total
variation inequality, and asked what they unlock.
**Outcome:** THM-3290, THM-3291, `05-knowledge/reference/CORE-PAPERS-GMC-TV.md`.
This reflection records provenance and reasoning, not truth; the theorem files
are the truth surface.

## 1. The first move was a change of category, not a computation

The supplied object was framed as a differential-operator statement:
`Delta^(6m)(P^m)=0` for `Delta=4d_x d_y+d_t^2`.  Verifying that is finite work
and I did it first (it checks, exactly, to `m=6`).  But the useful observation
came from asking what `Delta` *is*.

`rho=t^2+xy` is precisely the quadratic form whose Laplacian is `Delta` — under
`x=u+iv`, `y=u-iv` the operator is `d_u^2+d_v^2+d_t^2` and `rho=u^2+v^2+t^2`.
So on forms of degree `2k`,

```text
Delta^k f / (2^k k!) = E[f(G)],   G ~ N(0, I_3).
```

**A GVC statement about a power of the Laplacian in `n` variables is a GMC
statement in `n` variables.**  It is not a separate lane.  The supplied object
belongs to the repository's Gaussian-moments lane, next to THM-2022, THM-1490
and THM-1645 — and the repository already knew GMC is false for every `n>=3`.
So the object refutes nothing new.  Everything of value had to come from
*why* it works.

## 2. Why it works: two orders in collision

Radial-spherical split turns `L` into a spherical mean.  On `S^2`, **Archimedes'
theorem** says the uniform measure is the product of uniform `t` on `[-1,1]`
with uniform argument of `z=u+iv`.  The argument average is a *constant-term
extraction*, an exact annihilator.  Substituting `zbar=(1-t^2)/z`, the supplied
identity `xC=rho^3-t^2A^2` becomes, on `rho=1`,

```text
zC = 1 - t^2 A^2,     A = 1+z^2,
```

so `P=AC^2 = A(1-t^2A^2)^2/z^2` is a function of `w=z^2` and `t` alone.  With
`a=1+w`, the whole spherical average is one coefficient extraction followed by
one `t`-integral, and the `t`-integral is an antiderivative

```text
F_n(s)=int_0^s (1-u^2)^n du,     F_(2k)'  vanishes to order 2k at s=1.
```

The answer is a single generalized binomial:

```text
<x^(2delta) R_nu^k> = C(k-nu, k-delta) * 2^(2k)(2k)!/(4k+1)!!.
```

It vanishes because a polynomial prefactor of degree `k-nu` is annihilated by a
derivative of strictly higher order `k-delta`, the antiderivative being flat
there.  **The vanishing is a collision of a degree against a flatness order.**
It is not a cancellation among terms — the companion confirms directly that the
`z`-constant term is a *nonzero* polynomial in `t` whose `t`-average is zero.

That last sentence is the sentence I would keep if I could keep only one:
**the `t`-average manufactures a vanishing that no fixed slice has.**

## 3. What that buys

Three things, none of which were available from verification alone.

1. **Proof for all `m`.**  The supplied closed form
   `2^(8m+1)(6m+1)!(2m)!(12m+3)!!/(4m+1)!!` is now derived, not sampled, and
   the `m=0` exception is explained (`x^2` is harmonic because `x` is an
   isotropic linear form).
2. **A sharp threshold.**  `<R_nu^k>=0` exactly when `k>=nu`, with the
   sub-threshold value `(-1)^k C(nu-1,k) F_(2k)(1)`.  I predicted the threshold
   from the mechanism before computing it, and the direct polynomial algebra
   reproduced `-8/15`, `-16/15`, `128/315` on the nose.  That agreement is the
   strongest evidence the mechanism is the real one and not a re-description.
3. **An infinite family.**  `P_nu=R_nu^nu`, `Q_nu=x^(2nu)`,
   `D=Delta^((4nu+2)nu)` — GVC fails in three variables for
   `Delta^6, Delta^20, Delta^42, ...`.  The supplied object is `nu=1`.

## 4. The dimension boundary, read off the mechanism

THM-1645 already records that the GMC(2) *angular* layer is
Duistermaat–van der Kallen and the gap is purely radial.  The mechanism above
says what happens one dimension up, and it is not "more of the same":

- `n=2`: the sphere is a circle, the spherical mean **is** the `z`-constant
  term of a one-variable Laurent polynomial, and DvdK applies.  There is no
  integration left to be flat.
- `n=3`: the spherical mean is that same constant term **followed by an
  Archimedes `t`-integral**, and the antiderivative's flatness is exactly what
  produces the counterexample.

So the `2`/`3` boundary of GMC is the appearance of an Archimedes coordinate.
Concretely this is a no-go for a natural proof strategy: **you cannot prove
GMC(3) by slicing to DvdK**, because on every fixed slice `t=t_0` the DvdK
hypothesis *fails* — the constant term is nonzero — and only the average
vanishes.  Any future attempt in that direction should be checked against this
first.

## 5. One probe run and refuted before it was recorded

Having the threshold, I guessed that the minimal witness degree grows with
`nu`: `Q=x^(2delta)` fails for `delta<nu`, so perhaps no low-degree `Q` works.
A monomial sweep at `nu=2` (where `P_2` has degree 40) killed it in one run:
`x^2` does fail exactly as the formula predicts, but `y^2`, `xy` and `t^2` all
give a nonzero moment.  The minimal witness degree is `2` for both `nu=1` and
`nu=2`.

The lesson is narrow and worth stating: **the master formula is a statement
about the `x`-direction only.**  Its threshold indexes how far the `x`-shift
must go, not how large a general witness must be.  The caveat is now in
THM-3290 section 6 so the formula cannot be over-read later.  (Witnesses odd in
`t` always vanish, since `P_nu` is even in `t` — that part is free.)

## 6. The second supplied item is a contrast, not a bridge

The TV homogenization inequality `delta_N <= delta_I+delta_J-delta_I delta_J`
is, in the two-single-observation-block case, exactly a box constraint plus
AM-GM once you have the closed form `TV(Bin(2,x),Bin(2,y))=|x-y|(1+|x+y-1|)`.
Its equality locus is rigid and I classified it: same sign, **equal** gaps, and
every block pinned to one common face of the unit square — 62 predicted points
and 62 actual on the complete `17^4` grid, matching both ways.

The temptation is to pair it with THM-3290 as "two theorems about orbit
averages."  That is a rhyme.  The honest statement is a **type** distinction:

- `TV` is a *positive* functional; homogenizing can only contract it, and the
  contraction is exactly supermultiplicative in the affinity `1-delta`.
- The Gaussian moment is a *signed* functional; averaging can annihilate every
  moment while leaving a perturbation visible.

No map, no preserved predicate, no transfer — so it is filed as a contrast in
THM-3291 section 4 and nothing is claimed across it.  What does carry is the
warning itself, which is why it is proposed as a method card below.

## 7. A probe that did not run: the Factorial lane

The same radial split applies verbatim to the factorial functional
`L(x^alpha)=alpha!`: for `f` homogeneous, `L(f)=((n+d-1)!/(n-1)!) A(f)` with
`A` the *simplex* average, which is the repository's own THM-3018 framing.  So
the shape of the mechanism ports.  The engine does not.  On `S^2` the inner
average is a constant-term extraction — an exact annihilator supplied by a
compact group acting transitively on the fibres.  On the simplex the Dirichlet
decomposition gives an honest `[0,1]` integral with no group and no
annihilator, so there is nothing to make flat against.  **The obstruction to
porting Archimedes flatness to FC is the absence of a transitive compact
symmetry on the simplex fibres**, and it is recorded here so a later session
does not spend the same hour discovering it.

## 8. Candidate META-PATTERNS card (not promoted)

`META-PATTERNS.md` is at the shared startup byte budget; the card is parked
here.  It has evidence from two distinct threads.

> **Identify the functional before attacking the operator**
>
> **Trigger:** a claim of the form `D^m(P^m)=0` for a constant-coefficient
> operator `D`, or any "all moments vanish" hypothesis.
> **Action:** ask which quadratic form `D` is the Laplacian of, then split
> radial from spherical.  The statement usually becomes a moment statement in
> an existing lane, and the spherical mean often collapses to a coefficient
> extraction plus one interval integral.  Vanishing is then a *degree versus
> flatness-order* comparison, not a cancellation, and it generalizes by
> retuning either exponent.
> **Counterindication:** the collapse needs a compact group acting
> transitively on the angular fibres to supply an exact annihilator.  Simplex
> or orthant functionals (the Factorial lane) have no such group; the shape
> ports, the mechanism does not.  Also: a threshold in one distinguished
> direction is not a bound on general witnesses — sweep other monomials before
> claiming minimality.
> **Evidence:** THM-3290 (Gaussian/`S^2`: the `t`-average annihilates every
> moment; threshold predicted then confirmed) and THM-3291 with THM-1645 (the
> positive/signed split, and the `n=2` DvdK layer that has no `t` to average).

## 9. Honest remaining frontier

Nothing here touches `JC`, `JC(2)`, `DC(2)`, or THM-1435's VC-witness
dimension bracket, which is about the plain Laplacian and is untouched by
statements about its powers.  `GMC(2)` remains proved and is not in tension.
`FC(3)`/`SFC(3)` remain open and this session did not move them; it only
recorded why one attractive route is blocked.  The `nu=1` object's provenance
is unresolved and no priority claim attaches to it; the supplied arXiv
identifier belongs to an unrelated computational-geometry paper.
