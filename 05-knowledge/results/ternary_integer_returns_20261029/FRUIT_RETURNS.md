# A sharp thirteen-shift cover and a seventy-six-step coefficient budget

**PROVED with an exact algebraic boundary certificate:** a source-preserving
return controller on the real fruit curve and its three-isogeny quotient.
Every relevant real phase is covered; rational lifting remains a separate
arithmetic guard. No Collatz implication or literature-priority claim.

## 1. Recover the original object and the stronger incoming quotient

The inherited [fruit analysis](../prime_shells_20260921_elliptic.md) gives

\[
 E:y^2=x^3+109x^2+224x,\quad G=(-4,28),\quad T=(56,728),
\]

and the projective fruit coordinates
`[a:b:c]=[56-x+y:56-x-y:-56-12x]`. The positive locus `S` lies on the bounded
real component `C` and is exactly

\[
 x<-14/3,\qquad x^2+112x+784>0.                         \tag{11}
\]

The September proof already excludes global preservation of `S` by fixed
expanding affine elliptic maps. The replacement here is **translation followed
by an authenticated first return**, rather than another purported invariant
tripling tree. Large coordinates are not evidence of large rank: the inherited
79/80/81-digit solution is `9G` on this rank-one curve.

The incoming [Session97 quotient](../session97_overnight_arithmetic_carriers_20261029/fractions/README.md)
supplies the exact map

\[
 \phi(x,y)=\left(\frac{x(x-56)^2}{(x-4)^2},
 \frac{(x-56)(x^2+44x+224)y}{(x-4)^3}\right)
\]

to `E':Y^2=Z^3+421Z^2+43904Z`, with kernel `O,(4,±52)`. Put
`H=(-196,196)` on `E'`; its exact identities are `phi(G)=3H` and

\[
 nH\in\phi(E(\mathbb Q))\iff3\mid n.                   \tag{12}
\]

The quotient positive-real chamber `S'` on the bounded component `C'` is
`Z<-686/3`. On `C`, `phi^(-1)(S')=S`. One way to see the last assertion is
to normalize `a+b+c=1`: the quotient records `q=ab+ac+bc`, and
`abc=(6q-1)/7`. Since the three original coordinates are real, `q>1/6`
forces them all positive. The quotient itself does not guarantee that a
rational point reconstructs a rational triple.

## 2. The exact universal cover

Let `U(P)=P+2G` on `C`. On `C'` let `U'(Q)=Q+6H`. These are conjugate through
the isogeny in the sense `phi U=U' phi`, with three marked lifts retained.

**Return theorem.**

1. For every real `P in C`, at least one of `U(P),...,U^13(P)` belongs to `S`.
   The analogous statement holds for every real `Q in C'` and `S'`.
2. If the source already belongs to `S` (respectively `S'`), the first positive
   return time is always one of **4, 9, 13**.
3. Thirteen is sharp: `17G` first returns at `43G`, after thirteen translations.
   Other exact witnesses are `9G ->17G` at time four, and
   `21G+T ->39G+T` at time nine.
4. If time zero is allowed, every source enters by time twelve: apply the
   first assertion to `U^(-1)(P)`. The permitted times are then `0,...,12`.

This covers all real sources, not just the tested rational multiples. The
finite certificate establishing it is as follows.

The six boundary points of `S` are obtained by mapping
`[0:2+sqrt(3):1]`, `[0:2-sqrt(3):1]` and their cyclic coordinate permutations.
They lie in `E(Q(sqrt(3)))`. The two quotient boundaries are

\[
 (Z,Y)=(-686/3,\ \pm686\sqrt3/9).                      \tag{13}
\]

Translate each boundary backwards by the first thirteen shifts, including
time zero for the first-return calculation. There are 84 distinct endpoints
on `C`, or 28 on `C'`. Sort them in circular order by increasing `x` on the
upper half and decreasing `x` on the lower half. Translation preserves that
orientation, since it preserves `dx/(2y)`. The leftmost point of each bounded
component is inside its target region; its original boundaries therefore
alternate exit/entry. At the rational anchor `G` (respectively `phi(G)`),
evaluate all fourteen membership bits directly, then sweep the boundary list.

Every comparison uses the exact sign of `a+b sqrt(3)`, decided by rational
signs and `a^2-3b^2`; no numerical elliptic logarithm is used. On each open
cell the membership bits are constant. At every endpoint the checker separately
evaluates all fourteen translated points by the group law, verifying the
event sweep including the strict/open boundary convention. Thus all points
of the circle are covered. The minimum covering multiplicity is zero for
every horizon 1 through 12 and one at horizon 13. Restricting the same
arrangement to sources in the target gives exactly the first-return alphabet
`{4,9,13}`. The six-endpoint and two-endpoint calculations agree independently.

## 3. An effective strengthening of the nonlifting obstruction

Session97 already proves that infinitely many `35H+6kH` lie in the positive-real
chamber while none lifts rationally; its proof used circle density. The present
cover gives a **uniform gap at most thirteen**. More generally, for any odd
`n` not divisible by three, every thirteen consecutive forward translates
`nH+6kH` contain a positive-real point, and all remain nonlifting by (12).
Such hits have lower density at least `1/13` in the translation clock.

This is a stronger obstruction to dropping the lifting coordinate: even
uniform real entry, rational orientation, and bounded positive-return gaps
leave the arithmetic error intact. It is not a Sha example; the inherited
local obstruction already detects a failure of rational lifting.

## 4. Pay for a change of arithmetic phase: the exact budget 76

There is nevertheless a complete steering controller on the **named cyclic
subgroup**, provided changing the target by a declared translation is allowed.
Start at `nH` for any odd integer `n`, retaining its integer label.

1. Choose the unique `e in {0,2,4}` with `n+e=0 mod3`. Even correction preserves
   the real component. Set `m=(n+e)/3`, which is odd, and choose the named lift
   `P=mG` of the adjusted point `(n+e)H`.
2. Choose the first `j in {0,...,12}` for which `P+2jG in S`. It exists by the
   cover. Its image is `(n+e+6j)H`, and its fruit coordinates are positive
   rationals, hence scale to positive integers.

The total displacement in the original subgroup coefficient is

\[
 \boxed{0\le n_{\rm out}-n=e+6j\le4+6\cdot12=76.}       \tag{14}
\]

This bound is sharp for the specified rule: `n=53` selects `e=4`, lifts to
`19G`, waits twelve shifts, and outputs `43G`, whose quotient is `129H`.
The source coefficient has moved from 53 to 129, a difference of 76.

The representation `(n,e,j,chosen kernel lift)` is reversible: the output
coefficient recovers `n` by subtracting `e+6j`. Coordinate labels and torsion
owners remain attached. The rule does **not** assert that the original
nonlifting `nH` acquired a lift: it changes the point and charges that change.
Nor are these translations native Collatz operations. The number 76 here has
a specific derivation from a phase repair and a certified cover, rather than
being evidence for an identification with an unrelated 76-state carrier.

## 5. What actually transfers to the controller research

Source: a marked arithmetic state. Target: a state satisfying both a real
positivity test and a rational reconstruction test. Preserved: exact group
action, source label, real component, kernel owner, integer clock. Lost under
the quotient alone: the rational lift and cyclic coordinate owner. Required
sidecar: the lift plus charged phase correction and return time. Cheapest
decisive test: all algebraic boundary cells, rather than a long orbit sample.

The reusable move is to separate **phase correction**, **finite-cover return**,
and **source recovery**. Here all three can be certified. The analogous
Collatz gap is that its permitted native moves have not been shown to realize
a uniformly paid phase correction; one may not add an arbitrary translation
to the integer simply because it works on an elliptic group.

Reproduce with `python3 -B check_fruit_returns.py` and with `-O`, from this
directory. [fruit_checks.json](fruit_checks.json) pins both algebraic sweeps,
their endpoint digests, the three first-return witnesses and the sharp
coefficient-displacement witness. These are exact finite proof certificates,
not estimates of a recurrence from samples.
