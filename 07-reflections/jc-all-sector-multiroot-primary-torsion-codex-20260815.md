# All-sector primary torsion: local resonance needs global algebraization

## Status

This reflection records the mechanism behind provisional THM-3433.  The
theorem and its exact companion remain audit-required until a different agent
rederives the proof delta and the status is explicitly promoted.

## 1. Inheritance board

The closest proved mechanisms were:

- THM-3418: exact character-sector colimits for
  `P=ax+b+g(x)z^d`, with unquotiented nonwrap stages and quotient wrap stages;
- THM-3422: the complete one-root integral module, including the unique
  Prüfer break and its exact filtration;
- THM-3419: every generic character has rank `N=deg(rad(g))`, regardless of
  root multiplicity;
- THM-3427: the generic wrap defect retains a logarithmic fingerprint of the
  multiplicity vector without recovering the integral module.

The corrected near miss was “one locally resonant root gives one arm.”  It
fails because the formal root gauge may have fractional orders at the other
roots.  The smallest useful hostile is

```text
(d,sigma;e_1,e_2)=(4,2;3,1).
```

Root `1` is locally resonant since `4|2(3-1)`, but the simple second root has
nontrivial character monodromy since `4` does not divide `2*1`.  The double
local hostile `(2,1;3,3)` is even sharper: both roots are locally resonant,
yet neither arm algebraizes because each root blocks the other.

The least-used sidecar was transition injectivity.  It is what converts a
formal horizontal solution into an exact integral endpoint and what makes the
embedded root filtration strict.

The live concept board was:

```text
transition kernel | global power | root gauge | horizontal ODE |
primary endpoint | filtration strictness | Galois orbit
```

## 2. The transition kernel is a regime switch

In nonwrap character `sigma`, depth `k` uses `n=sigma+kd` and

```text
L_n(q)=d g q'/(an)-g'q/a.
```

The kernel equation has only the root-divisor solution

```text
q=C product_j (x-alpha_j)^(n e_j/d)
```

when the displayed orders are integral.  This notation suppresses the
irrelevant leading scalar and does not adjoin a root of it.  It is polynomial
exactly when

```text
d | sigma e_j                    for every root j.
```

This condition is independent of `k`.  Therefore there are exactly two
regimes:

```text
global power: every transition has a one-dimensional kernel;
otherwise:    every transition is injective.
```

The global-power regime initially looked dangerous for torsion because zero
transition directions are everywhere.  The opposite is true.  If

```text
H_k=product_j (x-alpha_j)^(k e_j+sigma e_j/d),
```

then

```text
L_n(H_k p)=nonzero_scalar * H_(k+1) p'.
```

Every polynomial dies after enough derivatives, so `H_k K[x]` has zero
colimit.  Quotienting it makes CRT available.  Each root channel becomes its
one-root diagram modulo the same dying subdiagram.  Both the transition and
the two-term Hamiltonian action intertwine with that root channel.  Since a
nonwrap character cannot satisfy both `d|sigma e_i` and
`d|sigma(e_i-1)`, every channel is a Laurent line.  Thus global transition
kernels produce a torsion-free `K[P]`-module sum rather than torsion.

This is a reusable reversal: a persistent family of finite-stage kernels can
be colimit-zero when the induced map on the kernel ideal is differentiation.

## 3. The horizontal equation sees every other root

In the injective regime, put `xi=(rho-b)/a`.  Moving `(P-rho)[p]_k` to stage
`k+1` gives the exact polynomial operator

```text
M_(n,xi)(p)
 =d g(x-xi)p'+(n+d)gp-n(x-xi)g'p.
```

Its unique formal logarithmic horizontal solution is

```text
p=product_j (x-alpha_j)^(n e_j/d)(x-xi)^(-(n+d)/d).
```

It is rational only when the combined root orders are integral.  If `xi` is
not a root, its residue there is nonintegral.  If
`xi=alpha_i`, its root orders are

```text
j!=i:  k e_j+sigma e_j/d,
i:     k(e_i-1)+sigma(e_i-1)/d-1.
```

The solution is polynomial exactly when

```text
e_i>1,
d | sigma(e_i-1),
d | sigma e_j for every j!=i.
```

This is the selected-root law.  The first congruence is the local one-root
break; all remaining congruences are the global algebraization debt.  The
`e_i>1` inequality cannot be dropped: at a simple root the formal exponent is
`-1` even though `d|sigma(e_i-1)` is automatic.

Two nonwrap roots cannot both be selected.  Selection at `i` forces
`d|sigma e_j`, while selection at `j` forces `d|sigma(e_j-1)`; subtracting
gives `d|sigma`, impossible for `sigma<d`.  Thus collision is prevented by
character arithmetic before the distinct boundary values are even used.

## 4. Why the one-root arm is exact, not merely visible

For a selected root, the other-root factor

```text
h_i u_i^k
 =product_(j!=i)(x-alpha_j)^(sigma e_j/d+k e_j)
```

is polynomial and satisfies the exact gauge identity

```text
L_n^g(h_i u_i^k p)
 =h_i u_i^(k+1)L_n^(gamma y_i^e_i)(p).
```

It embeds the entire THM-3422 one-root sector, hence at least one Prüfer arm.
The delicate question was whether a deeper element of that arm could be
represented at an earlier ambient stage outside the gauge.  That would make
the claimed finite thickness too small.

The valuation test closes the gap.  At another root `alpha_j`, write

```text
v=ord_(alpha_j)(p),
r=k e_j+sigma e_j/d.
```

If `v<r`, the leading coefficient of `L_n(p)` is proportional to
`dv-ne_j`, which vanishes only at `v=r`.  Its order is `e_j+v-1`, strictly
below the next gauge boundary `r+e_j`.  Hence a class can enter the gauge at
stage `k+1` only if it was already in the gauge at stage `k`.  The quotient
transitions are injective, so the embedding is filtration-strict.

The exact stage thickness is therefore inherited without distortion:

```text
k(e_i-1)+sigma(e_i-1)/d.
```

To exclude extra primary torsion, use the unique endpoint.  If
`lambda^m v=0`, then `lambda^(m-1)v` lies in the one-dimensional endpoint of
the visible Prüfer arm.  That arm is lambda-divisible, so subtract a matching
preimage and induct on `m`.  Every primary torsion element lies in the arm.
This exhaustion is performed after passing from the root-splitting field to
an algebraic closure, so every annihilator factors linearly.  Each element
and annihilator is defined over a finite subextension; the polynomial norm
and faithful flatness then descend the equality.  Merely calling the finite
root-splitting field a “split PID” would not rule out unrelated nonlinear
annihilator factors.

## 5. Character arithmetic and descent

For root `i`, selection says that `sigma mod d` annihilates

```text
(e_i-1,{e_j:j!=i}).
```

The number of selected characters is exactly

```text
gcd(d,e_i-1,{e_j:j!=i}).
```

Every selected character at that root has the same slope `e_i-1` but a
different intercept.  This answers the original regular-packet puzzle:
generic rank is uniformly `N` in every character, while integral torsion is
character-shifted and controlled by the whole multiplicity vector.

THM-3431 supplies a lawful but deliberately lossy comparison language.  The
first visible endpoint at root `i` has DeathBar

```text
[0,q_i),                    q_i=sigma(e_i-1)/d,
```

and depth `k` extends its right endpoint to `q_i+k(e_i-1)`.  Wrap uses
denominator level `q=k+1`, hence length `q(e_i-1)`.  The arms therefore give
a multiset of valuation bars.  This is related-only: THM-3431 proves that all
additive maps between the LRC and JC secondary coefficient objects vanish in
both directions, even when bar lengths agree.  The barcode forgets the site,
class, coefficients, and target predicate; it is not an LRC-to-JC map.

Nonwrap descent is unusually rigid.  The selected set is Galois-stable and
has at most one point, so a selected root must be rational.  A nonsplit orbit
cannot carry a hidden nonwrap packet.  For example, over `Q`,
`(x^2+1)^3` has two locally resonant geometric roots for `(d,sigma)=(2,1)`
but no arm; `x^3(x^2+1)^2` has exactly the rational-root arm.

Wrap is different.  Every repeated root is selected, so conjugate arms
package by irreducible factors and the descended torsion is the localization
quotient for `c=g/rad(g)`.  The wrap proof is repeated directly from the
THM-3418 CRT diagram and THM-3422; it does not use provisional THM-3430 as a
dependency.

## 6. What changed on the board

- Transition kernels became a dying quotient, not a torsion source.
- Local resonance became a necessary coordinate inside a global
  algebraization criterion.
- The horizontal ODE supplied both the endpoint classification and the
  uniqueness needed for primary exhaustion.
- Valuation strictness upgraded an embedded arm to an exact finite
  filtration theorem.
- Galois orbit size became a nonwrap torsion obstruction.
- Generic rank and torsion thickness are now cleanly separated: `N` measures
  punctures; the selected congruence vector measures integral boundary debt.

The theorem intentionally stops at primary torsion.  Outside the global-power
regime it does not identify the torsion-free quotient or split the complete
module, and it does not manufacture a coefficient correspondence to another
problem.  Torsion is still killed by generic localization, exactly as
MISTAKE-374 warns.  No polynomial mate or open `JC(2)` consequence follows.
