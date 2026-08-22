# Nonlinear wrap linearization: the degree disappears at the boundary

## Status

This reflection records the mechanism behind provisional THM-3430.  The
theorem remains audit-required until a different agent rederives its delta
from the pushed immutable checkpoint.

## 1. Inheritance board

The closest proved mechanisms were:

- THM-3418: every nonlinear monomial fiber has exact character-sector
  colimits, with a special quotient tower in the wrap sector;
- THM-3412: in its gradient-unimodular range, the linear-`z`
  principal-part tower has a complete free/torsion response and one Prüfer
  arm per repeated root;
- THM-3422: a one-root nonlinear sector is a weighted bilateral chain, and
  its wrap character always selects the repeated-root arm;
- THM-3427: generic ranks forget multiplicity, while the canonical wrap
  defect recovers its logarithmic divisor fingerprint.

The corrected near miss was the tempting rootwise decomposition of every
character.  THM-3422's `(d,sigma;e_1,e_2)=(4,2;3,1)` first-window hostile
shows that nonwrap periods couple the roots.  The least-used sidecar was the
actual wrap stage `K[x]/g^(k+1)`, rather than the generic Kummer curve.

The live concept board was:

```text
wrap quotient | d=1 principal parts | CRT root arms |
character monodromy | generic defect | finite torsion thickness
```

## 2. The cancellation that changes the problem

At wrap depth `k`, the incoming exponent is `n=(k+1)d`.  THM-3418's
transition is

```text
L_n=(d/(an))g partial_x-(1/a)g',
```

so

```text
L_((k+1)d)=[1/(a(k+1))]g partial_x-(1/a)g'.
```

The exponent `d` has disappeared.  The stage quotient, transition, and
Hamiltonian action are now byte-for-byte the `d=1` diagram.  This is much
stronger than a rank coincidence:

```text
wrap response for ax+b+g z^d
  = full response for ax+b+g w
```

as filtered modules after identifying the two abstract Hamiltonian
generators.  The source/target/map contract is unusually tight:

| field | content |
|---|---|
| source | nonlinear wrap colimit |
| target | linear principal-part response |
| map | identity on every `K[x]/g^(k+1)` stage |
| preserved | transitions, `P` action, free rank, primary support, filtration |
| lost | every nonwrap character |
| needed sidecar | the wrap quotient boundary |

The operation is therefore “restrict to the boundary character, then forget
the covering degree.”  Reversing the order is invalid.

## 3. Exact integral consequence

Put

```text
S=rad(g),       c=g/S,       N=deg(S),       E=deg(c).
```

The tempting direct appeal to THM-3412 does **not** work for arbitrary `g`.
For `P_1=ax+b+gw`, the gradient pair is `(a+g'w,g)`, and `a+g'w` is a unit
modulo `g` exactly when `g'` is nilpotent modulo `g`, equivalently when all
geometric roots of `g` are repeated.  The examples `g=x` and
`g=x^2(x-1)` expose the false shortcut.  They do not refute the wrap theorem;
they force a different proof.

The exact CRT diagram and the proved one-root theorem THM-3422 instead give,
after faithful splitting,

```text
W_d tensor K' ~= direct_sum_i
  (K'[P] direct-sum (Pr_(P-beta_i) if e_i>1)).
```

Descent is visible rather than implicit.  The global stage-zero classes

```text
[c x^j z^(d-1)],                 0<=j<N,
```

restrict to nonzero scalar multiples of
`alpha_i^j y_i^(e_i-1)` at the `i`th root.  Their coefficient matrix is a
diagonal matrix times a Vandermonde matrix, so they descend and split the
free rank `N`.  The stage-`k` torsion window splits as

```text
direct_sum_(i:e_i>1)
  K'[P]/(P-beta_i)^((k+1)(e_i-1)),
```

which descends by the PID elementary-divisor theorem to `K[x]/c^(k+1)`.
Choosing its cyclic generators compatibly and taking the colimit gives
`K[x,c^(-1)]/K[x]`.  Therefore

```text
W_d ~= K[P]^N direct-sum K[x,c^(-1)]/K[x].
```

The abstract split coordinates are noncanonical, but the torsion submodule
and its filtration are canonical.  Over a splitting field, a root of
multiplicity `e_i` contributes

```text
K[P,(P-beta_i)^(-1)]/K[P]
```

exactly when `e_i>1`.  At stage `k`, its finite block has length
`(k+1)(e_i-1)`.  Summing gives exact thickness `(k+1)E`.

This upgrades THM-3427's carrier bridge.  There, the transform

```text
c -> d rad(g) dlog(c)
```

compressed multiplicity into one generic defect.  The wrap colimit retains
the whole localization quotient and therefore the complete Prüfer arms.
Generic response sees the fingerprint; integral wrap response sees the
object that cast it.

## 4. Why CRT works here and nowhere automatically

After splitting `g=gamma product y_i^e_i`, define
`u_i=product_(j!=i)y_j^e_j`.  At stage `k`,

```text
K[x]/g^(k+1)
 =direct_sum_i u_i^(k+1)K[x]/g^(k+1).
```

Moreover

```text
L^g_(d(k+1))(u_i^(k+1)p)
 =u_i^(k+2)L^(gamma y_i^e_i)_(d(k+1))(p).
```

Thus CRT splits the entire directed diagram into the one-root diagrams of
THM-3422.  It is an independent explanation of the same theorem, and it
keeps Galois descent honest.

Nonwrap stages are unquotiented copies of `K[x]`; the CRT sum does not exist.
There is only a limited survivor.  If `g=gamma y^e u` and
`u^(sigma/d)` is polynomial, then

```text
u^(sigma/d+k)K[y] subset K[x]
```

is an injected one-root subdiagram.  When `e>1` and `d|sigma(e-1)`, its
endpoint is

```text
tau=[u^(sigma/d)y^(sigma(e-1)/d-1)z^(sigma-1)].
```

The positive control `(d,sigma;e_1,e_2)=(2,1;3,2)` carries such an embedded
arm.  The hostile `(4,2;3,1)` does not algebraize the other-root gauge.  This
explains both the existence of some multiroot nonwrap torsion and the failure
of a naive global root sum.

## 5. What changed on the concept board

- The wrap quotient and `d=1` principal parts became the same object.
- CRT is not merely a special-fiber dimension count; at wrap it is a
  diagram-level decomposition.
- THM-3427's defect is now understood as the logarithmic shadow of the exact
  carrier torsion.
- Character monodromy remains the gate for nonwrap algebraization.
- Finite thickness is no longer inferred from generic rank: it is exactly
  `(k+1)(deg(g)-deg(rad(g)))`.

The next safe question is whether the algebraized nonwrap root-gauge
subdiagram exhausts primary torsion.  Nothing here proves that.  A proof
would need a global exclusion for torsion outside those subdiagrams, not
another local formal gauge.  No Keller or open `JC(2)` consequence follows.
