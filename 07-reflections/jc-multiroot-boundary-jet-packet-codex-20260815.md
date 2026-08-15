# Boundary jets at a repeated root: the exact N versus N-1 packet

## Status

**HISTORICAL PROOF REFLECTION; PROMOTED TO THM-3436.**  The theorem file is
the current proof source.  This reflection preserves the derivation and exact
probe contract; its original candidate has now been independently audited.

## 1. Candidate statement

Let `K` have characteristic zero, let `g` be nonconstant, and

```text
P=ax+b+g(x)z^d,                    a!=0, d>=2.
```

After a finite normal splitting extension `K'/K`, write

```text
g=gamma product_j(x-alpha_j)^e_j,
beta_i=a alpha_i+b,
lambda=P-beta_i,
N=deg(rad(g)).
```

Base-change the target character sector before using the geometric boundary
coordinate:

```text
Cbar_(sigma-1)=C_(sigma-1) tensor_K K'.
```

Fix a **repeated** root `e_i>1`, target character `sigma-1` with
`1<=sigma<=d` (where `sigma=d` is wrap), and jet order `q>=1`.  Put

```text
R_q=K'[lambda]/lambda^q.
```

The candidate law is

```text
Cbar_(sigma-1)/(P-beta_i)^q Cbar_(sigma-1)
 ~= R_q^c,                                                (1)

c=N     if root i is THM-3433-selected in character sigma,
c=N-1   otherwise.                                       (2)
```

Selection means

```text
e_i>1,
d | sigma(e_i-1),
d | sigma e_j for every j!=i.                            (3)
```

Thus the first-window `N` versus `N-1` dimension repeats freely at every
nilpotent order; there are no hidden extension losses.  In particular the
`(P-beta_i)`-adic completion would be free of rank `N` in a selected
character and rank `N-1` otherwise.  This is a completion statement only.  It
does not yet assert an actual local direct-sum decomposition of the integral
module.

The connection contract is:

| item | exact content |
|---|---|
| source | the integral Hamiltonian sector modulo `(P-beta_i)^q` |
| target | a free packet over the principal Artin jet ring `R_q` |
| map | CRT into vertical/horizontal thickenings, then a formal horizontal gauge |
| preserved | character, jet order, repeated-root support, selected congruence, and module length |
| forgotten | the Prüfer arm itself, every `(P-beta_i)`-divisible channel, and nonsplit primary coordinates |
| required sidecar | the vertical finite-Neumann inverse and operator-level gauge, not fibre dimension alone |
| hostile | simple roots, where `lambda=yH` does not have comaximal factors; no claim is made there |

## 2. The exact thickened special fibre

Write `y=x-alpha_i` and

```text
g=y^e u,                       e=e_i>1,
lambda=a y+g z^d=yH,
H=a+y^(e-1)u z^d.                                      (4)
```

The factors `y` and `H` are comaximal because `H=a mod y`.  Therefore

```text
K'[x,z]/lambda^q
 ~=K'[x,z]/y^q direct-sum K'[x,z]/H^q.                  (5)
```

Both factors are stable under `D=Jac(P,-)`.  Indeed,

```text
D(lambda)=0,
D(y)=-d g z^(d-1) in (y^e),
D(H)=-H D(y)/y in (H).                                  (6)
```

Since `D(lambda^q A) subset lambda^q A`, quotienting the Hamiltonian
cokernel by `lambda^q` is exactly the cokernel of the induced derivation on
`A/lambda^q A`.  Thus `(5)` is a decomposition of the jet complex, not just
of its ambient ring.

## 3. The vertical factor is acyclic at every order

On `K'[y,z]/y^q`, split

```text
D=a partial_z+R.                                         (7)
```

Every term of `R` raises `y`-adic order by at least `e-1>=1`: the `g'`
coefficient has order `e-1`, while `g partial_y` has net order gain `e-1`.
Let `J` be coefficientwise integration in `z`, divided by `a`.  Then

```text
D J=1+R J,                                               (8)
```

and `RJ` is nilpotent modulo `y^q`.  The finite Neumann formula

```text
Q=J[1-RJ+(RJ)^2-... ]f                                   (9)
```

terminates and satisfies `D(Q)=f`.  It respects the character grading because
`J` raises the fibre exponent by one.  Hence every vertical target sector has
zero cokernel.  This remains true at `q=2`; no nilpotent vertical correction
creates a hidden class.

## 4. The horizontal operator is formally constant

Modulo `H^q`, the element `H` is nilpotent.  Equation `(4)` then makes `y`,
`u`, and `z` units.  Since `lambda=yH` and `y` is a unit, the ideals
`(H^q)` and `(lambda^q)` coincide on this factor.  Its weight decomposition
has coefficient ring

```text
B_q=R_q[x,rad(g)^(-1)].                                  (10)
```

For an input `p z^sigma`, direct substitution of
`g z^d=lambda-a y` gives the exact coefficient operator

```text
A_lambda(p)
 =a sigma p+(lambda-a y)[sigma(g'/g)p-dp'].              (11)
```

Put

```text
s=1-lambda/(a y),
r=s^(-sigma/d) mod lambda^q.                             (12)
```

The binomial series in `(12)` is finite modulo `lambda^q`, belongs to
`B_q^*`, and is defined in characteristic zero for every character, including
wrap `sigma=d`.  Direct differentiation gives the operator identity

```text
A_lambda(rp)=r s A_0(p) mod lambda^q.                    (13)
```

Both `r` on the domain and `rs` on the target are units.  Thus the entire
horizontal two-term complex is isomorphic to its special-fibre complex
tensored with `R_q`.  This is stronger than equality of fibre dimensions or
an Euler-characteristic calculation.  In particular, at `q=2`,

```text
r=1+[sigma/(ad)]lambda/y mod lambda^2,                   (14)
```

and `(13)` cancels the actual first extension term.

For wrap, `r=s^(-1)` and the same identity is literal.  Hence the wrap packet
has rank `N` at every repeated root and every jet order, agreeing with the
THM-3430 free rank without using its global direct-sum coordinates.

## 5. Special-fibre cohomology and freeness

At `lambda=0`, the horizontal Kummer component is

```text
z^d=-a/[y^(e_i-1) product_(j!=i)(x-alpha_j)^e_j].        (15)
```

The base is the affine line with `N` punctures.  Character `sigma` has
monodromy exponent vector

```text
sigma(e_i-1),              {sigma e_j:j!=i} mod d.       (16)
```

Its rank-one local system is trivial exactly under `(3)`.  The cochain model
of a wedge of `N` circles has one vertex and `N` edges.  Therefore

```text
(dim H^0,dim H^1)=
  (1,N)     in the selected/trivial case,
  (0,N-1)   otherwise.                                  (17)
```

Equation `(13)` now identifies the jet cokernel itself with

```text
H^1_special tensor_(K') R_q,                             (18)
```

which is free of the rank in `(2)`.  Freeness does not come from dimension
alone; it comes from the complex isomorphism.  Equivalently, over the
principal Artin local ring, the first fibre has `c` generators and total
length `cq`, forcing all `c` invariant factors to have maximal length `q`.

## 6. Wrap, nonsplit roots, and stopping boundaries

The formal gauge depends only on the chosen geometric root and is Galois
equivariant.  A nonsplit repeated-root orbit can therefore be checked after a
faithful splitting extension and then repackaged by its irreducible boundary
factor.  The exact probe includes two controls over `Q`:

```text
g=(x^2+1)^3:
  the two geometric roots are both unselected at (d,sigma)=(2,1);

g=x^3(x^2+1)^2:
  the rational multiplicity-three root is selected at (2,1),
  while the conjugate multiplicity-two roots are not.                   (19)
```

There is no separate collision case for the boundary parameters in this
family: if `alpha_i!=alpha_j`, then

```text
beta_i-beta_j=a(alpha_i-alpha_j)!=0
```

because `a!=0`.  Thus distinct geometric roots always give distinct
`P`-fibres, even after passing to a splitting field.  Any apparent collision
after descent is a Galois orbit of distinct parameters, not two local packets
supported at one geometric maximal ideal.

No claim is made for a simple root `e_i=1`.  Then `H mod y` is no longer the
constant `a`; the special-fibre components meet and the CRT/Neumann proof
above does not apply.  No claim is made that completion determines the full
local module, or that selected/unselected packets furnish canonical global
coordinates.  The jet quotient also kills every Prüfer arm because
multiplication by `lambda^q` is surjective on a Prüfer module.  Thus `(1)` is
compatible with THM-3433 but does not recover its torsion without the
DeathBar sidecar.

There is no LRC-to-JC map, polynomial mate, new Keller case, or open `JC(2)`
consequence.

## 7. Exact probe

The companion checks:

- `lambda=yH`, `D(lambda)=0`, and the exact `D(H)` identity on repeated-root
  profiles of multiplicities two and three;
- `D`-stability of both CRT ideals and `28` explicit finite-Neumann inverses
  through jet order four;
- `35` universal operator gauges through jet order five, including `10` wrap
  identities and `7` explicit `q=2` coefficients;
- `126576` special-fibre graph ranks for
  `2<=d<=10`, one to four roots, and multiplicities one through five;
- `506304` derived Artin-packet instances at jet orders `1,2,3,5`, split
  across `21686` selected and `104890` unselected geometric root-character
  profiles (bookkeeping evidence; the operator gauge proves freeness);
- the blocked `(4,2;3,1)` hostile, wrap controls, and the nonsplit profiles in
  `(19)`.

Reproduce with

```bash
python3 04-computation/jc_multiroot_boundary_jet_packet_probe_20260815.py
```

The finite grids audit the exact identities; they are not an extrapolated
cutoff.  THM-3436 records the completed independent audit of the CRT
derivation, vertical Neumann inverse, horizontal gauge, special local-system
rank, and faithful descent.
