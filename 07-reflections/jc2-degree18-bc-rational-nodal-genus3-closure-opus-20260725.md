---
source: codex-2026-07-25-degree18-bc-rational-closure
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Within the genuine nonsplit polynomial exact-square-prefix degree-eighteen
  branch of THM-2262/2297, the two rational roots of THM-2311's B--C
  ratio bank each give a geometrically irreducible trigonal curve with
  one ordinary node, ten simple finite branch values, three unramified
  points at infinity, and normalization genus three.  The genus/deck
  argument empties these two ratios.  Seven algebraic B--C ratios, the
  B--W/C--D/D--W banks, higher-support strata, split/even descent, JC(2),
  and DC(2) remain open.
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - jc2-degree18-bd-4075-smooth-cubic-closure-opus-20260725
  - jc2-degree18-bd-quartic-nodal-atlas-opus-20260725
script: 04-computation/jc2_degree18_bc_rational_branch_scout.py
output: 05-knowledge/results/jc2_degree18_bc_rational_branch_scout.out
script_sha256: 38413e6e16a2096ffb981da58d92606b1b18a4bf4cfc5fc6c96c29ca85f3a53d
output_sha256: d49a12cbbb563101510273a14878f96e8ae92c4ca23f3a0359ce9b18a4a79a43
hash_basis: LF-normalized working-tree bytes
---

# The two rational B--C ratios normalize to genus three

## 1. Scope and the rationalizing coordinate

On the `B`--`C` plane, THM-2311 defines

```text
t=C^2/B^3
```

and leaves nine nonzero ratios.  Two are rational:

```text
t=-2000/15309,                  t=-125/1134.        (1)
```

Because `B!=0`, make a constant weighted scaling over the algebraic
closure so that `B=1`.  Then `C^2=t`, and the invertible constant
coordinate

```text
x=Cy                                                (2)
```

removes the square root from the spectral equation.  Multiplying by
`t^3` gives the rational trigonal equation

```text
H_t(u,x)
 =-26040609t^3u^3
  +(49601160t^3+1607445t^2x^2)u^2
  +(-20995200t^3-2857680t^2x^2-138915tx^4)u
  -5598720t^3x+777600t^2x^2-435456t^2x^3
  +78120tx^4+1127x^6.                              (3)
```

The coefficient of `u^3` is a nonzero constant.  Thus (3) defines a
finite degree-three map to the `x`-line, and a factorization over
`Q(x)` would specialize without losing `u`-degree at every finite
`x`.

## 2. Irreducibility and geometric irreducibility

At `x=1`, exact reduction gives the following cubics, with coefficients
listed from `u^3` to `1`:

```text
t=-2000/15309:    (20,28,22,28) mod 37;

t=-125/1134:      (11,12,12,8)  mod 13.             (4)
```

Neither cubic has a root in its displayed finite field, so each is
irreducible.  Gauss's lemma and the constant nonzero leading
`u`-coefficient imply that (3) is irreducible in `Q(x)[u]`, hence in
`Q[u,x]`.

This is also geometric irreducibility, not merely rational
irreducibility.  The point

```text
(u,x)=(0,0)
```

is rational and smooth because

```text
H_x(0,0)=-5598720t^3!=0.                            (5)
```

If an irreducible rational curve split into several conjugate geometric
components, a rational point would lie on every conjugate component.
Their intersection would be singular.  Equation (5) excludes that
possibility.

## 3. Exact discriminant layers

Let

```text
Delta_t(x)=Disc_u H_t(u,x).                         (6)
```

The standard-library exact Euclidean algorithm gives, for both values
in (1),

```text
Delta_t=A_t(x)R_t(x)^2,

deg A_t=10,              deg R_t=1,

gcd(A_t,A_t')=gcd(A_t,R_t)=1.                       (7)
```

The primitive linear factors and their points are

```text
t=-2000/15309:
  R_t=200+189x,       (x_0,u_0)=(-200/189,20/189);

t=-125/1134:
  R_t=25+21x,         (x_0,u_0)=(-25/21,10/63).     (8)
```

At each point in (8),

```text
H=H_u=H_x=0,

H_uu H_xx-H_ux^2!=0,

H_uu!=0.                                            (9)
```

Thus it is one ordinary node.  The last inequality says neither
normalization branch is vertical over the `x`-line, so both are
unramified there.

Every root of the squarefree degree-ten factor `A_t` is a simple
discriminant value.  Since the leading `u`-coefficient is constant,
each gives one smooth simple ramification point of index two.  These
are all remaining finite branch values.

For reproducibility, the primitive ascending coefficient strings of
the two degree-ten factors have SHA-256 fingerprints

```text
6e173b70dd11450fbf17b6e3522362a5da3ed8a5e75e388386c62cd642a36327,

5165d043f0441a5f9b831706c96fc138ae88ad5df954856d84710b3134baee79.
                                                            (10)
```

## 4. Infinity and the exact genus

Put `v=u/x^2` and `r=1/x`.  At `r=0`, the leading equation is

```text
-26040609t^3v^3
 +1607445t^2v^2
 -138915tv
 +1127=0.                                          (11)
```

Its discriminant is nonzero at both ratios.  Hence there are three
distinct smooth points at infinity.  The implicit-function theorem
uses `r` as a local parameter on each branch, so the degree-three map
to the `x`-line is unramified there.

The normalization therefore has exactly ten units of ramification:
one from each root of `A_t`, none from the node branches, and none at
infinity.  Riemann--Hurwitz gives

```text
2g-2=3(-2)+10=4,

g=3.                                                (12)
```

## 5. Keller contradiction

A rational Keller trajectory gives a rational map from `P^1` to the
smooth projective normalization of (3).  Genus three forces that map to
be constant, so `u` and `x` are constant.

Because `C!=0` is a constant after weighted normalization, (2) makes
`y` constant.  If `y!=0`, THM-2262's first-flux square makes `T^2`,
then `T`, then the nonzero deck coordinate `q` constant, contradicting
the genuine nonsplit deck.  If `y=0`, THM-2262 Section 4 closes the
exceptional center with the retained third flux, Keller one-form,
rational-primitive lemma, and whole-polynomial Faber sidecar.

Returning through the constant weighted scaling preserves constancy and
the original first-flux equation.  Therefore no survivor in the scoped
branch has either ratio in (1).

## 6. Information ledger

```text
source:
  the two rational linear factors of THM-2311's B--C ratio bank;

map:
  weighted normalization B=1, rationalizing coordinate x=Cy, then
  the degree-three projection to the x-line;

preserved:
  the Keller trajectory, genuine deck, finite branch multiplicity, and
  the y=0 exceptional alternative;

destroyed:
  the chosen square root C of t, restored as a nonzero constant;

hostile boundary:
  each repeated discriminant value is a genuine node, not a ramification
  point; the ten simple values still force genus three;

next target:
  apply the same rationalizing-coordinate and simple-branch lower-bound
  test to the seven algebraic B--C ratios and the rational C--D/D--W
  points.                                             (13)
```

Run

```text
python3 04-computation/jc2_degree18_bc_rational_branch_scout.py
python3 -O 04-computation/jc2_degree18_bc_rational_branch_scout.py
```

Both executions must match

```text
05-knowledge/results/jc2_degree18_bc_rational_branch_scout.out
```

after LF normalization.
