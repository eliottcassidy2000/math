---
id: THM-3564
title: "Target-graph trisection degree-resonance irreducibility gate"
status: >
  PROVED + VERIFIED-EXACT.  Let phi be a nonzero polynomial in C[a,b] and
  d=deg_a(phi).  For the fixed THM-1300 Keller map, the pullback of the
  target coordinate graph c+phi(a,b)=0 is irreducible whenever d is not
  congruent to 1 modulo 3.  At a=infinity the three terms of the core cubic
  have a unique least valuation except possibly at the single resonant value
  k=(d+2)/3; this value is integral exactly for d congruent to 1 modulo 3.
  Resonance is necessary for reducibility, not sufficient.
source: kps-s188
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
related:
  - THM-1335-trisection-modulus-master-identities-trace-polynomiality
  - THM-3546-invariant-graph-keller-descent-criterion
  - THM-3560-jelonek-euler-gate-monomial-target-shear-no-go
companion: 04-computation/jacobian_target_graph_degree_resonance_irreducibility_kps_s188.py
output: 05-knowledge/results/jacobian_target_graph_degree_resonance_irreducibility_kps_s188.out
script_sha256: 6d5895a87268763412ebd87cba63003060ba89f28dfbbac47469dbbac5f3408c
output_sha256: a79f387fd26ccf756193296bb31388c3b1300e6575968546148b5955383a34f9
hash_basis: LF-normalized bytes
---

# THM-3564 -- target-graph trisection degree-resonance gate

**PROVED + VERIFIED-EXACT.**  In the complete family of triangular target
coordinates

```text
g_phi(a,b,c)=c+phi(a,b),                 0!=phi in C[a,b],              (1)
```

two thirds of the `a`-degree classes have irreducible pullback under the
fixed THM-1300 Keller map.  Reducibility can occur only in the sparse degree
ladder

```text
deg_a(phi)=1,4,7,10,... .                                             (2)
```

The theorem does not assert that a member of `(2)` is reducible.  It says
that `(2)` is the exact boundary at which the valuation obstruction stops.

All varieties and function fields below are over `C`.

## 1. Core cubic on a target graph

Write target coordinates as `(a,b,c)`.  THM-2473 gives the generic fibre's
depressed core cubic

```text
E(x)=Lx^3+(4-3bc)x-2c,

L=27a^2c^2-18abc+16a+b^3c-b^2.                         (3)
```

On the graph `T_phi=V(c+phi)` this becomes

```text
E_phi(x)=L_phi x^3+(4+3b phi)x+2phi,                    (4)

L_phi=27a^2 phi^2+18ab phi+16a-b^3 phi-b^2.            (5)
```

Put `K=C(b)` and regard `(4)` as a cubic over `K(a)`.  If
`d=deg_a(phi)`, then nonvanishing of the leading coefficient of `phi` gives

```text
                  deg_a L_phi   deg_a(4+3bphi)   deg_a(2phi)

d=0                     2                0                0,
d>=1                  2+2d               d                d.           (6)
```

There is no possible leading cancellation in `(5)`: its top term is
`27a^2 phi^2`.  In the `d=0` row, `phi!=0` is essential.  The excluded
polynomial `phi=0` gives the familiar reducible pullback
`F3=x(2-3xy-x^2z)` and is a sharp hostile to any unqualified statement.

## 2. The valuation trichotomy

Let `v` be the valuation of `K(a)` at `a=infinity`, normalized by
`v(a)=-1`.  Suppose that `r in K(a)` were a root of `(4)` and put
`k=v(r) in Z`.  A vanishing sum in a valued field must attain its least
valuation at least twice.

For `d=0`, the valuations of the three displayed terms in `(4)` are

```text
-2+3k,                    k,                    0.                      (7)
```

If `k<=0`, the first entry is strictly least; if `k>=1`, the last entry is
strictly least.  Thus no root exists.

Now let `d>=1`.  The three valuations are

```text
A=-(2+2d)+3k,             B=-d+k,               C=-d.                  (8)
```

If `k<=0`, then

```text
A-B=-(d+2)+2k<0,          A-C=-(d+2)+3k<0,       (9)
```

so `A` is uniquely least.  If `k>0`, then `B-C=k>0`; hence `B` cannot be
least.  The only remaining way to avoid a unique minimum is

```text
A=C    iff    3k=d+2.                                      (10)
```

The integer `k` in `(10)` exists exactly when

```text
d congruent to 1 modulo 3.                                 (11)
```

Consequently, when `d` is not `1 mod 3`, `(4)` has no root in `K(a)`.
A cubic over a field is reducible if and only if it has a root, so `E_phi`
is irreducible over `C(a,b)`.

## 3. From the generic cubic to the pullback surface

The graph `T_phi` is a smooth irreducible copy of `A2`.  Its pullback is

```text
X_phi=V(F3+phi(F1,F2)) subset A3.                         (12)
```

The map `F` is etale and quasi-finite.  Hence its base change
`X_phi -> T_phi` is etale and quasi-finite, so `X_phi` is smooth and reduced.
Every irreducible component of `X_phi` has dimension two, and quasi-finiteness
forces its image closure to have dimension two; therefore every component
dominates `T_phi`.  At the generic point of `T_phi`, THM-2473's eliminant
and recovery formulas identify the generic fibre algebra with the algebra
defined by the core cubic `(4)`.  Its irreducibility leaves exactly one
generic point.  Thus `X_phi` has exactly one irreducible component.

We have proved

```text
0!=phi in C[a,b],  deg_a(phi) not congruent to 1 mod 3
        ==> F3+phi(F1,F2) is irreducible in C[x,y,z].       (13)
```

In particular, on these two degree classes the complete pullback is the only
possible source hypersurface factor.  THM-3560's Euler gate can therefore
exclude coordinate factors merely by showing the complete pullback is not
`A2`.  The first genuinely unresolved factor-bearing target-graph cell is
the resonant mixed-quadratic row `deg_a(phi)=1`.

## 4. Scope and exact verification

Resonance in `(10)` is only a necessary condition for a rational root.  It
does not construct one, prove reducibility, or establish a planar Keller
descent.  Degrees `1 mod 3` require their own residual-coefficient equation.
Likewise, `(13)` concerns triangular target graphs; a general nonlinear
target coordinate need not have this form.

Reproduce the exact controls with

```bash
python3 04-computation/jacobian_target_graph_degree_resonance_irreducibility_kps_s188.py
python3 -O 04-computation/jacobian_target_graph_degree_resonance_irreducibility_kps_s188.py
```

The ordinary and optimized transcripts agree.  The companion independently
checks the coefficient degrees in `(6)` through `d=12`, checks the valuation
minimum on `0<=d<=30` and `-40<=k<=40`, and obtains exactly one exceptional
`k=(d+2)/3` for each resonant row.  These finite tables are hostile controls;
the all-degree conclusion is the two-case inequality `(9)--(11)`.

**QED.**
