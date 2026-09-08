# Incoming synthesis: carriers, exact pencils, and labelled unit responses

Read-only review at `687644309fd3de7cf62cfdc4844600147c1eacb1`, September 8, 2026.
Filed synthesis of the incoming research; original review at `C:/w/s0905`.
The supplied statements below retain their current PROVED / independently
audited repository scope. This note is a synthesis, not an independent
reaudit of those suppliers or a literature-priority claim. JC(2) is OPEN.

Five concepts organize the useful connections: the exact source ring,
the actual invariant and vector field, the radical differential, the
whole discriminant pencil, and labelled special-fibre principal parts.

## 1. Genus has closed the polynomial cusp carrier; its missing elliptic sidecar was a vector-field pole

**Sources:** `05-knowledge/results/continuing9_20260907_flow_carrier_genus.md`,
`continuing9_20260907_ylinear_genus.md`, and current
`planar_jc48_sep07_universal_genus.md`, `planar_jc48_sep08_cusp_ideal.md`.

**Source -> target -> map:** A putative nonzero rational scalar-time
Hamiltonian operation with invariant `S=c+Delta R(p,y)` maps to a rational
selfmap of each geometric generic invariant component. Retain commutation
with its actual vector field, finite relative constants, and literal
source-flow comparison.

**Preserved:** the invariant and commuting vector field. **Lost by genus
alone:** the marked affine divisor p=0. Its restored sidecar is the
nonempty vector-field pole divisor on any genus-one component. For
`S=Delta`, genus is exactly one, but
`delta_S(p)=-Delta*S_y/p` has an actual simple pole at generic p=0 points.
The incoming full cusp-ideal theorem uses this to close the elliptic case.

**Superseded next direction:** another coefficient family or a polynomial
y-degree extension inside `K+Delta K[p,y]` cannot produce rational scalar
time; the universal result already covers every such polynomial. Earlier
genus formulas remain correct and more quantitative.

**Decisive next test (OPEN):** a composition with genuinely different
invariants must retain an actual common preserved object, or be tested
directly in the original source ring. Individual nonrationality does not
exclude a rational composition; inverse cancellation is a compulsory
hostile. This review supplies no new composition obstruction.

## 2. Source-linear closure now feeds an exhaustive quartic search space

**Sources:** `continuing10_20260907_dg_linear_carrier.md`,
`planar_jc48_sep08_dg_quadratic.md`, `planar_jc48_sep08_quartic_boundary.md`.
The quadratic/cubic input is recovered canon THM-2071
`quadratic-fiber-square-parity-gate` and THM-2118
`all-degree-cubic-faber-boundary-flux-coprimality`, with full file paths in
the DG filtration note; these are not new low-degree planar results.

**Source -> target -> map:** Restrict a hypothetical pair of global
functions on W to `C[x,t]`; then use polynomial inversion of any
low-degree pencil member to express the original x in those global
functions. This contradicts its actual pole `x=1/r` along D.

**Preserved:** globality of both functions, constant Jacobian on U0, and
the entire constant output pencil. **Lost if one only quotes planar
inversion:** whether the inverse coordinate extends to W. The pole of x
is the decisive sidecar.

The current exhaustive filtration is
`L_n=span{x^a t^(n-j)(1+x^2 t)^j: 0<=a<=2n, 0<=j<=n}`,
with dimension `(2n+1)(n+1)`. Every nonconstant pencil member of source
degree <=3 is now excluded. Thus the continuing10 phrase “degree at least
two” is a correct local theorem boundary but a superseded search frontier.

**Decisive next test:** start in L4 (dimension45), using the paid global
square prefix `F=H^2+L`, `H in L2`, `L in L1` nonconstant. The approximate
root is global because its finite-pole proof accepts the actual second
chart Jacobian `lambda*r^2`. Neither H nor L inherits a mate or belongs
to the original output pencil. Discarding either distinction would
incorrectly close the quartic problem.

## 3. Leading exactness connects the old residue obstruction to an all-degree necessary filter; full pencils supply the quadratic converse

**Sources:** `planar_jc48_sep08_leading_exactness.md`,
`planar_jc48_sep08_boundary_exactness.md`,
`planar_jc48_sep08_exact_pencils.md`.

**Source -> target -> map:** For `F=A(x)t^n+...` and a rational mate,
formally solve `F(x,T)=v^-n` over `K=C(x)(a)`, `a^n=A`, and extract the
v^(n-1) coefficient of the mate. This proves exactness of `dx/a` in K.
For n=1 this recovers precisely continuing10's rational residue test.

**Preserved:** the primitive equation in its actual radical coefficient
field; no mate-degree bound. **Lost:** all lower coefficients. Pure
monomials have a converse; arbitrary polynomials do not. For a global
quartic square prefix the filter is `dx/sqrt(N)` where N is H's leading
octic section, and the full degree-eight classification applies in the
fixed affine x coordinate with the infinity weight retained.

There is a stronger, unsaturated bridge when the first function itself
is quadratic: `H=Nt^2+Pt+Q`, `D0=P^2-4NQ`. The map
`z=2Nt+P` gives the *whole* generic discriminant pencil
`z^2=D0+4cN` and relative equation `dG=-dx/z`. A rational mate exists iff
the independent span of `(D0,N)` is one of the six exact spaces
`span{1,x}`, `u^3 span{1,u}`, `u^5 span{1,u}`,
`u^7 span{1,u}`, `u^4 span{1,u^2}`, `u^5 span{1,u^3}`,
or, in the proportional case, the single N differential is exact.

**Decisive next test:** in an actual surviving nonconstant-L quartic,
extract the next formal coefficient on the retained radical field and
ask whether it forces L constant. This has already closed every
5+1+1+1 and every (4,3,1) placement, every (4,4), every single-root octic,
and (4,2,2). Those routes should not be rerun. Passing a single exact
leading N is insufficient; the full two-variable lower rows remain.

## 4. A concrete component-torsion bridge for continuing10's fourfold-root submersion

**Sources:** `continuing10_20260907_dg_linear_carrier.md`,
`planar_jc48_sep06_torsion.md`, `planar_jc48_sep08_unit_torsion.md`.
**Status of this bridge:** **PROVED + INDEPENDENTLY AUDITED** explicit analytic corollary; the
[independent referee](continuing10_20260908_fourfold_unit_bridge_audit.md)
checks the full source-ring and component argument. Exact order-three unit controls
already exist in the general torsion supplier. The point is to identify
the actual globally critical-free DG hostile with that component module.

Let `a*h!=0`,
`F=a(x-h)^4*t+a(x^2-4hx)+B0`, `f0=B0-3ah^2`, `g=F-f0`,
`z=x-h`, `K=z^3*t+z-2h`. Then

```
g=a*z*K,       G0=1/(3a*z^3),       D_F G0=1,
P3=g^3*G0=(a^2/3)*K^3 in C[x,t],    D_F P3=g^3.
```

F is smooth on A2 (indeed on all W); `C(x,t)=C(F)(x)` pays rational
constants `C(F)`. The two disjoint reduced components of g=0 are E1:z=0
and E2:K=0. G0 is regular on E2 and has the entire scalar principal part

```
pp_E1(G0)=-8a^2*h^3/(3g^3)-2a*h/g^2,      pp_E2(G0)=0.
```

There is no g^-1 term: subtracting the displayed expression from G0
leaves a rational function regular at z=0, because
`K^3+6hzK+8h^3` is divisible by z^3 and K(0)=-2h.
The complete torsion theorem therefore gives exact primary order3 for
`theta=[1]` in `C[x,t]/D_F C[x,t]`, and order `3+j` for its j-th
canonical connection derivative. Equivalently, g^3 has a polynomial
primitive, whereas g^2 does not: a fibre-constant correction cancelling
its pole on E1 necessarily introduces a pole on E2.

**Preserved:** fixed F, the exact polynomial ring and component labels.
**Lost by generic A1 geometry or rational integration:** unequal special
fibre principal parts. Passing modulo the diagonal removes exactly the
allowed rational repair H(F), rather than an arbitrary local repair.

**Connection to incoming work:** its boundary-linear family
`F=a*rb+c+h(b)`, `h in b^3 C[b]`, realizes exact unit order2 and is smooth
only on the affine source. Our existing fourfold hostile realizes order3
and is smooth on all W. Thus the same component mechanism survives
imposing the stronger global-submersion sidecar.

**Decisive next test:** require an explicit map of a moving-source
response into this fixed-F module that respects its derivation and
polynomial ring. Numerical similarity of response vectors is not such
a map. No such entry is supplied here.

## 5. Exact rational mates of high genus are compulsory hostiles to transferring the flow obstruction

**Sources:** `planar_jc48_sep08_genus_capacity.md`,
`planar_jc48_sep08_dg_genus.md`, `planar_jc48_sep08_binomial_exactness.md`.

**Source -> target -> map:** For a genuine global quadratic on
`W_m=(P1_x x P1_z) minus {z=x^m}`, the full discriminant pencil supplies
an exact relative differential and possibly a rational mate. The sharp
generic genus bound is m-1, attained for every m by the sparse pencil
`u^(2m+1)[a(c)u^(2m-1)+b(c)]`.

**Preserved:** the actual normalized fibre and its particular relative
differential. **Lost by “genus” alone:** which meromorphic differential
is being integrated and which operation is claimed. These explicit
genus>=2 rational mates for m>=3 do not contradict the cusp-carrier
nonrational scalar-time result. A primitive on a curve and a commuting
rational scalar-time selfmap are different predicates.

**Decisive next test:** any proposed transfer from the genus obstruction
to rational mate nonexistence must carry the vector field and its
commutation condition, then survive an extremal W3 example. Conversely,
successful rational integration still needs polynomial/global pole
repair; the recovered quadratic rigidity excludes those polynomial
mates. The current all-exponent binomial iff classifies the differential,
not arbitrary Hamiltonian flows or lower-row Keller pairs.

Reproduction for the new explicit bridge:

```
python -B 04-computation/continuing10_20260908_fourfold_unit_bridge.py
python -B -O 04-computation/continuing10_20260908_fourfold_unit_bridge.py
```

Filing changes status and routing only. Frozen exact producer and auditor source/output bytes are unchanged. The response module is C[x,t]; generic affine-source fibres are punctured lines, and become A1 on W after adding the boundary point.
