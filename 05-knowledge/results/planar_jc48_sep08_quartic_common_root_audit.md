# Independent audit: the quartic common-boundary-root obstruction

**Status: VERIFIED / independent analytic and exact-source audit PASS.**
Root referee, September 8, 2026. The audited theorem is in
[the common-root note](planar_jc48_sep08_quartic_common_root.md).
The producer is `orthogonal_returns`; root independently read the full
proof and source, reconstructed the argument below, and replayed normal
and optimized programs against the frozen output. `three_ray_geometry`
independently read the analytic proof as a second hostile reader and found
no gap. That second read does not substitute for the root checks here.

## 1. Precisely what passes

On the fixed smooth surface `W=(P1 x P1)\S`, `S={z=x²}`, take actual
global functions `H in L2`, `L in L1`, with `deg_t H=2`, and write
`H=N/s²`, `L=M/s`. The restrictions to `S` are sections of degrees eight
and four respectively. If those sections have no common zero anywhere
on the complete `S`, then `F=H²+L` has no rational mate in `C(x,t)`
with nonzero constant source Jacobian. Both charts and every normalized
component of a generic fibre are retained. No assertion that an arbitrary
Keller map extends to this particular surface is included.

This strengthens the earlier safe multiplicity-at-most-two criterion:
it permits every octic root partition, including the genuine
multiplicity-six higher-pole hostile. It does not assert that a common
zero is sufficient for a mate, or close the shared-zero cases.

## 2. Types, valuations, and local exhaustiveness

At a boundary root use coordinates `(s,w)` and write
`N=A(w)+sB(w)+s²C(s,w)`, with `m=ord A>=1` and `j=ord B`, including
`j=infinity`. The no-common-zero assumption makes `M` a unit there.
The compact fibre is exactly `E=N²+s³M-cs4=0`; up to a holomorphic
factor its relative differential is `s² dw/E_s`. That factor is a unit
at finite boundary points, and has an additional order-two zero at the
point meeting the closure of `D`. The latter can decrease pole budgets
and cannot create an omitted pole.

The local equation has Weierstrass degree two in `s` for `j=0` and
three for `j>0`, including `B=0`. This pays branch exhaustiveness,
not just the existence of the displayed branches.

For `m<3j`, the only Newton slope is `ord_w s=2m/3`; after clearing
its denominator the leading cubic has three distinct nonzero roots.
The higher term `s²C` in `N` has order strictly larger than `m`, and
`N_s` has order strictly larger than `m/3`. Thus `3Ms²` dominates
`E_s`, and the differential is a unit times `dw`. Puiseux ramification
adds the nonnegative order of the pullback of `dw`. This also proves
the full high-`j` and `B=0` tails, beyond the finite control universe.

For `m>3j` and `j>0`, the simple branch has `s`-order `2j`; using
`B²=-Ms` on its leading equation makes the leading derivative `Ms²`,
which is nonzero. The other two determinations have
`s~-(a/b)w^(m-j)`. Here `N_s~bw^j` and `N²=-s³(M-cs)` force
`ord N=3(m-j)/2`. Since `(m-j)/2>j`, `2NN_s` dominates the
derivative. The differential's exponent before ramification is
`(m-3j)/2>0`. These two determinations may belong to a single normalized
branch; the proof does not count them as two analytic branches when the
parity is odd. For `j=0` they exhaust the degree-two local equation;
the would-be third branch has nonzero limiting `s` and is absent from
this boundary germ.

The remaining slope is exactly `m=3j=3k`, `k>=1`. Setting
`s=w^(2k)Z` and dividing by `w^(6k)` gives
`Q=Q0(w,Z)-cw^(2k)Z4` with face `(a+bZ)²+M0 Z³`.
The terms from `s²C`, variation of `M`, and higher coefficients of
`A,B` all have positive `w` order in this equation. No face root is
zero. A triple root would imply incompatible identities
`9M0²Z0⁴=12M0²Z0⁴`; consequently the only exceptional face is one
double root and one simple root. Simple roots are regular branches.

## 3. The generic double-root split and actual pole budgets

At the double root the derivative `Q_ZZ` is a unit. The analytic critical
centre `z(w,c)` therefore exists, and for `R=Q(w,z(w,c),c)` the chain
rule gives the exact identity

    partial_c R = -w^(2k) z(w,c)^4.

Terms below order `2k` have no `c` dependence, and the coefficient of
order `2k` has nonzero linear slope `-Z0⁴`. Thus the generic split order
`lambda=ord_w R` is between one and `2k`; one does not need to assume
that a chosen coefficient perturbation is generic. There are only finitely
many boundary roots and finitely many additional excluded fibre values.

Writing `Q=R+v²` via the analytic square root of its quadratic Taylor
unit is a valid coordinate change with unit derivative in `Z`.
Accordingly `Q_Z` is a unit times `v` and the differential becomes
a unit times `dw/v`. For `lambda=2l`, two smooth branches each have
pole order `l`; for odd `lambda`, the normalization has `w=tau²`
and differential order `1-lambda`. These give total possible primitive
pole degree exactly `max(lambda-2,0)` in the unit-numerator chart,
and at most that in the chart with the extra zero. A simple pole with
nonzero residue forbids exactness directly; it is never silently treated
as the derivative of a pole of order zero.

The bound is at most `2k-2`. Since `sum m=8`, resonant multiplicities
can only be three and six. Multiplicity three contributes zero;
multiplicity six contributes at most two and can occur only once.
All other multiplicities contribute zero. This bounds the total possible
primitive pole degree over all compact normalized components by two.
It is stronger than necessary, because the contradiction only consumes
the budget on one particular component.

## 4. The degree consumer is on the same component

The restrictions of `H,L` to `D` have degrees at most two and one in
`b`. If `H|D²+L|D` were constant, the degree of a nonconstant square
could not cancel against the linear summand. Thus both restrictions
would be constant. In the complete infinity chart
`Ninf=r4q²N(1/r,1/q)`, `Minf=-r²qM(1/r,1/q)`, they are
`Ninf(0,q)/q²` and `Minf(0,q)/q`. Constancy forces both numerators
to vanish at `(r,q)=(0,0)`, the common point of `S` and the closure
of `D`. This contradicts the full-boundary no-common-zero hypothesis.

Therefore `F|D` is nonconstant and a generic fibre has an unramified
point `(0,b0)` on `D`. On its smooth normalized component, `r` is a
parameter, and `omega/dF` has order exactly two. A rational mate would
restrict to a meromorphic function with nonzero derivative equal to a
constant multiple of that differential, so its local degree there would
be three.

Choose the fibre also to avoid the finitely many critical values of the
surface morphism and the images of vertical pole components of the rational
mate. Its restriction to each compact component is then meromorphic.
Regularity of the differential in `W` prevents any pole of the restricted
function there, even if the original rational expression was not globally
regular on `W`. Every remaining pole is at the counted boundary. Thus
the degree of the map to `P1` on the component containing `(0,b0)`
is at most two and at least three, a contradiction. This argument uses
neither geometric irreducibility of the generic fibre nor an unproved
Hamiltonian constants-field identification. It also does not sum local
degrees at different `D` points with potentially different function values.

## 5. Exact replays and adversarial controls

The source has explicit always-active `RuntimeError` gates, so `-O`
removes no checks. Root read all of it and ran both modes independently;
each output matched the frozen 397 bytes with 568 gates.

| Artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| [Source](../../04-computation/planar_jc48_sep08_quartic_common_root.py) | 7441 | `df2f4ef7ae1fcca65aa9f99e9982ca26eb27768d2dad095a1402b054fbb1185b` |
| [Frozen output](planar_jc48_sep08_quartic_common_root.out), normal, optimized | 397 | `0240b334889275c2cf891fff1f995e7d4479cacab013e65000c18af8a73ebf19` |

The finite universe is declared: `m=1..8`, `j=0..9` and infinity,
all 22 partitions of eight, and every generic split order through `2k`
for `k=1,2`. Symbolic identities cover the full face discriminant,
the no-triple contradiction, and the critical-value derivative.
The analytic argument, not the sample monomial germs, covers arbitrary
higher coefficients and the infinite `j` tail.

The source reconstructs the genuine global multiplicity-six section,
both boundary charts, its residue-free double poles, and its nonconstant
`D` restriction. This hostile attains the complete possible budget of
two. Its affine critical points alone say nothing about rational mates;
the new degree contradiction is needed for that stronger conclusion.

For an independent sharp hypothesis control, take `H=t²`, `L=0`.
Then `F=t4` has rational mate `G=-x/(4t³)` with source Jacobian one.
It has common boundary roots and constant `F|D`, exactly the hypotheses
excluded from the theorem. In the other affine chart the mate is
`1/[4r7(1+r²b)³]`, so this is not a global regular mate on `W`.
The identity follows directly from `F_x=0`, `F_t=4t³`, and
`G_x=-1/(4t³)`. No finite data assert existence of a Keller counterexample.

A stronger shared-root control retains a nonconstant boundary restriction.
Put `h=x²+x4t=-b`, `H=h²`, `L=0`, `F=h4`, and
`G=1/(12x³h³)`. Then `J(h,x^-3)=3`, hence `J(F,G)=1`.
Here `N=x4z²`, `N|S=x8`, and the zero section `M` shares that
root. On the boundary chart `F=b4`, `G=-r³/(12b³)`. On each
generic component `h=constant!=0`, the rational conjugate has exactly
one pole of order three at `x=0` and local degree three at `D`.
Thus the common-root assumption cannot be dropped even after imposing
nonconstant `F|D`; the degree consumer is sharp. This is a rational
mate with poles, not a polynomial Keller pair. The second hostile
reader supplied this control and root independently checked its bracket.

The control also has a version with nonzero remainder section: keep
`h=x²+x4t`, take `H=h²`, `L=h`, `F=h4+h`, and
`G=1/[3x³(4h³+1)]`. The same identity `J(h,x^-3)=3`
gives `J(F,G)=1`. Now `N=x4z²`, `M=x²z` restrict to the
nonzero boundary forms `x8`, `x4`, sharing only the finite root zero.
On the other chart `F=b4-b`, `G=r³/[3(1-4b³)]`.
Every generic compact component again has degree three. This strengthens
the hostile without relying on a zero remainder section.

The bounded exact-statement recovery found the degree-of-primitive
mechanism already in
[THM-2071, quadratic-fiber-square-parity-gate](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md),
and recovered versions in
[THM-2214, nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten](../../01-canon/theorems/THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten.md),
§7.3, and
[THM-2723, split-exact-square-prefix-rational-primitive-pole-capacity](../../01-canon/theorems/THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity.md),
§4. That mechanism is credited. The present proof supplies the complete
DG octic local budget and the actual order-two zero on `D`. It does not
identify its section common roots with the different Faber-observable
common-root condition in THM-2129, nor drop the chosen-sheet entry of
THM-2778. The search found no pre-wave statement identical to the theorem;
it is not a literature priority determination.

Reproduction commands from this worktree root:

```sh
python3 -B 04-computation/planar_jc48_sep08_quartic_common_root.py
python3 -B -O 04-computation/planar_jc48_sep08_quartic_common_root.py
```

## 6. Promotion and boundary

**PASS** for the theorem as stated. The primary may be promoted to
`PROVED + INDEPENDENTLY AUDITED`; frozen computational wording that asks
for an audit is retained as production provenance. This note supplies it.
No frozen source or output was modified during this independent audit.
The next unproved case is an actual shared boundary zero, with the pair
of multiplicities and the full local equation retained. JC(2) remains open.
