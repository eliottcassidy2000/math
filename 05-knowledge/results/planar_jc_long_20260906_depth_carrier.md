# The maximal polynomial Hamiltonian carrier preserving actual depth

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
For polynomial generators in the cusp plane, preserving depth zero is
equivalent to preserving every actual source depth. The exact carrier is
larger than the universal affine-source carrier; their intersection is
computed below. No polynomial flow or Keller pair is asserted.

The [independent audit](planar_jc_long_20260906_depth_carrier_audit.md)
accepts all iff directions and uses generic quotient differentiation to
check the generator formulas, independently of literal source derivatives.
Its 56 live checks pass in byte-identical normal/optimized runs.

Use the characteristic-zero ring and sign convention of
[the completed Hamiltonian theorem](planar_jc_long_20260906_hamiltonian.md):

```
A=K[p,y], D=p^3-y^2,
B_2=K[x,u,p,y] subset K[x,t],
u=x^2t, p=t(1+x^2t), y=xtp,
P_d=sum_(a+b<=d) x^a u^b A,
delta={-,S}, {f,g}=f_xg_t-f_tg_x, S in A.
```

The actual module identification is inherited from
[THM-3989 / cusp-log-laurent-conductor-and-nondividing-depth-reduction](../../01-canon/theorems/THM-3989-cusp-log-laurent-conductor-and-nondividing-depth-reduction.md)
and [THM-4308 / source-normal-bracket-hasse-truncation-through-row-eight](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md).
The closest mechanism is the universal source-carrier theorem in
[planar_jc48_sep06_hamiltonian.md](planar_jc48_sep06_hamiltonian.md).
Its fixed-source control `S=D` is the canonical hostile to confusing the
two carrier conditions. The needed sidecar is the actual module on which
the derivation acts, before projecting to finite rows.

The three live concepts are preservation of `P_0`, preservation of the
whole filtration, and preservation of `-u/2+A`. Exact carrier/depth searches
found the incoming smaller source-carrier theorem but no maximal
depth-carrier classification. This is a repository check, not a literature
priority claim.

## 1. Exact equivalence and all generator images

For every `S in A`, the following are equivalent:

```
delta(P_0) subset P_0;
delta(P_d) subset P_d for every d>=0;
S=c+p^2 R,                    c in K, R in A.          (1)
```

**Necessity.** The inherited bracket `{p,y}=-D/p` gives
`delta p=-D S_y/p` and `delta y=D S_p/p`. Since `p` and `D` are coprime,
polynomiality of these two images forces `p|S_y` and `p|S_p`. Expand `S`
by powers of `p`. The first condition makes its constant coefficient a
scalar, using characteristic zero; the second kills the coefficient of
`p^1`. This is exactly `S=c+p^2R`.

**Sufficiency.** Literal differentiation, or the quotient identities
`x=yp/D`, `u=y^2/D`, gives

```
delta p=-pD R_y,
delta y=D(2R+pR_p),
delta u=(1+u)y(4R+2pR_p+3yR_y),
delta x=p(1+2u)(2R+pR_p)+py(2+3u)R_y.                 (2)
```

Thus `delta p,delta y in P_0` and `delta x,delta u in P_1`. The Leibniz
rule on each `x^a u^b f(p,y)` preserves the bound `a+b<=d`. This proves
all of (1), including that it is enough to test depth zero. These are
properties of the whole actual depth module at every `d`, not a finite
rank experiment.

## 2. Intersection with affine-source preservation

Within (1), the following are also equivalent:

```
delta(-u/2+H) in A for some H in A;
delta(-u/2+H) in A for every H in A;
D divides R.                                         (3)
```

Indeed `delta H in A` already, so either source condition is equivalent
to `delta u in A`. Put `E=2p partial_p+3y partial_y`. Formula (2) says

```
delta u=p^3 y (4+E)R / D.
```

The numerator factor `p^3y` is coprime to `D` in `A`. Hence the condition
is `D|(4+E)R`. The Euler derivation descends to the cusp ring
`A/(D)=K[z^2,z^3]`, acting on weight `n` as multiplication by `n`.
The operator `4+E` has eigenvalue `n+4` on every nonnegative weight.
Characteristic zero therefore makes it injective on the cusp ring.
Consequently `D|(4+E)R` if and only if `D|R`. This proves (3) and recovers
exactly the inherited smaller carrier `K+p^2D A`.

The control `S=p^2` separates the conditions. It preserves every `P_d`,
but `delta u=4(1+u)y` is not in `A`, so it preserves no affine source
`-u/2+H`. It also does not lower every positive depth: in the logarithmic
chart `p=s^2+tau`, `y=s(s^2+tau)`, `u=s^2/tau`, its image `delta u` has
leading term `4s^5/tau`. In the literal source chart `delta x` begins at
order `t`, so it is not in `D B_2`, whose elements start at order at least
`t^3`. Thus the stronger uniform convergence mechanism of the completed
Hamiltonian theorem requires its smaller carrier; mere depth preservation
does not supply that mechanism.

Conversely `S=D` satisfies the incoming fixed-source test at `H=0`, but
`delta p=2yD/p` is not in `A`. This is an actual fixed-source generator
outside the depth-preserving carrier, so the two original conditions are
not interchangeable.

## 3. No nonconstant carrier is locally nilpotent

Every nonconstant generator in (1) has a non-locally-nilpotent derivation
on `K[x,t]`. This extends the incoming factorial-closure argument from
`K+D A` to `K+p^2 A`.

Suppose instead that `delta` were locally nilpotent. In a characteristic-zero
domain, `deg_delta(fg)=deg_delta f+deg_delta g` for nonzero polynomials;
this follows from the unique highest nonzero term of the Leibniz formula.
Its kernel is consequently factorially closed. The invariant nonzero
polynomial

```
S-c=p^2 R=t^2(1+x^2t)^2 R(p,y)
```

would make both `t` and `1+x^2t` invariant. Subtracting one and applying
factorial closure again makes `x` invariant. Thus `delta` would vanish on
both polynomial source coordinates, forcing its Hamiltonian to be constant.
The cusp-plane substitution is injective by its rational inverse, giving
a contradiction. This excludes polynomial additive-group actions with
this infinitesimal generator. It does not decide whether an individual
scalar specialization of a separately convergent exponential is polynomial.

## 4. Exact controls and reproduction

The standalone [source](../../04-computation/planar_jc_long_20260906_depth_carrier.py)
uses literal source differentiation and the full Laurent depth criterion
without importing a repository producer. The [frozen output](planar_jc_long_20260906_depth_carrier.out)
passes **1,589 always-active exact gates** with byte-identical normal and
optimized output. The finite universe has eight carriers, each acting on
all 60 generator monomials with `a+b<=3` and `c+e<=2`; 28 generator
monomials of total `(p,y)` degree at most six for the necessity check;
the same 28 cusp eigenvalue and intersection controls; and the two
separating generators. Universal claims rest on the proofs above.

```bash
python3 -B 04-computation/planar_jc_long_20260906_depth_carrier.py > /tmp/depth_carrier_normal.out
python3 -B -O 04-computation/planar_jc_long_20260906_depth_carrier.py > /tmp/depth_carrier_optimized.out
cmp /tmp/depth_carrier_normal.out /tmp/depth_carrier_optimized.out
cmp /tmp/depth_carrier_normal.out 05-knowledge/results/planar_jc_long_20260906_depth_carrier.out
```

Raw LF-byte SHA-256 manifest:

```
source fc8f000dd985c7b87350327228b1b208b9a1e6d547ad3d083d42facc131f5ba7
output 3899e36bad96bc3878d3ee5b44c24b88e5f5c24568d0d400962ea9a812b22c60
```

No theorem identifier is allocated. The independent analytic and generic
quotient-rule audit is accepted and linked at the start of this note.
