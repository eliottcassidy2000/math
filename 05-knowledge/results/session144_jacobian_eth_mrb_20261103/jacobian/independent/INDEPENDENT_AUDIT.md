# Session 144: independent Jacobian classification audit

Status: PROVED (explicit algebraic arguments below) and FINITE-EXACT (the companion program). No global classification of Keller maps, formal proof-assistant build, novelty claim, or publication credit is asserted. Source canon was read through `origin/main` in the read-only source checkout, not edited. In particular, no central mistakes ledger was read or changed.

## Accepted checkpoint 1: composition tower and exact scope correction

Let the source coordinates be `(x,y,z)` and let `F` be the displayed THM-1310 map. Reorder its target coordinates to obtain `G=(F3,F2,F1)`, so source and target weights are both `(1,-1,-2)` and `det DG=2`. Put

- `v=xy`, `t=x^2 z`, `l=2-3v-t`;
- `K=(1+v)^2 t+v^2(4+3v)`;
- `A=v+3K`, `B=(1+v)K`;
- `H=(l A,l^2 B)`.

Then `G=(x l,A/x,B/x^2)` in these reordered coordinates, the displayed divisions being polynomial after substitution. Direct expansion gives `det DH=2l^2`, `H(0)=0`, and `l(0)=2`.

For `j>=0`, define `f_j=l o H^j`; for `m>=1`, define `l_m=product_{0<=j<m}f_j`. Composition gives the exact formulas

- the quotient map of `G^m` is `H^m`;
- its orbit scale is `l_m`;
- `det D(G^m)=2^m`;
- `det D(H^m)=2^m l_m^2`.

The generic degree is `3^m`: the fixed-map degree is three (the complete complex-fiber proof below is an independent certificate), and degrees of dominant generically finite self-maps multiply. This is generic field degree, not ordinary polynomial degree.

### The higher-degree scale already exists

`l_2=l(2-3l A-l^2 B)` has total degree seven, with highest homogeneous part `v^3(t+3v)^4`. Hence THM-4727's unrestricted sentence asking whether solutions with `deg l>=2` exist is false. The same composition tower gives arbitrarily high degree scales. This is an explicit family of solutions of the twisted equation, not a classification of all of its solutions.

### Disjoint smooth divisors

For `i<j`, both coordinates of `H^(i+1)` are divisible by `f_i`. Every coordinate of every iterate `H^r` has zero constant term. Consequently `f_j-2` is divisible by `f_i`. Thus `(f_i,f_j)=(1)` over characteristic zero: their zero sets are disjoint.

At a zero of `f_j`, none of `f_0,...,f_(j-1)` vanish, so `det D(H^j)=2^j product_{i<j}f_i^2` is nonzero. Since `dl=(-3,-1)` is nonzero, `df_j=dl D(H^j)` is nonzero. Each divisor `V(f_j)` is therefore smooth and reduced. Smoothness does **not** prove irreducibility; this audit makes no irreducibility assertion.

The identity `3B-(1+v)A=-v(1+v)` shows that `A=B=0` implies `v=0` or `v=-1`. The latter is excluded by `A(-1,t)=2`, while the former forces `t=0`. Hence `H^-1(0)={0} union V(l)`. Induction gives

`(H^m)^-1(0)={0} union union_{0<=i<m}V(f_i)`.

These are fibers of the quotient map `H^m`, not fibers of the etale map `G^m`. The positive-dimensional quotient fiber does not contradict finite fibers of `G^m`: orbit information has been discarded.

### Degree recurrence

The degrees of the two components of `H` are `(4,6)`. If input polynomials have degrees `a<b`, the unique highest-degree monomials after substitution are `-3V^2 T^2` in `H1` and `V^3 T^3` in `H2`. Thus the next degrees are `(2(a+b),3(a+b))`, with no possible top-degree cancellation. It follows that

`deg(H1^j)=4*5^(j-1)` and `deg(H2^j)=6*5^(j-1)` for `j>=1`.

Therefore `deg f_0=1`, `deg f_j=6*5^(j-1)` for `j>=1`, and

`deg l_m=(3*5^(m-1)-1)/2`.

### Scope and equivalence requirements

The tower disproves the unrestricted higher-scale existence question. It does not supply primitive or indecomposable solutions outside known generators. Such an open question must first specify allowable compositions and an equivalence relation. Generic degree distinguishes the iterates under independent polynomial-automorphism composition on the left and right, but the degree of the displayed scale is coordinate dependent and is not asserted invariant under arbitrary such automorphisms. Moreover, unrestricted left-right composition is not a congruence for a monoid of endomorphisms, so quotienting by it does not automatically give a well-defined composition monoid.

An additional necessary condition in the general torus lemma is stronger than equality of weight sums. Since all source and target weights are nonzero, equivariance forces `F(0)=0`; an invertible derivative at the origin intertwines the two torus representations. Their entire weight multisets must agree. Equal sums alone do not imply existence: `(1,-1,-3)` and `(1,-2,-2)` both have negative-weight sum four but cannot be the source and target weights of an equivariant Keller map. The determinant formula itself is correct as a rational-function identity.

Calling `{x=0}` a "fixed plane" is informal and potentially misleading: it is an invariant hyperplane, not pointwise fixed by the torus.

## Accepted checkpoint 2: complete complex fibers of the fixed map

This classifies fibers of this particular map only, not all dimension-three counterexamples.

Write target coordinates `(a,b,g)` and set

`L=27a^2 g^2-18abg+16a+b^3g-b^2`,
`Q=27ag^2-9bg+8`, `k=4-3bg`,
`N(x)=Lx^3+kx-2g`.

Direct expansion gives `Q^2=k^3+27g^2 L`.

### Correct exact fiber equations

For `x!=0`, let `u=1+xy`. Reconstruct

`y=(u-1)/x`, `z=(x(5-3u)-g)/x^3`.

The first two coordinates of `F=(a,b,g)` are equivalent to

`E1=a x^3+g u^3-xu(u+1)=0`,
`E2=b x^2-x(2+4u)+3g u^2=0`.

These equations are exact even at `u=0`. Their resultant in `u` is `g x^3 N(x)`. For `g!=0`, their degree-one subresultant is

`g x (D u+E)`, with `D=kx-3g` and `E=(9ag-b)x^2+2x`.

The discarded relation in THM-1310's conic pair satisfies `u E2-3E1=-x C1`, where `C1=3a x^2-u(bx-(u-1))`. Thus `(C1,E1)` is equivalent to `(E1,E2)` only after imposing `u!=0`, or after separately restoring the wall equation. At `u=0`, exact fibers require `a=0` and `bx=2`. A minimal spurious conic solution is target `(0,0,0)`, `x=1,u=0`: reconstruction gives `(1,-1,5)`, whose image is `(0,2,0)`.

### Case g=0

The point `x=0` contributes exactly `(0,b,a-4b^2)`. For `x!=0`, `E2` gives `u=(bx-2)/4`, and `E1` becomes `Lx^2=-4`. Thus there are two further points when `L!=0` and none when `L=0`.

### Case g!=0, L!=0

No root of `N` is zero. The resultant criterion applies without loss, since the leading coefficients in `u` of `E1,E2` are nonzero constants `g,3g`.

If `Q!=0`, the discriminant `-4LQ^2` is nonzero, so `N` has three distinct roots. At any root, `D=0` would force `E=0`, and substituting `x=3g/k` into `E` would force `Q=0`. Thus each root gives exactly one `u`, hence one source point.

If `Q=0`, then `k!=0` by `Q^2=k^3+27g^2L`. The cubic factors as a nonzero constant times `(x-3g/k)^2(x+6g/k)`. At the repeated root `D=E=0`, so the quadratic `E2` divides `E1`. Its discriminant is `108g^2/k`, which is nonzero, giving two distinct values of `u`. The remaining root has `D!=0` and gives one point. Therefore there are still three distinct source points, despite only two distinct x coordinates.

### Case g!=0, L=0

If `Q!=0`, then `k!=0`, and `N` has the unique root `x=2g/k`. Here `D=-g!=0`, giving exactly one source point.

If `Q=0`, then `k=0` and `N=-2g` has no root. Therefore the fiber is empty. The missing curve is exactly

`g!=0, b=4/(3g), a=4/(27g^2)`.

### Complete result and nonproperness

Over the complex numbers,

- `L!=0`: exactly three distinct source points;
- `L=0,Q!=0`: exactly one source point;
- `L=Q=0`: no source point.

Thus the generic field degree is three. The nonproper value set is exactly `V(L)`: off it, the three local inverse branches around the three distinct fiber points account for every nearby fiber, so preimages of a small neighborhood cannot escape. At a point of `V(L)`, approach through generic targets. If all three source points stayed bounded, a convergent subsequence would give three local branches through the limiting fiber; local invertibility prevents two from merging, contradicting its size zero or one. Hence there is an escaping sequence.

This argument establishes exact nonproperness with a direct topological certificate rather than inferring it merely from a vanishing leading coefficient.

### Trace law and fold warning

The sum of the x coordinates over the three actual source points is zero **off `V(L)`**, counting each source point once, even when two have the same x coordinate. It is not an all-target statement. Target `(0,1,1)` has the singleton preimage `(2,-1/2,9/8)`, so its fiber trace is two. At target `(1/27,1,1)`, the x coordinates are `(3,3,-6)`; forgetting multiplicity gives the different support sum `3-6=-3`.

The divisor `V(Q)` is a failure of the chosen x coordinate to separate source points, not ramification of the Keller map. Outside `V(L)`, its fibers remain three distinct points.

## Other definite correction: local versus global deformations

Noninjectivity is locally open in a fixed finite-dimensional coefficient family of Keller maps near a collision: use two disjoint local inverse neighborhoods. It does not hold for every member of an arbitrary Keller deformation. The polynomial family

`F_h(X)=h^-1 F(hX)` for `h!=0`, `F_0=DF(0)`

extends polynomially because `F(0)=0`. It has determinant `-2` for every parameter, is noninjective for every `h!=0`, and has the invertible linear endpoint `F_0(x,y,z)=(z,y,2x)`. This gives an explicit counterexample to THM-1310's global deformation wording.

## Reproduction and trust boundary

Run `python verify_audit.py --output audit_checks.json` and `python -O verify_audit.py --output audit_checks.optimized.json`, then compare the files. Every obligation is evaluated by explicit checks, not Python assertions. The script includes seven rejected corruptions and an exhaustive independent source enumeration versus fiber reconstruction for all targets over `F_5,F_7,F_11` (1,799 targets total). Finite-field tests certify the reconstruction algorithm in those exact universes, not the all-complex classification; the latter is proved above. The ordinary and optimized outputs must be identical.
