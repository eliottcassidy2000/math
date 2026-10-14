# All-degree accessory clutches: a compact-fibre time-form obstruction

**Date:** 2026-10-14. **Status:** PROVED by a complete Newton–Puiseux puncture ledger; independently checked. **Novelty unconfirmed.** This excludes all polynomial clutches in the stated fixed source architecture. It does not establish that every such clutch has a critical point, and does not resolve the planar Jacobian conjecture.

## 1. Statement and inherited inputs

Work over an algebraically closed field of characteristic zero. Suppose

- `V,A in k[x]`, with `deg V=16`, `deg A=8`;
- at every root `alpha` of `V`, its multiplicity is either one or at least three, and `A` has a simple zero there.

For any polynomial `B`, define

`P_B(x,z)=(V(x)z^2+B(x)z)^2+A(x)z+x`.

### Theorem 1.1: every polynomial clutch is excluded as a Keller component

There is no polynomial `H(x,z)` for which `Jac_(x,z)(P_B,H)` is a nonzero constant.

More precisely:

1. If `B` and `V` have a common zero `alpha`, then `P_B` has the explicit critical point `(alpha,-1/A'(alpha))`.
2. If `gcd(B,V)=1`, then on every connected component of a sufficiently general smooth fibre of `P_B`, its Hamiltonian time differential extends to a nonzero holomorphic differential on the smooth projective completion. Such a differential cannot be the derivative of a rational function, so no Jacobian mate exists.

The two fixed accessory pairs in THM-3225 satisfy these hypotheses: their degrees are `(16,8)`; the unique `S` root is simple; the `T` roots have multiplicities among `3,4,5,6`; and `A` has a simple zero at every boundary root. These are trusted inputs from [THM-3225](../../01-canon/theorems/THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape.md), Sections 1, 3, and 4. No accessory coefficients are reclassified here.

The response identity `2VA'-AV'=2V` is not required by the proof once the listed degrees and local orders are available. The multiplicity-two exclusion is likewise unnecessary for the local calculation, but the statement retains the actual inherited boundary pattern.

## 2. The exact differential obstruction

Fix `P=P_B`. On a smooth fibre `C`, its Hamiltonian derivation is

`D_P=P_x partial_z-P_z partial_x`.

The tangent time form is

`omega=-dx/P_z=dz/P_x`.

These expressions agree where both denominators are nonzero, because `dP=0` on the fibre. Since the fibre is smooth, they cover it and define a regular, nowhere-zero differential on its affine part; equivalently, `omega(D_P)=1`.

If `Jac(P,H)=j != 0`, then on `C`

`d(H|C)=j omega`.

Suppose `omega` extends holomorphically to the smooth projective completion of a connected component of `C`. The restriction of `H` is a rational function there. A pole of order `d>0` of a rational function gives a pole of order `d+1` in its differential in characteristic zero. Thus a rational function with holomorphic differential has no poles, and is constant on the projective curve. Its differential is zero, contradicting `j omega != 0`.

No generic irreducibility hypothesis is needed. Applying this argument to one connected component already suffices; the puncture calculation below works on all of them.

The time-form/period method is classical. For a primary treatment, see Bustinduy, Giraldo, and Muciño-Raymundo, *Jacobian mates for non-singular polynomial maps in C^n with one-dimensional fibers*, Journal of Singularities 9 (2014), 27–42, [DOI 10.5427/jsing.2014.9b](https://journalofsingularities.com/volume9/article2.html), Section 2 and Remark 4. The new claim here is the explicit all-degree regularity calculation for this source architecture, not a new general exactness principle.

## 3. Generic fibre and exhaustive list of possible punctures

Assume `gcd(B,V)=1`, so `B` is nonzero. Set `n=deg B`.

Choose a fibre value `c` outside the following finite set:

- the critical values of `P`;
- every root `alpha` of `V`;
- for each simple root of `V`, the number `alpha-A'(alpha)B(alpha)/V'(alpha)`;
- in the degree-nine resonant case described in Section 5, the finite limit of `x-AB/V` at infinity.

The critical values form a finite set: on each positive-dimensional component of the critical locus the differential of `P` vanishes, so `P` is constant in characteristic zero, and there are only finitely many components. The remaining exclusions are explicitly finite. Such a `c` exists, and its fibre is smooth.

On the projective normalization of a fibre component, a missing affine point either has `x=infinity`, or has finite `x=alpha` and a pole in `z`. If `V(alpha) != 0`, the degree-four term `V^2z^4` prevents a solution with finite `x` and a pole in `z`. Hence the only possible finite-`x` punctures lie over roots of `V`.

There are no vertical fibre components: at a fixed `x=alpha`, either `V(alpha) != 0` gives a nonconstant quartic in `z`, or `B(alpha) != 0` gives a nonconstant quadratic. Thus `x` is nonconstant on every component. The fibre equation has degree four in `z` over `k(x)`. In each local calculation below we construct four Puiseux roots with simple initial equations. The Newton–Puiseux/Hensel theorem therefore supplies the corresponding full branches, and the count of four proves that no additional branch has been missed. Puiseux roots may belong to a single ramified place; no assertion that four roots mean four distinct punctures is used.

## 4. Every finite boundary puncture is harmless for the time form

Let `t=x-alpha` at a root of `V`, and write

`V=v t^m+...`, `A=a t+...`, `B=b+...`,

with `v,a,b != 0`.

### 4.1 Two finite branches

Setting `t=0` in the fibre equation gives

`b^2 z^2+alpha-c=0`.

Our choice `c != alpha` makes its two roots simple and nonzero. They give two finite affine branches with `P_z != 0` at `t=0`.

### 4.2 Two escaping branches

Use the exact rational shift

`z=-B/V+h`, `E=x-AB/V-c`.

Then the equation and its `z` derivative become

`(Vh^2-Bh)^2+A h+E=0`,

`P_z=2(Vh^2-Bh)(2Vh-B)+A`.                         (1)

If `m>=3`, then

`E ~ -(ab/v)t^(1-m)`.

The other two roots have

`h ~ eta t^(-(m-1)/2)`, with `b^2 eta^2=ab/v`.

Both initial roots are nonzero and simple. The term `Vh` tends to zero, `A h` is of strictly higher order than `E`, and (1) gives

`P_z ~ 2b^2 eta t^(-(m-1)/2)`.

After the uniform ramified substitution `t=s^2`, the form `omega` has order `m` in `s`, so it is holomorphic. If the actual branch is unramified, holomorphicity on this cover still implies holomorphicity on the branch.

If `m=1`, then `E` tends to

`E_0=alpha-ab/v-c`,

which is nonzero by the generic-value exclusion. The escaping branches have `h` tending to the two nonzero simple roots of

`b^2 h^2+E_0=0`.

Equation (1) gives `P_z -> 2b^2 h != 0`. Thus `omega` is holomorphic here too.

This accounts for all four local roots at every finite boundary point. In particular, no pole of the time form is being hidden at a high-multiplicity `T` root or at the simple `S` root.

## 5. All infinity branches, including the degree-nine resonance

Write the leading coefficients of `V,A,B` as `v,a,b`, respectively. We describe growth as `x` tends to infinity. Every initial coefficient displayed below is nonzero.

### 5.1 Degrees `0<=n<=8`

All four roots have

`z ~ eta x^(-31/4)`, with `v^2 eta^4+1=0`.

Indeed, `Vz^2` dominates `Bz`, while `(Vz^2)^2` balances `x`; `Az` is lower order. The four roots of the initial quartic are simple. Differentiation gives

`P_z ~ 4v^2 eta^3 x^(35/4)`.

With `x=s^(-4)`, `omega` has order `30` in `s`. It is holomorphic, with positive vanishing order.

### 5.2 Degrees `n>=9`: two small roots

Two roots have

`z ~ eta x^(1/2-n)`, with `b^2 eta^2+1=0`.

Now `Bz` dominates `Vz^2` and its square balances `x`; `Az` is lower order. The initial roots are simple, and

`P_z ~ 2b^2 eta x^(n+1/2)`.

With `x=s^(-2)`, the order of `omega` is `2n-2>0`.

### 5.3 Degrees `n>=9`: two roots near `-B/V`

Use the same shift and `E` as in (1). The growth degree `e` of `E` is exactly

- `e=n-8` if `n>=10`, with leading coefficient `epsilon=-ab/v`;
- `e=1` if `n=9` and `ab/v != 1`, with leading coefficient `epsilon=1-ab/v`;
- `e=0` if `n=9` and `ab/v=1`, with leading coefficient `epsilon=lim_(x->infinity)(x-AB/V)-c != 0`.

The last inequality is precisely the extra generic-fibre-value exclusion. It handles the genuine degree-nine leading cancellation instead of silently assuming it away.

The remaining two roots have

`h ~ eta x^(e/2-n)`, with `b^2 eta^2+epsilon=0`.

In every case, `Vh/B` tends to zero and `A h` has smaller growth than `E`. Explicitly the positive gaps on the common cover `x=s^(-2)` are

`4n-e-32 > 0`, `2n+e-16 > 0`.

Thus the displayed two simple quadratic roots are valid initial forms, and (1) gives

`P_z ~ 2b^2 eta x^(n+e/2)`.

On the same cover, `omega` has order `2n+e-3>0`. In particular, on the degree-nine resonance wall it has order `15` on this cover, rather than acquiring a pole.

These two roots are distinct from the small roots: their main term is `z=-B/V`, of growth `x^(n-16)`, whereas the small roots have growth `x^(1/2-n)`. The inequality `2n>16+1/2` holds for every integer `n>=9`. We have therefore accounted for all four roots at infinity in every degree.

## 6. Conclusion and exact scope

Every puncture has now been checked. The regular, nowhere-zero affine time form extends holomorphically to each smooth projective fibre component. It remains nonzero there. Section 2 rules out a polynomial constant-nonzero-Jacobian mate.

If instead `B(alpha)=V(alpha)=0`, then the square term and its first derivatives vanish on that fibre, while `A(alpha)=0`; hence

`P_z(alpha,z)=0`, `P_x(alpha,z)=A'(alpha)z+1`.

The stated finite critical point follows, and again no mate can exist. This proves Theorem 1.1 for every polynomial `B`.

The earlier critical-resultant hypothesis H-J is strictly stronger: it asks whether every boundary-avoiding clutch retains a finite critical point. That remains OPEN. The present theorem bypasses that question and excludes polynomial mates directly. The toy exception `V=x,A=2x,B=x/2-4` explains the distinction, but does not satisfy the degree hypotheses of this theorem.

No genus value, connectedness claim, general planar-JC conclusion, or literature priority claim is needed or made.

## 7. Exact companion

`all_degree_accessory_time_form_obstruction_20261014.py` checks the exact shifted fibre and derivative identities, all leading polynomial models, the simple-root conditions, and the positive valuation gaps on every degree sector. Infinite sectors are checked as affine expressions in nonnegative integer parameters, rather than inferred from a finite degree census. It also verifies the Broughton pullback and generic-fibre time-form identities in the separate synthesis note. The mathematical completeness argument is the four-root Newton–Puiseux proof above.

Normal and optimized Python must produce identical JSON. This worker does not modify a live checkout or push to a remote.
