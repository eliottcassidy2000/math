# Balanced-degree clutches have a compact-fibre obstruction

**Date:** 2026-10-14. **Status:** PROVED; exact symbolic valuation certificate. **Novelty unconfirmed.** This generalizes the degree `(16,8)` checkpoint without modifying it. It concerns the displayed clutch architecture, not arbitrary planar polynomial maps.

## 1. Two precise statements

Let `k` be algebraically closed of characteristic zero, let `A,V in k[x]` satisfy

`deg A=a>=1`, `deg V=2a`,

and set

`P_B=(Vz^2+Bz)^2+Az+x`.

### Theorem 1.1: the boundary-unit sector

For every polynomial `B` with `gcd(B,V)=1`, the polynomial `P_B` has no polynomial mate `H` satisfying `Jac(P_B,H)=j != 0`.

There is no assumption here on the multiplicities of the roots of `V` or on the orders of `A` at those roots. The response identity is not needed.

### Corollary 1.2: all clutches under the simple-zero boundary hypothesis

If additionally `A` has a simple zero at every root of `V`, then the same no-mate conclusion holds for **every** polynomial `B`, including those meeting the boundary.

Indeed, if `B(alpha)=V(alpha)=0`, then `P_z(alpha,z)=0` and `P_x(alpha,z)=A'(alpha)z+1`; thus `(alpha,-1/A'(alpha))` is a critical point. Otherwise Theorem 1.1 applies. Multiplicity two in `V` is permitted.

The condition `a>=1` cannot be dropped wholesale: if `V,A,B` are constants, then `P_x=1` and `H=z` is a polynomial Jacobian mate.

The fixed accessory pairs of THM-3225 have `a=8` and satisfy the corollary's boundary hypothesis. The earlier toy exception has `(deg V,deg A)=(1,1)` and lies outside this balanced-degree hypothesis; its no-mate proof instead uses a logarithmic residue.

## 2. Why regularity on the completed fibre suffices

Take a sufficiently general smooth fibre `C={P_B=c}`. Its Hamiltonian time form is

`omega=-dx/P_z=dz/P_x`.

It is a regular, nowhere-zero differential on the affine curve. A constant-Jacobian mate would imply `d(H|C)=j omega`. If `omega` is holomorphic on the smooth projective completion of a connected component, then `H|C` has no poles, because differentiation increases the order of any pole in characteristic zero. Hence `H|C` is constant, contradicting the nonzero differential.

This argument needs neither a genus calculation nor generic irreducibility. The classical time-form/period method is discussed by Bustinduy, Giraldo, and Muciño-Raymundo, [Journal of Singularities 9 (2014), 27–42](https://journalofsingularities.com/volume9/article2.html), Section 2 and Remark 4. The assertion proved here is the following complete regularity ledger for all balanced degrees.

## 3. Finite boundary: arbitrary local orders

A puncture with finite `x=alpha` must have `V(alpha)=0`; otherwise the nonzero quartic coefficient forbids a pole in `z`. Write `t=x-alpha` and

`V~v t^m`, `A~d t^r`, `B~b`,

where `m>=1`, `r>=0`, and `v,d,b != 0`.

The two finite roots at `t=0` satisfy

`b^2z^2+A(alpha)z+alpha-c=0`.

Choose `c` so this quadratic has distinct roots. Its discriminant is `A(alpha)^2-4b^2(alpha-c)`, and the `z` derivative is nonzero at each root.

For the two escaping roots, put

`z=-B/V+h`, `E=x-AB/V-c`.

The exact equations are

`(Vh^2-Bh)^2+A h+E=0`,

`P_z=2(Vh^2-Bh)(2Vh-B)+A`.                        (1)

If `r<m`, put `d_0=m-r>0`. Then `E~-(db/v)t^(-d_0)`, and the escaping roots satisfy

`h~eta t^(-d_0/2)`, `b^2 eta^2=db/v`.

These are nonzero simple initial roots. The correction `Vh/B` tends to zero; the relative order of `A h` above the leading terms is `(m+r)/2>0`. Thus

`P_z~2b^2 eta t^(-d_0/2)`.

On `t=s^2`, the order of `omega` is `d_0+1>0`, so it is holomorphic on the actual branch as well.

If `r>=m`, then `E` has a finite limit. Exclude the single value of `c` making that limit zero. The two escaping roots have nonzero constant limits satisfying `b^2h^2+E_0=0`, and `P_z` has nonzero constant limit. The time form is again regular.

This includes `r=0`, `r=m`, multiplicity two, and arbitrarily high boundary multiplicities. There is no missing local-order case.

## 4. Infinity: the general balanced-degree ledger

Let `n=deg B`, and let the leading coefficients of `V,A,B` be `v,d,b`.

### 4.1 Low clutch degrees `0<=n<=a`

All four roots have

`z~eta x^(1/4-a)`, `v^2 eta^4+1=0`.

The term `Vz^2` dominates `Bz`; its square balances `x`; `Az` grows only as `x^(1/4)`. The four initial roots are simple, and

`P_z~4v^2 eta^3 x^(a+3/4)`.

On `x=s^(-4)`, the order of `omega` is `4a-2>=2`.

### 4.2 High clutch degrees `n>=a+1`: two small roots

Two roots satisfy

`z~eta x^(1/2-n)`, `b^2 eta^2+1=0`,

`P_z~2b^2 eta x^(n+1/2)`.

The relative gaps are positive because `n-a>=1`. On `x=s^(-2)`, the order of `omega` is `2n-2>=2`.

### 4.3 High clutch degrees: two roots near `-B/V`

For the shift in (1), the growth degree `e` of `E` is

- `e=n-a` for `n>=a+2`;
- `e=1` for `n=a+1` and `db/v != 1`;
- `e=0` for `n=a+1` and `db/v=1`, provided `c != lim_(x->infinity)(x-AB/V)`.

In the last case the limit exists; the leading degree-one terms have canceled. This explicitly includes the only infinity resonance.

The escaping roots satisfy

`h~eta x^(e/2-n)`, `b^2 eta^2+epsilon=0`,

where `epsilon` is the nonzero leading coefficient of `E`. The relevant strict inequalities are

`2a+e/2-2n<0`, `a-n-e/2<0`.

They show that `Vh/B` tends to zero and `A h` is lower order. Therefore

`P_z~2b^2 eta x^(n+e/2)`.

On `x=s^(-2)`, the time form has order `2n+e-3>=1`. This guarantees holomorphicity. At the smallest resonant edge `a=1,n=2,e=0`, the two actual escaping punctures are unramified and the order is zero at each; the form is regular there, not necessarily vanishing. The ramified-cover order is one. No strict-positive-order claim on those actual punctures is needed.

## 5. Completeness and choice of the fibre value

There are finitely many excluded values: the critical values; the finite-boundary quadratic-discriminant values; the finite limits in Section 3 when `r>=m`; and the one resonant infinity value if applicable. Critical values are finite because `P_B` is constant on every positive-dimensional component of its critical locus in characteristic zero. Choose `c` outside this set.

At every local place we have constructed four distinct Puiseux determinations, each arising from a simple initial root. Newton–Puiseux/Hensel lifting supplies full branches. The fibre equation is quartic in `z`, so these are exhaustive. Several determinations can represent one ramified normalization point; the proof uses orders on a ramified cover and does not mistake them for separate punctures. A pole of a differential remains a pole under a finite ramified cover, so regularity on the stated covers suffices.

There are no vertical fibre components: `gcd(B,V)=1` ensures that at any fixed `x`, the fibre equation remains a nonconstant quadratic or quartic in `z`. Thus the local `x`-branch ledger accounts for every component.

All finite punctures and all infinity punctures have now been checked. Section 2 proves Theorem 1.1, and the explicit critical point proves Corollary 1.2.

## 6. Scope and reproducibility

This is a no-Jacobian-mate theorem for a specified polynomial architecture. It does not prove that every boundary-unit clutch has a finite critical point. The stronger fixed-accessory critical-point hypothesis H-J remains open. It makes no general planar-JC or literature-priority claim.

The exact companion checks the symbolic domination gaps on all infinite degree and boundary-order sectors, including the smallest resonant edge. It is a certificate for the displayed arithmetic inequalities; the four-root proof above supplies mathematical completeness. Normal and optimized runs must match byte-for-byte.
