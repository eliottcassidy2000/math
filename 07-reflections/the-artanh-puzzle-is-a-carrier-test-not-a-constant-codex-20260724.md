---
source: codex-2026-07-24-artanh-multifrontier-session
status: >
  SYNTHESIS WITH PROVED/FINITE-EXACT ANCHORS AND SCOPED NO-GO RESULTS.
  THM-2143 is the proved LRC Gibbs/artanh carrier. THM-2146 is the
  finite-exact defect-sensitive union base. The GMC and planar-Jacobian
  sections diagnose why the same numerical margin does not transfer; they
  do not claim new closure theorems.
tags: [artanh, lrc14, gibbs, overlap-pgf, gmc2, dvdk, jacobian, square-face, no-go]
related: [THM-2143, THM-2146, THM-2111, THM-2127, THM-2129, THM-2132, THM-2136]
---

# The artanh puzzle is a carrier test, not a magic constant

The supplied fragment proves, exactly,

```text
(2457/6592) log(8847357/2974400)-log(1285/896)>1/25. (1)
```

Its durable instruction is not “look for `1/25` everywhere.”  It is:

> first produce a target-preserving positive ratio, then use the Cayley
> coordinate and the odd artanh series to certify its sign without floats.

This distinction separates one real LRC transfer from two tempting but
invalid cross-frontier transfers.

## 1. LRC: the missing carrier is the danger-depth PGF

Let `H(t)` count dangerous runners and let

```text
Z(z)=integral z^H,                  0<z<1.             (2)
```

THM-2143 proves

```text
measure{H=0}>0
 iff there exists rational z in (0,1) with Z(z)>z.     (3)
```

The forward implication uses `z->0`; the reverse is the pointwise cover law
`H>=1 => z^H<=z`.  Once (3) produces two positive rationals,

```text
log(Z/z)=2 artanh((Z-z)/(Z+z)).                        (4)
```

Thus the odd-power sandwich in the puzzle is a legal exact wrapper.  It does
not supply (3), and it forgets isolated height-exactly-`h` witnesses.  The
sidecar for the tight boundary remains the labelled phase-height/owner
carrier.

The two displayed external temperatures are not universal LRC fugacities:
exact `h=3/41` endpoint sweeps make both fail on a fixed defect-seven hostile
battery, although `z=1/10` succeeds there.  This is not a refutation of a
derived two-temperature argument.  It is a type check: the row-to-ratio map
must be proved.

The same principle sharpens the defect-seven analytic lane.  Charging six
residual bands separately destroys their overlaps.  THM-2146 computes the
exact small-core union mask and leaves one aggregate covariance against the
seven-speed body.  At defects `7,...,11` its independent base improves the
old `5/41` factor by

```text
153101/69300, 1649/900, 19/12, 427/330, 11/10.        (5)
```

The PGF and the union mask are the same operation at different resolutions:
both retain overlap multiplicity before scalarization.

## 2. GMC(2): Cayley linearizes return order, but scale kills a margin

For the normalized root-product series `Y` in the additive proof of
THM-2111, put

```text
W=(Y-1)/(Y+1).
```

The formal identity has the desired shape:

```text
2 artanh(W)=log Y=sum_(m>=1) CT(f^m)t^m/m.             (6)
```

If `m_*` is the first nonzero return, then

```text
ord_t(W)=m_*,
[t^m_*]W=CT(f^m_*)/(2m_*).                            (7)
```

This is a clean coordinate reformulation of the existing return theorem,
not a new DvdK-free proof.  The hostile operation is scalar dilation:
`f->lambda f` preserves the entire zero/nonzero return pattern while sending
`W_f(t)` to `W_f(lambda t)`.  Every fixed-temperature Archimedean margin can
therefore be driven to zero.

What survives is a normed tail estimate.  If `C=sum|c_q|`,
`r=C|t|<1`, and the first `L` returns vanish, then

```text
|sum_(m>=1) CT(f^m)t^m/m|
 <=r^(L+1)/((L+1)(1-r)).                              (8)
```

So a log margin can force a low return only after adding coefficient norm
and parameter scale.  THM-2111 already gives the stronger
coefficient-independent compound-degree bound.  The puzzle contributes a
coordinate, not the missing algebraic seed or good-prime sidecar.

## 3. Planar Jacobian: the logarithm collapses to the old time form

On the quartic-square carrier write

```text
P=H^2+R,          S^2=-R,          u=H+S.
```

Then `P=(H+S)(H-S)` and

```text
{P,u}=2{H,S}u.
```

If a rational mate satisfies `{P,Q}=1`, its differential on the generic
fibre is

```text
dQ=du/{P,u}=-dH/{H,R}=-dx/P_z=dz/P_x.                (9)
```

The branch logarithm

```text
log((H+S)/(H-S))=2 artanh(S/H)                        (10)
```

therefore differentiates to the standard Hamiltonian time form already used
by the THM-2127 residue method.  The involution `S->-S` reverses (10) but
leaves (9) invariant.  A numerical inequality between branch logarithms
cannot descend to the base without extra orientation data.

The infinity calculation explains the existing wall.  In the strict
noncancelling branch

```text
H=z^2+p/2,     R=qz+delta,
d=deg p,       rho=max(deg q+d/2,deg delta)<2d,
```

one has

```text
z~c x^(d/2),   H~h x^(rho/2),
-dx/P_z ~ (N/(4ch)) s^(N((d+rho)/2-1)-1) ds
                                      for x=s^(-N).   (11)
```

Hence:

- `d+rho>2`: the time form is holomorphic at the square-split infinity
  branches, so the residue test is automatically silent in the genuinely
  high-degree lane;
- `d+rho=2`: the remaining integer cases already have an elementary
  critical locus or critical point;
- `d+rho<2`: leading order includes tame controls and cannot decide
  invertibility;
- cancellation in `q sqrt(-p/2)+delta` creates a secondary Newton edge,
  returning exactly to the THM-2132--2136 coarsened-power/Hermite problem.

Thus the missing planar coordinate is global de Rham exactness or
cross-factor Newton compatibility, not a sharper artanh estimate.

## 4. Updated concept board

| lane | source | target map | preserves | loses | next sidecar |
|---|---|---|---|---|---|
| LRC strict branch | danger-depth law | `p -> Z(z)` | positive zero-depth mass | tight isolated points, owners | phase-height equality carrier |
| LRC defect `>=7` | six residual combs | exact union mask | all overlap cancellation | boundary ownership | labelled `6+7` covariance/relation spectrum |
| GMC return | root product | Cayley/log series | first return order | scale-normalized margin | coefficient norm or compound degree |
| JC square face | split factors | logarithmic derivative | Hamiltonian time form | branch sign and global exactness | secondary Newton/de Rham data |

The cross-domain pattern is now precise: artanh is the unique linearizer of a
positive ratio of two linear quantities, but the ratio must already be a
faithful carrier.  Exact numeration cannot repair a missing carrier.
