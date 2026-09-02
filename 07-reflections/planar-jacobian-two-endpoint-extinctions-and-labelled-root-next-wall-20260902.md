# Planar Jacobian endpoint extinctions and the labelled-root next wall

**Status:** research synthesis.  The proved content is exactly
[THM-4334](../01-canon/theorems/THM-4334-beta-zero-exact-weight-twelve-endpoint-wall-extinction.md)
and
[THM-4337](../01-canon/theorems/THM-4337-zeta-zero-exact-weight-twelve-endpoint-wall-extinction.md),
both relative to the inherited exact-weight-twelve seam.  This reflection is
not a proof dependency.  `JC(2)`, `DC(2)`, and exact-seam entry remain open.

## Inheritance pass and portfolio

The anchor was the `Z=0` endpoint residual left by THM-4327.  The closest
proved mechanism was its lower-face/good-differential/proper-flat interface.
The canonical hostile was MISTAKE-531: at `Lambda=U+W=0`, twelve ordinary
top-edge nodes coalesce into one `A_23` contact, so a simple-root response
packet cannot be transported.  The corrected near miss was MISTAKE-540:
cyclic face genus has to be computed in the actual invariant function field,
not in a convenient but potentially nonbirational cover.  The least-used
sidecar was `K=[y^2]H`, whose row can create a rational face and whose local
term partly enters the balanced leading coefficient.

The live concept board was

```text
endpoint owner | supporting plane | invariant face field | edge orbit
graph genus | A23 splitter | optional K face | labelled Laurent root
```

The niche was the distinction between an owner disappearing at a toric
endpoint and two surviving edge roots colliding in the torus.  The wildcard
was the owner-address permutation in
[THM-4335](../01-canon/theorems/THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal.md):
it did not prove a Jacobian claim, but it supplied the correct labelled-root
organization for the next wall.  The concurrent THM-4336 rank-four appender
result did not preserve a comparably precise predicate here and was not used.

## Anchor outcome

THM-4334 closes

```text
Z=beta_11=0,                 U*W*zeta_3!=0,
```

including `Lambda=0`.  Its new face is a connected cyclic-nine genus-three
curve, the central face has genus three, graph `b1=11`, and the two positive
orders are `34,27`.  At the `A_23` contact, `zeta_3!=0` forces the shortened
critical ladder to split through `C_2` or `C_3` before the normalized `t^6`
terms.

THM-4337 closes the complementary large zeta-owner stratum

```text
Z=zeta_3=0,                  U*beta_11!=0,
W,K,Lambda arbitrary.
```

Its exact four-stratum atlas is:

| Regime | Faces | Positive genera | Graph `b1` | Orders |
|---|---|---:|---:|---|
| `W!=0,K=0` | `M,V7` | `3+3` | `11` | `63,70` |
| `W!=0,K!=0` | `M,V7,Kf` | `3+3+0` | `11` | `63,70,196` |
| `W=0,K=0` | `M,V9` | `6` | `11` | `27,28` |
| `W=0,K!=0` | `M,V9,Kf` | `6+0` | `11` | `27,28,84` |

Every row sums to arithmetic genus seventeen.  A broad inherited-deletion
over-atlas and an import-free 293-plane reconstruction recover the same four
complexes.  The latter also rebuilds every outer and internal edge from
coefficient-labelled exponent dictionaries.

The actual face fields matter:

```text
C:   y^2=W*x*(x^6-U),                                  genus 3;
V7:  P^7=(beta^2/W^3)*x^3/(1-x)^2,                     genus 3;
V9:  S^3=(1-U*P^6)/(beta*P^4),                         genus 6;
Kf:  P=(1-KV^2)/(beta V^3), V=SP,                      genus 0.
```

Their complete edge schemes show that only the `W!=0` top edge can repeat,
and its discriminant is `Lambda^2`.  At `Lambda=0`, put `s=v(sigma)>0` and
`nu=v(b)>s`.  Deep repetition gives

```text
c_1=alpha_11+beta_11=0,
C_1=alpha_11/c^2=-beta_11/c^2!=0.
```

There are two forced splitters:

```text
d_b=6(nu-s),                       d_1=s+nu.
```

The first is `b^12q` for `1<nu/s<7/5`, they tie at `7/5`, and `C_1t` is
first above `7/5`.  Their Keller orders `6s+2nu` and `(5s+9nu)/2` are
positive.  At equality the face `X(C_1-X^5/c)` has derivative `-5C_1` at
every nonzero root.  Both gaps precede every normalized `t^6` term.

The `K` timing was a useful firewall.  Its local contribution is

```text
K*t^6*r^2  ->  K*y*r^2=K*y+2K*t^6*y^2+K*t^12*y^3.
```

Thus `Ky` is absorbed into the balanced coefficient `c`, while the remainder
starts at `t^6`; neither can erase the earlier splitter.

The exact `Z=0,U!=0` residual is now

```text
beta_11=0 and W*zeta_3=0.
```

This is an endpoint reduction, not seam entry and not `JC(2)`.

## Niche signal: owner addresses become labelled Laurent roots

The incoming LRC work separates a tail label from the sheet address it owns.
The useful Jacobian analogue appears after `beta_11=0,W!=0`, where the new
edge polynomial is

```text
A(P)=K+Theta*P+xi_10*P^2+W*P^3.
```

The connection has the following typed form.

- **Source:** LRC tail labels with componentwise sheet-owner addresses.
- **Target:** roots of `A` in the labelled three-point finite etale algebra
  over the coefficient stratum.
- **Map:** tail-to-address incidence becomes coefficient/root incidence.
- **Preserved predicate:** distinct active owners correspond to distinct
  nonzero torus roots after Laurent saturation.
- **Destroyed information:** the scalar discriminant forgets which root
  moves to `0` or infinity and the lower hull forgets root collision order.
- **Needed sidecar:** labelled root idempotents, Laurent saturation by
  `K*W`, and local intersection multiplicities.
- **Cheapest decisive tests:** one simple cubic, one interior double-root
  cubic, and one root-exit cubic.

The symmetric collision invariant is

```text
Disc(A)=xi_10^2*Theta^2-4W*Theta^3-4xi_10^3*K
        -27W^2*K^2+18W*xi_10*Theta*K.
```

It organizes three different geometric events:

```text
K*W!=0, Disc(A)!=0:  three distinct torus attachments;
K*W!=0, Disc(A)=0:   an interior contact;
K=0 or W=0:          an owner exits through a toric endpoint.
```

These cannot be merged into the single slogan “the discriminant vanishes.”
The hostile controls are

```text
P^3+P^2+P+2,                 Disc=-83, simple;
(P-1)^2(P+1),                interior double root;
P(P^2+P+1),                  root exit after Laurent saturation.
```

This is the sharp next object: a labelled Laurent-root stratification is more
informative than immediately enumerating a larger unlabelled plane atlas.

## Wildcard and generated tasks

The other underused object is the critical-value/conductor sidecar at the
`U=0` repeated `(2,5)` and `(2,3)` cusps.  Supporting planes and positive
primary orders do not control those contacts by themselves.  The next board
is therefore:

1. **Anchor:** over `beta_11=0,W!=0`, construct the Laurent-saturated cubic
   root algebra, stratify `Disc(A)`, and compute the local tail at a generic
   interior double root.  Determine whether the first nonzero normal Hasse
   coefficient always splits it with positive good-form order.
2. **Niche:** over `beta_11=W=0`, rebuild the honest lower atlas before
   specializing edges.  Separate the zeta-nonzero quadratic wall from the
   beta/zeta double-owner intersection; use
   `U+alpha_11 X+xi_10 X^2=(X+1)^2` as the first hostile.
3. **Wildcard:** at each `U=0` repeated cusp, recenter at the unique critical
   section and compute the invariant critical-value series, conductor, and
   Kahler/dualizing discrepancy.  Stop if the proposed coefficient is not
   invariant under the remaining source-coordinate gauge.
4. **Entry lane:** search for a valuation or degree mechanism forcing a
   hypothetical counterexample into one of the extinct exact-seam strata.
   Endpoint extinction without this arrival statement cannot settle `JC(2)`.

For every lane the cheapest negative test comes first: a repeated nonzero
edge root, a root at a toric endpoint, a vanished proposed splitter, or a
coordinate change that alters it.  A positive signal is promoted only after
the face field, edge orbit, lost information, and proper-flat consumer have
all been identified.
