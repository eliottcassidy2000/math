# A degree-one derivative gcd does not see hidden DVR thickness

**Status: accepted hostile audit of
[THM-3538](../01-canon/theorems/THM-3538-fixed-keller-newest-prime-prescribed-coordinate-index-criterion.md).**
The theorem remains proved.  This reflection records the sidecar that makes
its equivalence `(16)` valid and prevents exporting it as a general
special-fibre maximality test.

## The portable implication that fails

Let `R=k[[pi]]` with `char(k)!=2`, and let `t^2=pi`.  The integral closure of
`R` in the tame quadratic extension is `k[[t]]`.  In

```text
B=k[[t]] x R,                  theta=(pi t,1),
```

the observation has polynomial

```text
m_theta(T)=(T^2-pi^3)(T-1).
```

Its reduction is `T^2(T-1)`, so

```text
gcd(mbar_theta,mbar_theta')=T.
```

The visible special fibre is exactly one doubled shadow.  Nevertheless,
`R[pi t]` has basis `(1,pi t)` inside the normalized basis `(1,t)`, hence
`length_R(B/R[theta])=1`; the cross-block resultant `1-pi^3` is a unit.  The
extra index is the invisible square factor in
`Disc(T^2-pi^3)=4pi^3`, not an additional reduced root collision.  Thus a
degree-one derivative gcd does not imply local maximality for a general DVR
order: the first failed implication is `(d)=>(b)`, and consequently
`(d)=>(a)`.

## Why the fixed Keller section is different

THM-3538 does not apply `(16)` to an arbitrary doubled block.  Degree-one
ancestry identifies its distinguished section `q_0` with the generic point
of the irreducible divisor `V(L)`.  At the exact point

```text
q_*=(2/27,1,1) in V(L),
(h_y(q_*),h_z(q_*),h_u(q_*))=(1,49/32,1/2).
```

Therefore none of the three square-factor functions vanishes identically on
`V(L)`, so each `h_theta(q_0)` is a unit at the generic newest-prime DVR.
The hostile mechanism `h_theta(q_0) in (pi)` is unavailable.  Once that
internal thickness is excluded, every remaining nonunit internal factor or
cross-block resultant creates an additional or thicker derivative-gcd
factor, which is exactly the implication used in `(16)`.

## Boundary retained

| item | exact conclusion |
|---|---|
| observable | derivative gcd of the reduced total coordinate polynomial |
| information it loses | positive valuation of the square factor inside an already doubled block |
| required Keller sidecar | `q_0` is the generic `V(L)` section and `h_theta(q_0)` is a unit |
| hostile control | `(T^2-pi^3)(T-1)` |
| theorem consequence | equivalence `(16)` remains valid for the fixed newest-prime model |
| non-consequence | no carrier-unit or index-zero assertion for any `n>=5` |

The exact positive rows remain only `n=2,3,4`.  This audit neither weakens
THM-3538 nor supplies an all-level prescribed-coordinate maximality theorem.
