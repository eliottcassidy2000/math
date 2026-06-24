# LRC14 Kernel Homotopy Boundary Defects

The new proof pass is aimed at a small but load-bearing gap between two live
routes.  HYP-2981 wants interval-certified Fejer packets.  HYP-2982/HYP-2983
want smoothing, large-sieve, and Kaczynski-style boundary language.  The bridge
is not "choose a better kernel."  The bridge is:

```text
kernel change = labelled homotopy or named defect.
```

If a row has a strict safe component `(a,b)`, a symmetric kernel supported
inside radius `<(b-a)/2` can sit at the midpoint without crossing a danger
boundary.  That support radius is a certificate field.  If a row is zero-open,
the support radius is zero and any positive smoothing crosses the boundary.
Then the deformation has emitted a boundary defect, not a scalar estimate.

## Exact Readout

The S164 script reuses the exact Haar/Baire component engine and the taut bridge
audit.  The named packets separate cleanly:

```text
AP, GW 12->24: safe_mu=0, eps<0, six taut endpoint transfers.
near/K33 12->36: safe_mu=1/1260, eps<1/5040.
P10+GW: safe_mu=1/980, eps<1/3920.
covering 12->168: safe_mu=263/30030, eps<23/23520.
```

The AP/GW taut transfers are exactly the six denominator-14 unit points, with
owner pair sums `0 mod 14`.  That is the defect payload.  It is not optional
metadata.

The one-swap sanity scan through replacement `<=160` is a useful calibration:
`1910` primitive rows are open-stable and exactly one is zero-open, the
Goddyn-Wong move `12->24`.  The first positive escape is `12->36`, with
`safe_mu=1/1260`.

## Proof Use

The practical rule for future analytic smoothing work is:

```text
declare the packet,
declare the kernel family,
declare the support or interval certificate,
declare the boundary defects created by the homotopy.
```

Without those fields, smoothing is just another scalar quotient.  With them, it
can become a controlled bridge between Fejer certificates, Ramanujan exact
period packets, Kaczynski approach classes, and state-lift residuals.

Tournament Analysis supports the ordering.  Proof carriers beat raw kernels:

```text
packet Fejer certificate
> open component certificate
> boundary defect atom
> Kaczynski approach class
> Ramanujan exact-period packet
> analytic smoothing kernel
> raw kernel scalar.
```

The next theorem target is familywise: lift the support-radius and defect
ledger from named packets to HYP-2963 packet families.  A smoothing proof should
then say that every deformation either preserves an open-component/Fejer
certificate or routes the emitted defect to AP/GW, K33/state-lift,
Ramanujan/exact-period, or an existing interval certificate.
