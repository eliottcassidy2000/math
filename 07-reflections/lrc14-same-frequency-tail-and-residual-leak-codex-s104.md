# LRC same-frequency tail and the residual leak

S104 grew out of the KPS S31k additive-energy coefficient.  The natural
question was whether the proved `m=1` positivity is just the first harmonic of
a convergent packet.  It is.

For fixed residue `r mod 7`, the same-frequency coefficient satisfies
`Gamma_k(m)=C_{k,r}/m^4`.  All six residue constants are positive for
`k=8..13`, and the absolute tail after `H=12` is already around `10^-6`.  So
there is a genuine positive additive-energy Fourier tail.

The surprise is that the tail is too positive.  On the AP rows it predicts
about twice the actual `p0-p0_decorr` deviation.  The missing term is not a
small leftover.  It is a negative labelled correction that grows with
AP-likeness.

This turns the proof target around.  The right object is

```text
R_sf(E) = p0(E) - p0_decorr(k) - Gamma_k^sf A*(E).
```

AP has a very negative `R_sf`.  Moving away from AP leaks some of that
negative correction back.  The exact bounded scans suggest the leak is always
smaller than the same-frequency energy advantage AP loses:

```text
R_sf(E)-R_sf(AP) <= Gamma_k^sf (A*(AP)-A*(E)).
```

For `k=8` the worst ratio in the anchored `max<=14` bank is only `0.469`.  For
`k=9` it jumps to `0.933`, and the row is not arbitrary:

```text
(0,2,4,6,7,8,10,12,14)
```

That is an even AP plus the midpoint bridge.  This is exactly the same
scaling/tiling signal that mac-mini S39 isolated: the tight LRC structures are
scaling-invariant exact tilings, while raw additive energy is
translation-invariant and therefore too coarse.

The tournament analogy is a cut/cycle split rather than a runner tournament.
The positive same-frequency additive-energy packet is the cut-like two-body
carrier.  The residual is the cycle/current correction: hidden folds,
support-cycle terms, and in the repeated-packet branch the octahedral face-curl
from HYP-2887.  A scalar energy proof would forget the labels that make the
negative correction appear.

So the session did not finish LRC14, but it sharpened the open constant.  The
live inequality is now a residual-leak bound, with an explicit worst family to
attack first.

Post-pull KPS S31l makes the same point from the next moment layer.  The
`s=2` same-frequency packet is positive, but `s>=3` moment coefficients are
mixed and often negative.  That means the residual is not an error term to take
in absolute value.  It is the signed convexity content of the problem.  The
right tournament analogy is therefore not term-by-term cycle maximization; it
is the Jensen proof of the H-maximizer, where the whole signed functional is
controlled by the extremal score/difference profile.
