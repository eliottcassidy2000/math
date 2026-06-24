# LRC14 admissible smoothing dispatcher

The useful move in this pass was to stop asking for a universal analytic
estimate and instead ask which smoothing clock is allowed to inspect which
packet.

HYP-2982 says the large-sieve/Selberg normalizer `sum mu^2/phi` is real but
squarefree-blind.  HYP-2983 says Goldbach-style proof architecture has the
right shape only if smoothing, exceptional sets, and exponential sums remain
labelled.  HYP-2985 fuses those into a dispatcher:

```text
AP/GW boundary        -> Kaczynski/endpoint boundary clock
K33 / covering        -> Fejer interval clock
q=27 petals/splices   -> Fejer plus Ramanujan prime-power clock
late denominator wall -> large-sieve normalizer with exact-period side channel
true-wide generic     -> Kaczynski/Abel off-resonance clock
true-wide resonant    -> Freiman finite atlas or state lift
```

The four clocks worth keeping separate are endpoint-owner time, exact-period
denominator time, smoothing/certificate time, and far-approach boundary time.
This is a better answer to "does the orbit hit the safe box even once" than a
single global reset clock: the global reset is often infinite or too large,
but these four clocks are local proof coordinates.

The creative proof target is now:

```text
admissible-smoothing lemma:
  every HYP-2963 primitive packet has an allowed local smoothing policy;
  if no local policy applies, the obstruction labels produce a THM-572
  state-lift obligation.
```

This should make the Fejer, Ramanujan, large-sieve, and Kaczynski agents stop
competing for one scalar endpoint.  They are different clocks in one proof
sheaf.  The proof should glue their certificates, not average them.

The wild edge is the "explicit explicit formula" version of LRC14:

```text
safe/noncover packet =
  endpoint boundary main term
  + Ramanujan exact-period correction
  + Fejer/Toeplitz dual certificate
  + Kaczynski approach defect
  + finite state-lift residual.
```

Each summand is only meaningful with its address coordinate attached.  That is
the common lesson across the recent divisor, smoothing, and boundary passes:
quotients are allowed to forget only after the forgotten thing has a route home.
