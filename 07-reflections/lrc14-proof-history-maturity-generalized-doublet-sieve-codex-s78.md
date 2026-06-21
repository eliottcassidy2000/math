# LRC14 proof history matured into addressed residual families

codex-2026-06-22-S78, HYP-2810.

The repo's LRC14 work has been through several plausible centers of gravity.
The consecutive/AP maximizer route was useful but too rigid: it was false at
the k=12 cap layer and unnecessary for the direct cap inequality.  The uniform
decorrelation-error route was also too scalar: large error and small room live
in different row families, so `decorr_sup + err_sup` is the wrong currency.
The next durable layer was the regime split: bounded rows, single-far rows, and
genuine-wide rows.  THM-563 made the single-far branch exact by signed endpoint
periodicity, which forced the remaining proof obligation into genuine-wide
rows with at least two far elements.

Post-fetch HYP-2811/HYP-2812 changes the proof order again.  The clean route is
now the gK8 missed-sector moment

```text
L_yK8 = 10q0 + q3 + 10q6,
```

not a direct `p0` extremality statement.  The new concentration-extremality
claim is that `L_yK8` is maximized by bounded/concentrated rows and strictly
drops under wide decorrelation.  This explains why the k=12 odd-bridge breaker
can beat `Q(11)` in `p0` while remaining comfortably below the bounded maximum
in `L_yK8`: gK8 charges the extreme missed-count atoms `q0` and `q6`, and
spreading pushes mass into the middle.  The exact S78 gK8 overlay on the
span-18 genuine-wide domain found zero `L_yK8` violations, with the same r=2
leaders at k=10,11,12.

The current mature address looks like

```text
(missed-sector distribution q, gK8 moment, bounded base B,
 far pair {M,M+g}, bridge parity, full-set gcd, remove-one scale witnesses,
 frozen slow coordinate y, finite bridge window).
```

After rebasing over HYP-2814/HYP-2815/HYP-2816, the concentration mechanism is
sharper: the wide move should be viewed first as `q6` suppression.  On the
all-missed atom, a far speed is almost an independent clock, so it asks for the
new speed to land in the one surviving sector, costing roughly a `1/7` factor
per far coordinate.  At small boundary speeds this ratio is not literally
`1/7`, but the new work says it remains a strict contraction.  The proof
target is therefore a generated-profile majorization lemma: the admissible
gain in `q0` from spreading cannot compensate for the gK8-weighted loss in
`q6`, with endpoint-period/R-tail machinery handling the low-`f` boundary.

This address explains why mac-mini S7's k=12 over-`Q(11)` obstruction did not
break the cap route.  It is not a new unstructured regime; it is a generalized
doublet with an even-AP base and two odd bridges.  S78's exact span audit
supports this: through `max(E)<=20`, the k=12 global maximizer and all seven
over-`Q` top rows remain `r=2`, while the first `r=3` row is already below
`Q(11)`.

The 2026 public below-frontier proof by Sungkawichai-Trakulthongchai gives a
second frame for the same lesson.  Their lifting/backward-projection sieve
works by shrinking the improper ansatz modulo primes after factoring obvious
symmetries.  For the fourteen-runner frontier, the right connection may be:

```text
public improper ansatz survivor
  -> our addressed generalized-doublet / bridge profile
  -> frozen-room plus P/R tail
  -> finite cap atlas.
```

So the bold finish is not "make Weyl error uniformly small."  It is:

1. Prove gK8 concentration extremality: wide decorrelation cannot increase
   `10q0+q3+10q6` beyond the bounded maximum, ideally by the new `q6`
   contraction/Krawtchouk majorization route.
2. Use the bounded gK8 certificate to compare that maximum with `10cap_k`.
3. If the concentration proof is only asymptotic, use the HYP-2811/HYP-2808
   R-tail route; the Mordell-Tornheim constant now has closed form
   `12*zeta(3)`.
4. Keep the generalized-doublet and finite bridge atlas as the resonant
   exception layer, especially for even-AP plus odd-bridge families.
5. Compare any remaining finite addressed families with the prime-sieve
   improper survivors for `k=13`.

The strategic correction is sharp: retain addresses until the last possible
step.  Every time this project scalarized too early, it either proved a false
stronger claim or hid the finite family carrying the actual obstruction.
