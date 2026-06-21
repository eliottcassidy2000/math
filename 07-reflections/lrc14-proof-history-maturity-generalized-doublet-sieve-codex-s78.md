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

The current mature address looks like

```text
(bounded base B, far pair {M,M+g}, bridge parity, full-set gcd,
 remove-one scale witnesses, frozen slow coordinate y, finite bridge window).
```

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

1. Prove cap-dangerous genuine-wide rows have exactly two far elements, or
   prove all `r>=3` rows have a survival/live-depth margin.
2. Prove arbitrary-gap generalized doublets by a frozen-room plus P/R-tail
   theorem, with full-set primitivity and remove-one witnesses in scope.
3. Build a finite bridge atlas for the low-`M,g` rows, especially even-AP plus
   odd-bridge families.
4. Compare the resulting finite addressed families with the prime-sieve
   improper survivors for `k=13`.

The strategic correction is sharp: retain addresses until the last possible
step.  Every time this project scalarized too early, it either proved a false
stronger claim or hid the finite family carrying the actual obstruction.
