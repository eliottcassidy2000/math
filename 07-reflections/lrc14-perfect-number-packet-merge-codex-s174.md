# LRC14 Perfect-Number Packet Merge - S174

The prior perfect-number work is useful here only if it stays typed.  The
even-perfect chain is the exact `n=2` unit-excess control:

```text
a=2^(r-1), q=2a-1 prime, sigma(aq)/(aq)=2
```

The LRC14 analogue uses the same power-of-two spine but changes the denominator
rule to `q=14a-1`.  When `q` is prime, it is uniformly deficient:

```text
sigma(aq)/(aq)=14(2a-1)/(14a-1)
2 - sigma(aq)/(aq)=12/(14a-1)
```

The S174 computation makes the guardrail concrete.  Prime `q14` rows in the
bounded scan are `(a,q)=(1,13),(16,223),(256,3583)`, and each has the predicted
defect.  Composite `q14` rows flip the behavior: every composite `q14` row in
the scan is abundant, starting with `a=2,q=27=3^3` and defect `-2/9`.

That means the proof object cannot retain only the product `a*q`, the power
address `a=2^k`, or automatic membership.  It must keep at least:

```text
exact M/Farey address
unit-excess apex n
prime/composite q flag
factorization
abundancy defect
Kpq/product-incidence route
automaton transition state
```

Assumption challenge: tournament vertices were not chosen as runners.  I
considered runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,
divisor atoms, and proof obligations.  The chosen vertices are proof carriers
and packet side channels because the preserved predicate is a labelled LRC14
packet verdict.  The quotient destroys runner identity and scalar product
comparisons on purpose.

The practical next pull is to add the S174 packet fields to HYP-2963 sidecars
and rerun route-purity on `q=14p-1` rows with the HYP-3012 automaton labels
attached.  The first shared warning is `p=2,q=27=3^3`: it appears in the
Fermat-Catalan / perfect-power ledger and in the composite-q abundancy flip.
