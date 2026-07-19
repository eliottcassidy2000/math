# CORRECTED: n=12 equality rigidity is one input, not the compact or LRC(14) residual itself

*boxeph-2026-07-18-S114, corrected by codex-2026-07-18-S74 after THM-1099,
THM-1143, and THM-1149.  Historical filename retained for links.  No global
LRC(14), compact-floor, crown-collapse, or twelve-speed equality proof is
claimed.*

## The conflation

The original reflection identified three statements:

```text
(R12)  |C|=12 and M(C)=1/13  =>  C=d{1,...,12};

(CORE) primitive Cover14 V with M(V)<1/13
       => some deletion is d{1,...,12};

(LRC14) every thirteen-speed family has M(V)>=1/14.
```

They are not known to be equivalent.

- `(R12)` classifies a deletion **if that deletion is already known tight**.
- `(CORE)` additionally needs tight-deletion extraction.  THM-1149 proves an
  exact alternative: every deletion may instead be loose, producing an
  all-loose essential crown.
- `(LRC14)` has threshold `1/14`; the compact `1/13` floor is a stronger
  sufficient residual.  THM-1099 proves only the forward implication chain.

No reverse embedding of an arbitrary tight twelve-set into a primitive
compact Cover14 strict row was provided either.  Thus the S114 equivalence
claim was unsupported in both directions.

## The corrected composition

The proof route now has four separately typed arrows:

```text
primitive + Cover14 + rho<13 + M(V)<1/13
  => some deletion C has M(C)=1/13                  (A)
  => C=d[12]                                        (B)
  => 13d divides the deleted speed                  (C)
  => contradiction                                 (D).
```

- `(A)` is **crown collapse**, open.
- `(B)` is the repo's open twelve-speed equality/rigidity probe, HYP-4382.
- `(C)` is THM-1149's proved Farey-toothpick regeneration blocker.
- `(D)` is proved: primitivity forces `d=1`; the only possible 14-carrier is
  the extra speed, hence it is divisible by 182 and has
  `rho>=182/12=91/6>13`.

The striking point is that `(C)+(D)` needs only a tight classified deletion,
primitivity, one 14-carrier, and `rho<91/6`.  The exact remaining difficulty
is upstream, in `(A)` and `(B)`.

## Why the maximum-deletion formulation is too narrow

The original statement fixed `C=V\{v_max}`.  Essential crowns have thirteen
owners, and no theorem says the maximum deletion must be the tight one.
Even flexible deletion does not solve the issue: the exact row

```text
{1,2,3,5,7,8,9,10,11,12,17,19,104}
```

is primitive, has `rho<13` and `M=8/105<1/13`, and all thirteen deletions
are loose.  It covers `2,...,13` but misses `14`, so it is not an `INVcov`
counterexample.  It proves that Cover14 must enter the extraction step, not
merely be checked after a local residue chart has been classified.

The lift `7 -> 112=7+105` supplies Cover14 while preserving that old chart;
the global maximum jumps to `3/20`.  This is why the integer-lift sidecar in
THM-1099 is load-bearing.

## What remains valuable from S114

The twelve-speed equality problem is height-unbounded and remains a serious
input.  THM-1143 makes its shallow object precise: a labelled mechanical
`A_12` root word acting on a finite mass-eleven simplex.  THM-1150 proves the
all-height local owner-stalk lowering law and isolates PEHD13.  Neither result
closes the deep `s>=2` branch.

The older diagnosis that soft scalar tools do not force AP structure also
survives, but it must be aimed at the correct arrow.  For `(B)`, one seeks a
global equality characterization of a tight twelve-speed object.  For `(A)`,
one seeks a third-order, lift-sensitive crown-collapse theorem.  They are
different obligations and may require different invariants.

## Fano/chi_7 and tournament reframing

The exact crown balance is

```text
mu(N=1)=sum_(k>=3)(k-2)mu(N=k).
```

This says the obstruction first appears at triple incidence.  A Fano or
`chi_7` quotient should therefore use triple-overlap obligations,
lift-compatible chambers, or root-labelled cuts as vertices.  Bare runners,
private-stalk midpoint order, and pair-overlap totals all forget the
co-location which pays for private mass.

The useful tournament remains telemetry: orient proof obligations by first
wall or owner cost, record score histograms, SCCs, edge flips, and the tie
Hamiltonian path, but retain the hyperedge sidecar.  On the exact all-loose
row the private-stalk tournament is transitive with 94 singleton SCCs and one
Hamiltonian path; that simplicity is evidence that pair orientation has
forgotten the hard object.

## Honest frontier

The sharpened frontier is not one allegedly definitive equivalence.  It is a
two-input composition:

1. close shallow and deep twelve-speed equality rigidity;
2. prove Cover14 crown collapse, or directly contradict the all-loose crown
   using the edge/curtain word and shared lift congruences.

THM-1149 finishes the regenerated AP branch once those inputs meet.  LRC(14)
remains open.

Cross-links: THM-1099, THM-1143, THM-1149, THM-1150, HYP-4382,
HYP-7665, HYP-7675, and MISTAKE-170.
