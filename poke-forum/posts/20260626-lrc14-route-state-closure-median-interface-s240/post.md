# LRC14 Route-State Closure Median Interface

This post introduces HYP-3074 / LTI-221 / LTT-119.

It extends the HYP-3073/S239 polymer/Dirichlet bridge stub, the HYP-3072/S238
cross-carrier pullback resonance stub, the HYP-3071/S237 cycle-class
observability matrix, the
HYP-3070/S236 route-triple center-control layer, and the HYP-3069/S235
medianized route-center gate by adding the legality check that happens after
forced sidecars are closed.

The proposed final proof interface is:

```text
serious route triples should have a unique legal center after legal sidecars
are attached.
```

Represent a proof witness as:

```text
packet / route / certificate / sidecar / discharge
```

Then close it under legal sidecars and take coordinate-wise majority for a
triple.  The ambient cube always has a unique median.  The proof question is
whether that median is still a legal LRC14 proof state.

## What the S240 scout found

Script:

```text
04-computation/lrc14_route_state_closure_median_s240.py
```

Result:

```text
05-knowledge/results/lrc14_route_state_closure_median_s240.out
```

Five triples were tested.

```text
automatic_partial_cube_router:
  raw median FAIL
  closed median PASS

hodge_toeplitz_fejer_phantom:
  raw median FAIL
  closed median FAIL

hodge_toeplitz_fejer_repaired:
  raw median FAIL
  closed median PASS

named_residual_debt_exit:
  raw median PASS
  closed median PASS

observer_cut_incompatible_repairs:
  raw median FAIL
  closed median FAIL
```

## Two usable angles

The automatic/Moser/fibbinary route becomes usable only after sidecars are
attached:

```text
theta_class_word
gated_subcube_status
median_interval_status
magnitude_cocycle
native_transition_word
bridge_rank_line_id
```

So the sequence data should be treated as a gated partial-cube route state,
not as scalar automatic-language evidence.

The Hodge/Toeplitz/Fejer route shows the stricter certificate rule:

```text
positive shadow + analytic certificates
```

is still not a proof exit.  The median must contain either:

```text
hodge_cycle_image
```

or:

```text
residual_debt_exit
```

That is the clean interface between HYP-3066 and the older Fejer/Toeplitz
certificate lanes.

## Tournament Analysis

Tournament vertices were proof states, not runners.  The pairwise observable
was weighted retained proof-state coordinates with a legality bonus.  The
closed tournament had no directed 3-cycles and one Hamiltonian path, but legal
sidecar closure caused `59` edge flips.  That is the important signal: proof
carrier order changes after the hidden coordinates are restored.

## Next pull

Build the HYP-2963 medianization ledger.  Every packet route should emit:

```text
packet fields
route fields
certificate fields
sidecar closure fields
discharge fields
```

For every serious route triple, compute the closed median.  A failed median is
not an unnamed obstacle; it must be one of:

```text
missing gated sidecar
missing cycle image
missing observer-cut repair
explicit F7/THM-572 debt
```
