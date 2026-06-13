---
source: codex-2026-06-01-S548
status: formalization session; theorem extracted from S546/S546b
tags: [lonely-runner, formalization, endpoint-core, p-adic, zero-branch, theorem, entropy]
---

# Formalizing the zero-branch star: why one p-adic branch cannot be the core

**Prompt (user):** spend a long session formalizing recent novel work, and let
anything encountered spur mathematical investigation and discovery.

The recent work around HYP-2036, HYP-2037, and HYP-2038 has two voices:

1. p-adic/product-tree branches, which say where rational sieve witnesses live;
2. entropy/order-parameter probes, which say global spread and criticality are
   the right large-scale coordinates.

The formal question is where these meet.  S546b had computed that covered
prime-power zero branches kill unit witnesses but have no local endpoint core.
S548 turns that into THM-391.

## The Formal Object

Fix `n >= 2`, `2 <= q <= n`, a set of nonzero q-grid centers

```text
C subset {1/q, ..., (q-1)/q},
```

and speeds `S` all divisible by `q`.  The zero-branch star is the interval
family

```text
I(c,s) = (c - 1/(n s), c + 1/(n s)).
```

In the LRC prime-power case, `q=p^d`, `C` is the unit set `(u,q)=1`, and these
are exactly the local danger intervals around the THM-369 unit witnesses `u/q`
after the zero branch is covered.

## Theorem Extracted

THM-391:

```text
Every q-grid zero-branch star has empty strict endpoint-protection core.
Peel layers are explicit: speeds peel in increasing speed order, with
|C| * multiplicity(s) intervals removed at the layer for speed s.
```

The proof is short.  Each radius is at most `1/(nq) <= 1/(2q)`.  Different
q-grid centers are separated by at least `1/q`, so intervals from different
centers do not strictly protect one another's endpoints.  At a fixed center the
intervals are nested; a largest active interval has no larger interval to
protect its endpoints.  Therefore any nonempty proposed core contains an
interval with an unprotected endpoint.

The discovery is the q-agnostic form: primality is unused.  Prime powers matter
because they label p-adic tree branches; the local star geometry is universal.

## Verification

`lrc_zero_branch_star_theorem_s548.py` checks the theorem and layer formula over
bounded exact q-grid stars:

```text
bounded exact cases checked: 3255
max intervals in bounded check: 21
all bounded cores empty and all peel layers match the formula
```

Selected examples include `n=18,q=9`, `n=18,q=18`, `n=14,q=7`, and a non-prime
`n=12,q=6` star.  All peel exactly by the formula.

## What This Clarifies

HYP-2036's product-tree scan remains important, but now its local negative part
is proved.  A branch with `z_q>0` kills the obvious rational witness.  It does
not create a local obstruction core.

HYP-2037/HYP-2038 now have a cleaner job description.  Entropy and
box-dimension see global spread/criticality.  THM-391 says that a proof cannot
collapse that global signal onto a bare zero branch.  It must keep at least one
export label:

- descendant endpoint,
- event owner,
- compactified critical wall,
- cross-prime product-tree coordinate,
- or Gabor zero-column harmonic label.

So the next nontrivial tournament/trienerment should not have vertices `q`.
It should have vertices such as:

```text
(q, endpoint descendant),
(q, owner speed, endpoint sign),
(p^a, r^b, coupled wall cell),
(q, unit u, Gabor harmonic m).
```

## New Mathematical Direction

THM-391 suggests a two-stage proof grammar:

```text
sieve branch covered
  => unit witness killed
  => zero-branch star peels by THM-391
  => endpoint debt exported to descendant labels
  => only labelled descendant cores can obstruct LRC.
```

This is a useful firewall.  It prevents future searches from mistaking
zero-branch occupancy, p-adic entropy, or full-support zero-flow abundance for a
local endpoint obstruction.  Those are real signals, but the proof core must
live one layer lower.
