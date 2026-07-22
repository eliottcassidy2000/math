---
id: THM-2075
title: "Safe-child homeomorphism and dyadic wall-word conjugacy"
status: >
  PROVED from THM-2073's unique safe-child law. Doubling is a homeomorphism
  between consecutive closed safe sets in the dyadic tower. Component count,
  Euler characteristic, and endpoint count are invariant; every component
  has one constant binary address, its length and total measure halve exactly
  at each level, and no odd guard can create an endpoint over the interior of
  the quotient safe set. Iteratively, every endpoint is inherited from a
  terminal-core speed. The address-decorated component wall word is therefore
  conjugate to the terminal word. Ambient cyclic order is not preserved if
  addresses are discarded. This is a carrier theorem, not LRC(14).
source: codex-2026-07-21-LRC-safe-child-homeomorphism
depends_on:
  - THM-2073
related:
  - THM-2072
  - THM-775
  - THM-2047
---

# THM-2075 -- safe-child homeomorphism

For `0<delta<1/2` and a finite set `Q` of positive integers, write

```text
G_Q(delta)={sigma in R/Z:||q sigma||>=delta for every q in Q}. (1)
```

The application below uses `delta=1/14`; the elementary topological lemma
does not depend on that value.

## 1. One safe child is a homeomorphism

Let `h` be odd and put

```text
P=2Q union {h}.                                          (2)
```

Assume the exact safe-child condition

```text
for every sigma in G_Q(delta), exactly one of
L_0(sigma)=sigma/2,    L_1(sigma)=(sigma+1)/2
belongs to G_P(delta).                                  (3)
```

Then doubling restricts to a homeomorphism

```text
D:G_P(delta) -> G_Q(delta),    D(theta)=2theta.          (4)
```

### Proof

If `theta in G_P(delta)`, then for every `q in Q`,

```text
||qD(theta)||=||2q theta||>=delta,
```

so (4) is well-defined. Every `sigma in G_Q(delta)` has a preimage by the
existence half of (3), and it has only one by uniqueness. Thus `D` is a
continuous bijection. Both safe sets are closed subsets of the compact
circle, so the domain is compact and the codomain is Hausdorff. A continuous
bijection between them is a homeomorphism. QED.

This recovers THM-2072's antipodal exclusion in a stronger form: `G_P` cannot
contain both `theta` and `theta+1/2`, because they would have the same image
under the injective map (4).

## 2. Components, Euler characteristic, and exact metric scaling

Each set in (1) is a finite union of closed intervals and isolated points.
It never equals the whole circle because `0` is not safe. Let `I` be one of
its connected components. The inverse section

```text
s=D^(-1):G_Q(delta)->G_P(delta)                           (5)
```

has a constant sheet bit on `I`: there is a unique

```text
epsilon_I in {0,1}
```

such that

```text
s(sigma)=(sigma+epsilon_I)/2    for every sigma in I.    (6)
```

Indeed, after lifting the proper circle arc `I` to a real interval, its two
inverse sheets are disjoint. The continuous image `s(I)` is connected and
cannot jump between them. This also covers a singleton component.

Consequently:

```text
#components(G_P)=#components(G_Q),
chi(G_P)=chi(G_Q),                                       (7)

length(s(I))=length(I)/2,
measure(G_P)=measure(G_Q)/2.                             (8)
```

Here `chi` is ordinary Euler characteristic; every nonempty component is
contractible, so it is simply the number of components. Equations (7)--(8)
are exact, including isolated points, which contribute zero length and one
component.

The component endpoints are inherited. More precisely,

```text
theta in boundary(G_P)  iff  D(theta) in boundary(G_Q).  (9)
```

The forward direction follows because if `D(theta)` were in the relative
interior of a safe interval, formula (6) on a small neighborhood would put
`theta` in the interior of its selected child interval. The reverse direction
follows in the same way from the local affine homeomorphism (6).

Therefore the odd guard cannot be the sole creator of a new endpoint above
the interior of `G_Q`. If `theta` is a boundary point of `G_P`, then
`sigma=D(theta)` is a boundary point of `G_Q`, so some `q in Q` satisfies

```text
||q sigma||=delta.
```

The even speed `2q in P` then satisfies

```text
||2q theta||=delta.                                      (10)
```

The guard `h` may tie at the endpoint, but an inherited even-core owner is
always present.

## 3. Iteration along the THM-2073 tower

Retain the tower

```text
C=Q_0,
Q_i=2Q_(i+1) union {h_i},    0<=i<r,                    (11)
```

from THM-2073, at `delta=1/14`. Its unique safe-child theorem supplies (3)
at every level. Hence

```text
D^r:G_C -> G_(Q_r),    theta |-> 2^r theta              (12)
```

is a homeomorphism. For every connected component `I` of `G_(Q_r)`, there is
a unique address

```text
a_I in {0,1,...,2^r-1}                                  (13)
```

such that the inverse of (12) on `I` is

```text
s_I(sigma)=(sigma+a_I)/2^r.                              (14)
```

The binary digits of `a_I` are precisely the successive safe-child choices.
They are constant on the whole component, not merely pointwise labels.

The complete iterated consequences are

```text
#components(G_C)=#components(G_(Q_r)),
chi(G_C)=chi(G_(Q_r)),                                   (15)

length(s_I(I))=2^(-r)length(I),
measure(G_C)=2^(-r)measure(G_(Q_r)).                     (16)
```

Moreover every endpoint `theta` of `G_C` descends to an endpoint
`sigma=2^r theta` of `G_(Q_r)`. If `q in Q_r` owns that terminal endpoint,
then the actual speed `2^r q in C` owns the original endpoint:

```text
||q sigma||=1/14  implies  ||2^r q theta||=1/14.         (17)
```

Thus none of the odd guards `h_0,...,h_(r-1)` is needed as a sole endpoint
owner anywhere in the tower. Their exact role is sheet selection.

## 4. The address-decorated wall word

For a terminal component, record its two endpoint owner sets, its length, and
its weak/strict convention, and attach the address `a_I` from (13). Equations
(14) and (17) transport its geometry and all inherited terminal-owner labels
to `G_C`:

```text
terminal component I
  -> child component (I+a_I)/2^r,
endpoint owner q -> endpoint owner 2^r q,
length ell -> ell/2^r.                                  (18)
```

This is an exact conjugacy of the **address-decorated component wall word**.
It retains precisely the data needed to reconstruct the safe-set geometry at
the top of the tower and a nonempty inherited owner set at every endpoint.
An odd guard may tie an inherited endpoint; recovering the complete tie set
requires evaluating the known guards and is a separate label sidecar.

There is one important guardrail. Different terminal components may have
different addresses, so the map (14) can permute their ambient cyclic order.
If the addresses are forgotten, the unadorned circular wall word need not be
order-conjugate. Component incidence, owner labels, and metric scaling
survive; raw cyclic order survives only together with (13).

## 5. Consequences for the LRC(14) residual

THM-2073 reduced every non-hereditarily-primitive seam core to a terminal
hereditarily primitive `Q_r`. The present theorem shows that this reduction
loses no topology and introduces no new endpoint obligations:

```text
all safe-set components and endpoints already exist at the terminal core;
the dyadic guards only choose binary lifts and shrink lengths.              (19)
```

This changes the live terminal question. A proof need not classify fresh
guard-created wall combinatorics at every depth. It may classify terminal
component/endpoint words once, then test whether any assignment of the finite
address bits can meet the original two-tail folded-owner condition.

The theorem does not perform that final test. In particular, homeomorphic
safe sets can occupy different circle positions after their address choices,
and those positions matter to the odd tails. Euler characteristic or measure
alone cannot close the seam.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption is that odd guards should be treated as new wall
vertices. Equation (10) proves the opposite: in a valid safe-child tower they
are sheet-choice bits, while terminal-core speeds remain the endpoint-owner
vertices. The quotient to component count preserves Euler data but destroys
addresses and cyclic placement; the lossless carrier is the terminal
component set with endpoint labels and the address function `I |-> a_I`.

One can orient two components by their numerical addresses and use terminal
component order as the tie Hamiltonian path. This produces a transitive
tournament with score histogram `(0,1,...,k-1)`, no directed cycles, and only
singleton SCCs. Those fingerprints add nothing: the LRC predicate depends on
the affine placements (14) and endpoint owners (17), not pairwise address
order. The address-decorated wall word, rather than a tournament, is the
faithful object. QED.
