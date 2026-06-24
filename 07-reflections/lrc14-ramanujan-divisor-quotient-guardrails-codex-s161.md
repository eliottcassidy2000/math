# LRC14 Ramanujan-Divisor Quotient Guardrails

Session: codex-2026-06-24-S161
Linked hypothesis: HYP-2978
Child route: HYP-2979

## Core Lesson

The reusable abstraction is not "find the right quotient."  It is:

```text
Before quotienting, state the contract of forgetting.
```

A quotient is allowed to forget a coordinate only when that coordinate is:

```text
invariant under the proof step,
reconstructible from retained data,
annihilated by a proved dual or orthogonality relation, or
parked in a labelled residual bucket.
```

This is the same lesson across the user's current reservoirs.

```text
irreducibility:        coefficient/factor quotients must remember convolution lifts
unital designs:        incidence quotients must remember pair-ownership
Faulhaber moments:     even/odd collapse must remember which moments drive balance
Pollock defects:       simplex-count quotients must remember dimension and defect degree
unit-distance carriers:norm quotients must remember layer/owner information
tiling and solids:     volume quotients must remember boundary and half-copy seams
multiplicative funcs:  product quotients must remember shared prime support
LRC14:                 runner quotients must remember endpoint, q/Farey, Haar, and state labels
```

The divisor-function neighborhood made the rule sharper.  `sigma_k`, `tau`,
`phi`, `J_k`, and `psi` are not merely scalar features.  They are ledgers for
which divisor fibers, unit fibers, squarefree supports, and prime-power layers
survive a quotient.  Dirichlet convolution is the honest multiplication law
because it keeps the hidden factor split `ab=n`.  Ramanujan sums are similarly
useful because `c_q(n)` is not just a residue statistic; it is the primitive
unit trace of an exact-period layer.

## What The Computation Says

The S161 script audited named LRC14 rows using scalar channels:

```text
qdiv
open/zero-open state
mod-14 residue multiset
c_14 Ramanujan profile
unit counts at 14, 27, 41
lcm divisor scalars
guarded packet signature
```

Every scalar channel except the deliberately labelled guarded packet signature
had route-mixing collisions.  That is the useful negative result.  AP/GW
boundary atoms, q-witness rows, K33/state-lift rows, petals, and covering rows
can share the same coarse divisor or primitive-shell data.  A proof step using
one of those channels alone would silently identify rows that need different
exit theorems.

So HYP-2978 should be read as a quotient admissibility theorem target rather
than a direct LRC14 proof:

```text
If Q is a quotient used in the LRC14 proof, then every fiber of Q must carry a
constant next-step certificate, or Q must attach the labels needed to split the
fiber into certificate-constant residual buckets.
```

This is almost tautological, but it is powerful because it gives a falsifier:
find two named rows with the same Q-value and different proof route.  The S161
audit does exactly that for multiple tempting scalar quotients.

## Tournament View

The tournament vertices in this pass were not runners.  They were quotient
channels:

```text
raw_divisor_counts
squarefree_psi_support
totient_jordan_unit_capacity
gcd_strata
ramanujan_primitive_shell
exact_period_packet
labelled_lrc_packet_sheaf
```

The pairwise observable was payload retention: orient `A -> B` when `A` retains
strictly more declared proof payload than `B`.  This yields a transitive
tournament with the Hamiltonian path:

```text
labelled_lrc_packet_sheaf >
exact_period_packet >
ramanujan_primitive_shell >
gcd_strata >
totient_jordan_unit_capacity >
squarefree_psi_support >
raw_divisor_counts
```

The challenged assumption was that tournament vertices have to be runners or
arcs.  For quotient safety, the more natural vertices are proof obligations,
divisor fibers, unit twists, primitive Fourier modes, endpoint-owner labels,
C27/K33 state packets, and dual certificates.  The LRC predicate preserved by
this quotient-channel tournament is not "which runner wins?" but "does this
compressed object still force the next implication without mixing exit routes?"

## Cross-Domain Payoffs

Irreducibility suggests the convolution-lift guardrail.  A polynomial
coefficient vector is reducible if it has a hidden factor-coefficient grid whose
diagonal sums equal the visible coefficients.  Quotienting to visible
coefficients alone is admissible only if the missing grid is either impossible
or explicitly represented as a residual certificate.

Unital designs suggest incidence ownership.  A block quotient that keeps only
point counts forgets the key property: every pair is owned by exactly one
block.  This mirrors LRC endpoint ownership and C27/K33 packet ownership.

Faulhaber moments suggest parity contracts.  In the anchor equation, only odd
moments drive the balance after midpointing.  Forgetting even moments is honest
there because symmetry annihilates them; that is an example of a legal
forgetting operation.

Pollock and tiling analogies suggest defect ledgers.  A volume count without
boundary defect degree is like an lcm scalar without endpoint labels: it
detects capacity but not obstruction.

Multiplicative functions give the cleanest arithmetic warning.  Multiplicative
laws are simple only on coprime inputs.  Non-coprime product identities need
Busche-Ramanujan-style correction terms.  In LRC language, any product-like
quotient that forgets shared prime support must supply a correction term.

Regular A-function Ramanujan sums add one more guardrail: even "which divisors
are allowed" is data.  A quotient must declare the divisor subset before it can
use a convolution or primitive-shell average.

## Most Promising Next Move

The next high-leverage test is to remove the route label from the guarded
packet signature and add only mathematically justified labels:

```text
q/Farey threshold
Haar/Baire open-state label
Ramanujan exact-period profile for selected q
unit-capacity profile from phi/J_k/psi
endpoint-owner multiset
C27/K33 state-lift marker
Fourier/Toeplitz or danger-count dual certificate marker
```

Then audit HYP-2963 packet families for route-mixing.  If route-mixing remains
zero without the route label itself, the resulting packet is a plausible finite
endpoint for an LRC14 proof interface.  If route-mixing remains, the first
collision tells us exactly which label the quotient still illegally forgets.

That is the proof style now suggested by this session: not a miracle scalar
invariant, but a labelled packet theorem whose labels are minimized by repeated
collision falsification.
