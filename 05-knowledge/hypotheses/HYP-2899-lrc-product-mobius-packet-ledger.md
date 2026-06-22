---
id: HYP-2899
status: STRUCTURAL PROOF-TARGET / product-lattice ledger, not a proof of LRC14
source: codex-2026-06-22-S110
tags: [lrc14, euler-totient, coprime-density, mobius, multiplicative-functions, half-tiling, product-lattice, tournament-analysis]
related:
  - HYP-2628
  - HYP-2856
  - HYP-2886
  - HYP-2898
  - HYP-2896
  - HYP-+2897
  - THM-442
  - THM-549
  - THM-550
  - THM-551
  - HYP-2689
  - HYP-2890
  - OPEN-Q-108
results:
  - 04-computation/lrc_totient_tiling_recursion_s110.py
  - 05-knowledge/results/lrc_totient_tiling_recursion_s110.out
---

# HYP-2899: product-Mobius packet ledger for LRC14

The owner asked to connect coprime density, Euler-totient/multiplicative
functions, and the three tournament tiling recurrences:

```text
full:       A+B+C-D-E-F+G
even half:  A+B-C
odd half:   A+B-C+D-E-F+G
```

S110's answer is that these are not separate coincidences.  They are the same
instruction at different address levels:

```text
retain primitive packet labels before scalarizing.
```

On the denominator side, the packet labels are divisor/exact-period labels.  On
the tiling/far side, the packet labels are Boolean subset labels.

## Exact audit

Script:

```text
04-computation/lrc_totient_tiling_recursion_s110.py
```

Stored output:

```text
05-knowledge/results/lrc_totient_tiling_recursion_s110.out
```

### Divisor Mobius side

The Euler-copy rule

```text
sum_{d|n} c(d) = n
```

has the unique solution

```text
c(n) = phi(n) = sum_{d|n} mu(d) n/d.
```

This is the exact-period packet law used in HYP-2628 and HYP-2886: residues
`a/q` split by reduced denominator `d|q`, with packet size `phi(d)`.  The
coprime-density floor is then forced by totient sums:

```text
sum_{b<=Q} phi(b)      ~ 3Q^2/pi^2
sum_{b<=Q} phi(b)/b    ~ 6Q/pi^2
```

S110 checks the convergence numerically.  For example, at `Q=210`,

```text
sum_phi/Q^2      = 0.304172   vs 3/pi^2 = 0.303964
sum(phi/b)/Q     = 0.608143   vs 6/pi^2 = 0.607927
```

This is the HYP-2856 Farey floor in multiplicative-function form.

### Boolean Mobius side

The three tiling recurrences are Boolean Mobius kernels:

```text
full B3:
  A+B+C-D-E-F+G       signs +++---+

even half B2:
  A+B-C               signs ++-

odd half prompt order:
  A+B-C+D-E-F+G       signs ++-+--+
```

The odd half sign string has two useful addresses.  Incoming S31q observes
that in prompt order over `A..G = 1..7`, the signs `++-+--+` are the Legendre
`chi_7` split with the zero/triple slot forced positive: positive on quadratic
residues `{1,2,4}` and on slot `7`, negative on nonresidues `{3,5,6}`.  In the
half-tiling address, however, the three geometric corner generators are `A,B,D`,
with pair overlaps
`C,E,F`; in that address the same odd half recurrence becomes

```text
A+B+D-C-E-F+G         signs +++---+.
```

So the odd half-tiling is simultaneously a `chi_7` numeric sign word and a
`B3` geometric inclusion-exclusion word.  The product-ledger rule is to keep
both addresses until the proof decides which quotient is being used.  The real
tiling difference is the pushforward to size/depth offsets:

```text
full staircase:
  offsets n-1,n-2,n-3      ->  3, -3, +1

even half:
  offsets n-1,n-2          ->  2, -1

odd half:
  offsets n-1,n-2,n-3,n-4  ->  2,  0, -2, +1
```

The zero at `n-2` is not missing structure.  It is geometric cancellation:
the prompt's `-C` and `+D` live at the same size but occupy different overlap
slots.

S110 verifies:

```text
full Delta^3 recurrence:                 true for n=5..20
even half recurrence:                    true on even n
odd half recurrence:                     true on odd n
global half recurrence (x-1)^3(x+1):     true
```

## Product-lattice proof target

The LRC14 witness/cap proof should keep both packet labels:

```text
Div(D) x B_r
```

where `Div(D)` is the divisor/exact-period denominator lattice and `B_r` is the
Boolean lattice of far runners, cover arcs, or tiling generators.

At the scalar level, the exact-period main term looks like

```text
N(S,D) main term = (6/7)^13 * phi(D).
```

At the Boolean three-generator level, independent danger probability `x=1/7`
has packet mass

```text
B3 nonempty IE = 1-(6/7)^3 = 127/343.
```

These two facts should be multiplied only after the corresponding residual
labels have been retained.  Scalarizing early to `phi(D)` alone loses the
CRT/multiplicativity defect; scalarizing early to a seven-letter sign alone
loses the denominator/fiber where the exact packet lives.

The proposed OPEN-Q-108 subtarget is therefore:

```text
For covering rows, delete divisibility-killed denominators, then bound the
residual packet error on Div(D) x B_r below the product-Mobius main ledger.
```

This refines HYP-2886's exact-period packet target and HYP-2890's support-six
residual target.  Incoherent high-denominator/product-lattice packets should
route to the `zeta(2)` / L2 / Parseval floor; coherent low-rank packets should
route to finite AP/Freiman/Clebsch-Bruhat/octahedral atlases.

Incoming HYP-2898/S111 reinforces the same guardrail from a different lab:
small even `q=8,10,12,14` cases preserve Bonferroni floors, bounded p0 caps,
AP-facing Fejer/difference-profile structure, and labelled residual leaks,
while scalar additive-energy monotonicity fails.  That fits this ledger: the
positive Mobius/totient principal is real, but the proof cannot discard the
Legendre/Eisenstein/Fejer labels before bounding the residual.

Incoming mac-mini S44 gives the complementary resonance-killing version:
a speed `s` kills every primitive Farey point of denominator `b` exactly when
`b|s`, so it kills a whole `phi(b)` packet at once.  The small denominator
survival lattice through `14` has size `Phi(14)=sum_{b<=14}phi(b)=64`, making
the covering hard core a totient-weighted overdetermined killing problem rather
than a list of independent scalar denominators.

## Link to the latest proof/disproof dialectic

Incoming mac-mini S43 / HYP-+2897 says the proof/disproof split is:

```text
non-covering rows:
  q-witnessed by THM-523.

covering rows:
  overdetermined and empirically loose, but barely-safe seeds need
  equidistribution rather than a naive union-bound measure inequality.
```

S109/HYP-2896 closed the one-tail branch over `C={1,...,11,13}`.  S110 adds the
packet interpretation of why that branch escapes:

```text
w=84m kills the small q-witness denominators,
but the binding witness denominator D=84m+5 is coprime to 2,3,7.
```

Thus q-killing by divisibility does not remove the coprime packet floor.  It
pushes the witness to a new unit denominator packet.  S110's table shows
`gcd(84m+5,84)=1` for `m<=12` and, in fact, identically for all `m`.

## Tournament Analysis

Vertices are proof carriers:

```text
divisor_mobius_phi
exact_period_unit_packets
full_B3_tiling_recursion
even_B2_half_recursion
odd_B3_half_recursion
product_lattice_packet_ledger
scalar_projection_guardrail
```

Observable: how much primitive packet/address data survives before
scalarization.  Switch: carriers preserving both divisor and Boolean Mobius
labels rank above carriers that collapse them.  Tie Hamiltonian path:

```text
divisor_mobius_phi
  > exact_period_unit_packets
  > full_B3_tiling_recursion
  > even_B2_half_recursion
  > odd_B3_half_recursion
  > product_lattice_packet_ledger
  > scalar_projection_guardrail
```

Challenged assumption: the three tiling recurrences are unrelated formulas.
They are all Mobius kernels; the differences come from which quotient pushes
the kernel to sizes/depths.  Preserved predicate: primitive packet capacity and
signed cancellation.  Destroyed predicate: the raw geometry of individual
runner arcs, which must be reintroduced only inside finite coherent atlases.


## VERIFIED SHARPENING (S45): the prime-covering reduction pins the hard case at 30030|v
- 64% of random 13-sets have a surviving prime <=13 => M>=1/13>1/14 (EASY). A counterexample must be
  PRIME-COVERING: a runner divisible by EACH of {2,3,5,7,11,13} (6 constraints) AND kill b=14.
- The radical handle, verified: smallest prime not dividing v = lcm(2..X) is 7,11,13,17 for X=5,7,11,13.
  So a committed speed gives a prime-witness <=13 (M>=1/13) UNLESS v is divisible by ALL primes <=13
  (i.e. 30030 | v). **The unique hard committed-speed case is 30030 | v** (the lcm/2310-radical family) --
  exactly where the prime-witness fails (next prime 17 > 14) and equidistribution (Node 3) is required.
- Clean roadmap: surviving-prime (easy, 64%) -> prime-covering -> {bounded: Node 2 three-gap/AP-majorization;
  unbounded with 30030|committed-speed: Node 3 equidistribution/Weyl}. The two nodes are exhaustive.