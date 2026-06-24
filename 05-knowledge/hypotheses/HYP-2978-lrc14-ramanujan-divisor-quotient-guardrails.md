---
id: HYP-2978
title: LRC14 Ramanujan-divisor quotient guardrails
status: EVIDENCE / quotient-admissibility proof guardrail with S162 interval scaffold; not an LRC14 proof
source: codex-2026-06-24-S161
related:
  - HYP-2979
  - HYP-2977
  - HYP-2976
  - HYP-2975
  - HYP-2974
  - HYP-2981
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2963
  - HYP-2956
  - HYP-2946
  - HYP-2938
  - HYP-2887
  - HYP-2886
  - HYP-2264
  - THM-406
  - THM-572
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc14_ramanujan_divisor_named_channel_guardrails_codex_s161.py
  - 05-knowledge/results/lrc14_ramanujan_divisor_named_channel_guardrails_codex_s161.out
  - 04-computation/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.py
  - 05-knowledge/results/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.out
  - 07-reflections/lrc14-ramanujan-divisor-quotient-guardrails-codex-s161.md
  - 04-computation/lrc14_packet_fejer_interval_scaffold_codex_s162.py
  - 05-knowledge/results/lrc14_packet_fejer_interval_scaffold_codex_s162.out
  - 07-reflections/lrc14-robin-robbins-fejer-interval-scaffold-codex-s162.md
---

# HYP-2978: LRC14 Ramanujan-Divisor Quotient Guardrails

Core rule:

```text
A quotient Q is admissible for an LRC14 proof step P only if P is constant
on every Q-fiber, or Q carries a named certificate for the labels forgotten
inside nonconstant fibers.
```

Equivalently, each quotient must declare:

```text
preserved predicate
forgotten labels
compensating transform or side-channel
defect certificate when a forgotten label is load-bearing
```

This is not a proof of LRC14.  It is a proof-safety theorem target: before a
scalar, tournament, divisor, moment, Ramanujan, or packet quotient can rule out
a strict counterexample, it must declare exactly which LRC predicate survives
the quotient and which lost labels are reattached.

## External Seed

The divisor-function neighborhood gives the right arithmetic warning.  The
divisor functions `sigma_k(n)` are multiplicative, `sigma = Id * 1` under
Dirichlet convolution, and admit Ramanujan-sum expansions.  Ramanujan sums

```text
c_q(n) = sum_{(a,q)=1} exp(2*pi*i*a*n/q)
       = sum_{d | gcd(q,n)} d*mu(q/d)
```

are simultaneously divisor/Mobius data and primitive-root phase data.  Thus a
scalar divisor quotient is already hiding a cyclotomic packet ledger.  For
LRC14 this is exactly the missing bridge: `tau`, `sigma`, `omega`, perfect-number
product loads, Jordan/totient unit capacity, squarefree support, and unitary
divisor counts are useful only if they retain or recover exact-period phase,
endpoint owner, and route labels.

The immediate LRC14 hook is HYP-2974's divisor-curried Fourier coefficient:

```text
hat F_S(k) = sum_{v in S, v|k} sin(pi*(k/v)/7)/(pi*(k/v)).
```

Mode `k` sees both the divisor fiber `v|k` and the quotient `k/v`.  Ramanujan
sums retain primitive residue-period data after averaging over units instead of
collapsing to a bare divisor count.

HYP-2979 is the companion retained-packet route: exact-period Ramanujan
projectors for q-ladders, endpoint sums, and primitive unit phase packets.  It
uses `c_14` as a primitive-unit trace with four values keyed by `gcd(14,n)`:
`6,-6,-1,1`.  HYP-2978 supplies the admissibility guardrail for that route: the
projector is not enough merely because it is arithmetic; it must be homogeneous
for the proof predicate or carry endpoint/safe-measure/K33 certificates.

## Named-Channel Audit

Artifact preserved from the concurrent S161 remote checkpoint:

```text
04-computation/lrc14_ramanujan_divisor_named_channel_guardrails_codex_s161.py
05-knowledge/results/lrc14_ramanujan_divisor_named_channel_guardrails_codex_s161.out
```

This audit checks arithmetic identities through `n<=80`, `q<=24`:

```text
tau=sigma0, sum phi, phi=mu*id, psi=id*|mu|,
sum J2, and the Ramanujan c_q Mobius/gcd formula
```

with no mismatches.  It then audits named rows: AP, GW, residue liar `12->26`,
near/K33 `12->36`, petals, P10 splices, and covering repairs.  The warning is
sharp: AP, GW, the residue liar, K33, petals, and P10+K33 all share the same
coarse lcm divisor scalars; several also share the same `c_14` profile.  Those
channels mix AP/GW boundary, q-witness, K33, petal, and covering proof routes.

Named-channel collision counts:

```text
qdiv_only route-mixing collisions                 1
open_state_only route-mixing collisions           1
mod14_residue_multiset route-mixing collisions    1
ramanujan_14_profile route-mixing collisions      1
unit_counts_14_27_41 route-mixing collisions      2
divisor_lcm_scalars route-mixing collisions       1
guarded_packet_signature route-mixing collisions  0
```

Tournament Analysis over quotient channels is transitive:

```text
labelled_lrc_packet_sheaf
  > exact_period_packet
  > ramanujan_primitive_shell
  > gcd_strata
  > totient_jordan_unit_capacity
  > squarefree_psi_support
  > raw_divisor_counts
```

with `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}` and no directed 3-cycles.

## Web pass and source trail

The divisor-function page makes the useful warning explicit.  `sigma_a(n)` is
both a divisor-fiber count/sum and a convolutional coefficient: its Dirichlet
series is `zeta(s) zeta(s-a)`, its Lambert series re-sums the same divisor
fibers into harmonic/modular coefficients, and for `k>0` it has a Ramanujan
expansion

```text
sigma_k(n) = zeta(k+1) n^k sum_m c_m(n)/m^(k+1).
```

So `sigma`, `tau`, and `psi` are not proof predicates by themselves.  They are
pushforwards of divisor fibers, useful only if the kernel of the pushforward is
declared.

One-hop reading sharpened the algebra:

- Dirichlet convolution says arithmetic functions form a local Dirichlet ring;
  multiplicative functions form a subgroup under convolution, and identities
  such as `d=1*1`, `sigma_k=Id_k*1`, `phi*1=Id`, and `mu*1=epsilon` are packet
  laws rather than loose analogies.
- Jordan's totient `J_k` counts primitive `k`-tuple capacity and satisfies
  `sum_{d|n} J_k(d)=n^k`, so it is the higher-rank exact-period version of the
  `phi` packet used in HYP-2899.
- The divisor summatory function counts lattice points under `jk<=x`; the
  divisor problem is exactly a boundary-defect problem after quotienting the
  product lattice by the hyperbola.
- Ramanujan sums retain the primitive unit shell and have orthogonality over
  lcm periods.  Murty's survey records the Carmichael-style shifted
  orthogonality `M(c_r(n)c_s(n+h))=c_r(h)` when `r=s` and `0` otherwise, which
  directly suggests shifted packet autocorrelations for LRC endpoint sums.
- The supercharacter reading of Ramanujan sums says the same thing in
  representation language: quotient by unit-group orbits is legitimate only
  because the orbit character is retained.

External pages read in this pass:

```text
https://en.wikipedia.org/wiki/Divisor_function
https://en.wikipedia.org/wiki/Multiplicative_function
https://en.wikipedia.org/wiki/Dirichlet_convolution
https://en.wikipedia.org/wiki/Ramanujan_sum
https://en.wikipedia.org/wiki/Jordan%27s_totient_function
https://en.wikipedia.org/wiki/Divisor_summatory_function
https://en.wikipedia.org/wiki/Lambert_series
https://en.wikipedia.org/wiki/Eisenstein_series
https://arxiv.org/abs/1201.1060
https://mast.queensu.ca/~murty/HRJ-2013.pdf
```

## S162 Addendum: Robin/Robbins and Packet-Anchored Fejer Intervals

The requested divisor-function reading adds a useful distinction.  Robin's
number-theory theorem is a scalar extremal inequality for `sigma(n)/n`; it is
useful here mainly as a warning that a divisor average is not a proof object
unless its prime-power and phase packet can be reconstructed.  Robbins'
graph-theory theorem is the better structural analogy: a strong orientation
exists exactly when no bridge prevents reachability.  In LRC14 terms, a quotient
is only admissible when it keeps the certificate graph's bridges: endpoint
owner, exact q-phase, Fejer degree/center, labelled packet route, and
K33/state-lift debt.

The one-hop divisor-function trail reinforces the same guardrail:

```text
sigma_k            multiplicative, but prime-power exponents are load-bearing
Dirichlet convolution  sigma = Id * 1 and Id = sigma * mu need an inverse channel
Ramanujan sums     c_q(n) is the primitive-root trace / Mobius-inverted shell
divisor summatory  divisor counts live on a hyperbola/lattice carrier
superabundant/Robin extremality  scalar maxima still need labelled divisor packets
```

S162 uses that lesson on HYP-2974's Fejer certificates.  It computes rational
interval enclosures for selected Fejer quadratic forms and records the labelled
packet fiber

```text
P(S) = (route, family, q_class, packet_route, state_lift, q_threshold).
```

All five selected rows have interval upper endpoint `<0`: `12->36` at degree
`159`, `P10+GW` at degree `280`, `12->168` at degree `63`,
`two drop(12,13)->add(14,29)` at degree `41`, and `single swap 6->63` at
degree `266`.  This is the first concrete bridge from floating Fejer evidence
to packet-keyed interval certificates.  It remains a scaffold until the `pi`
and trigonometric interval backend is formalized and the full HYP-2963 bank is
grouped by packet fiber.

## Guardrail theorem target

The abstract lesson across the repo can now be stated as a quotient-kernel
rule.  Let `P` be the predicate needed by the next LRC implication, and let
`pi:X->Q` be a proposed quotient.  The quotient is admissible only if every
fiber of `pi` satisfies at least one of:

```text
P is constant on the fiber;
the forgotten coordinate is reconstructible from retained packet data;
the forgotten coordinate is killed by a proved orthogonality/dual certificate;
the fiber is assigned to an explicit residual bucket with endpoint/state labels.
```

This is the common shape behind the older lessons:

- irreducibility/coefficient quotients need Newton slope, fixed-divisor, or
  residue-depth side channels before sign/count shadows are trusted;
- unital prompts split geometric unitals, algebraic identity preservation, and
  unit groups instead of using the same word for incompatible kernels;
- Faulhaber moment shadows preserve positivity but not compatibility/integrality
  packets;
- Pollock dyadic residues fail singly, while carry-pair defect lifts recover the
  load-bearing predicate;
- unit-distance scalar edge counts need spine/bulk/direction/ear labels;
- tiling and regular-solid analogies work only when curvature, norm recursion,
  and local vertex-figure defects are retained;
- HYP-2899's `Div(D) x B_r` product lattice is the exact denominator/far-packet
  version of the same rule.

Thus the fundamental object is not a connection between two scalars.  It is a
connection whose kernel has been named and discharged.

## Tournament-adjacent proof angle

Use Tournament Analysis on quotient channels or proof obligations, not on
runners.  The pairwise observable is:

```text
Does channel A separate every route collision that channel B mixes, while
preserving the LRC predicate: strict witness, AP/GW boundary equality, or named
state-lift debt?
```

Ties are broken by the Hamiltonian path

```text
exact-period packet
> endpoint-owner label
> dual certificate
> residual bucket
```

The S161 named-row audit already gives the first transitive channel tournament.
The proof-facing version should be stronger:

```text
Every post-THM-571 Moon-core packet maps into the labelled LRC packet sheaf.
If a coarser quotient is used, its kernel must be contained in the AP/GW
equality atoms plus packets killed by Ramanujan orthogonality, Toeplitz PSD,
moment duals, taut-current potentials, or HYP-2908/THM-572 state lifts.
```

A hypothetical LRC14 counterexample is then not just "not seen" by current
features.  It must lie in the intersection of the kernels of q/Farey witness
ladders, Ramanujan primitive shells, Haar-open fronts, Toeplitz PSD tests,
danger-count and multiplicity moment duals, endpoint taut currents, C27/K33
owner labels, and covering boundary-moment packets.  HYP-2978's job is to make
that intersection exact and then empty it, or name the first genuinely new
residual.

Remaining work:

1. Extend the collision audit from named rows to packet families in the
   HYP-2963 bank.
2. Turn the quotient-kernel rule into a theorem-facing lemma over the current
   labelled packet sheaf.
3. Test HYP-2979's shifted Ramanujan autocorrelations against HYP-2973 count
   moments and HYP-2974 Toeplitz/Fejer failures.
4. Stress-test this admissibility criterion against new packet-family mixed
   fibers, using the one-swap addendum below as the first broad test case.

## One-Swap Fiber Audit

Second artifact from the continued S161 pass:

```text
04-computation/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.py
05-knowledge/results/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.out
```

This audit uses exact rational danger-interval union arithmetic on:

- `10` named rows: AP, GW `12->24`, K33 near `12->36`, petal rows, two-swap
  splices, and covering rows `6->98`, `12->84`, `12->168`;
- the one-swap AP neighborhood through `add<=220`, deduplicated, for `2694`
  rows total.

Proof-route target:

```text
boundary-zero, or qdiv=<d> / qdiv>14 plus exact safe-measure bucket
```

Named packet highlights:

```text
AP                         qdiv=14 safe=0          c14(v)=-6 c14(sum)= 6 c14(diff)=-36
GW 12->24                  qdiv=14 safe=0          c14(v)=-6 c14(sum)= 6 c14(diff)=-29
near 12->36                qdiv=14 safe=1/1260     c14(v)=-6 c14(sum)= 6 c14(diff)=-29
petal 10->20               qdiv=14 safe=1/980      c14(v)=-6 c14(sum)= 6 c14(diff)=-29
covering 12->84            qdiv=15 safe=563/105105 c14(v)= 1 c14(sum)=-1 c14(diff)=-36
```

The important signal is negative: `c_14(v_i+v_j)` catches the AP/GW/K33/petal
zero-credit trace (`c14(sum)=6`) and separates covering repairs (`-1`), but it
does not by itself distinguish AP/GW equality from positive K33/petal exits.
Endpoint/safe-measure labels or a K33 state-lift certificate must be reattached.

Each quotient was bucketed and checked for mixed proof-route fibers.

```text
quotient              classes  bad_fibers  bad_pair_collisions  largest_fiber
scalar_divisor           2403         138                  239              5
unitary_divisor          2677          12                   18              3
qcover                     37          10                76948            628
ramanujan_speed           164          75                72586            395
ramanujan_pair           1564         265                 2291             38
exact_period_packet      2491          14                   14              4
endpoint_measure         2112           0                    0             76
full_row                 2694           0                    0              1
```

Representative scalar-divisor mixed fiber:

```text
swap 5->173   route=qdiv=14:open         safe=50569/1731730
swap 6->118   route=qdiv=14:small<=1/100 safe=5353/708708
swap 10->122  route=qdiv=10:open         safe=83711/3846843
swap 11->179  route=qdiv=11:open         safe=232049/4560920
swap 13->181  route=qdiv=13:open         safe=966677/35121240
```

All five have the same scalar signature:

```text
sum tau=37, sum omega=15, sum bigOmega=20, sum sigma=309,
gcd14 counts=((1,6),(2,6),(7,1)).
```

Thus scalar divisor data is not an admissible proof carrier for the qdiv/safe
route predicate.  It can be a feature, heuristic, or irreducibility ledger, but
not a final quotient unless a refinement recovers the lost phase/route labels.

Representative exact-period packet mixed fiber:

```text
AP          route=boundary-zero          safe=0
swap 12->96 route=qdiv=14:small<=1/100   safe=5219/840840
```

This is the sharper warning: even qcover + Ramanujan speed/pair packets +
residue multiset can still forget the open-vs-boundary predicate.  HYP-2979 is
therefore useful only after endpoint owner labels, exact safe measure, or
K33/state-lift debt are included.

## Fiber Tournament

Vertices were quotient/proof carriers, not runners.  Pairwise observable:
smaller `(bad_pair_collisions, bad_fibers, classes)` wins; after equal proof
admissibility, the more compact quotient wins.

```text
score_hist={3:1, 4:1, 0:1, 1:1, 2:1, 5:1, 7:1, 6:1}
directed_3_cycles=0
SCC_sizes=(1,1,1,1,1,1,1,1)
Hamiltonian_path_count=1
tie_Hamiltonian_path=
  endpoint_measure > full_row > exact_period_packet > unitary_divisor >
  scalar_divisor > ramanujan_pair > ramanujan_speed > qcover
```

The ranking is not saying endpoint-measure proves LRC14.  It says that for the
particular route predicate used in the audit, `qdiv + exact safe-measure bucket`
is the most economical homogeneous quotient.  To become a theorem, that
endpoint-measure quotient must be lifted back to structural labels: endpoint
owner pairs, AP/GW zero-credit current, exact-period Ramanujan packets, and
HYP-2908/THM-572 K33 state-lift debt.

## Guardrail Across Prior Themes

The same rule explains several older failures and useful analogies:

- irreducibility: degree or factor count is not enough; the irreducible packet
  and field/monodromy labels matter;
- unital designs: point/line counts are not enough; incidence ownership matters;
- Faulhaber/moment methods: a scalar moment can be positive while support
  compatibility is lost;
- Pollock defects: polygonal/tetrahedral counts forget the degree/dimension
  carrier unless the defect label is retained;
- unit-distance carriers: Euclidean norm-1 kissing layers and field-norm unit
  layers are different quotients;
- tiling/solid analogies: cell counts are not proof carriers until the boundary
  unit and gluing labels are declared;
- perfect/Farey product: `ab=|E(K_{a,b})|` is a scalar load, not a Kuratowski
  obstruction unless graph incidence/minor labels are retained.

## Theorem Target

The proof-facing statement should be:

```text
Let P be any LRC14 route predicate used to discharge a residual packet
(weak witness, strict open mass, boundary equality, qdiv route, Toeplitz
dual failure, K33 state-lift debt, etc.).  A quotient Q may be used as a
proof carrier for P only if every Q-fiber is P-homogeneous, or if Q is
paired with a refinement certificate R whose labels make P recoverable.
```

For the Ramanujan/divisor lane, the live target is:

```text
Every non-AP/GW LRC14 residual either has a positive endpoint/safe-measure
certificate visible after exact-period Ramanujan refinement, or it carries
a named K33/HYP-2908/THM-572 state-lift debt.  Scalar multiplicative
quotients alone cannot certify the route.
```

Next concrete tasks:

1. Add endpoint-owner labels to the HYP-2979 exact-period projector and rerun
   the fiber test against the HYP-2963 bank.
2. Compare `c_14(v_i+v_j)` directly against HYP-2970 endpoint credit
   `K=14(rm-sn)+r+s`.
3. Test shifted Carmichael/Ramanujan autocorrelation of danger multiplicity
   against HYP-2973 danger-count moments and HYP-2974 Toeplitz PSD failures.
4. Treat multiplicative functions as irreducibility ledgers, not proof-ending
   scalars, unless the lost-label certificate is explicit.
5. Turn regular A-functions into declared divisor-subset quotient metadata:
   even the allowed divisor set is load-bearing before convolution or primitive
   shell averaging.
6. Test divisor-summatory hyperbolic-simplex boundary defects against the
   Pollock, tiling, and solid-carrier analogies.
7. Try Euler pentagonal sigma recurrences as sparse boundary-moment ledgers in
   the HYP-2974/HYP-2977 harmonic stack.
8. Use Busche-Ramanujan non-coprime correction terms as the model for any
   product quotient that tries to forget shared prime support.
