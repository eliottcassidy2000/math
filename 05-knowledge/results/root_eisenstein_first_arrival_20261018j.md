# Eisenstein source charges and first arrival energy

**Status: PROVED elementary reductions and finite-horizon mass bound;
FINITE-EXACT censuses; OPEN uniform quadratic bound and universal entry.**
Session `clock_duality-20261018j`. Proofs are self-audited, with independent
integer implementations for the numerical Gram checks. No theorem ID or
priority claim is assigned.

The first-arrival problem has a useful two-coordinate lift. Three cyclic
rotations of the actual source weights form one Eisenstein integer. Its
norm retains cancellation without choosing one residue phase. A uniform
bound allowing a factor `3/2` of source norm would already imply
`E_V(N) <= 15N` for `N >= 2^14`, inside the existing large-source budget.
The norm bound is **OPEN**. Its stronger, contractive version passes the
declared finite tests, including every arrival prefix in the smaller test.

## 1. Inheritance and live concepts

- Closest proved mechanism: [THM-4690, first-arrival ladder and clock Gram](../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md),
  building on [boundary transport R5 and R6](root_sparse_boundary_module_20261017i.md).
  The certified forward-closed frontier has 627 seeds and descending ranks.
- Hostile: arbitrary signed transport is not an unweighted contraction;
  first-arrival cross-port orthogonality is refuted in
  [HYP-9325](../hypotheses/HYP-9325-one-tick-law-of-first-arrivals-into-the-627-seed-frontier.md).
  A new explicit hostile below also excludes matrix contraction on the
  three-phase source space.
- Corrected near miss: moving mass must remain in finite-horizon ledgers;
  ROOT-supported current assigns zero to sources outside the ROOT basin.
  The whole-support correction in the inherited R1--R4 is retained.
- Least-used sidecar: the two unused residue phases, with the *same*
  original sources, terminal ports and integer edge clocks.

Concept board: marked first-arrival clock polynomials; the Eisenstein
charge lattice; ordered graph boundary tensors; coefficient-copy
identifications; the golden three-port defect. The anchor is first-arrival
energy, the niche is geometric gluing, and the wildcard is the exact
star-defect/dynamical-parameter family in the companion note.

The incoming [graph/Pythagorean synthesis](graph_pythagorean_20261018_synthesis.md)
keeps both legs of a boundary-state triangle. Here that advice becomes an
explicit identity for the actual energy and its unused phase coordinate,
rather than an analogy between scalar sizes.

## 2. Exact source lift

Let `S={s>0: s odd, s != 5 mod 8}`, `B(s)=4s+1`, and use the native
transport `P` and edge clock `A=a+2j` of THM-4690. The clock here is that
canonical integer edge clock; for an uncompressed source the endpoint
braid correction of R5 is still required. No modular clock is substituted.

For odd `m<4N` define three weights, indexed by `r in Z/3`:

```text
K_N^(r)(m) = kappa(m+r)  if m<2N,
             lambda(m+r) otherwise,
kappa = (0,-3,3),        lambda = (-1,-1,2).
c_N^(r)(s) = sum_{j: B^j(s)<4N} K_N^(r)(B^j(s)).
```

The row `r=0` is exactly the original source, not a replacement. The other
rows rotate weights only; routes, first arrivals and clocks do not change.

Put `omega^2+omega+1=0`, with positive imaginary part, and define

```text
A_N(s) = sum_{j: B^j(s)<4N} omega^(B^j(s)) times
          (omega-omega^2 if B^j(s)<2N, otherwise omega).
```

**Proposition E1 (PROVED).** `A_N(s)` belongs to `Z[omega]`, and

```text
c_N^(r)(s) = 2 Re(omega^r A_N(s)),
sum_r c_N^(r)(s)=0,
sum_r c_N^(r)(s)^2 = 6 Norm(A_N(s)),
Norm(A_N(s)) belongs to {0,1,3}.
```

*Proof.* Both weight triples are the real parts stated above. Along a
B-chain, `B^j(s)=s+j mod 3`. If its last point is in `[2N,4N)`, there is
only one such point and

```text
(omega-omega^2) sum_{j<L} omega^j + omega^(L+1) = omega.
```

Thus its charge is the unit `omega^(s+1)`, independently of L. If all L
points are below `2N`, the charge vanishes when `3|L`; otherwise its norm
is 3. The norm/three-phase identity follows by summing a full period. QED.

In coordinates `A=a+b omega`, the three rows are

```text
(2a-b, -a-b, -a+2b),      Norm(A)=a^2-ab+b^2.
```

These coordinates are recovered integrally from the rows:
`a=(c^(0)-c^(1))/3`, `b=-(c^(0)+2c^(1))/3`.

## 3. First-arrival polynomials and the sufficient norm target

Keep the 627-seed certified set `V`. For each source s that enters V,
write `u(s), t(s), dep(s)` for its entry seed, integer clock and native
depth. Define the finite-horizon, port-labelled polynomial

```text
Q_(u,N,d)(z) = sum_{s: u(s)=u, dep(s)<=d} A_N(s) z^t(s).
```

Sources not yet entering are retained as moving sources. This definition
requires no universal-entry assumption and has finite support for each N.
For coefficients one may also retain a separate formal source label;
aggregating it is an explicitly declared operation.

For an uncompressed source `m=B^j0(s)`, also retain m, j0 and the terminal
braid jlast. Its physical shortcut clock after a nonempty native path is
`t+2j0-2jlast`; at native depth zero it is zero. Restoring these labels
changes the formal clock polynomial but not its value at -1. The Grams
in the receipts use the seed-source clock convention of THM-4690.

Set `W_u=Q_(u,N,d)(-1)` and `C_u^(r)=2 Re(omega^r W_u)`. Then

```text
T_(N,d) := sum_{u,r} (C_u^(r))^2 = 6 sum_u Norm(W_u),
E_(N,d) := sum_u (C_u^(0))^2 <= (2/3) T_(N,d).
```

**Proposition E2 (PROVED).** If `W_u=a_u+b_u omega`, the sharper identity is

```text
E_(N,d) + 3 sum_u b_u^2 = (2/3) T_(N,d).                 (1)
```

The missing leg is nonnegative and explicit. Taking the trace forgets the
orientation of this two-dimensional phase vector, but (1) records exactly
what is lost. The inequality is sharp when every `b_u=0`.

Let `T_src(N)=sum_{s,r} c_N^(r)(s)^2`. The proposed target is

```text
(ET)  T_(N,d) <= (3/2) T_src(N)
      for every N>=4096 and every finite native depth d.  OPEN
```

The stronger `(ET1)` replaces `3/2` by `1`. These are source-specific
statements. They are not operator contraction claims on arbitrary inputs.
The final-arrival version is a weaker sufficient target; it is not
silently inferred from, or substituted for, the finite-horizon version.

**Proposition E3 (PROVED implication and source bound).** Put
`J(N)=1+ceil(log_4(12N+1))`. Then

```text
T_src(N) <= (102/7) N + 72 J(N) + 6,
T_src(N) = (102/7) N + O(log N).
```

Consequently `(ET)` implies

```text
E_V(N) <= (102/7) N + 72 J(N) + 6 <= 15N  (N>=2^14).   (2)
```

*Proof of the source estimate.* For each L, let `e_L` count seed chains
whose point at level L lies in `[2N,4N)`. Let `i_L` count chains with L
points below `2N` and the next point above `4N`. The corresponding source
intervals have lengths respectively `2N/4^L` and `(4N+1)/4^L`. Seeds occupy
three of eight residues, so the count differs from `3/8` times length by
at most 3. All nonempty levels have `L<J`. By E1,

```text
T_src/6 = sum_(L>=0) e_L + 3 sum_(L>=1, 3 does not divide L) i_L.
sum_(L>=0) 4^-L=4/3,
sum_(L>=1, 3 does not divide L) 4^-L=20/63.
```

The leading term is `17N/7`; endpoint errors are at most `12J`, and the
extra `+1` in `4N+1` contributes less than 1. Truncation tails are bounded,
giving also the asymptotic equality. Multiply by 6.
For the last inequality, `J<=log_2 N+4`; at `N=2^14`,
`72 log_2 N+294 <= (3/7)N`, and the right-minus-left derivative is positive
beyond this point. E2 and ET give (2).

To pass to final arrivals, for fixed N the set of sources that ever hit V
is finite, so their arrival depths have a finite maximum even if other
sources never enter. Taking d beyond that maximum proves the stated
bound without asserting entry of the other sources. QED.

The existing R6 current inequality gives `R(N)^2<=627 E_V(N)`. Bound (2)
lies below the `C=15.2` calibration in HYP-9325's stated density-exponent
implication chain. This note proves the reduction, not ET, that external
chain anew, or convergence of every positive Collatz orbit.

## 4. A genuinely uniform but weaker certificate

Define the unit-step norm

```text
ell(a+b omega)=max(|a|,|b|,|a-b|).
```

**Proposition E4 (PROVED).** At every finite horizon, the total ell-mass
of all frozen and moving charges is at most the initial mass, and

```text
sum_s ell(A_N(s)) <= (41/21)N + 9J(N)+1.
```

*Proof.* Multiplication by the transport sign is an isometry; merging is
subadditive. A source of norm 1 has ell-mass 1; a source of norm 3 has
ell-mass 2. Repeat the `e_L+2i_L` count from E3. QED.

Every zero-sum multiset of the six unit directions decomposes into three
types of opposite pair and two alternating equilateral triples. Indeed,
remove opposite pairs first. Remaining support has at most three
directions; a positive zero relation then requires the alternating
triple, with equal multiplicities. This supplies a finite cancellation
alphabet with source provenance. For a nonzero sum, a pair at angle
120 degrees can be replaced by its unit sum; this is another orientation
of the same three-term relation.

This is an unconditional finite-horizon mass certificate, including the
moving complement. It does **not** bound the quadratic energy by O(N):
mass could concentrate at one terminal. The precise unpaid strengthening
is ET. A cancellation at z=-1 must keep its clock polynomial; it certifies
divisibility by z+1 in `Z[omega][z]`, not equality at every frequency.

**Clock-cost lemma (PROVED).** A parity cancellation need not have small
integer-clock cost. For an odd lag L and nonzero Eisenstein charge a, put
`Q(z)=a z^t(1+z^L)`, t>=0. Then Q(-1)=0, but the unique quotient
`R=Q/(z+1)` has L coefficients of size ell(a), and

```text
sum_j ell([z^j]R) = L ell(a),
R(-1)=Q'(-1)=(-1)^t L a.
```

Indeed `(1+z^L)/(1+z)=1-z+...+z^(L-1)`. Consequently exact phase
cancellation alone cannot pay a bound on lag or on the quotient mass.
This hostile is for general packets; it does not assert that arbitrary
such pairs occur as first arrivals of this specific source family.

The incoming [Lucas-block CM detector note](lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md)
distinguishes a period's vanishing from its first-carry defect, with
weight and split/inert information retained. That is useful methodological
inspiration here: inspect the first surviving layer after a value has
cancelled. The concrete ROOT layer is Q'(-1), not a modular-form carry;
no map between those arithmetic predicates or theorem transfer is claimed.

## 5. Exact tests and hostile

The C program builds every first-arrival route with a 10,000-step cap and
64-bit arithmetic. Every orbit multiplication is overflow-guarded; for
its allowed N<2^26 the Gram entries are bounded by (6N)^2<2^63. The cap is never reached. It updates the three
source rows at `4N-3`, `4N-1`, and `2N-1` with exact integer Gram matrices.
The independent Python program sums complete B-chains and aggregates
routes directly; it does not use the incremental source update.

**FINITE-EXACT:**

- Every `4096<=N<2^24`: final `T_V<=T_src`, and the clock-diagonal trace
  is at most `T_src`. Thus all three declared final-arrival predictions
  pass. Maximum `T_V/T_src` in this range is below 0.990249, at N=5997.
- Every `4096<=N<2^20`, every native arrival prefix d: `T_(N,d)<=T_src`.
  Maximum ratio is below 0.991829, at `(N,d)=(5997,66)`. The largest source
  depth in this universe is 190. Later prefixes are unchanged.
- Dyadic N from 8 through `2^21`: every first-arrival prefix is saved.
- The 19 common dyadic points independently match all 27 entries of the
  source, arrival and clock-diagonal Gram matrices (513 checks).
  The two C modes also match all 20 common dyadic Gram records exactly.
- Source shell and exact norm identities: 74,321 seed cases, for all
  N=1..200 and N=1024,4096,7910,16384.
- A third implementation follows the literal uncompressed odd inputs at
  N=1204 and N=5997, independently checking 108,202 integer-clock
  identities. At N=5997 it gives source trace 87,372, final trace 86,520,
  and maximum prefix trace 86,658 at depth 66. At N=1204 these are
  17,538, 20,076, and 20,154 at depth 45. No compressed-clock equality
  is assumed without the braid boundary correction.

The proposed large-source cutoff matters: at N=1204 the final trace is
more than 1.14 times source trace. Small-source contraction is false.

**REFUTED stronger claim:** source Gram minus arrival Gram is positive
semidefinite. At N=524288, for `v=(-1,3,-2)`,

```text
v^T (G_arrival-G_source) v = 2,763,090 > 0.
```

Thus cancellation can coexist with expansion in another phase direction.
The repaired target keeps only the trace on this particular source
family. Its truth at unbounded N and d remains open.

## 6. Reproduction and next decisive test

Scripts and receipts live in `04-computation/experiments/`:

```bash
python3 frontier_phase_20261018j.py --max-power 21 --out frontier_interfaces_20261018j/root_phase_prefixes.json
cc -O3 -Wall -Wextra frontier_phase_20261018j.c -o /tmp/frontier_phase_20261018j
/tmp/frontier_phase_20261018j first_arrival_20261018/V_closed_Y1000.txt 24 > frontier_interfaces_20261018j/root_phase_everyN.json
/tmp/frontier_phase_20261018j first_arrival_20261018/V_closed_Y1000.txt 20 prefix > frontier_interfaces_20261018j/root_phase_prefix_everyN.json
python3 ternary_interfaces_20261018j.py --out frontier_interfaces_20261018j/ternary_controls.json
python3 frontier_interface_audit_20261018j.py --out frontier_interfaces_20261018j/independent_audit.json
```

The first final-arrival predictions were written before the scan in the
tool sequence and committed with its output, not independently committed
as a preregistration. Prefix testing was an exploratory continuation.

The [manifest](../../04-computation/experiments/frontier_interfaces_20261018j/manifest.json)
records the scripts, inherited inputs, commands and receipts.

The next proof task is a source-specific bound on concentration of the
surviving Eisenstein charges. The six-direction cancellation alphabet
gives exact local certificates, but an arbitrary local contraction rule
cannot replace that task: the displayed Gram hostile already forbids it.
Keep the entry port and full clock polynomial while grouping charges;
test any proposed potential against the N=5997 prefix and the anisotropic
N=524288 hostile before extrapolating.
