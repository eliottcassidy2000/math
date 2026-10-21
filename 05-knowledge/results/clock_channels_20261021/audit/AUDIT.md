# Independent audit: finite static observers and integer clocks

Date: 2026-10-21. Result: **PASS, no mathematical correction required.**

Audited frozen proof: `FINITE_OBSERVERS_INTEGER_CLOCKS.md`, SHA256
`a20ec35bcd408192091c93181d8d38c5ce8c7672982a162ac50bef93facc6304`.
The six entries of the producer's manifest all match. Its normal and optimized
reruns both reproduce the retained 263-check receipt byte for byte.

The independent checker has **46,083 successful checks** in both normal Python
and `python3 -O`; these receipts are byte-identical. The count includes repeated
integer ancestry and trajectory controls, not 46,083 distinct mathematical
theorems. No producer module is imported.

## 1. Exact source, map, and clock audit

The independent source generator uses the recurrence

    s_0 = 3,  s_(j+1) = 64 s_j + 35,

and then checks `9 s_j + 5 = 32*64^j`. This separates source construction from
the producer's closed-form implementation. For `j=0,...,96`, the audit checks
positive integral sources congruent to 3 modulo 8 and directly reconstructs the
primitive B-ancestry. There are exactly `3j+1` removals, each of the form
`upper = 4*lower+1`, ending at the actual integer 1. The primitive braid clock is
therefore `1+2(3j+1)=6j+3`.

The algebraic all-j proof is valid: `s_j` is integral because `64=1 (mod 9)`;
its reduction modulo 8 is 3; and

    (3s_j+1)/2 = (4^(3j+2)-1)/3 = B^(3j+1)(1).

Each nonterminal B-iterate of 1 is 5 modulo 8, so the stripping count is exact.
The actual odd word is `(1,6j+4)`, since the first value is odd and the next
numerator is exactly `4^(3j+2)`. Thus literal shortcut ROOT cost is `6j+5`,
and ordinary cost is `6j+7`. The braid cost is not substituted for either.

The independent checker separately generates the full ordinary trajectory and
the shortcut trajectory, then compresses the ordinary trajectory and compares
it against the independently generated shortcut trajectory. It recovers the odd
word from actual runs of even ordinary steps. All three costs agree with the
proof, including the producer's 31 retained rows.

The printed decimal example is parsed from the proof, rather than supplied as
an independent test constant: `s_0=3`, `s_6=244335917283`, primitive costs 3 and
39, and shortcut ROOT costs 5 and 41 all replay exactly.

## 2. First descents and stopped frontiers

The four printed shortcut values are exact. For `s>=3`, their differences from
s are, respectively,

    (s+1)/2, (5s+5)/4, (s+5)/8, (5-7s)/16.

The first three are positive, the fourth negative. Hence first literal descent
occurs at shortcut clock 4 for every family member. The first odd checkpoint
is larger than the source, while the next odd checkpoint is 1. First odd descent
therefore occurs only at the unbounded shortcut ROOT clock. The audit also
recovers ordinary first-descent clock 6, although the proof does not need that
extra observation.

For a fixed finite frontier V containing 1, a source already in V has stopped
depth and clock zero. Otherwise its first primitive edge reaches 1 in V,
so first-arrival depth is 1 and clock `6j+3`. Since `s_j` increases without bound,
all sufficiently large sources exceed `max(V)`. This proves the eventual claim
for every such V, without needing the membership list of the inherited
627-point frontier. The audit exercises three explicit finite frontiers,
including nonconsecutive family members, and confirms both stopped cases.
It does **not** claim to enumerate the inherited frontier.

## 3. Modular periods, exact dyadic threshold, and full clocks

For positive L,

    s_(j+L)-s_j = 32*64^j*(64^L-1)/9.

The quotient `(64^L-1)/9` is an odd integer. Thus the exact 2-adic valuation of
the difference is `6j+5`. If L is a multiple of the order of 64 modulo `9 M_o`,
its odd part is divisible by `M_o`, and congruence modulo `M=2^a M_o` holds
precisely when `6j+5>=a`. The proof only requires sufficiency; the checker also
tests failure below this sharp threshold.

Independent direct multiplicative-order searches confirm orders 3, 2, and 1
modulo 27, 5, and 7. Consequently L=6 works modulo 945, and the dyadic condition
for 840 holds even at j=0. All actual P, shortcut ROOT, and ordinary ROOT costs
of `j+6` and j differ by 36. P costs agree modulo 12 but differ as integers.

The audit uses nine odd modulus parts and twelve dyadic exponents, reaching
`a=200`, and checks 342 source-pair threshold cases. It separately checks Euler's
period choice. For each modulus and each of eight finite clock moduli q, it
uses `lcm(L,q/gcd(q,6))`; the full clock difference stays positive while becoming
zero modulo q. These finite controls agree with the all-modulus proof.

## 4. Factor criterion and information-loss scope

The factor criterion is correct: an observed next state exists as a function of
the current observation exactly when observed successors are constant on every
observer fibre. Retaining the integer edge cost requires the pair consisting
of observed successor and cost to be constant there. Exhaustive finite models
on one through four states independently verify both formulations.

The small companion citations have the claimed type: the actual sources 63
and 273 share the 15/35 register, but P gives respectively `(95,1)` and `(51,4)`
with different observed successors. These are unstopped single-edge braid
costs, not fully odd costs or stopped first-arrival costs. The actual odd steps
3 to 5 and 7 to 11 break the idoneal-flag factor. Independent complete reduced-
form enumeration of the small orders confirms the flags used here. No global
idoneal-completeness theorem is required or audited in this packet.

For uniform j on `0,...,J-1` and any deterministic q-valued static observer O,
`H(j|O)=H(j)-H(O)>=log2(J)-log2(q)` is valid. The injectivity of `6j+3` rules out
exact reconstruction from any finite static state on the entire infinite
family. It does not establish a density assertion for all integer sources.

The proof explicitly excludes the invalid stronger conclusion. Finite alphabets
with unbounded words, finite controllers with unbounded counters, and finite
parameterized rule schemes are not ruled out. Nor does the argument prove that
every convergence proof must retain each exact clock. The restriction is to a
claim that the compressed state itself preserves the specified source-and-clock
information. The fixed-order genus-group comparison respects its different
consumer and does not transfer algebraic losslessness to a timed dynamical map.

## 5. Source verification and reproducibility

The inherited native primitive convention and the distinction between shortcut
and ordinary clocks were read directly from these repository sources, whose
hashes match the packet's provenance:

- `05-knowledge/results/native_cancellation_105_20261018g.md`:
  `91b52f8e7bd7332ff1edcbd1dfdd0b43afca15cf17c7e451b32f3d5621c0b839`.
- `05-knowledge/results/guarded_children_105_20261018g.md`:
  `75a59e714f2aba35ebf6b9dcaff3e52e06d520ff42542bc06be35bf3f25c575b`.

The cited idoneal companion proof has the recorded SHA256
`5b7790cc4adf658ab14cf893941fdbe0bd7a2d2c04aed135141cc105c6c1c229`.
Its cited explicit source examples were additionally replayed here. Classical
genus theory and the companion's global analytic background are outside this
independent clock audit; the principal grounded family does not depend on them.

Run the checker with a path to the frozen producer packet:

    python3 independent_clock_channels.py --packet ../session18-clock-channels
    python3 -O independent_clock_channels.py --packet ../session18-clock-channels

The same relative location is the default. The checker uses only the Python
standard library and reads no repository files outside the selected packet.
The separate source checks described above are preserved in
`AUDIT_PROVENANCE.json`. Neither mathematical proof nor infinite-family claims
are inferred solely from finite successful replays.

No repository or remote files were modified during this audit.
