# One uniform digit per transitive suffix chart

2026-11-03. **PROVED:** a fixed transitive suffix chart has at most one
clock-uniform affine ternary digit, even if each proposed output uses an
arbitrary new port set and its own fixed integer clock shift.
**FINITE-EXACT:** the independent controls and bounded odd-mask scan below.
**OPEN:** classification of the possible digits over all port masks,
source-compatible two-digit closure outside this chart, and the tournament
Hamilton-path spectrum conjecture.

The inherited mechanism is Session182's
[terminal-port criterion and binary port encoding](../session182_physical_returns_20261104/induction/PROOF.md),
read at commit `b486271d24`, Sections 2, 4, and 6. Its TT3-to-TT5 uniform
+4 family is the positive control. The canonical hostile to extending that
one success is its own output mask: it permits digit -28 in this chart,
not another +4. The least-used sidecar is the joint port correction j;
the marginal f alone determines it here only because the entire suffix is
transitive. This note makes no priority claim beyond this scoped deduction.

## 1. Actual charts and the preserved count

Let TT_m have vertices 0,...,m-1 and arcs i -> k for i < k. Let W be a
port subset containing m-1. Form M_t(W), for t >= 1, from a transitive
prefix TT_t, this suffix, and a new root r, with

    TT_t -> TT_m,   W -> r,   r -> TT_t and TT_m \ W.

These are actual strong tournaments: the suffix's terminal component is
the singleton {m-1}, which is a port. Conversely, within this attachment
model strongness requires that terminal vertex to be a port. Physical
names and any native-operation certificates are additional data; neither
is reconstructed from the count formulas.

An ordered two-slot cover puts every suffix vertex in exactly one of two
directed path slots; either slot may be empty. Let f count covers whose
first slot is nonempty and ends in W. Let j count covers whose two slots
are nonempty, the first ending in W and the second starting in W. Then

    f = sum_(i in W) 2^i,
    H_t(f) = H(M_t(W)) = 1 - j(f) + 2^t f.              (1)

The formula for f follows by fixing the first slot's endpoint i: vertices
above i must lie in slot two, vertex i in slot one, and each earlier vertex
can independently choose a slot. Since m-1 is a port, f is positive and
m = bitlength(f). Thus f specifies both the suffix order and W uniquely.
There is no unrecorded freedom to append terminal nonports while remaining
strong in this chart.

The count identity in (1) is the specialization of the inherited exact
attachment formula H(M(A;B,W)) = H(A)(H(B)-j)+P_A(2)f: both transitive
Hamilton counts equal one, and P_TT_t(2)=2^t because every vertex may choose
one of the two ordered slots. Thus its dependency is the cited attachment
proof, not a new bijection asserted here. The independent checker also
counts complete Hamilton paths in small actual attachments, rather than
using (1) as its only control.

## 2. The exact joint correction and binary-shift invariance

Write b_i = 1 if i is a port, and 0 otherwise. Directly counting the two
possible endpoint orders for every unordered pair of distinct ports gives

    j(f) = sum_(i<k, b_i=b_k=1) (2^(k-i-1) + [k=i+1]).  (2)

Here the first slot ending at k and the second starting at i allow each
vertex strictly between i and k to choose either slot, giving 2^(k-i-1).
If the first slot ends at i and the second starts at k, all intermediate
vertices would have no legal slot; this is possible exactly when k=i+1,
with one cover. Equal endpoints are impossible because the slots are
disjoint. This proves (2), including adjacent ports and empty-slot
exclusion, without appealing to a numerical recurrence.

Consequently, for every a >= 0,

    j(2^a f) = j(f).                                    (3)

Geometrically this prepends a nonport transitive block to the suffix.
Only differences of port positions occur in (2). All summands in (2) are
even, so j(f) is even. These are all-order statements, with no bound on
the number of ports, suffix order, or leading nonport block.

One useful independent evaluator scans b_0,b_1,... with initial g=j=e=0:

    j_new = j + b_i*(g+e),
    g_new = 2*g+b_i,
    e_new = b_i.                                        (4)

The old g is the bit-reversal value of the previously scanned bits and e
is their final bit. Equation (4) is the inherited ordered-component law;
the checker compares it with both (2) and literal ordered subset covers.

## 3. Uniqueness, including independent output clock shifts

**Theorem.** Fix a positive integer port mask f. Suppose a second fixed
positive mask f', an integer k, and an integer d satisfy

    H_(t+k)(f') = 3 H_t(f) + d                           (5)

for every sufficiently large integer t (so both prefix lengths are at
least one). Then necessarily

    2^k f' = 3f,
    d = 3j(f) - j(3f) - 2.                             (6)

Conversely, whenever f' is a positive integer satisfying the first
equation, (5) holds at every legal clock with d as in (6). In particular,
from one fixed source chart there is at most one uniform digit, independent
of the output clock shift k.

Proof. Substitute (1) in (5). Evaluating at any two distinct legal clock
values cancels the constant terms and gives 2^k f'=3f. Equating the constant
terms then gives d=3j(f)-j(f')-2. If k >= 0, the first equation and (3)
give j(f')=j(2^k f')=j(3f). If k < 0, they give f'=2^(-k)3f and the same
conclusion. This proves necessity. Substitution proves sufficiency.
In fact equality at any two distinct legal clock values already forces
the all-clock relation; an infinite testing assumption is unnecessary. QED.

The two proposed output letters may choose entirely different suffix
orders, port sets, and fixed clock shifts. Their slopes still force their
masks to be binary shifts of 3f, and hence the same j. Therefore +2 and +4
cannot both be uniform letters from the same source while both targets
stay inside these TT-prefix/TT-suffix charts. Increasing the number of
ports without bound, freely changing those ports, or changing either
output clock by a fixed amount does not evade this obstruction.

## 4. Exact examples and the remaining escape routes

The forced digit d(f)=3j(f)-j(3f)-2 is unchanged by f -> 2^a f. For the
small positive controls:

| f | j(f) | j(3f) | forced digit |
|---|---:|---:|---:|
| 1 | 0 | 2 | -4 |
| 3 | 2 | 4 | 0 |
| 5 | 2 | 14 | -10 |
| 7 | 6 | 12 | 4 |
| 11 | 8 | 16 | 6 |
| 21 | 12 | 62 | -28 |

Thus the inherited f=7 to f'=21 return is an actual uniform +4 family,
H_t=-5+7*2^t -> -11+21*2^t. Its target chart cannot repeat this uniform +4
operation. The f=11 control prevents mistaking the obstruction for a bound
d <= 4. A small literal construction may still succeed at just one prefix
length; that does not satisfy the theorem's two-clock hypothesis.

A both-digit uniform grammar must leave the declared transitive charts
for at least one output, or give up using a fixed output chart/clock shift
for all input clock values. Nontransitive suffix components, nontransitive
prefixes, output port masks or structures chosen as a function of t, and
other attachment models are outside this theorem. In a larger carrier the
marginal slope f and joint correction j can be separate coordinates; this
proof does not impose their TT-specific relation on such states. It also
does not obstruct isolated source-compatible returns at a prescribed t,
or prove anything about which native Hamilton-complement guards exist.

This is an operation-language boundary, not a global nonexistence result.
It narrows which charts could satisfy the six-entry accepted-source
[coverage criterion](coverage/COVERAGE.md). Neither that conditional
criterion nor the exact multiport forest carrier supplies the remaining
two-digit closure by itself.

There is a concrete sharpness witness in the companion
[two-digit entry construction](transport/TWO_DIGIT_ENTRY.md). From the
same TT3/all-port source, keep the first two vertices as ports and replace
its terminal singleton by a directed triangle with that old terminal as
its sole port. The resulting suffix is TT2 -> C3, with signature

    (h,p,f,g,j,e,s) = (3,48,21,39,16,1,3).

Its terminal component contains a port, so the full attachment is strong,
and its count is -13+21*2^t = 3*(-5+7*2^t)+2. Together with the inherited
TT5 alternating-port +4 target, this supplies both uniform letters at that
one entry. Exactly one suffix leaves the TT chart, as the theorem requires.
The companion proof supplies the literal two-birth construction and its
own independent controls. This sharpness witness gives no closure proof
from either output and no native Hamilton-complement legality assertion.

## 5. Reproducible finite controls

Run from the repository root:

    python 05-knowledge/results/multiport_defect_grammar_20261103/check_transitive_chart.py
    python -O 05-knowledge/results/multiport_defect_grammar_20261103/check_transitive_chart.py

Both modes use explicit exceptions, not Python assertions, and write the
same deterministic `transitive_chart.json` receipt. The checker includes:

- Every mask 1 <= f <= 255, evaluated independently by all ordered TT
  subset covers, the pair formula, and the bit recurrence.
- Fixed output clock shifts -4 <= k <= 4 whenever 3f/2^k is integral,
  checking the derived return at two distinct legal clocks and testing
  literal shift invariance.
- Actual full attachments for every mask 1 <= f <= 31 at prefix lengths
  1, 2, and 3: directed subset-DP Hamilton counts and independent strongness.
- All 32,768 positive odd masks f < 2^16, comparing the pair and bit
  evaluators for both f and 3f before computing the forced digit.

In that finite odd-mask universe there are no digit-2 masks, and the only
digit-4 mask is f=7. By the proved shift invariance this classifies masks
with odd part below 2^16, but **does not** classify every positive mask.
No all-order absence of digit 2 or uniqueness of odd f=7 for digit 4 is
claimed. The finite controls audit the formulas; the all-order uniqueness
theorem rests on the proofs above, not on the scan.
