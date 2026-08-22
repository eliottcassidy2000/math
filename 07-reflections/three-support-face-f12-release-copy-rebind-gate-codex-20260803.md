# A released face tag does not rebind, but one explicit state copy closes the reachable fork

**Status: NON-CANONICAL FINITE-EXACT PARTIAL SCOUT.**  This note continues
the
[three-face all-choice closure](three-support-face-all-choice-reset-fork-closure-codex-20260803.md).
Exact evidence is the
[release/copy gate scout](../04-computation/fc3_three_support_face_f12_release_copy_rebind_gate_scout_20260803.py)
and its
[stored output](../05-knowledge/results/fc3_three_support_face_f12_release_copy_rebind_gate_scout_20260803.out).
No theorem ID is claimed.  The positive copy result is conditional on a new
operation that is absent from the inherited carrier; it proves no `FC(3)`,
`SFC(3)`, GMC, positivity statement, or shared physical chronology.

## Inheritance pass and portfolio

The closest proved mechanism is
[THM-3286](../01-canon/theorems/THM-3286-three-face-availability-helly-defect-and-binary-origin-width.md),
which gives the static colouring

```text
F12 | (F13,F23).
```

The canonical hostiles are the four Boolean-fork blockers in the preceding
all-choice closure, with smallest member

```text
C=(3,4,5,6,7,8).
```

The corrected near miss is the scheduling statement that `F12` resets no
later than the shared fork.  It is an exact inequality, but it does not turn
THM-3286's face-origin label into a reusable token.  In
[THM-3278](../01-canon/theorems/THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary.md)
the origin bit is a persistent face tag and a static cut coordinate, not an
occupancy token, state register, or physical transition.

The least-used sidecar was the exact terminal object of the `F12` lane.  The
portfolio was:

- **Anchor:** decide whether an inherited `F12` reset can lawfully rebind at
  the four shared blockers;
- **Niche:** classify the blocker reached from every one of the `94` starts,
  independently of row choices; and
- **Wildcard:** add one explicit diagonal copy of the current shared state and
  test whether that single missing operation removes every later response
  obstruction.

The live concept board contained the persistent face tag, the terminal
`F12` state, the shared-lane ancestry, active-face deletion, reset potential,
and a possible state diagonal.  The applicable method cards were **Audit
sections, not only fibres** and **Close a candidate section under its next
native operation before trusting its scalar shadow**.

## What the inherited carrier actually releases

An inherited node is `(A,n)`, where `A` is the active face block and `n` is
one raw physical multiset.  An edge exists only when every face in `A` admits
the same literal pair

```text
(row r, one-pole successor n').
```

The target active set deletes exactly those faces for which `n'` is the
pinned reset.  It never acquires a face.

The full all-choice closure has `69` row-labelled active-set changes.  Every
one is

```text
({F12},n) --r--> (empty,Q12),
Q12=(1,2,2,3,4,5).
```

Only seven of these terminal edges use rows `{16,17}`.  There are

```text
pair active-set changes = 0,
face-acquisition edges  = 0.
```

The exact cross-product of the `69` release edges with the four reachable
blockers has `276` candidates and no inherited rebind.  At every blocker the
shared decorated fibre is already empty; none equals `Q12`; and no terminal
edge targets it.

There is a subtle positive-looking control.  `Q12` is separately a legal
decision state on both `F13` and `F23`, at reset distance six on either face.
Retagging it is nevertheless not a repair: it would acquire a new face and
replace the current shared ancestry by the unrelated `F12` ancestry.  Syntax
and facewise typing do not supply the missing map.

Thus the current-carrier obstruction is typed:

```text
deleting the last F12 obligation
  != producing a reusable state register at the current pair state.
```

Strictly speaking, no origin object becomes free in the proved model.  There
is only an empty active set at `Q12`.  Treating that event as a reusable
resource is already a new controller sidecar.

## The shared fork has a row-independent normal form

The all-choice computation sharpens the preceding four-blocker census.  From
each of the `94` starts `n`, every unrestricted shared path ends at one and
the same blocker, and the same endpoint law holds when every edge is
restricted to rows `{16,17}`:

```text
B(n)=C union (support(n) intersect {1,2}).              (1)
```

The four endpoint counts are

| blocker | starts |
|---|---:|
| `C union {1,2}` | 47 |
| `C union {1}` | 16 |
| `C union {2}` | 24 |
| `C` | 7 |

This is not merely a graph digest.  Its mechanism is coordinatewise.  The
two resets ask for opposite low-root occupancies:

```text
              root 1   root 2   root 3   root 4
F13 reset        1        0        2        1
F23 reset        0        1        1        2.
```

For roots `1` and `2`, both reset directions remove copies while the
multiplicity exceeds one.  At multiplicity one, one face wants to keep the
root while the other wants to delete it, so the common cone freezes at the
presence bit.  At multiplicity zero the converse add/keep conflict freezes
absence.  The inherited triple carrier starts with at most one copy of root
`3` and at most one of root `4`; common reset directions fill these and roots
`5,...,8` to one.  This proves the normal form in (1) from the directed-cone
geometry and explains why changing response rows cannot change the endpoint.

The unrestricted number of row-labelled shared paths from one start to its
unique blocker ranges from `4,698` to `159,909,458,497,526`; with rows
`{16,17}` it ranges from `6` to `95,694`.  These multiplicities change, but
the endpoint does not.

The separate-lane timing comparison remains exact.  Every `F12` path has
length `d_12(n,Q12)` and every pair path has length

```text
(d_13(n,Q13)+d_23(n,Q23)-4)/2.
```

Their difference has histogram

```text
0:7, 1:24, 2:16, 3:16, 4:31.
```

This compares two reset-directed DAG lengths from the same start.  It does
not interleave their states or declare simultaneous physical time.

## The smallest explicit copy extension tested

Now add one new operation, separately from the inherited audit.  A combined
record carries

```text
an F12 completion certificate ending at Q12,
the current shared ancestry ending at b=B(n), and
two origin-state registers after the split.
```

At a reachable blocker `b`, define

```text
COPY_b:
  (released certificate; shared {F13,F23} state b)
    |--row=None-->
  ((F13,b),(F23,b)),                                   (2)
```

with either assignment of the two register names.  The old shared register
retains the corresponding restriction of the shared ancestry.  The other
register carries an explicit `COPY_FROM_SHARED_ANCESTRY` link.  Equation (2)
copies `b`, not the terminal state `Q12`; has no response row; and changes
the pair reset potential from `4` to `4`.  It is a diagonal on controller
states, not a one-pole physical edge.

After (2), every move is again an inherited, row-labelled, reset-directed
one-pole move on one singleton face.  Exhausting all choices gives:

| blocker | `F13` distance / paths | `F23` distance / paths | `{16,17}` path counts |
|---|---:|---:|---:|
| `C union {1,2}` | `2 / 198` | `2 / 230` | `2,2` |
| `C union {1}` | `1 / 9` | `3 / 6800` | `1,6` |
| `C union {2}` | `3 / 6084` | `1 / 10` | `6,1` |
| `C` | `2 / 171` | `2 / 200` | `2,2` |

No singleton continuation has a blocker.  Both register assignments work at
all four states.  Because (1) gives one endpoint per start and the release
inequality holds on all starts, the explicit extension succeeds on all `94`
start/blocker pairs.  Every post-copy completion uses exactly four singleton
one-pole moves in total.  The stronger control also passes: rows `{16,17}`
alone suffice before and after the copy on this generated closure.

At the minimal blocker the canonical paths can be seen directly:

```text
F13: C --16,+1--> C union {1} --16,+3--> Q13,
F23: C --16,+2--> C union {2} --16,+4--> Q23.          (3)
```

The copy is therefore the only non-inherited event in (3).

## Hostile controls and failure boundary

The full-domain `{16,17}` hostile remains load-bearing.  Across all `6,668`
active coloured contexts there are `38` canonical-pair blockers.  Granting
COPY only at the four generated blockers leaves `34`.  Even granting the
same local copy continuation at all seven unrestricted common-row blockers
leaves `31` row-specific canonical-pair blockers.  Thus the four-state copy
gate does not promote `{16,17}` to a global three-face controller.

The response-positive two-pole diagonal is also unchanged.  At `C`, its only
jointly typed targets are obtained by `(+1,+2)` and `(+1,+4)`.  Both preserve
the distance vector `(2,2)`, land on unrestricted blockers, and have no
common first one-pole step.  They are not edges of (2) or of any post-copy
path.  The new extension copies a state and then follows two honest
single-face descents; it never reinterprets the old macro as chronology.

The strongest surviving statement is conditional:

```text
if the factorial/Gaussian carrier supplies the diagonal (2) with its
ancestry semantics, then response selection creates no further obstruction
on the complete 94-start generated closure.                         (4)
```

What remains open is exactly the antecedent.  Neither THM-3286 nor the
product-Gamma response bank contains a comultiplication of tagged controller
states.  A future positive argument must derive (2) from an actual
coefficient, polarization, or proof-obligation composition law and show
that the target predicate survives the duplication.  Merely calling a
persistent face tag a recyclable token would repeat the original type error.

The next cheap hostile is the `31` row-specific full-domain blockers: test a
general copy-on-block rule there, but keep it separate from the four-state
extension.  The next substantive positive test is algebraic rather than
graph-theoretic: construct or refute a diagonal map on face-tagged
product-Gamma response obligations.  The LRC warning that copies are not
common atoms is useful procedural hostility here, not a mathematical bridge.

## Reproduction

Run

```text
python3 04-computation/fc3_three_support_face_f12_release_copy_rebind_gate_scout_20260803.py
python3 -O 04-computation/fc3_three_support_face_f12_release_copy_rebind_gate_scout_20260803.py
```

and compare LF-normalized bytes with the stored output.  The source and
output SHA-256 hashes are respectively

```text
90387770e45f2b560cdebc547ac295cb682935445fcbbbf0a7484f6a0a6e7200
6ba33a9f578a08a63b06d81b994818dc9c18f64e3d4cfce971ef5b8eee3b7357
```

The script reconstructs all `8,397` face/state responses, closes every
relevant row-labelled DAG, contains no assertion node or floating literal,
and pins its inherited source, transcript, and theorem dependencies.
