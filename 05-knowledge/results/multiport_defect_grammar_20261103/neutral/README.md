# The failed inner collar cannot be neutralized within three suffix reversals

2026-11-03. **FINITE-EXACT**, independently replayed. This is a complete
declared local menu, not an all-order impossibility theorem.

The actual H1537279 grandchild in the inherited
[`zero_ready_strong_switch` packet](../../zero_ready_strong_switch_20261103/README.md)
has outer root39, sole inport40, transitive prefix50->51->55->56 and a
strong14-vertex suffix B. Its suffix has H(B)=80319 and no zero-ready
vertex. The minimum defect E_r(B,2)-H(B) is2960 at physical32. This is the
specific hostile to pretending the inner zero-ready collar survives every
return.

We ask whether a **count-preserving reset** can restore any inner zero-ready
mark, allowing:

- zero through three arbitrary reversals among B's91 edges;
- every actual B vertex as the new outer inport;
- every transitive outer-prefix length at least4;
- any zero-ready vertex in the resulting B, not just32.

The suffix must remain strong. The exact target equation is

    H(B') + 2^s E_w(B',2) = 1537279,   s>=4.

All old names remain present; a new outer inport would require its two actual
root incidences to change, in addition to the budget of suffix reversals.
The arithmetic search permits that extra freedom, so its negative conclusion
does not rely on an uncompiled native guard. No proposed reset is claimed
without an actual witness.

## Complete result

|Exactly this many suffix reversals|Examined|Strong B'|B' with some zero-ready mark|
|---:|---:|---:|---:|
|0|1|1|0|
|1|91|89|2|
|2|4095|3938|130|
|3|121485|115389|4192|
|Total|125672|119417|4324|

None of the4324 eligible suffixes, at any actual outer inport and any allowed
prefix length, solves the target equation. The nearest numerical distance
within that eligible menu is80; proximity is not a realization.

A preliminary two-reversal menu had nearest distance1446. Forcing reversal
of the specific damaging pair{22,32}, with at most two other reversals,
checks4096 graphs and also finds no reset. These earlier menus are recorded
for provenance; the complete three-reversal result supersedes their scope.

## Why the filter and exponent cutoff are exact

For a strong tournament B', the inherited zero-ready classification says r
is zero-ready iff it has one in-neighbour w and w is not in the first strong
component of B'-r. A tournament's SCCs form a directed chain, so the latter
condition is equivalent to w failing to reach every vertex of B'-r.
The discovery program uses this graph criterion before expensive counting.
Its proof is in the inherited
[ZERO_READY.md](../../triangle_return_170807_20261103/ZERO_READY.md).

For each selected suffix, direct endpoint DP independently confirms exactly
the same set of zero-ready roots. Since E_w(B',2)>0, the target expression
strictly increases with s. Testing can stop at the first overshoot. In this
universe s>=21 already exceeds the target even with h=F=1, so the program's
integer safety bound cannot omit a solution.

## Independent audit and controls

`search.cpp` uses the inherited forward subset recurrence from the previous
packet. `audit.py` independently regenerates every edit set, uses Warshall
closure for the SCC filter, and sends all4324 survivors to a backward deletion
recurrence. It agrees on every count, every zero-ready root, and the absence
of a solution. The order14 arithmetic has the a priori bound
P_B(2)<=15*14!<2^41; the independent counter additionally checks overflow.

Positive control: the actual preceding H512425 state's suffix has
H=E_32(2)=46073. Hostile control: the actual target suffix has
H=80319,E_32(2)=83279, and every other root has an even larger defect.

Run from the repository root:

```text
python 05-knowledge/results/multiport_defect_grammar_20261103/neutral/scan.py --depth 3
python 05-knowledge/results/multiport_defect_grammar_20261103/neutral/audit.py
python -O 05-knowledge/results/multiport_defect_grammar_20261103/neutral/audit.py
```

The mathematical stopping reason is now sharper: the damaged collar is not
repaired for free by another three-suffix-edge neutral reset, even with a new
outer mark and longer prefix. A grammar must retain its nonzero defect, use
a larger repair, change more of the attachment, or admit another structural
chart. This says nothing against an unrestricted such repair or global
spectrum completeness.
