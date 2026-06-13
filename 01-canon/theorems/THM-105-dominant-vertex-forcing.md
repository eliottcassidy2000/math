# THM-105: Dominant Vertex Forcing

**Status:** PROVED (algebraic)
**Filed by:** kind-pasteur-2026-03-08-S43

## Statement

If a directed 3-cycle C = {a,b,c} in a tournament T is **dominated** — meaning
some external vertex d satisfies either d→a, d→b, d→c (d dominates all) or
a→d, b→d, c→d (d is dominated by all) — then for any TT-cocycle z ∈ Z_1:

  sum_C(z) = 0

where sum_C(z) = z(a,b) + z(b,c) + z(c,a).

## Proof

### Case 1: d dominates all of C (d→a, d→b, d→c)

The cycle is a→b→c→a. Since d→a and d→b, we have the directed 2-path d→a→b.
Since also d→b, this is a transitive triple (d,a,b).
By the TT-cocycle condition: z(d,a) + z(a,b) = z(d,b).

Similarly d→b→c with d→c gives TT(d,b,c): z(d,b) + z(b,c) = z(d,c).
And d→a with a→c... but we need d→c→a? No: d→c and c→a, so d→c→a is a
2-path, and d→a makes it TT(d,c,a): z(d,c) + z(c,a) = z(d,a).

Now sum the three TT equations:
  z(d,a) + z(a,b) = z(d,b)
  z(d,b) + z(b,c) = z(d,c)
  z(d,c) + z(c,a) = z(d,a)

Adding all three:
  [z(d,a) + z(d,b) + z(d,c)] + [z(a,b) + z(b,c) + z(c,a)] = [z(d,b) + z(d,c) + z(d,a)]

The bracketed terms on both sides cancel, leaving:
  z(a,b) + z(b,c) + z(c,a) = 0

### Case 2: d is dominated by all of C (a→d, b→d, c→d)

The cycle is a→b→c→a. Since a→b→d with a→d, TT(a,b,d): z(a,b) + z(b,d) = z(a,d).
Since b→c→d with b→d, TT(b,c,d): z(b,c) + z(c,d) = z(b,d).
Since c→a→d with c→d, TT(c,a,d): z(c,a) + z(a,d) = z(c,d).

Adding all three:
  [z(a,d) + z(b,d) + z(c,d)] + [z(a,b) + z(b,c) + z(c,a)] = [z(a,d) + z(b,d) + z(c,d)]

Again: z(a,b) + z(b,c) + z(c,a) = 0. QED.

## Consequence

Dominated cycles contribute nothing to b_1. Only **free** (non-dominated) cycles
can carry nonzero TT-cocycle sums, so b_1 depends only on free components of
the 3-cycle adjacency graph.

## See Also
- THM-104 (cycle sum equality)
- THM-106 (free cycle bridge theorem)
- THM-107 (at most 1 free component)
