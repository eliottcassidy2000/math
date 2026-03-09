# THM-106: Free Cycle Bridge Theorem

**Status:** PROVED (algebraic)
**Filed by:** kind-pasteur-2026-03-08-S43

## Statement

Let C = (a→b→c→a) be a **free** 3-cycle in a tournament T, and let v be any
external vertex (v ∉ {a,b,c}). Then v creates a "bridge" 3-cycle B_v that
shares a directed edge with C.

Precisely: B_v is a directed 3-cycle on {v} ∪ {two vertices from C} that
shares one of the three directed edges {a→b, b→c, c→a} with C.

## Proof

Since C is free, vertex v neither dominates all of {a,b,c} nor is dominated
by all of {a,b,c}. Therefore:
- v does NOT beat all three: at least one of a,b,c beats v.
- v does NOT lose to all three: v beats at least one of a,b,c.

Let out(v) = {x ∈ {a,b,c} : v→x} and in(v) = {x ∈ {a,b,c} : x→v}.
Then |out(v)| ∈ {1, 2} (not 0 or 3 by freeness).

### Case |out(v)| = 1, |in(v)| = 2

Say v→x, y→v, w→v where {x,y,w} = {a,b,c}.
The predecessor of x in C beats v (since x has one predecessor in C and
we need to find a shared edge).

Subcases by which vertex v beats:
- **v→a:** Then b→v and c→v. Path: c→a (cycle edge) with v→a and c→v gives
  c→v→a. But we need c→v (yes) and v→a (yes) and check c→a: yes (cycle edge).
  So (c,v,a) is a 2-path c→v→a with c→a: TT triple. But we need a 3-cycle.
  Check: v→a→b→...→v? v→a (yes), a→b (yes), b→v (yes). So (v,a,b) is a
  3-cycle v→a→b→v. This shares directed edge a→b with C.

- **v→b:** Then a→v and c→v. Check: v→b→c→v? v→b (yes), b→c (yes), c→v (yes).
  3-cycle (v,b,c) sharing edge b→c with C.

- **v→c:** Then a→v and b→v. Check: v→c→a→v? v→c (yes), c→a (yes), a→v (yes).
  3-cycle (v,c,a) sharing edge c→a with C.

### Case |out(v)| = 2, |in(v)| = 1

Say v→x, v→y, w→v where w is the one vertex beating v.
The successor of w in C is beaten by v.

Subcases:
- **a→v, v→b, v→c:** Check a→v→b with a→b (cycle edge): (a,v,b) shares edge
  a→b? Need 3-cycle. a→v→b→...→a? a→v (yes), v→b (yes), b→a? Only if b→a.
  But a→b in the cycle. Try: v→c→a→v? v→c (yes), c→a (yes), a→v (yes).
  3-cycle (v,c,a) sharing c→a with C.

- **b→v, v→a, v→c:** v→c→a with... Try: v→a→b→v? v→a (yes), a→b (yes),
  b→v (yes). 3-cycle (v,a,b) sharing a→b with C.

- **c→v, v→a, v→b:** v→b→c→v? v→b (yes), b→c (yes), c→v (yes).
  3-cycle (v,b,c) sharing b→c with C.

In every case, v creates a 3-cycle sharing a directed edge with C. QED.

## Bridge Edge Pattern

The specific shared edge depends on v's orientation:
| v beats | v loses to | Bridge shares |
|---------|-----------|---------------|
| a only  | b, c      | a→b          |
| b only  | a, c      | b→c          |
| c only  | a, b      | c→a          |
| b, c    | a         | c→a          |
| a, c    | b         | a→b          |
| a, b    | c         | b→c          |

Note: beating 1 vertex or beating 2 vertices that are "opposite" in the cycle
produce the SAME shared edge.

## Consequence

Every external vertex is connected to a free cycle's component in the 3-cycle
adjacency graph. This is the key lemma for THM-107 (at most 1 free component).

## See Also
- THM-105 (dominant vertex forcing)
- THM-107 (at most 1 free component, uses this)
