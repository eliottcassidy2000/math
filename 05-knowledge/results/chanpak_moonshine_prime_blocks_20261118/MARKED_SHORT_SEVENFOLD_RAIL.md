# The short sevenfold return is physical as a marked arc load

The lower count q_i in PAIRED_SEVENFOLD_SUM.md follows q→8q−3. We now
realize that recurrence at every height as the exact number of
Hamilton paths using a specified arc, not as the entire Hamilton
count. The upper witness 7q_i remains an ordinary tournament count
at every height by the earlier two-port construction.

Use the inherited retained-root binary word law (ANS_INVERSE.md).
A rooted tournament T with root a keeps h=H(T), e=number of full
paths ending at a, d=H(T−a), and the weighted start count
t=Σ_(a∈S) S_a(T[S])H(T−S), with H(empty)=1. Letter 0 appends a
sink; letter 1 appends a sink and reverses just the arc from
a to it. Both keep the root. Their physical updates are

    0: (h,e,t,d) -> (h,0,2t,d);
    1: (h,e,t,d) -> (h+t−e,d,2t−d,d).

Start with the arc 0→1 rooted at 0, state (1,0,2,1). Appending
011 twice gives state (87,1,101,1). Starting from there,
each additional 011 sends

    (h,1,t,1) -> (h+6t−2,1,8t−3,1).

Thus after i additional blocks t_i=q_i=(704·8^i+3)/7 and

    7h_i = 6q_i + 3 + 4i.

One further letter 1 creates a candidate copy of the arc load:
the new vertex w beats only a. Hamilton paths using w→a are
in bijection with terms of t_i. The suffix starts at a and
uses S; all omitted vertices form the prefix immediately
before w, which beats none of those vertices. This last prefix
is allowed to be empty. The other paths contribute h_i−1,
so the entire output has h_i−1+q_i, not q_i. This describes
exactly what a reader would need to remove or redirect.

## A scoped test of readers

For the first three physical states (orders 8,11,14), the replay
tries every orientation of the arcs to one new vertex and checks
the exact count using named endpoints of two-path covers. None
of these 2^8+2^11+2^14 attachments hits the corresponding q_i.
The count formula for a proposed apex with outgoing set U is:
full paths starting/ending at the apex, plus each ordered
nonempty two-cover of the old tournament with its first end
outside U and second start in U. This is the inherited apex
quadratic. Counts are obtained via integer subset DP. The
script also verifies the marked-arc counts by counting with
and without the specified arc.

This is not an all-height impossibility for other readers.
It separates getting the p+7p−3 recurrence physically into a
marked coordinate from converting that coordinate to total
H. The finite failure indicates that simply appending one
reader vertex to these particular sources is insufficient,
even though independent unmarked lower witnesses for those
targets are known. Changing old contacts, rerouting the
marked coordinate, and using a larger reader remain open.

Run python3 physical_short_rail_replay.py.
