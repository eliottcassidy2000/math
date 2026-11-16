# Prime mantissas and the blocked repeatable port

The sequence 27→53→105→209→417 is A_e=1+13·2^e for e≥1.
One step is h→2h−1, very different from ordinal multiplication by
3. [Earlier exact work](../inverse_fences_landing_20261116/INVERSE.md)
gave 27→53 by a marked append, ruled out every literal induced
53→105, and gave a cost-optimal repair with one old arc reversal.

## Continued literal witnesses

Start at H27 with rows (62,124,120,112,32,64,17), where bit j
of row i indicates i→j. When adding one vertex, give its outgoing
mask in the old vertex labels and complement its contacts in
the other rows. Apply the following steps to a *single* evolving
witness, retaining all vertices. A single reversal suffices at
each step after the first, through 833; the displayed move to
1665 uses two. Only for 53→105 is a reversal proved necessary
over all choices of tournament witness.

| New H | Old arc to reverse | New outgoing mask |
|--:|:--|--:|
| 53 | none | 16 |
| 105 | 0–1 | 181 |
| 209 | 2–3 | 208 |
| 417 | 3–8 | 656 |
| 833 | 4–7 | 692 |
| 1665 | 0–6 and 2–5 | 4084 |

A different H105 seed, rows (6,12,24,49,99,71,15,383,127),
continues to 209,417,833,1665 by successive arcs
0–1, 1–3, 1–10, 3–5 and masks 31,2,1088,3087, respectively.
[The replay](check_ray.py) recomputes all counts by subset
path DP. These are exact finite trajectories and not a
uniform rule for selecting an arc on later heights.

## No finite one-apex fixed-run library for any prime mantissa

For a one-apex tape use the state (h,A,p) from the earlier note.
Let L,R be fixed binary words, with A the binary value of L, p its
last bit (0 if empty), and B=Σ_i(1−R_i)2^i in tape order, zero-based.
For s≥1,
`H(L0^sR)=F+G·2^s` with
`F=H(L)+H(R)−1−A+p` and `G=A(1+B)`.
The inserted zeros contribute A(2^s−1)+p, and the old left
vertices contribute A·2^s·B to the zeros of R. Other increments
are just the two standalone counts minus their duplicated baseline.

To hit infinitely many points of 1+q·2^e, with q a fixed odd
prime, G>0 is necessary. The equality forces 2^(e−s) to
converge to G/q, hence be eventually constant by discreteness
of the powers of two in a compact positive interval. Thus F=1,
and exactly one of A and 1+B has odd part q, while the other
is a power of two.

If A is a power of two, then L is 10^u after ignoring leading
zeros, and H(L)−A+p=1. So F=1 forces the right word to be
acyclic, hence to have the form 0^v1^w: any 1 before a 0 makes
a directed triangle. This makes 1+B a power of two as well,
a contradiction. In the other case 1+B is a power of two, which
likewise gives H(R)=1. Deleting trailing zeros from L leaves
H(L)−A+p unchanged by the recurrence, and deleting leading
zeros has no effect. After these deletions L ends in 1 and
has binary value q>1. The stronger inequality H+p≤A+1 is preserved by both
letter updates after the first 1, and the final
1 makes the inequality strict enough to give H≤q−2: just
before the final digit A_0≥1, H≤A_0+1<2A_0+1=q and both H,q
are odd. But F=1 now requires H=q. Contradiction. The converse
tournament under reversing the tape gives the reverse-complement
word, so a pump of ones cannot help either. A finite library
therefore has only finitely many successes. Arbitrary words with
target-changing internal structure are not covered by this theorem.

## The pair 27 and 63

The 27 witness made from three separate triangles attains
H=3^τ and has minimum positive induced gap 18 by the
[factored-gap formula](../triangular_backbone_chanpak_20261116/FACTORED_GAPS.md).
Yet it cannot sit as an induced subtournament of ANY H63 witness.
Every H63 witness has all nontrivial count in one strong component:
a nontrivial factoring of 63 would need H=7 or H=21, both impossible.
Busch's strong order bound puts this H63 component at order at
most 8 (63³<5⁸). All three old directed triangles must lie in it,
but they demand nine distinct vertices. This obstruction is specific
to the factored H27 witness; it says nothing about a different H27 core.

## Further experiment, not an all-order bound

At H209, the order-8 atlas and all new-vertex masks yield 474
order-9 H209 hits (repetitions allowed), from 1926 source rows
with count at most 209. None are ready for indefinitely repeating
the fixed-root near-sink append. Neutral singleton padding makes
this exclude all such ready roots up to order 9, not more.
In the narrower order-10 search, take either of the two explicit
H105 rows in near209.cpp, optionally reverse one old arc, then
insert a new vertex in all ways. There are 13 hits of H209;
again none are ready. These are controls on the next construction,
not a proof that H209 has no ready root at all.
