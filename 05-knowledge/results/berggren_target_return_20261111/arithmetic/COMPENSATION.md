# The known H13 actuator cannot fit the growing-defect overshoot

**PROVED IN NOTE, conditional only on the sibling lane's displayed source
formulas.** The H13 actuator calculation and baseline obstruction are exact.
**FINITE-EXACT:** direct subset-DP verification of the base two graphs and
32 dyadic formula controls. This rejects the specified actuator family as a
literal installation; other source switches and cheaper actuators remain open.

The repeatable-return lane supplies actual indegree-one sources with

    h=891D+215, F=1083D+262, b=F-h=192D+47,
    D=2^k, k>=1.

Its ordinary two-birth step has count

    h+2F = 3057D+739,

whereas the desired +2 target is

    3h+2 = 2673D+647.

The debt is therefore

    w=2(b-1)=384D+92.

The inherited H13 pure-H switch width is `w(A)=3A-(A mod2)+2`. Since this
debt is2 modulo6, its unique address in that menu is the even integer

    A=(w-2)/3=128D+30.

For every k>=1 the binary recurrence gives

    g(128D+30)=68D+17.

For example, strip the final zero, then four terminal ones, from
`2^(k+7)+30`; the remaining power of two evaluates to `2^(k+2)+1`.
The first zero contributes `64D+16`, giving the displayed sum.

The lower and upper H13 sources at the smallest right address B=0 have

    H-=6900D+1649,
    H+=7284D+1741,
    (E,U,d)=(7,13312D+3165,13).

Their count difference is exactly w. Every B>0 strictly increases the lower
baseline: the increment is

    13(g(B)-1)+(104A+32)B+6(B mod2) > 0.

Thus the smallest available actuator already exceeds the overshoot by

    H--(h+2F)=3843D+910 > 0.

At the first D=2 this is the actual pure-H pair

    (15449,7,29789,13), (16309,7,29789,13),

whose gap860 is the correct debt, but whose lower baseline15449 exceeds
both6853 and5993. The packet's independent full subset dynamic program
checks these two literal 15-vertex sources. Their retained root has
indegree k+6 in the general family, so it is not an indegree-one source.

This excludes more than whole-source equality. Hamilton path count is
monotone under taking induced subtournaments: given any directed Hamilton
path of T-v, insert v before the first vertex if possible, after the last
if possible, or at a consecutive transition `u->v->w`. At least one such
position exists. Choose a deterministic position; deleting v recovers the
old path, so the resulting map is injective. Iterate vertex deletions.

Consequently no member of this H13 actuator family can occur as an induced
subtournament of the overshoot source, since its H is already greater.
Adding vertices to such an actuator cannot produce that source either.
The same conclusion holds for the smaller desired target. It remains
possible that a different actuator has the same gap at a smaller baseline,
or that a graph rewrite uses a different context with a separately proved
response. The present count difference alone does not supply either map.

The generic distinction is now concrete: a selectable even response width
is insufficient for debt payment. A literal local replacement also needs
a compatible baseline, actual embedded source, retained root, and matching
context response. Here the baseline itself rules out the known installation
before any larger embedding search is justified.
