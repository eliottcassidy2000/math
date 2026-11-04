# Source-specific continuations and an empty ready fibre

Use H for ordinary directed Hamilton vertex words. A b-ready marked source
has E_r(2)-H=b, where E_r(2) counts ordered spanning two-path covers with
the first path nonempty and ending at r; the other may be empty. The
inherited reset adds x beaten only by r, then y beaten only by x, and gives
H'=3H+2b. These are actual incidences, not formal scalar commands.

## R1. Complete1-ready fibres at the four requested values

**PROVED reduction + FINITE-EXACT complete grammar.** The inherited
[READY §1](../two_entry_funnels_20261103/READY.md) identifies every strong
1-ready source as a newborn over a strong unique-start parent. If the
parent has h paths and its marked two-bin count is F, the child has
H=h+F-1. Cutting any Hamilton word immediately after the parent mark
injects h words into these covers, so F>=h and h<=(H+1)/2<=28471.

[Busch, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf)
gives the sharp strong minimum f(21)=46875, and f is increasing. Thus the
parent has at most20 vertices. This published theorem supplies the cutoff;
the unique-start rank/unrank and graded path-bin recurrence are inherited
from [FUNNELS](../polarized_path_calculus_20261103/FUNNELS.md) and
[the previous compiler](../return18979_palindrome_20261103/funnel_bank.cpp).

The complete parent grammar has u3=1,u4=3,u_n=3u_(n-1)-u_(n-2), giving
24,157,816 classes through20. `one_start_bank.json` gives precisely:

|Child H|Parent order and rank|Parent h,F|Child order|
|---:|---|---|---:|
|56933|15,10293|4783,52151|16|
|56933|15,13140|2961,53973|16|
|56935|none|none|none|
|56939|16,11421|11573,45367|17|
|56941|15,10996|5557,51385|16|

Independent forward endpoint DP checks every grammar class through9,
first/middle/last rank at10..20, and all positive ready targets:642 parent
queries plus12 positive targets. The exhaustive compiler uses checked
wide intermediates; absence is a finite computation, not a guessed pattern.
[OBSTRUCTION O1](../ready_return_pentagon_20261103/OBSTRUCTION.md) extends
the empty56935 fibre to weak sources: for H>2, weak1-ready sources are
strong1-ready cores followed by transitive tails, which preserve H.

Therefore **no H-preserving return by any operations can make56935
1-ready**. The value56935 itself is already realized by the saved source.
This does not exclude170807, nor a3H+2 operation using a different input law.
The two-entry bank here is a sufficient bounded bank (4005 core/port objects
through10 and tails to parent20), not a complete2-ready enumeration.

## R2. Seven literal moving-mark programs

**FINITE-EXACT.** Sources are the four `next_rows` in
[the previous returns](../return18979_palindrome_20261103/returns.json).
Each has an old live mark x and a latest vertex y beaten only by x. Move
the live mark to y, or, when one more retained vertex is required, create
a new vertex beaten only by x and mark it. Keep x as a named in-port and6
as the spare. Alignment is a search heuristic; it is not a permitted reset.
The actual program is the recorded sequence of operations below.

|H|offset|old mark -> new mark|permanent birth|old deletions|native complements|child H|
|---:|---:|---|---|---|---:|---:|
|56933|2|51 ->52|no|4,5|28|170801|
|56933|4|51 ->52|no|5,22|20|170803|
|56935|4|50 ->51|no|none|34|170809|
|56939|2|57 ->59|yes|none|26|170819|
|56939|4|57 ->58|no|13|26|170821|
|56941|2|56 ->58|yes|none|32|170825|
|56941|4|56 ->58|yes|none|28|170827|

Each old-edge edit uses the inherited one-helper compiler: choose a legal
Hamilton path with a->b consecutive, neither b nor the live mark last;
birth a helper beaten by every old vertex except a; complement off the
path with the helper last; re-mark to the path with helper between a,b;
complement again; delete the helper. This flips only a-b and keeps the
live mark's actual entry/exit frame common between the two paths.
Every path is checked in its **successive actual relation**.

At most one temporary helper is present at a time. Every completed edge
repair is strong. Intermediate counts need not equal H. Old-name deletions
are explicit; these are not fixed-order metagraph paths using only complement
and re-marking. After repair the exact profile satisfies E_r(2)-H=b;
two actual births give the listed child. Source, target and child relations
and all physical names are stored in `moving_returns.json`.

`audit_moving.json` independently replays all194 native guards, all births
and deletions, and verifies forward endpoint counts. It retains every
ordered two-bin pair response on the three physical ports at targets and
children, alongside their **full relation**. Pairwise responses alone are
not claimed sufficient for later joint operations. F49/F343 codecs and five
Q4 presentations are recalculated from each actual golden vector; ancestral
clock receipts remain frozen. No field equality authorizes a source reset.

The12-vertex hostile has equal H462915 and inherited automorphism orders
81,27, but the same named cut256 gives1061775,1395279. The literal count
distinction is independently replayed; group orders are inherited.

## R3. Why the easiest escapes do not iterate

**PROVED uniform identity.** Let T be any b-ready source with h paths.
Form its standard two-birth child U with H=3h+2b, then replace its newest
vertex y by a uniform transitive block of m nonempty clones. For every m>=1,

    H(U_m)=(h+b)2^m+h.                                      (R3)

Proof: y has only one external in-neighbour x. Its clones can occupy at
most two path runs. A single run contracts to a Hamilton path of U, giving
H possibilities. With two runs, the first must start the whole path. The
middle segment consists of a path in T ending at r, then x; after the
second run comes the remaining, possibly empty, path of T. These are
exactly E_r(2)=h+b ordered covers. Assigning the m ordered clones onto
two nonempty runs gives2^m-2 choices. Hence H+(h+b)(2^m-2), proving(R3).

The exhaustive menu consisting of a first indegree-one or outdegree-one
birth at any old vertex, followed by an arbitrary nontrivial birth cut,
has precisely three hits on each of the four sources. All three are the
same TT3 replacement with three placements of the old named y. Independent
matrix checks are in `audit_clones.json`. One unrestricted birth has no hit.

For fixed b, another3H(U_m)+2b within this same clone family would require

    2^m' = 3*2^m+2.

At m1 this gives m'=3. For m>=2 the right side has2-adic valuation1 but
is larger than2, so it is not a power of two. Thus these escapes provide
a uniform **one-time** continuation, not a reusable ternary return.

The boundary agrees with incoming
[Session175's fixed-library obstruction](../session175_tournament_spectrum_20261104/clone_obstruction/PROOF.md):
lengthening finitely many transitive inflation templates cannot cover a
positive proportion of the prime lane. Its stronger counting theorem is
inherited, not inferred merely from this one family.
