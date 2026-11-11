# An exact ordinal source compiler, paid balanced returns, and the unit budget

**PROVED:** the seven-coordinate ordinal law, its arbitrary-depth compiler,
positive-defect budget, unmarked-catalyst classification, and balanced return
families below. **INHERITED:** the selected-tail formula and the Session261
six-vertex module trade. **FINITE-EXACT:** twelve small literal compiler controls.
**OPEN:** uniform unit rearming, the missing 9-grid levels, and Chan–Pak.
No new numerical spectrum coverage or historical-priority claim is made.

## 1. Inheritance and the actual constructive problem

Closest mechanisms: Session261 `constructors/MARKED_RESELECTION.md`,
`MODULE_REBASE.md`, and `NINE_OFFSET_EXTENSION.md`; and
`multiplicand_saturation_chanpak_20261111/chanpak/OFFSET_TRANSPORT.md`.
The hostile is the Fermat fan: the correct offset alone does not keep the
odd part of U equal to 9. The corrected near miss is treating equality of
current responses as equality of rooted continuation menus. The least-used
sidecars here are BOTH endpoint responses and ordinal-component order.

We accept the current universal domination-defect inequality as inherited.
None of our operations prescribes its defect independently of H: every
state below is emitted as a literal graph, and its base H is retained.
Our live board is: base cost / offset / start response / end response /
ordinal composition / hidden rooted menu. Anchor: uniform source supply.
Niche: a genuinely closed small interface. Wildcard: palindromic balance.

## 2. A closed seven-coordinate interface

For a nonempty marked tournament (T,O), let H count Hamilton paths, Q count
ordered spanning two-path covers with either slot allowed empty, S and E
count Hamilton paths starting and ending in O, U count two-path covers
whose second path starts in O (first may be empty), V count those whose
first path ends in O (second may be empty), and W require both conditions.
Set D=H-W. All seven coordinates are nonnegative integers, even though D
can be negative. All unselected sources have S=E=U=V=W=0.

For A⇒B, with inherited marks in both blocks, the exact law is

    H = h_A h_B,                  Q = q_A q_B,
    S = s_A h_B,                  E = h_A e_B,
    U = u_A q_B + h_A u_B,
    V = v_A h_B + q_A v_B,
    W = w_A h_B + h_A w_B + u_A v_B + e_A s_B.       (1)

Proof: every path traverses its A portion before its B portion. Independent
ordered two-slot covers give the Q product, including empty portions.
For U, its second path starts either in A, giving u_A q_B, or in B,
forcing all A into the first path, giving h_A u_B. Converse gives V.
For W classify the first end and second start by their blocks: (A,A),
(B,B), (B,A), (A,B), respectively giving the four displayed terms.
These cases are disjoint and exhaustive; notably the last term must not
be dropped. It is the boundary coupling of the two endpoint statistics.

Thus seven numbers suffice under **every future ordinal word**, at every
length. This is a genuinely uniform restricted-grammar controller, unlike
agreement through three sampled clone lengths. It does NOT suffice for
arbitrary module substitution: that still requires the richer cover tensor.
Associativity follows either from the literal operation or its bijection.

Concrete initial profiles (H,Q,S,E,U,V,W):

    unmarked singleton: (1,2,0,0,0,0,0)
    marked singleton:   (1,2,1,1,1,1,0)
    unmarked C3:        (3,12,0,0,0,0,0)
    Session261 unit32:  (9,84,4,4,32,32,8)
    Session261 offset9: (9,84,1,1,9,15,0).

The last two literal rows/marks are (2,4,9,19,7), {0,3}, and
(2,28,25,17,1), {0}. These are actual initial cases, not abstract states.
`ordinal_interfaces.py` emits the graph and computes the profile for any
nonempty word of actual marked atoms.

## 3. A positive budget that survives arbitrary word length

Normalize each block's lowercase responses by its H, writing
q_i=Q_i/H_i, u_i=U_i/H_i, etc. Expanding (1) gives

    W/H = sum_i w_i
        + sum_(i=1)^(k-1) e_i s_(i+1)
        + sum_(j<i) u_j v_i product_(j<l<i) q_l.       (2)

A direct proof labels the blocks containing the two required endpoints.
If both lie in one block, all other blocks use one path. If the first end
precedes the second start, their blocks MUST be adjacent, since no path
could cover an intervening block. If the first end follows the second
start, intervening blocks may split freely across the two paths. This gives
exactly (2), with empty products equal to 1.

**Positive-defect budget.** An ordinal word has D>0 iff the right-hand side
of (2) is strictly less than 1. In particular it cannot contain two marked
unit-seed blocks each with H>=3: each contributes w_i=1-1/H_i>=2/3 before
any cross term is paid. This remains true with arbitrarily many unmarked
catalysts between them. No local search cutoff is involved.

Example: two copies of unit32 produce H=81, W=1184, D=-1103, U=V=2976.
Already the two within-block terms exceed the budget; the cross terms
make the failure much worse. Positive source supply is not closed under
blind ordinal multiplication of successful unit seeds.

## 4. Exact classification of unmarked ordinal balancing catalysts

Let B carry all marks and put T=L⇒B⇒R, with L,R unmarked, allowing empty
sides with H=Q=1. Formula (1) gives

    D(T)=H(L)H(R)D(B),
    U(T)=H(L)Q(R)U(B),
    V(T)=Q(L)H(R)V(B).                               (3)

If D(B)=D(T)=1, then H(L)=H(R)=1. A tournament has a unique Hamilton
path iff it is transitive: a nontrivial strong component has a Hamilton
cycle and hence at least three Hamilton paths, while strong-component
counts multiply. Thus L,R are transitive, including their empty cases.
If their orders are l,r, balancing U(T)=V(T) is exactly

    U(B)/V(B)=2^(l-r).                               (4)

Consequently a positive-response unit seed can be balanced by **arbitrary
unmarked ordinal catalysts iff U/V is an integral power of 2**. This is
both necessary and sufficient, and gives the actual initial padding sizes.

Applied to the inherited clone family U=16*2^m,
V=4*3^m+9*2^m+2: at m=1 the ratio is 1. For m>=2, V has valuation 1
and odd part greater than 1, whereas U is a power of 2. Thus NO unmarked
ordinal catalyst can restore balance while retaining unit defect. This
strictly generalizes merely testing source/sink padding. It does not rule
out marked catalysts, nonordinal attachments, or internal rebase operations.

## 5. Positive all-order balanced source supply with paid offset

Let A and B be ANY nonempty unmarked tournaments with equal H=h and Q=q.
Then the single-mark sandwich A⇒{z}⇒B has

    H=D=h^2,      W=0,      U=V=hq.                  (5)

This is a balanced physical source, so full converse preserves its exact
(D,U), and every nonempty selected tail gives h^2+2^t hq. No d is chosen
independently of base cost. Source returns can be abundant at paid offset
although the same mechanism cannot preserve D=1 unless h=1.

For arbitrary k>=0 take left and right words of k atoms chosen independently
from the inherited pair

    A=(2,4,9,19,39,15),  B=(2,4,57,3,11,27).

Their complete cover profiles agree, with H=17 and C2=190, hence Q=224.
Each sandwich has the exact uniform family

    H=D=289^k,       U=V=3808^k,       W=0.            (6)

There are 2k labelled module positions independently available for the
physical exchange A↔B. Each exchange preserves every current marked module
context, by the inherited complete-profile theorem. Yet the rooted clone
menus differ: A's vertex1 has H(m)=9*2^m-1, B's vertex0 has H(m)=8*2^m+1.
Thus every k supplies actual interchangeable modules with different next
operations, not just copied scalar labels. The initial case k=0 is the
marked singleton with D=U=V=1; k=1 has 13 vertices and D289,U=V3808.
A nonempty tail constructs 289^k+2^t3808^k for every k>=0,t>=1.

This return is reusable under arbitrary ordinal contexts and fresh tail
lengths. It is NOT a claim that after arbitrary rooted inflation the same
trade can still be applied. Restoring that stronger guard remains open.
Balanced transitive padding on both sides multiplies U,V by 2^b without
changing D, giving another independent physical parameter b>=0.

A smaller initial sandwich C3⇒{z}⇒C3 has D9,U=V36. It supplies the entire
9-grid row s=3 using tail length t=r+1, including r=1. It is less compact
than Session261's asymmetric source U9, but has a new exact balance guard.
Among equal-H single-mark sandwiches, this is the ONLY level s>=1 with
D=1+2^s: (h-1)(h+1)=2^s forces h=3,s=3. Thus this symmetry mechanism cannot
silently be extrapolated to the missing 9-grid levels.

## 6. Sharp constructive reduction and source connections

Given any finite or computably supplied alphabet of actual marked atoms,
(1), (2), and the tail law give an exact word-level construction problem:
find a word and t>=1 with

    H(word)-W(word)+2^t U(word)=N.

Every accepted word compiles to a tournament with exactly N paths; every
word in THIS grammar is covered by the law. There is no claim the word
search always terminates or covers all N. The universal missing gate is
now explicit: supply an alphabet/word-selection rule hitting the required
numbers while respecting the positive budget and the future return guard.
The concrete alphabet above initializes the controller, but is not asserted
to meet that gate. Target anchoring, cost, and initial cases are separate
obligations, not consequences of the scalar recurrence.

Recorded random seed 262913 sampled three full-tree Markdown paths from
canon/results; paths and pool size are in `random_sources.json`.

1. Session214 `surgery/PROOF.md`: SCC product and its 7477 hostile suggest
   retaining ordinal order and keeping the grammar-specific completeness
   claim separate from unrestricted coverage. Our map is SCC product →
   selected two-slot product; H multiplicativity survives, scalar closure
   does not; the seven coordinates are the missing sidecar.
2. `two_entry_funnels_20261103/TWO_ENTRY.md`: its first-entry decomposition
   separates a controlled entry count from later SCC factors. Applying the
   same separation to TWO path slots gives precisely U=u_A q_B+h_A u_B,
   rather than incorrectly multiplying U by H. Marks and empty slots are
   retained; unique-cycle recognition is not imported.
3. THM-4930 `eternal-colours-which-odd-path-counts-are-ever-self-converse.md`:
   palindromic condensation physically pays a square factor in H. Here the
   single central mark converts that square into D=h^2 and supplies exact
   endpoint balance. The transfer preserves the anti-isomorphism when
   B=A^op; equal H,Q is enough for response balance without an actual
   anti-isomorphism. This distinction allows the hidden-menu A/B trade.

## Reproduce

    python3 -B 05-knowledge/results/session262_interface_supply_20261111/constructor/ordinal_interfaces.py

The receipt has twelve literal small controls, including the negative-defect
and balanced initial cases. The all-order statements follow from the
bijections and inherited trade, not extrapolation of those controls.

The same receipt also checks the actual balanced-return initial cases k=0,1
and their t=1 tail outputs, respectively 3 and 7905. Emit an arbitrary member
without exponentially recounting its whole graph:

    python3 -B 05-knowledge/results/session262_interface_supply_20261111/constructor/ordinal_interfaces.py --k 4 --tail 3 --padding 2

The Python API `balanced_trade_source(k, choices, padding)` accepts all 2k
independent binary A/B choices and returns literal rows, the central mark,
and the exact compiled seven-coordinate profile. It rejects invalid lengths
and tail parameters rather than introducing absent-module extrapolations.
