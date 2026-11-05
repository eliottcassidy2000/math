# Binary stationary blocks and logarithmic witness length

Session 195, 2026-11-05. **PROVED all-order structural bounds for the rooted grammar.** This is separate from the frozen two-gadget finite-spectrum packet. It does not prove that every admissible odd integer is attained.

## 1. Setup

Use the four physical coordinates `(h,e,R,d)` and operations from the [returning-interface theorem](../rooted-interfaces/PROOF.md), [combined grammar](../combined-growth/PROOF.md), and [two-gadget extension](../gadget-extension/PROOF.md):

    S(h,e,R,d) = (h-e+R, h-e, h+R, h),
    K(h,e,R,d) = (h-e+R, d, 2R-d, d),
    D(h,e,R,d) = (h, 0, 2R, d),
    M_G(h,e,R,d) = (gh-se+LR, Eh, Qh+BR, dG h)

for non-singleton rooted gadgets G. S is the root-moving singleton step; K retains the old root after adding a vertex beating only that root; D adds a universal sink. For G, `B=P_(G-b)(2)` and `dG=H(G-b)`, with both path slots allowed empty.

All physical states satisfy `1<=d<=h<=R` and `0<=e<=d`: deleting the final root from a path ending there injects those paths into the Hamiltonian paths of the root-deleted tournament. A non-singleton gadget has

    B >= 2dG >= 2.                                      (1)

Indeed the two covers `(P,empty)` and `(empty,P)` for each Hamiltonian path P of the nonempty graph G-b are distinct. Also `h' >= R` for S, K, and every M_G. The root-moving transfer gives `R' >= 2R` for M_G by (1).

The only initial state is `x0=(1,1,1,1)`. Consecutive initial K steps fix this state and can be removed from a construction word without changing any later coordinate or final count. The literal graphs may change, but the shorter word is another physical witness with the same count.

## 2. A new invariant

**Theorem 1.** Every reachable state other than x0 satisfies

    R >= 2d.                                            (2)

**Proof.** After S, `R'=h+R>=2h=2d'`. After any non-singleton M_G,

    R'=Qh+BR >= 2dG R >= 2dG h = 2d'.

D sends a state satisfying (2) to one with the same d and larger R. K sends such a state to `R'=2R-d>=3d>=2d`. At x0, K fixes the state, while D gives R'=2 and d'=1. These cases establish the invariant by induction. QED.

The invariant is stronger than the general physical inequality R>=d, and it uses this precise grammar rather than a claim about all rooted tournaments.

**Corollary 1.1 (strict Hamilton-count contraction in reverse).** At every noninitial reachable state, every non-D step satisfies

    2h' >= 2h+R >= 3h.                                 (2a)

Indeed `e<=d<=R/2`. For S and K, `h'=h+R-e>=h+R/2`. For a general gadget,

    h'=gh-se+LR >= gh+(L-s/2)R >= h+R/2,

because `L>=g>=s`, `g>=1`, and hence `L-s/2>=1/2`. Thus an inverse non-D step contracts h by a factor at most 2/3, except when its predecessor is x0. That exception is necessary: S and some transitive gadgets applied at x0 preserve h=1. Only D can preserve h once the state has left x0.

## 3. Every pair of steps doubles R

Let N denote any non-singleton root-moving gadget. At a noninitial reachable input, two successive operations always satisfy

    R_after_two >= 2 R_before.                          (3)

The following table gives lower bounds, with the row operation first and the column operation second. All h,R,d refer to the input before either step.

| First / second | S | K | D | N |
|---|---:|---:|---:|---:|
| S | 2R+h | 2R+h | 2R+2h | 2R+2h |
| K | 3R-d | 4R-3d | 4R-2d | 4R-2d |
| D | 2R+h | 4R-d | 4R | 4R |
| N | 3R | 3R | 4R | 4R |

For example, S gives `R1=R+h,d1=h,h1>=R`, so following it with S gives at least `2R+h`, and following it with K gives exactly `2R+h`. After K, `R1=2R-d,h1>=R,d1=d`. After N, `R1>=2R,h1>=R`, and Theorem 1 gives `d1<=R1/2`, so a following K gives at least `3R1/2>=3R`.

Every table entry is at least 2R because `R>=2d` and h is positive. At x0, after removing the initial fixed-state K loops, the first operation is S, D, or N. Those three table rows remain valid and at least 2R at x0 as well. Moreover any such first step has `R1>=2`.

Therefore any reduced construction word satisfies, at each step j,

    R_j >= 2^ceil(j/2).                                 (4)

This is an all-order inequality for every word over an arbitrary bank of actual nonempty rooted gadgets, not a numerical fit to small states.

## 4. Logarithmic witness theorem

**Theorem 2.** If a positive integer H is attained by this grammar, it has a physical construction with at most

    2 floor(log2 H) + 1                                 (5)

elementary steps. If the largest permitted gadget has k vertices, it has a witness tournament with at most

    1 + k (2 floor(log2 H) + 1)                         (6)

vertices. In the two-gadget extension, k=6 and the bound is `12 floor(log2 H)+7`.

**Proof.** Choose any word attaining H. Remove all initial K loops at x0. Remove its terminal D steps, which do not change h. If no steps remain, H=1 and the singleton suffices. Otherwise the last step is not D. Every preterminal R is at most H: a subsequent non-D step has output h at least that R, D steps only increase R, and h never decreases. In particular `R_(t-1)<=H` for a reduced word of length t. Equation (4) gives

    2^ceil((t-1)/2) <= H,

hence (5). Every elementary step adds at most k vertices, proving (6). QED.

This supplies a true logarithmic depth bound for all counts in the grammar and a finite witness-size certificate. It does **not** imply that an arbitrary target H belongs to the grammar. The earlier count-bound BFS remains preferable to brute-force enumeration of every word, since it merges repeated states.

## 5. Canonical stationary blocks

A block using only K and D keeps d fixed. For a physical seed d is odd, since it is the Hamiltonian-path count of a tournament; equivalently, oddness of d is directly preserved from the singleton by the actual finite bank's odd deletion counts.

Each K output has odd R and e=d; each D output has even R and e=0. Consequently the last letter of a nonempty stationary block is determined by R parity:

    even R: last letter D, predecessor R = R/2;
    odd R:  last letter K, predecessor R = (R+d)/2.       (7)

For R>d these inverses strictly reduce R. In terms of q=R-d, the inverse K halves q, while the inverse D sends q to `(q-d)/2`, which is smaller than q/2. Thus a given final R and fixed seed R determine at most one stationary word, apart from removable initial K loops when the seed has R=d. Recover the word by (7), stop when the seed R is reached, then replay it on the seed's full state to check the remaining coordinates. If an inverse is nonintegral, crosses below the seed R, or reaches an invalid predecessor, the target is not in that stationary orbit. This is an exact logarithmic-length normal form.

If a physical state has R=d, then h=d. Applying the valid K operation gives `R'=d` and `h'=2d-e`; since the output must still satisfy h'<=R', e=d. Thus the only physical R=d state is `(d,d,d,d)`, and K really fixes it. There is no hidden nontrivial zero-growth orbit.

For the singleton seed d=1, put q=R-1. Appending K appends binary digit 0 to q; appending D appends binary digit 1. The complete stationary state language is

    x(q) = (f(q), 1-(q mod 2), q+1, 1),
    f(0)=1,
    f(2n+1)=f(n),
    f(2n)=f(n)+n+(n mod 2).                             (8)

The leading zero K loops may be ignored. Equation (8) follows directly from the physical K,D updates and gives a canonical physical graph for every nonnegative q.

For count membership in this stationary subclass, strip terminal D letters. A nontrivial remaining word ends in K, so q=2n and

    H=f(2n)>=n+1, hence q<=2H-2.

Therefore H is attained in the stationary singleton subclass if and only if H=1 or `f(q)=H` for some even q in `2..2H-2`. This is another exact all-order decision rule. This does not establish eventual coverage; for example 7 and 13 fail this finite criterion.

## 6. How this changes the remaining problem

The returning grammar can now be analyzed as root-moving transitions separated by canonical binary stationary blocks. For a non-singleton last gadget its output determines

    h_previous = d_output/dG,
    R_previous = (R_output-Q h_previous)/B.

If s>0, it also determines

    e_previous = (g h_previous+L R_previous-h_output)/s.

The old deletion coordinate is forgotten by that step and still has to be supplied by an actually reachable predecessor; when s=0 the endpoint coordinate also remains to be resolved. For a singleton S step,

    h_previous=d_output,
    e_previous=d_output-e_output,
    R_previous=R_output-d_output,

with the necessary identity `h_output=e_output+R_output-d_output`.

These are exact inverse constraints, not a proof that arbitrary numerically feasible predecessors are physical. A successful all-target induction must still select reachable initial interfaces for all desired counts. The logarithmic bound and binary normal form sharpen that obligation, but no finite attractor or cofinite coverage theorem is claimed here.
