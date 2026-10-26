# Root 27, a single ternary tree, and the forced carry leaf

2026-10-26. **PROVED:** the elementary cross-section, tree/carry and
pumped-family statements for the explicit arithmetic macro map.
**FINITE-EXACT:** the bounded checker controls and identified hostile
examples. **OPEN:** the root-27 inequality and blank-state nonhalting.
Native-machine identification is inherited as **CITED**; no machine
transition clock or formalization is asserted here.

## Inheritance and the precise target

The closest proved mechanism is [the previous ancestry/debt construction](../../lca_ports_reverse_20261025/arithmetic/ANCESTRY.md):
for odd `x>1`,

    U_-(x)=3^(k-1)(3x-1)/2^k,  k=v_2(3x-1),
    a'=a+2k-3.

Its integer root label is computable, `root_-(27)=27`, and the remaining
target is

    v_2(y-1) <= 6+2K-3J  whenever root_-(y)=27.     (A)

`J` and `K` belong to the unique integer ancestry from 27. The bound
keeps chronological guards, not just endpoint drift. The canonical
hostile remains an ordinary integer reset with arbitrarily large new
depth; [the density-zero result](../../lca_ports_reverse_20261025/arithmetic/DENSITY.md)
does not exclude a single reset on the marked ray.

The recovered older ingredients are [THM-2228, Mahler carry tails and integral stabilization](../../../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md),
and [the October 5 credit grammar and its chiral tree encoding, sections 7–8](../../trees5_tournaments4_converse_20261005.md).
The former retains the exact integer in an affine carry cylinder; the
latter already has ordered ternary trees for abstract `+2,-1` credit.
This note supplies their exact fixed-source arithmetic interface. The
least-used sidecar is the *position* of the first mismatching carry bit.

Board: two odd cross-sections / native guards / tree shape / exact affine
carry / first mismatch / rational periodic fuel. Anchor: (A). Niche:
forced leaf expansions. Wildcard: the old denominator-19 carrier.

## 1. Change cross-section: the marked seed becomes 7

The inherited abstract scalar map is `x -> floor(3x/2)`. Setting `z=x-1`
gives the integer conjugate

    T(z)=(3z+c)/2=ceil(3z/2),  c=z mod2,
    a'=a+3c-1.                                    (B1)

The cited machine parameter sector is `b>=0`, hence `z=b+3>=3`.
Statements and controls at `z=1,2` use the explicitly extended arithmetic
macro system; all blank-source uses at 7, 26 and 39 stay in the cited sector.

At `a=0,c=0` the native macro is blocked and the process halts; at
`a=0,c=1` the macro is legal and increases credit to 2. Define the
number of pending slots `q=a+1`. A legal `c=1` increases `q` by 2,
while `c=0` decreases it by 1. Reaching `q=0` is the **formal execution
of the first blocked symbol**, not an executed native macro.

The blank parameter state `(a,b)=(0,4)` corresponds to `x=8`, hence
`(q,z)=(1,7)`. Its first four carry bits and pending slots are

    z: 7 -> 11 -> 17 -> 26 -> 39,
    c:   1     1     1     0,
    q: 1 ->  3 ->  5 ->  7 ->  6.                 (B2)

The earlier first odd-`x` source is `(a,x)=(6,27)`, at the third cut;
one more legal macro gives `(a,z)=(5,39)`. This explains the initial
macro minimum 5 after that odd source. It is not a claim that 5 is a
proved minimum over the infinite future.

Now induce `T` on **odd z**, rather than odd x. Put

    ell=v_2(3z+1),
    V(z)=3^(ell-1)(3z+1)/2^ell,
    q'=q+3-ell.                                   (B3)

One `c=1` is followed by `ell-1` zero carries. The whole return is
executable precisely when `ell<=q+2`, equivalently `q'>=1`. If
`ell>=q+3`, a blocked zero occurs after exactly `q+2` executed macros
in this run; the formal return `(B3)` must not be mistaken for a reached
native state. The source guard is `q>=1`, allowing `a=0` because z is odd.

Every positive odd input satisfies `V(z)>z`. This is pointwise growth,
not order preservation: `V(5)=27 > V(7)=11`. Its integer inverse is
unique: for `y=3^r m`, `3` prime to `m`, a predecessor exists iff

    m=(-1)^(r+1) (mod3),
    z=(2^(r+1)m-1)/3,    ell=r+1.                  (B4)

The proof is the same direct valuation calculation with the opposite
sign; its legal positive inverse is smaller than y. In particular 7 is
a root for V, because `7=1 (mod3)` fails the required residue `-1`.
Write this different root function `root_+`; do not identify it with
`root_-` by subtracting 1 on the odd sections.

For the unguarded V-ray from 7, let `J` be its number of returns and
`K_J` the sum of its `ell` values. Then

    q_J=1+3J-K_J.

These J and K count V-returns, not the U_--returns in A.

Consequently the original target is equivalently

    K_J <= 3J for every J on the V-ray from 7.      (B5)

The first three V letters are `(1,1,2)` and reach `(q,z)=(6,39)`.
They are already certified legal. This is an exact change of section
with its credit and clock retained; it does not turn V into ordinary
accelerated Collatz. In particular, the original slope `3/2` in credit
and the Collatz multiplier threshold `log_2(3)` are different predicates
for different maps. V itself grows at every step.

## 2. Exact tree language, including the blocked last symbol

Let the preorder word of an ordered full ternary tree be generated by

    Tree = 0  |  1 Tree Tree Tree.                 (C1)

A `0` is a leaf and consumes one pending slot. A `1` is an internal
vertex, replacing one pending slot by three. Thus a word is one tree
exactly when its pending count, starting at 1, stays positive at every
proper prefix and reaches zero at the end. A forest of q ordered trees
has the same condition starting at q.

For a fixed positive integer source z and native credit `a=q-1`,

    the native process halts
      iff its carry itinerary has a finite prefix of q trees.       (C2)

The forest's **last 0 is the blocked symbol**. A forest word of length
`n` therefore records exactly `n-1` executed macros before the halt.
Completing an earlier component tree is not a native halt: other
pending trees remain. Arithmetic is transported through their ordered
composition, not restarted independently at each subtree.

For the blank source, (C2) is the question whether the itinerary of 7
begins with one finite tree. Its forced prefix `1110` leaves six slots,
so it is equivalently whether the itinerary of 39 begins with a forest
of six trees. A six-tree forest with E internal nodes has length `3E+6`;
the corresponding complete tree from 7 has `E+3` internal nodes, length
`3E+10`, and `3E+9` executed macros before the halt. These are macro
counts, never native Turing-machine transition counts.

## 3. Arithmetic survives as an ordered affine carry

For a binary word `w=c_0...c_(n-1)`, set

    C_w=sum_(i<n)c_i*2^i*3^(n-1-i),
    Q_z(w)=3^n*z+C_w.

The exact carry-cylinder theorem, independently proved by iterating B1,
is

    w is the first n carry bits of z
      iff Q_z(w)=0 (mod2^n).                       (C3)

The congruence also forces each intermediate parity: backwards modulo
2,3 being odd gives the first bit, and division by 2 continues the
induction. The intermediate integer states are exactly those of T.
For concatenation, with lengths `|u|=m`, `|v|=n`,

    C_(uv)=3^n*C_u+2^m*C_v.                        (C4)

Thus `(length,carry)` is an associative ordered composition law. The
single-tree fixed-source target becomes the explicit arithmetic claim

    v_2(7*3^n+C_w) < n for every finite tree word w. (C5)

Equivalently, use `39` and all six-tree forest words. Neither inequality
is proved here. This replaces an unspecified orbit by a precise class
of finite arithmetic certificates, retaining the hard fixed-source filter.

Counts and cycle-lemma rotations cannot replace this carry. For example
the two valid six-tree forests `100000000` and `010000000` have identical
size and one internal vertex, but their source residues modulo 512 are
341 and 398. More generally left rotation obeys

    2*C_(rot w)=3*C_w-c_0*(3^n-2^n).               (C6)

Rotation changes the source congruence, even when a ballot count or a
rotation multiplicity survives. The checker verifies (C6) and the
literal source distinction; no arithmetic quotient by rotation is used.

## 4. A leaf expansion has exactly one chance to improve precision

Replace a leaf at position p in `w=u0v` by `1000`, with `|u|=p` and
`|v|=s`. This preserves the number of component trees and increases
length by 3. Formula C4 gives the exact carry update

    Q_z(u1000v)
      =27*Q_z(u0v)+2^p*(3^(s+3)-38*C_v).           (L1)

The parenthesized factor is odd. If `h=v_2(Q_z(w))<|w|` is the first
mismatch with the actual source, then

    p != h  =>  v_2(Q_z(u1000v))=min(p,h),
    p = h   =>  v_2(Q_z(u1000v))>=h+1.              (L2)

The second case exists only when position h is a leaf. An expansion
before h destroys already correct precision; one after h leaves the
first mismatch unchanged. **Only the first mismatching leaf can improve
precision.** This is an exact obstruction to treating tree expansion
choices as interchangeable arithmetic repairs.

There is a canonical frontier construction. Start with q leaves, `0^q`.
When its first mismatch lies before its end, the prefix already agrees
with the source and the unprocessed suffix consists entirely of leaves.
Expand its first mismatching leaf; the matched prefix increases and the
remaining unprocessed suffix is again all leaves. Induction proves this
description, until a complete matching forest is found, if ever.

At a frontier `w=pref 0^q`, where `h=|pref|` and the actual state after
pref is odd z, one has

    Q_source(w)=2^h*3^q*z.

The active expansion changes this to

    2^h*3^(q+2)*(3z+1).

Hence precision increases by `ell=v_2(3z+1)`, the frontier length
increases by 3, and the new number of pending leaves is `q+3-ell`.
When this is positive the new odd state is exactly V(z) from B3.
When it is nonpositive, the expanded forest already matches and gives
a halt certificate. Thus the forced arithmetic leaf exploration is
**exactly the odd-z return dynamics**, with no extra freedom hidden in
the tree grammar. The remaining invariant question is still (B5), now
expressed as growth of matching precision versus growth of frontier length.

## 5. Every fixed pumped shape has a finite rational fuel certificate

Let P be a nonempty tree word with its last leaf removed. Then its credit
change is zero, all its partial changes are nonnegative, its length L is
a positive multiple of 3, and its first bit is 1. Every `P^e0`, `e>=0`,
is a full tree word. Write

    A=3^L, B=2^L, C=C_P, N_z=(A-B)*z+C.

Here C and `A-B` are odd. For every positive integer source z,

    P^e0 is realized at z  iff  v_2(N_z)=L*e.       (P1)

Proof: concatenation gives
`C_(P^e)=C*(A^e-B^e)/(A-B)`. After appending 0,

    Q_z(P^e0)=3/(A-B)*(A^e*N_z-C*B^e).

If `v_2(N_z)!=Le`, the valuation of this difference is the smaller
one, hence below the required `Le+1`. If they agree, the two odd
normalized terms cancel, giving valuation at least `Le+1`, exactly the
required cylinder test. Since z is positive and C positive, `N_z>0`:
the fuel is finite and the rational periodic source
`-C/(3^L-2^L)` is never the ordinary integer z.

For the right comb `P=100`, this specializes to

    C_((100)^e0)=27*(27^e-8^e)/19,
    (100)^e0 is realized at z iff v_2(19z+9)=3e.     (P2)

The rational source `-9/19` is exactly THM-2228's denominator-19
periodic carrier. For z=7 the valuation is 1, so **every right comb is
excluded**, not merely those below a cutoff. This exclusion uses only
the first mismatching bit and should not be mistaken for substantial
evidence about arbitrary tree shapes.

More generally, the checker takes every nonempty tree through five
internal vertices, removes its last leaf, and obtains 344 fixed blocks P.
For source 7, their fuel values have histogram
`{1:183,2:111,3:10,4:31,5:8,6:1}`; none is divisible by its block length.
Thus this explicitly bounded bank certifies 344 **infinite pumped
families**. The bank computation is FINITE-EXACT; its all-e consequence
is the proved formula P1. Arbitrary trees still have unbounded, nonrepeating
shape and ordered carries. In particular, selecting the entire unknown
tree minus its last leaf as P merely restates C5.

The failure boundary is concrete. Source `z=15` also has
`v_2(19z+9)=1`, so every right comb is excluded. Nevertheless its exact
first halting tree is

    1111000010000010100000,

with seven internal vertices and length 22. The native abstract state
starts at `(a,z)=(0,15)` and halts at `(0,78894)` after 21 executed
macros; the final displayed 0 is blocked. Thus the superficially similar
family `x=2^h,a=0` cannot have universal survival: `h=4` is a counterexample,
whereas the blank source is `h=3`. Both the comb shortcut and the bare
binary shape of the seed lose a decisive arithmetic coordinate.

## 6. A tempting two-prime charge fails on the actual marked ray

For the original U_- section let `d=v_2(x-1)` and `e=v_3(x-1)`. A
valuation-1 step sends

    (a,d,e) -> (a-1,d-1,e+1),

so `a+e` and `a-d` are invariant through a consecutive ones-run. A
reset `k>=2` has output divisible by 3, hence new `e=0`; the change of
the charge `a+e` is exactly `2k-3-e`.

This charge is not monotone even on root_-(27). At its odd-return index
19 the exact state is `x=57971215,a=21,k=2,e=4`. Its charge falls from
25 to 22. The checker reconstructs this finite hostile from 27 rather
than substituting an unrelated root. The survivor is the exact telescoping
law: ones-runs transfer a unit of 2-adic depth to 3-adic depth, while
resets may erase accumulated 3-adic depth. A useful rank would have to
pay for that reset loss; this local charge does not do so.

## 7. Scope and reproduction

The source/target map is the **actual macro itinerary** to an ordered
tree frontier with its exact `(length,carry,source)` data. It preserves
native halt and every finite guard. Forgetting the source or carry loses
arithmetic realization; forgetting the last blocked symbol loses one
macro in the clock. A fixed repeating branch has rational finite fuel,
but the unrestricted frontier remains unbounded.

For fair independent carry bits, the Kraft mass t of one finite tree
obeys `t=(1+t^3)/2`. Finite-height iteration from zero gives its least
nonnegative solution `t=(sqrt(5)-1)/2`. Six trees therefore have mass
`t^6`, recovering the inherited Haar formula after the forced prefix.
This counts cylinders. It assigns no probability to the fixed source 7,
39 or 27 and supplies no pointwise inequality.

Run `python -B check.py` and `python -O -B check.py` in this directory.
The checker uses exact integers, independent literal macro replay,
bounded complete tree and forest universes, leaf-expansion controls,
rotation hostiles, pumped-family certificates and the actual-ray charge
failure. No longer orbit census is used as a substitute for an invariant.
The saved `normal.json` and `optimized.json` are byte-identical `PASS`
receipts, SHA-256
`90c7c8c4751871104028ab613b2d57f7627479586dd61bf36688c09db805baba`.

The strongest next formal question is C5, or equivalently B5: can the
forced frontier's exact carry prevent its matching precision from catching
its length? The present result makes that comparison exact, proves a
source-preserving test for every fixed pumped family, certifies the bounded
bank's exclusions, and identifies why shape counts or local charge
monotonicity do not finish it.
