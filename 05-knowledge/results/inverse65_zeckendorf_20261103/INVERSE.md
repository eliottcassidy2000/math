# Exact inversion and the cost of retaining a source

Status: **PROVED** for the fixed 65-symbol bank, with exact finite matrix
certificates. No all-primes surjectivity theorem is asserted.

## I1. Objects and conventions

An atom has labelled source rows and the Hamilton spine `0,...,n-1`.
Concatenate atoms as ordered SCCs, toggle every nonadjacent spine pair, and
reverse the labels to put the processing port at the newest vertex. This is
the literal `actual_target` operation in `bank.py`. All orientations remain
Boolean. A word in this packet is processed in reverse source-SCC order.

The inherited three-coordinate observer is `(h,x,y)`, with physical coordinates

    H=h, A=h+x-y, B=y-x, C=x+y.

Here B counts Hamilton paths ending at the port; C counts paths starting there
plus the marked two-path covers of
[G1-G3, golden_hamilton_transfer/RECURRENCE.md](../golden_hamilton_transfer_20261103/RECURRENCE.md).
This is a faithful observer for future H under same-end joins, not a faithful
encoding of every adjacency or arbitrary re-marking operation. Literal word
and labels remain sidecars. Its formal root is `(1,0,0)`; an actual singleton
is `(1,0,1)`. The scalar inverse returns the latter for H=1.

Every nonempty reachable state lies in

    h>=1 odd, 0<=x<=y<=h+x, h<=x+y, x,y<=16h.             (I1)

The inequalities A,B>=0 are physical. Cutting a Hamilton path immediately
before the port injects it into C's objects, proving C>=H. The bound 16h and
monotonicity of h follow by exact matrix inequalities on the four extreme
rays of `h>=0, 0<=x,y<=16h`; `verify.py` rechecks all 65 matrices.

## I2. A multiplicatively increasing norm

Put `u=(6,9,5)^t` and `N(v)=v u`. Every matrix P in the bank satisfies

    5 P u >= 9 u                                         (I2)

entrywise. This finite certificate is reproduced by `branch_probe.py`.
Consequently a word of length ell satisfies

    6*(9/5)^ell <= N(v) <= 230h,
    ell <= log_(9/5)(230h/6).                            (I3)

There are at most six actual vertices per symbol. This is a certificate-size
bound, not a polynomial bound in the bit length of h. The finite cone permits
O(h^3) memoized states; see the independently obtained inverse and fibre
polynomials in [incoming THEORY I3](../golden_inverse_defects_20261103/THEORY.md).

Every P is nonsingular. Given a state w, each candidate last letter determines
one rational predecessor `v=w adj(P)/det(P)`. Retain it only when it is integral,
obeys (I1) or is the formal root, satisfies `h(v)<=h(w)`, and satisfies
`9N(v)<=5N(w)`. Recurse, accepting exactly at the root.

**Soundness:** a returned word multiplies from the root to w and reconstructs
an actual Boolean tournament. **Completeness:** every actual word supplies
its own retained last-letter predecessor. **Termination:** N strictly decreases
on every reverse step and is a positive integer on all admitted states.
Thus `decode_state` is an exact membership-and-witness inverse for this bank.
It selects a representative; it does not recover an unknown original word
from the three coordinates alone.

For scalar h, bounded-h forward closure is finite by (I1), and monotonicity
means it discards no potential witness. Invert any state of first coordinate h.
Alternatively enumerate the finite cone fibre and apply the inverse. This is
a uniform decision algorithm, not proof of universal positive answers.

**FINITE-EXACT:** all 747 reachable states with h<=151 decode; all 942 absent
states in the tested cone universe `h<=15, x,y<=2h` reject. The state
`(37,18,19)` rejects although H=37 has another witness. Exact graphs with
H=499,1009,19861,21001 are independently counted by subset DP.

## I3. Why a simple digit induction cannot retain the same source

The singleton a has

    a: (h,x,y) -> (h+2x,y,h+x),
    H_a=h+2x, H_aa=h+2(x+y), H_aaa=3h+4x+2y>=5h.         (I4)

Every other letter satisfies `H_new>=5h` on the physical cone. To check this
linear inequality it suffices to check rays `(2,1,1),(1,0,1),(0,1,1)`;
the first two normalize h, and the last is the recession ray. These rays
generate `0<=x<=y<=h+x, h<=x+y`. `branch_probe.py` checks all letters.

Thus any continuation whose result is strictly between h and 5h is either
a or aa. For h>=3 the three values `3h,3h+2,3h+4` all lie in this interval.
At most two can be continuations of one retained state. Similarly at most
one of `2h-1,2h+1` is possible, since both lie below 3h and aa has H>=3h.
Repeated finite-seed prefix induction has at most two children per source
where the full ternary digit tree requires three; a fixed finite seed bank
cannot remove the exponential branch deficit. This is scoped to same-end
append-only operations, without resetting to an unrelated representative.

The obstruction is not an objection to the existence of each desired scalar
somewhere in the family. It specifically identifies the failed implication
from scalar reachability to source-compatible recursive reachability.

The concurrent [THEORY I4](../golden_inverse_defects_20261103/THEORY.md) sharpens
this: no state supports both positive offsets3h+2 and3h+4. The only possible
cross-cases would require x=h+1, x+y=h+2 (hence y=1), or x=h+2, x+y=h+1
(hence y=-1). Both violate the physical cone. We credit this stronger result.

## I4. Exact positive ternary offsets from two singletons

For every non-a P and t=0,1,2, the three physical cone rays certify

    (C-H)(v P a^t) >= H(v P a^t)/96.

The minimum 1/96 occurs for r16 followed by a. For aaa, (I4) gives
`C_new-H_new=x+y` and `H_new<=6(x+y)`. Every word of at least four symbols
has a nonempty prefix followed by one of these suffix types. Therefore

    C-H >= H/96  for every such word.                    (I5)

To get `H_aa=3H+2` or `3H+4` one needs C-H=1 or 2. Equation (I5) bounds
H by192 for long words. Exhaust all reachable states up to192 and every
word of length at most3 (278915 words). The complete union contains exactly
41 distinct qualifying states, with largest H=155; the receipts include a
word for each. Therefore neither positive ternary offset can be supplied by
aa once H>155. This does not exclude offsets obtained by a in other states.

## I5. Complete H=55 and hostile controls

Monotonicity and (I3) justify enumerating all words of H<=55. There are334.
Exactly six have H=55, all with seven vertices:

| Processing word | (h,x,y) | H after a |
|---|---|---:|
| a a r26 | (55,48,57) |151|
| a a r34 | (55,43,50) |141|
| r17 a a | (55,20,43) |95|
| r26 a a | (55,20,43) |95|
| a a q1 a | (55,20,41) |95|
| a q4 a a | (55,20,43) |95|

All 768 one-vertex cuts fail to give111, independently counted on the literal
graphs. These are six presentations, not a claim of six isomorphism classes.
The inherited source rows `[118,108,24,113,98,68,4]` have H=55 and their cut96
does give111. Hence that source is not isomorphic to any of these six.
The bank also attains111 directly at order7; that is not an attachment of
one of its H=55 sources.

The 81/27 control is replayed literally: same67 cyclic triples, common
H=462915, automorphism orders81 and27. After the same sink and port8 operation
the counts are2390013 and1754919. A scalar, cyclic-triple list, or unmarked
quotient cannot be substituted for the required attachment state.
