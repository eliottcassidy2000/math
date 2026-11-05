# Reflection modulo four and the complete self-converse spectrum through73

Separate Session197 strengthening of the unchanged through41 packet. **Every even-order self-converse tournament has H=1 modulo4.** The complete self-converse Hamilton-count spectrum through73, at arbitrary order, is

    {1,3,5,9,11,13,15,17,25,27,29,31,33,37,41,
     45,49,51,57,61,65,67,69,71,73}.                  (1)

Thus the ordinary prime targets19,23,43,47,53,59, the7p target35, and the prime-induction seed63 are all absent from the self-converse class. The additional omitted composite values39 and55 also have unrestricted witnesses. This is a model restriction, not a new unrestricted spectrum omission or a literature-priority claim.

## 1. The all-order congruence

Let tau be an anti-involution of a self-converse tournament. Its existence at every order and the empty/singleton conventions are proved in `CONVERSE_MASS.md`. Acting on a directed odd cycle means applying tau and reversing its cyclic order, so the result is again a directed cycle. This is an involution on cycles counted modulo rotation.

From the inherited odd-cycle collection formula, collections with at least two cycles have weight divisible by4. Hence

    H(T) = 1+2*(number of directed odd cycles) mod4.

All cycles not individually invariant under tau occur in pairs. If I_tau counts the invariant directed odd cycles, then universally

    H(T) = 1+2 I_tau mod4.                           (2)

An anti-involution fixes at most one vertex. At even order it fixes none, since an involution's fixed-point count has the parity of the vertex count. An invariant odd cycle would have an odd-sized tau-invariant vertex set, impossible for a fixed-point-free involution. Therefore I_tau=0, proving H(T)=1 mod4.

At odd order tau has one fixed vertex, and any invariant odd cycle must contain it. No claim that every odd cycle through that vertex is invariant is needed. The directed triangle has H3, showing that the even-order assumption matters. The self-converse strong four-tournament[2,8,3,5] has H5, showing that1 modulo8 would be false.

For an even-order self-converse7p tournament, p must consequently be3 modulo4. This is compatible with the14-vertex example whose quotient2163467 is3 modulo4, and does not apply to the odd-order11-vertex example with quotient8069=1 modulo4.

## 2. Complete bounded support at arbitrary order

Use the exact SCC reduction proved in `CONVERSE_SPECTRUM.md`:

    H(T)=d^2 h,

where d is an unrestricted tournament count and h is a strong self-converse count, or1 when there is no central component.

Busch's externally checked strong-order theorem, with the same precise source as `PRIME_CUBES.md`, gives5^(n-1)<=H^3. The integer inequality

    5^8=390625 > 73^3=389017

therefore bounds every central strong component of count at most73 by order8. This is not an order bound on an arbitrary weak self-converse tournament.

The complete normalized self-converse universe has2^floor((n-1)^2/4) presentations. Exhausting it through8 gives exactly (1) below the cutoff. In particular the order8 contribution through73 is

    {1,5,9,17,25,29,33,37,41,45,49,57,65,73},

all1 modulo4 as the theorem requires. Both endpoint and startpoint subset recurrences independently count every one of the4096 order8 presentations; the prior through7 independent controls remain unchanged.

Since d is odd and d^2<=73, it belongs to{1,3,5,7}. Multiplying entries of (1) by1,9,25 or49 creates no new value at most73. Allowing7 here is only a harmless upper-bound enlargement; no existence of an H7 input is assumed. Every value of (1) already has an actual self-converse witness in the finite universe. This proves completeness at every order.

The mod4 argument gives a smaller independent reduction for the omitted3-mod4 targets43,47,59 and63: a central strong witness must have odd order, hence at most7 under this cutoff. The finite set through7 omits them. For63, the only potentially relevant nontrivial square factor is9, which would demand a central H7 witness; none exists in that bounded self-converse bank. The full square-closure calculation treats this automatically.

## 3. Explicit unrestricted escape constructions

The earlier through41 packet gives actual witnesses for19,23,35 and39. Six further examples are two-apex tournaments: a->b, increasing chain, and a letter c has bit0 for a beating the chain vertex and bit1 for b beating it.

    H43, word(3,1,2,0):
      [14,20,56,50,33,3]
    H47, word(1,1,1,2,2):
      [30,96,122,114,98,65,1]
    H53, word(1,0,2,1,0):
      [38,16,122,115,97,66,3]
    H55, word(1,1,2,3,0):
      [46,48,122,114,97,64,3]
    H59, word(1,0,2,0,0):
      [6,16,122,115,97,67,3]
    H63, word(1,0,0,1,2,2):
      [38,192,250,243,227,194,129,1].

All displayed counts are independently checked by subset DP and literal Hamilton-word enumeration, with physical word-to-graph equality. These give explicit escapes from the self-converse restriction; they are not claims of new numerical records or of an infinite prime-generating family.

The point is structural: even the complete self-converse construction universe cannot supply the desired unrestricted prime/7p induction, and it cannot even supply the finite seed63, no matter how much singleton padding or how many paired weak components are added.

## 4. Reproduction

Run `python3 verify_even_converse.py` and `python3 -O verify_even_converse.py`, comparing receipts and witness files. The script checks both independent Hamilton recurrences on every normalized graph through8, the even-order congruence on that whole universe, the general invariant-cycle congruence on the smaller named control universe, the exact strong-order bound and SCC square closure, all new physical escapes, and one self-converse witness for each entry of (1).

All earlier packet files are unchanged. The universal congruence is proved by the cycle involution, and arbitrary-order completeness by the strong-order bound plus SCC reduction, not by a finite extrapolation.
