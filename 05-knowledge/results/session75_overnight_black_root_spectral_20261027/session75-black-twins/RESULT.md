# Same Hamilton count and every named cut, different black distance

**Session75 bounded FINITE-EXACT addendum; independent external audit pending.**
Exactly two fixed order10 tournaments are compared. No search range is enlarged.
All earlier frozen milestone bytes remain unchanged.

## 1. Exact twin distinction

The Phase1 two-R5 one-back-arc tournament has rows

    A=(966,1004,1016,1009,995,193,384,768,544,96).

The Phase2 rooted-converse E10 has rows

    E=(1008,1001,995,997,494,448,640,768,576,48).

Their **named** score vector is identical:

    (6,7,7,7,7,3,2,2,2,2).

Consequently their entire named cut-capacity functions agree, not merely their
minimum cuts: c(S)=sum_(v in S)d(v)−binom(|S|,2). The checker also directly
counts all1,022 cuts from each adjacency matrix and verifies this equality.
Both have unique minimum shore {5,6,7,8,9}, λ=1 and δ=2.

| Observable | Phase1 two-R5 A | Rooted-converse E10 |
|---|---:|---:|
|Directed Hamilton words H|657|657|
|Automorphism order|1|9|
|tc=H/Aut|657|73|
|Blue Hamilton words|33|99|
|Black Hamilton words|624|558|
|All anti-automorphisms|1|9|
|Involutive anti-automorphisms|1|3|
|Minimum-cut-decreasing words|9|9|
|Blue / black decreasing words|3 / 6|9 / 0|
|Black distance to nonstrong|1|2|
|Black distance to singleton + strong core|2|2|

Thus **(H, the full named cut function) does not determine black distance to
nonstrong**. The failure is not a numerical inconsistency: cut capacity and
Hamilton-word totals omit the incidence between legal crossing ports and
positional-reversal anti-automorphisms. Automorphism order and colour counts
already distinguish these two objects.

All counts here concern actual labelled Hamilton words, not isomorphism classes
of words, not complement-line multiplicities, and not the quotient count tc.

## 2. Exact distance certificates

Every Hamilton word of each source is exhausted and classified using the full
anti-automorphism predicate. In A, the black word

    (6,7,8,9,5,0,2,3,4,1)

reverses to a nonstrong tournament, giving distance1. In E, all558 black words
reverse to strong tournaments, excluding distance1 completely. Its inherited
two-black core route supplies the upper bound2. Both sources have δ=2 and
literal two-black routes to a singleton and strong core, proving core distance2.
`normal.json` retains named R/K source and target matrices, failed anti-equations,
actual intermediate strongness, and the resulting core certificates.

The checker independently enumerates all degree-compatible automorphisms and
anti-automorphisms; it checks ordinary isomorphisms and confirms A and E are
nonisomorphic. No score equality is used as a substitute for path colour.

## 3. Incoming THM-4811 is the same E10, with concurrent provenance

The incoming source is commit
`87d61f67673a6b10cbafff1bf9e1377f3c41e56c`, file
`01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md`.
Its C1/C2 E10 is E(C3;K1), built by `exceptional_family.py` from
D=C3[K1,C3,K1], with the reversed bridge at D's first singleton. Directly
transcribing that pinned recipe gives

    (974,1012,1016,1010,993,513,288,96,160,448).

An explicit orientation-preserving isomorphism from our E labels to those
incoming labels is

    (4,1,3,2,0,9,6,8,7,5).

Every ordered pair is checked. Incoming H=657, Aut=9, blue=99, black=558 and
nine all-blue descent words therefore agree exactly with our E10 receipt.
They were never counts for the distinct two-R5 witness.

**Credit:** this reconciles concurrent independent construction and computation.
Session75 does not claim first discovery of E10 or the rooted-converse
obstruction. The incoming statement's complete order10 census and its other
claims are outside this two-object addendum's audit scope; they are not
silently inherited as verified here.

## 4. The full odd-cycle-collection profiles distinguish the twins

Define F(x3,x5,x7,x9) by summing over unordered collections of vertex-disjoint
directed odd cycles. A cycle is counted as a directed cyclic word modulo
rotation, not reversal; unused vertices are allowed. Its variable records its
length. The exact polynomials are

    F_A = 1 + 14x3 + 30x5 + 52x7 + 30x9
            + 35x3² + 48x3 x5 + 6x3 x7 + 4x5² + 4x3³,

    F_E = 1 + 14x3 + 54x5 + 90x7 + 36x9
            + 22x3² + 36x3 x5 + 9x5².

These are checked twice: subset/end dynamic programming rooted at each cycle's
least label, and independent literal cycle-word enumeration followed by
unordered disjoint-packing recursion. Every literal directed odd-cycle word is
retained in the receipts. The independent Hamilton enumeration agrees with
F_A(2,2,2,2)=F_E(2,2,2,2)=657.

Even recording only the number of cycle components distinguishes the pair:

    Q_A(t)=1+126t+93t²+4t³,
    Q_E(t)=1+194t+67t²,
    Q_E(t)−Q_A(t)=2t(2−t)(2t+17).

Evaluation at t=2 loses this distinction. In particular A has four collections
of three disjoint directed triangles, whereas E has none, despite their equal
individual triangle count14. This is a separation for this pair, **not** a
claim that the full cycle polynomial generally determines black distance.

## 5. Reproduction

    python3 -B check_twins.py > normal.json
    python3 -B -O check_twins.py > optimized.json
    cmp normal.json optimized.json

The script is self-contained, stdlib-only, and uses explicit exceptions rather
than assertions. Normal and optimized receipts agree byte for byte. Pinned
source bytes and exact hashes are bound in SOURCE_BINDINGS.json and SHA256SUMS.
