# Independent audit of H4 meridians of type (5)(2)

**Status: PROVED / INDEPENDENT AUDIT PASS.**
The complete analytic proof, unfiltered pair universe, frozen source,
and independent normal/optimized replays of
[the primary](planar_jc48_sep08_h4_mixed52.md) are accepted. A separate
referee implementation also reconstructed every typed overlap/joint-orbit
row by a different enumeration order and a different orbit algorithm.
No correction was needed. Producer files were not edited; root owns
the primary's promotion and integration.

## 1. Exact accepted theorem and finite-to-unbounded bridge

Let `a,b,c,d` act on any finite set, each with exactly one 5-cycle,
one disjoint transposition, and otherwise fixed points. Under

    aba=bab, bcb=cbc, cdcdc=dcdcd,
    ac=ca, ad=da, bd=db,

the theorem concludes `a=b=c=d`. Their common permutation has
distinct nontrivial orbits of sizes five and two, so the generated
action is intransitive even without any fixed padding. Equality
does not mean that the common permutation is the identity.

The proof does not put all four generators on fourteen labels.
Only a pair is reduced: each member moves seven labels, so its
support union has size at most fourteen. Restricting to that union,
conjugating the first member to

    sigma=(1 2 3 4 5)(6 7),

and adding fixed padding preserves every tested word and every
typed support. All other labels of the original action are fixed
by this pair and impose no additional condition. Thus the complete
fourteen-label partner inventory proves all-pair facts in arbitrary
ambient size. The later common-support argument is a separate
unbounded set argument; no full-tuple extrapolation occurs.

## 2. Full pair universe and literal predicates

The source constructs each partner uniquely by choosing its five-set,
fixing that set's smallest element as the first cyclic entry, taking
all 24 orders of the other elements, and choosing an unordered pair
outside it. The count is

    binom(14,5)*24*binom(9,2)
      =binom(14,7)*504
      =1,729,728.

Fixing a first cyclic entry removes precisely cyclic rotation,
without identifying inverse orientations. Each transposition occurs
once as an unordered two-set, and disjointness from the five-set
is part of the declared cycle type rather than a mathematical
filter. There are no support, parity, braid, or orbit filters
before testing the words.

For every raw partner the producer compares tuple composition with
direct nested action on all fourteen letters. The three predicates
are exactly the ordinary relation, the fifth relation, and
commutation. The construction guarantees a bijection; admitted
partners additionally have their literal cycle type and canonical
length gauge rechecked. Joint orbits are computed only for admitted
partners, which is sufficient for every claimed table. Forward
generator edges generate the same finite orbits as the full group.

The fixed-width raw digest contains fourteen image bytes and then
three Boolean bytes for every partner, including every rejected
one. Thus the frozen replay checks the complete declared universe,
not merely a digest of surviving entries.

## 3. Independent reconstruction of every typed/orbit row

I wrote a separate referee enumerator without importing any producer
function. It chooses the transposition first, then the five-set
from the remaining twelve labels, and puts the least five-set
element last in the cyclic order. The raw count is independently

    binom(14,2)*binom(12,5)*24=1,729,728.

It checks the words by direct letter action. Instead of the
producer's forward orbit traversal, it forms a union-find structure
on undirected edges of both generators; its component sizes are
the joint orbit sizes. It recovers every entry below, with the
typed coordinates ordered by cycle lengths five and two.

| Relation | Typed entries | Joint nontrivial orbit sizes | Count |
| --- | --- | --- | ---: |
| Ordinary | (3,2,2,0) | (7) | 10 |
| Ordinary | (4,0,0,1) | (3,6) | 420 |
| Ordinary | (4,0,0,2) | (2,6) | 35 |
| Ordinary | (5,0,0,1) | (3,5) | 84 |
| Ordinary | (5,0,0,2) | (2,5) | 6 |
| Fifth | (2,1,1,0) | (10) | 2,100 |
| Fifth | (4,0,0,2) | (2,6) | 105 |
| Fifth | (4,1,1,1) | (7) | 10 |
| Fifth | (5,0,0,2) | (2,5) | 6 |
| Commuting | (0,0,0,0) | (2,2,5,5) | 504 |
| Commuting | (0,0,0,2) | (2,5,5) | 504 |
| Commuting | (5,0,0,0) | (2,2,5) | 84 |
| Commuting | (5,0,0,2) | (2,5) | 4 |

The respective totals are **555, 2,221, 1,096**. The independently
serialized complete typed/orbit table has the same digest as the
producer:

    cb623d547de2fe30b050373c9bb5c2acc215156dcc1563451cb189d46ca6bc81.

In particular, ordinary-related five-blocks intersect in at least
three labels, and fifth-related five-blocks intersect in at least
two. The full table includes cyclic orientations that share supports
without agreeing as permutations; none are collapsed to a support
representative before testing the relation.

## 4. Commuting blocks and the uniform majority proof

The elementary commuting lemma is valid independently of the census.
A permutation commuting with `x` sends each `x`-orbit to an orbit
of the same size. There is exactly one five-point orbit and one
two-point orbit. Both are preserved. On the five-set the restriction
is a power of the original 5-cycle; primality of five makes that
restriction either trivial or a full 5-cycle. On the two-set the
restriction is identity or its transposition. Therefore typed
blocks do not cross, and each pair of same-length blocks is equal
or disjoint.

The commuting total also has a direct centralizer count. If the
new five-block is outside the original support, the transposition
is either the original one or occupies the remaining two fixed
labels; each case gives `binom(7,5)*24=504`. A common five-block
has four nonidentity powers, and an outside transposition has
`binom(7,2)` choices, giving 84. Common five- and two-blocks
give four. This independently gives 1,096 and all typed cells.

Let `A5,B5,C5,D5` denote the four long supports. If `C5` were
disjoint from `A5`, then

    |C5 intersect B5| <= 5-|A5 intersect B5| <=2,

contradicting the ordinary `b,c` overlap of at least three.
Since `a,c` commute, their equal-or-disjoint alternative gives
`C5=A5`. Commutation of `a,d`, followed by the positive fifth
overlap of `c,d`, then gives `D5=A5`. Finally `b,d` commute,
and the positive ordinary overlap of `a,b` gives `B5=A5`.
Thus all four preserve the same five-set and its complement.

Only now is restriction legitimate. On the common five-set,
`c,d` are powers of `a`, so they commute. Their odd fifth
relation forces equality there. On the complementary invariant
set, each is exactly one transposition. The product of two
transpositions has order one, two, or three; the fifth braid
relation for involutions is equivalent to `(cd)^5=1`, so only
order one is possible. Their restrictions on the complement
are equal too, hence `c=d` globally. The commuting ordinary
relations successively give `b=c` and `a=b`.

Every implication is reversible where cancellation is used and
applies to the actual permutations. The proof does not take
powers of arbitrary Artin-related permutations and assume that
the relation survives. It does not use a single-cycle theorem,
an ambient tuple bound, or transitivity as a hidden hypothesis.

The primary's unused centralizer sidecar also checks out. Its
ordinary joint-orbit sizes are distinct, so a common centralizer
preserves those nontrivial orbits individually. A centralizer
element's nontrivial cyclic action on a transitive orbit is
semiregular. An element with just one 5-cycle and one transposition
can act nontrivially in this way only on an orbit of size five
or two. The source correctly leaves singleton orbit permutations
unrestricted. This sidecar is not needed for the equality proof.

## 5. Hostiles and the actual geometric consumer

The two literal seven-label hostiles were checked with their original
orientations:

    tau3=(1 6 4 7 2)(3 5),
    tau5=(1 4 3 6 2)(5 7).

They satisfy the stated ordinary and fifth relations respectively,
but their extracted transpositions do not. Seven is the minimum
possible support for a permutation of this type, so the ambient
minimality claim is justified. These controls are individual pairs,
not full H4 representations.

The further partner `(1 2 4 5 3)(6 7)` has both supports equal
to those of `sigma` and satisfies the fifth relation, but is
neither equal to nor commuting with `sigma`. This directly checks
the distinction between support coincidence and cyclic-order
coincidence. The common invariant support plus the global
commuting sidecar is essential in Section 4.

The equal quadruple gives the positive control and has the two
distinct nontrivial orbits five and two. For the actual geometric
consumer, the independently proved marked supplier must provide
four positive meridian permutations satisfying these six literal
relations, all of this declared cycle type, and generating the
transitive action. Those are the precise hypotheses under which
the equality theorem excludes the type. The primary does not
realize every abstract tuple as a Keller monodromy, lose the
marked access conditions, or assert a new numerical degree
floor without auditing the other cycle types.

## 6. Frozen replays and pins

The complete source was read, including all relation predicates,
canonical cycle and orbit routines, raw-stream order, terminal
restriction controls, and hostiles. All checks raise explicit
exceptions and remain active under optimization. The source uses
no imported mathematical implementation and no implicit admissible
pair filter.

Independent full commands were

```sh
python3 -B 04-computation/planar_jc48_sep08_h4_mixed52.py > /tmp/h4-mixed52-audit-normal.out
python3 -B -O 04-computation/planar_jc48_sep08_h4_mixed52.py > /tmp/h4-mixed52-audit-optimized.out
python3 -B /tmp/h4_mixed52_independent.py > /tmp/h4-mixed52-independent.out
```

Both producer replays completed with **1,749,134 always-active
gates** and are byte-for-byte identical to the frozen output.
The separate short-block-first enumerator recovers the complete
unfiltered universe, all tables and totals with the distinct
union-find orbit implementation. The final primary and source
were reread and their pins checked without modification.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Primary before promotion | 13,395 | `8db372de7ee41223a45b913dbf7c0da7ece8676593cbf1a6677bf61b9c7a3edd` |
| Frozen producer source | 8,729 | `bb0bb94252ab3e818e26172d2bc224ec1039964fd0e9f290a9cb6057ada301b9` |
| Frozen output and each independent replay | 1,155 | `6ccb824360d4c5416b88d1fead54503b58e00adea3e412ee0efbccb8c1b5e175` |
| Independent short-first/union-find implementation | 1,985 | `2cbfb36126387bbf27472e4adf0cb9f39e00163ad9d722c155b7a1c40b5ff2d7` |
| Independent table output | 575 | `102102b9b499b3c877f5c0932c492e89f6605e5d78b609263686d5d5c3ec7f08` |

The complete producer raw-stream digest is
`f8dac4bf55a921332df0c02c5b262c067ac9747e89c807d6d937e7596533d801`.
The shared complete typed-table digest is
`cb623d547de2fe30b050373c9bb5c2acc215156dcc1563451cb189d46ca6bc81`.
No analytic or source repair remains. This audit accepts the uniform
equality theorem and its explicitly conditioned actual consumer.
