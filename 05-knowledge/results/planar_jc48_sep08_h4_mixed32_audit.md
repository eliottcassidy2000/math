# Independent audit of the H4 mixed (3)(2) equality theorem

**Status: INDEPENDENT ANALYTIC / FINITE-EXACT / SOURCE AUDIT PASS.**
The complete theorem and source were read; fresh normal and optimized
replays agree byte-for-byte with the frozen output. A separate literal
letter-action census, without importing any producer function, recovers
the entire pair classification. No mathematical or source correction is
requested.

Audited primary:
[planar_jc48_sep08_h4_mixed32.md](planar_jc48_sep08_h4_mixed32.md).
Audited source and output:
[source](../../04-computation/planar_jc48_sep08_h4_mixed32.py),
[output](planar_jc48_sep08_h4_mixed32.out).

## 1. Exact accepted statement and the finite-to-unbounded bridge

The accepted theorem concerns permutations `a,b,c,d` of an arbitrary
finite set. Each has one three-cycle and one disjoint transposition,
with arbitrary additional fixed letters. The exact hypotheses are

    aba=bab,  bcb=cbc,  cdcdc=dcdcd,
    [a,c]=[a,d]=[b,d]=1.

The conclusion `a=b=c=d` is accepted. Their generated action is
therefore never transitive, because a common permutation of this type
has separate nontrivial orbits of sizes three and two. There is no bound
on the ambient degree in this conclusion.

The finite computation is a valid uniform supplier because it concerns
**pairs**, not arbitrary four-generator tuples. The moved supports of
two permutations of this type have union size at most ten. Their union
is invariant under both permutations, and outside it both are identity.
Restriction to the union therefore preserves each literal braid identity
in both directions. Add jointly fixed letters if necessary, then
simultaneously relabel so that the first permutation is `(123)(45)` on
ten letters. This does not change cycle factors or any labelled
intersection cardinality.

Every possible partner is obtained by choosing a three-element support,
one of its two cyclic orientations, and a disjoint two-element support.
The complete universe has

    2 binom(10,3) binom(7,2)=5040

distinct partners. The source enumerates all of them before applying the
two braid predicates. There is no preliminary support or transitivity
filter in that pair universe. Its later ten-letter full-tuple check is
only supplemental; the proof does not assume that an arbitrary H4 tuple
has total moved support at most ten.

## 2. Independent pair-table verification and its precise consequences

For a type-(3)(2) permutation `s`, the canonical factors

    P_s=s^4,   Q_s=s^3

are respectively its original three-cycle and transposition. They have
disjoint supports, commute, and multiply back to `s`. This is checked
literally for every partner in the producer. It follows analytically
from their orders three and two.

I independently constructed all 5,040 partners as dictionaries of their
action on ten letters, evaluated the alternating words directly on each
letter, and formed the four support cells from the chosen disjoint
supports. No producer composition, power, braid, cell, or orbit function
was imported. The two resulting tables are exactly:

| Relation | Cell `(P/P,P/Q,Q/P,Q/Q)` | Count |
| --- | --- | ---: |
| 3 | `(3,0,0,2)` | 1 |
| 3 | `(3,0,0,1)` | 10 |
| 3 | `(2,0,0,2)` | 15 |
| 3 | `(2,0,0,1)` | 120 |
| 5 | `(3,0,0,2)` | 1 |
| 5 | `(2,1,1,1)` | 6 |
| 5 | `(2,1,1,0)` | 30 |
| 5 | `(1,0,0,2)` | 60 |

Thus all 146 ordinary partners have both cross-length intersections
empty. The union of the two three-cycle supports and the union of the
two transposition supports are disjoint invariant subsets for the pair.
Restriction of the original ordinary relation therefore pays ordinary
braiding separately on the two factor pairs. This is a legitimate
restriction after the supports have been separated, not an assumption
that taking powers preserves relations.

For length five the exact conditional statement is also accepted:
if either cross-length cell is empty, both are empty and the two
transposition supports coincide. A transposition is determined by its
two-element support, so the canonical transpositions themselves coincide.
The 36 mixed-cell partners are retained in the table and do not satisfy
that hypothesis. Only 61 of the 97 fifth-edge partners preserve the
fifth braid after cubing; all 146 ordinary partners preserve the third
braid. The source's separate factor checks agree with these counts.

The producer's centralizer orbit verification is complete. Its
generators give `C3 x C2 x S5`: the two unequal nontrivial cycle lengths
are preserved, and all five fixed labels can be permuted. Positive
iteration of its finite-order generators also generates their inverses.
The queue closure visits all elements of each orbit. Disjointness,
membership in the entire raw universe, predicate invariance, and
support-cell invariance are all checked. The 64 orbits have total size
5,040 and independently reconstruct the weighted accepted table.
This orbit check does not thin the raw universe before its verification.

## 3. The ordinary three-cycle triple lemma is complete

The small structural lemma used by the producer is valid: if single
three-cycles `p,r` commute and a single three-cycle `q` ordinarily
braids with both, then `p=r`.

First, commuting single cycles of the same length have equal or disjoint
supports. If these two three-cycle supports were disjoint, the ordinary
odd-run bound would force the three-letter support of `q` to meet each
in at least two letters, which is impossible. The bound itself follows
by applying `pqp=qpq` to two consecutive points of the `q`-cycle both
fixed by `p`: the two resulting points would be `q(x)` and `q^2(x)`,
contradicting that `q` moves `x`.

On a common three-element support, either `r=p` or `r=p^-1`. In the
second alternative, a `q` with that same support cannot braid with both:
all the permutations on that support in question commute, and a commuting
odd braid forces equality. The only other possibility is that `q`
shares exactly two support letters with `p`, with one external letter.
Conjugating by a power of `p` and relabelling that external letter fixes
the canonical possibilities

    p=(123),   q=(124) or (142).

Only `q=(142)` ordinarily braids with `p`. For that choice `p^-1`
does not: on letter 2 the two words `p^-1 q p^-1` and `q p^-1 q`
give 4 and 3, respectively. This independently confirms the producer's
two-orientation control. Every possible support/orientation case has
been covered; no larger ambient enumeration is being inferred from the
four-letter example.

## 4. The common transposition and equality reduction are valid

Apply the separated ordinary pair relations to `(a,b)` and `(b,c)`.
Because `a` commutes with `c`, their three-cycle factors commute as
well. Section 3 gives

    P_a=P_c.

Also `[a,d]=1` implies `[P_a,Q_d]=1`. A three-cycle commuting with a
transposition has disjoint moved support: the transposition's
two-element support would otherwise contain a full moved orbit of
length three. Consequently `P_c` and `Q_d` are disjoint. The actual
fifth relation `(c,d)` now satisfies the conditional table hypothesis,
and hence

    Q_c=Q_d.

Next, `[b,d]=1` and the ordinary relation between `b,c` make `Q_b`
both commute and ordinarily braid with this common transposition.
Thus `Q_b=Q_c`. Finally `[a,c]=1` and the ordinary relation `(a,b)`
give `Q_a=Q_b`. All four canonical transpositions are equal to one
fixed transposition `Q`.

This already gives an invariant proper two-element subset in the
original ambient set, so transitivity is impossible even before
proving equality of the three-cycles. For the stronger conclusion,
each `P_i` commutes with `Q`, and `s_i=P_i Q`. Every odd alternating
word therefore has the same final factor `Q` on both sides, which
cancels. All three commuting relations also descend. The four
three-cycle factors satisfy the full original H4 relations.

The dependency
[H4 single-cycle theorem](planar_jc48_sep08_h4_single_cycle.md)
is PROVED and independently audited on disk. I reread its odd-run,
support-coincidence, and centralizer equality proof. It applies with
cycle length three and arbitrary fixed padding here, so it gives
`P_a=P_b=P_c=P_d`. Its hypothesis is invoked only after the common
transposition has been established. Multiplying back by `Q` proves
the asserted equality on every original ambient label.

No step assumes that cubing is a homomorphism, that the generators are
involutions, or that an H4 Artin image is a Coxeter image.

## 5. Hostiles, source audit, and geometric limits

The minimal hostile `(123)(45),(124)(35)` passes the fifth braid, but
its cubes are `(45),(35)`. Their product has order three, so they fail
the fifth braid. Five is the least possible moved-support size for
this mixed type. The example refutes the proposed local power-transfer
map; it does not refute the repaired four-generator theorem, since it
does not supply the remaining H4 relations.

Taking all four generators equal to `(123)(45)` passes every relation
with arbitrary fixed padding. The producer retains this positive
equality control at ambient degrees 5, 6, 11 and 25. It remains
intransitive even with no fixed padding because the three-cycle and
transposition orbits are separate.

The complete producer was read, including the literal cycle constructors,
unique-partner check, both word evaluation paths, odd conjugator check,
all raw support cells, complete centralizer closure, orientation control,
supplementary full tuple bank, and semantic serialization. The bounded
tuple filters are exactly necessary relations: fixing `b=sigma`,
`a,c` are among its 146 ordinary partners and `d` among its 62
commuting partners, followed by the remaining two commutators and the
fifth relation. Exactly one tuple remains. This agrees with the theorem
but is not used to pay arbitrary ambient degree.

The actual geometric consumer is appropriately separated. The current
[mixed-cusp braid supplier](planar_jc48_sep08_mixed_cusp_braid.md) is
PROVED and independently audited, and gives generating positive
meridians with chain `(z,y,w,x)` and labels `3,3,5`; its complement
may have additional relations. The
[intrinsic mixed boundary family](planar_jc48_sep08_three_cusp_boundary.md)
is also promoted and carries the stated connected good-locus transport.
I checked these current statements and their exact matching of the H4
labels; this audit does not rerun their already audited geometric tube
certificates. A connected complement cover has transitive monodromy,
so the present equality/intransitivity theorem excludes its meridian
type `(3)(2)` once that supplier is used. It does not construct a Keller
map, identify the full complement with Artin H4, or exclude other mixed
cycle types such as `(4)(2)` or `(3)(3)`.

## 6. Fresh replays, pins, and independent reconstruction

Fresh commands from the worktree were

```bash
python3 04-computation/planar_jc48_sep08_h4_mixed32.py
python3 -O 04-computation/planar_jc48_sep08_h4_mixed32.py
```

Normal, optimized and frozen bytes agree: **51,345 always-active gates**,
3,607 output bytes. Every producer check raises explicitly under failure
and remains present under optimization.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Source | 8661 | `9f83b3345e0d39da9f7abf1bede93b45deecd3234501d5331ebf28a4b640f0f6` |
| Output | 3607 | `399694228f01eb1409c1476d2a42e25eda8eae32811b3367d9a65fb7f5da44ef` |

The output also reproduces complete raw-pair semantic hash
`9d33769e2dd066673ed4c828be3d280694ec702f1d364581f00111ded1559031`
and complete orbit hash
`554bcdaaf35b135645f4c5e4b215e3fca68e56f76b7c1b70c99305972589fadd`.

For reproducibility, the separate table reconstruction used only this
elementary scheme: choose disjoint supports of sizes three and two,
construct the two possible oriented dictionary actions, and apply the
odd palindromic alternating words to each of the ten letters. It formed
support cells directly from the chosen supports and counted accepted
partners. All 5,040 were visited and both complete tables in Section 2
were recovered. This path does not use powers to recover supports and
does not call the producer's permutation, braid, or cell functions.

**Final acceptance: PASS.** The proof genuinely removes an ambient
degree cutoff using the support size of each local pair. The complete
finite pair theorem, the common-transposition argument and the already
proved single-cycle theorem together establish equality uniformly.
All producer files remain unchanged by this audit.
