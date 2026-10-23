# Full-orbit incidence excludes center orders 18, 20, 24 and 28

October 22, 2026. **PROVED**, from the explicitly inherited Hamiltonicity,
pancyclicity and half-supported-prefix inputs below: no strong self-converse
tournament of black path-orbit deficit six has order **18, 20, 24 or 28**.
Any such tournament of order **21 must be vertex-transitive**. The
non-vertex-transitive cycle type $q=3,m\ge5$ is also excluded. These are
structural arguments, not completed frame censuses.

The incoming six-order frontier was $15,18,20,21,24,28$. This packet alone
leaves $15,21$; the separately owned order-15 packet supplies its own closure.
Nothing here enumerates or excludes the vertex-transitive order-21 lane.

**Audited joint close-out:** the independent
[order-15 theorem](../nontransitive_fifteen_20261022/README.md) and
[complete VT21 companion](../nontransitive_fifteen_20261022/VT21.md),
together with the non-VT21 theorem here, close both remaining positions.
Thus all six orders are excluded. The division of proof obligations remains
explicit: four orders and the non-VT reduction here, two separate complete
finite necessary universes in the companion packet. No metagraph wiring is
inferred from this classification.

## 1. Source, inherited mechanisms and missing coordinates

The source is a literal strong self-converse tournament $T$. Write
$A=\operatorname{Aut}(T)$, let $H$ count directed Hamilton words, put
$h=H/|A|$, and let $b$ count the $A$-orbits of words whose positional reversal
is an anti-automorphism. The deficit $\Delta=h-b$ counts **black Hamilton-word
orbits**, not vertices, arcs, or unlabelled tournaments. A claimed lower bound
must retain actual Hamilton words on this same source.

The closest proved mechanisms are the proper-shore deletion, odd-shore
balance, and half-supported-prefix lemmas in
[the audited odd-shore proof, §§1–3](../odd_shore_deficit_centers_20261021/session22-centers/PROOF.md).
The last lemma says: a directed $r$-cycle $S$, with $r$ a power of two,
$n\ge2r$, strong complement, and a terminal with $t$ exits from $S$ gives
at least $2(t-1)$ black orbits on a self-converse source. Entry configurations
are transported by an actual anti-automorphism. The proof retains both the
prefix paths and their converse/reversal mates; it does not merely double a
count without checking disjointness. Moon's strong-tournament Hamiltonicity
and pancyclicity are the **CITED/inherited** inputs used there and here.

The current cycle-index reduction is
[the combined cutoff](../deficit_six_centers_20261021b/COMBINED_CUTOFF.md).
For a suitable marked Hamilton cycle with successor $R$, let

$$m=|A\cap\langle R\rangle|,\qquad q=n/m.$$

At the remaining orders one may choose $2\le q\le8$, with odd
$m\in\{1,3,5,7\}$. Its $q$ rotation shores are the orbits of $R^q$.
Each is regular, and strong when $m\ge3$. Literal marked-cycle arcs give a
perfect matching from each shore to the next, including the wrap.

The full-orbit refinement, and the warning that anti-isomorphisms may permute
whole orbits rather than fix each shore, are inherited from
[center transitivity](../center_transitivity_20261021c/README.md).
The complete prior finite exclusions include
[order 12](../deficit_twelve_20261021c/README.md),
[order 14](../deficit_fourteen_20261021c/README.md), and
[orders 25/35](../center_transitivity_20261021c/ORDER35.md).
[VT15](../center_transitivity_20261021c/VT15.md) excludes only the
vertex-transitive part of order 15; that scope is retained here.

**Hostile boundary:** the known six-vertex center really has $\Delta=6$.
The two-family bound below is only four when $q=3$ or $4$; a third family
is essential in those cases. The four-prefix construction needs $m\ge5$,
and cannot be applied to a three-vertex shore. **Corrected near miss:**
regularity or equal scores do not certify a full automorphism orbit.
The least-used sidecar is the *entire full-orbit label word*, including its
run boundaries and the fixed action of the anti-coset on those labels.

The live objects are full vertex orbits, rotation shores, cyclic necklaces,
Hamilton block paths, removable prefixes, and literal complement targets.
Their precise map is

$$\text{Hamilton word on }T\longmapsto\text{its full }A\text{-orbit label word}.$$

Different label words certify different path orbits. Equal label words do
**not** certify equal orbits. Edges and Hamiltonicity are lost in this map,
so every finite control retains the source, paths, anti-isomorphism and
targets. No new intrinsic tournament is invented on the label words.

## 2. Odd full orbits and the anti-involution pairing rule

Every tournament automorphism group has odd order: an involution would
transpose a pair and reverse its arc. Thus each full vertex orbit has odd
size. A full orbit is a union of rotation shores, so its shore multiplicity
is odd too.

For a self-converse source the automorphisms and anti-automorphisms form a
group of order $2|A|$. Cauchy's theorem supplies an involution $\alpha$;
it lies in the anti-coset. An anti-automorphism fixes at most one vertex,
since two fixed vertices would force their arc to reverse. Consequently:

* at even $n$, $\alpha$ has no fixed vertex and pairs **all** full odd
  orbits with equal-size partners;
* at odd $n$, $\alpha$ has exactly one fixed vertex and stabilizes exactly
  one full odd orbit; every other full orbit has an equal-size partner.

A stabilized odd orbit would itself contain a fixed vertex, which explains
the second assertion. These rules are necessary, not a construction of
anti-isomorphisms from a size partition.

In particular, the possible shore-multiplicity partitions at even $q$ are

| $q$ | Paired odd full-orbit multiplicities |
|---:|---|
| 4 | $1+1+1+1$ |
| 6 | $1+1+1+1+1+1$, or $3+3$ |
| 8 | $1+1+1+1+1+1+1+1$, or $3+3+1+1$ |

For example, $3+1$ at $q=4$ would stabilize two unequal odd orbits, which
is impossible at even order. This is a full-orbit statement, not a score
partition heuristic.

## 3. A singleton full orbit separates two path families

Let $w=(w_0,\ldots,w_{q-1})$ record the full $A$-orbit containing each
rotation shore. Suppose one label occurs exactly once, and $m\ge3$.
Then

$$\boxed{\Delta\ge2\bigl(q-\gcd(2,q)\bigr).} \tag{1}$$

**Native path construction.** For every starting shore $s$ there is a
Hamilton word that traverses each entire shore, in cyclic shore order.
Start a Hamilton path in the first strong shore. From its last vertex use
the actual marked matching into the next shore, and choose a Hamilton path
there starting at the supplied endpoint (cut a Hamilton cycle before it).
Continue. No final return edge is needed. These paths retain the literal
matching endpoints; the construction does not assume uniform cross-arcs.

The $q$ chosen cycle cuts have label words
$(\operatorname{rot}_s w)^m$. The $q$ whole-shore paths have label words
$\operatorname{expand}_m(\operatorname{rot}_s w)$. The singleton label makes
$w$ aperiodic, so each family contains $q$ distinct label words and hence
$q$ distinct full path orbits.

Every anti-isomorphism induces the same permutation $\tau$ of the full
$A$-orbit labels: the quotient of any two anti-isomorphisms is in $A$.
If two starting positions $s,t$ in either family are blue, applying the
fixed $\tau$ to the reversed label words gives

$$\operatorname{rot}_{2(s-t)}w=w.$$

Aperiodicity forces $2(s-t)=0\pmod q$. Thus at most $\gcd(2,q)$ of
either family's orbits can be blue. Finally, the singleton full orbit is
one run of length $m$ in every block path, whereas its vertices are separated
by $q-1$ other vertices in every cut path. The two families are therefore
disjoint even as full path orbits. This proves (1).

The two-family construction and its prime-$q$ necklace version were
independently developed in the concurrent order-15 lane; this packet uses
the singleton label to state the arbitrary-$q$ form and the even-order
consequences. It requires no independence assumption on the cross-edge
orientations and no enumeration of all Hamilton words.

## 4. Three families when every shore is a full orbit

Assume now that every rotation shore is a full $A$-orbit, $q\ge3$, and
odd $m\ge5$. Then

$$\boxed{\Delta\ge2\bigl(q-\gcd(2,q)\bigr)+4.} \tag{2}$$

The two families from §3 contribute the first term. Choose a directed
four-cycle $S$ in one shore $A_0$. It is a proper subset. Deleting $S$
leaves a strong tournament by the inherited matching-chain argument.
Against another shore $B$, every vertex of $A_0$ has the same number $k$
of outgoing arcs, so $\#(S\to B)=4k$. All $m$ vertices of $B$ cannot
have exactly two incoming arcs from $S$, since $4k=2m$ contradicts odd $m$.
There are therefore at least three exits or at least three entries at some
terminal. The inherited half-supported four-prefix lemma supplies four
black path orbits, including distinct mates.

Each chosen prefix path begins with four vertices of one full orbit,
then leaves that orbit, and later returns to its remaining $m-4>0$ vertices.
It therefore has both a run of length at least four and a full orbit split
over at least two runs. Reversal and the anti-isomorphism preserve these
properties, with orbit labels permuted. Such paths and their mates cannot
belong to the cut family (all runs have length one) or the block family
(each full orbit occurs in exactly one run). All contributions in (2)
are thus disjoint.

For $q=3$ or $4$, (2) is at least eight. This is the extra four-orbit
contribution missing from the two-family estimate alone.

## 5. Two equal full orbits: an eight-prefix obstruction

Suppose $T$ has exactly two full vertex orbits $A,B$, of a common odd size
$r\ge9$. Then

$$\boxed{\Delta\ge8.} \tag{3}$$

Each induced orbit is regular and strong. Full-orbit transitivity makes
the outgoing degree from $A$ to $B$ constant, say $k$; the incoming degree
from $A$ is constant on $B$ and equals the same $k$ by double counting.
Strongness of $T$ gives $1\le k\le r-1$. Thus every surviving vertex of
$A$ has arcs in both directions to the intact strong $B$, and deleting
**any proper subset** of $A$ leaves a strong tournament.

Choose an internal directed eight-cycle $S\subset A$. Its complement is
strong. Its cross incidence is $8k$, which cannot equal $4r$ because $r$
is odd. Some terminal therefore has at least five exits or entries.
Since $n=2r\ge18>16$, the two eight-prefix supports in the inherited
mate argument are disjoint, and that lemma gives $2(5-1)=8$ black orbits.
This proof does not require either full orbit to be cyclic or to coincide
with a rotation shore.

## 6. Exact consequences for the remaining orders

The inherited necessary marked-cycle rows give:

| Order | Possible $(q,m)$ | New consequence |
|---:|---|---|
| 18 | $(6,3)$ | §2 gives six singleton shores or two full orbits of size 9. The first violates (1); the second violates (3). **Excluded.** |
| 20 | $(4,5)$ | §2 makes all four shores full orbits; (2) gives at least eight. **Excluded.** |
| 24 | $(8,3)$ | Both §2 partitions have a singleton shore orbit; (1) gives at least twelve. **Excluded.** |
| 28 | $(4,7)$ | Same four-full-orbit argument as order 20. **Excluded.** |
| 21 | $(3,7)$ or $(7,3)$ | In the first type, non-VT forces three singleton full orbits and violates (2). In the second, every non-VT odd partition of seven contains a singleton and violates (1). Hence **necessarily VT**. |

The $q=3,m=5$ non-VT part of order 15 is likewise excluded by (2).
The $q=5,m=3$ part is handled by the concurrent prime-necklace proof, and
the previously frozen VT15 bank handles the VT part. Those latter inputs
are not silently replaced by our seven controls.

## 7. Exact controls, consumer and limits

[check.py](check.py) and [inputs.json](inputs.json) retain seven literal
sources, at $(q,m)=(3,5),(3,7),(4,5),(4,7),(6,3),(8,3),(2,9)$.
They are **controls of the proofs, not all sources of those orders**.
Every source has a marked Hamilton cycle, an actual anti-involution, and a
rotation transitive on each shore. Distinct full scores certify that the
shore partition in these controls is exactly the full automorphism-orbit
partition. The two-nine-orbit control exercises (3); its existence does not
claim that a deficit-six source can have the same marked-cycle index.

The bank retains every cut and whole-shore Hamilton word, the directed
prefix cycle, its literal exits, common suffix, and source anti. The consumer
reconstructs the mates, authenticates full ordered adjacency keys, checks
blackness directly, and verifies every edge of each complement-except-path
target. Different keys are exactly different full automorphism path orbits
on the retained source. Targets are not inferred from a scalar orbit count.
All singleton necklaces coming from odd set partitions through eight shores
are checked against every possible anti-label action arising from a reversed
rotation. The parity partition tables are independently regenerated.

The producer's capped Hamilton search is available only under `--discover`.
Default, `--write` and `--verify` modes consume the supplied words without
search. Type aliases, wrong sources/anti-actions, duplicate cycle vertices,
wrong terminals, altered paths and altered targets are rejected. The actual
six-center is exhaustively checked as a positive $\Delta=6$ boundary.
No conclusion follows from merely failing to find a word.

```text
python -B 05-knowledge/results/center_incidence_20261022/check.py --discover
python -B 05-knowledge/results/center_incidence_20261022/check.py --write
python -B -O 05-knowledge/results/center_incidence_20261022/check.py --verify
```

Normal replay: **19,331 always-active checks**, LF receipt SHA256
`13effe30119b6397fb58f3e247aa9e53f1bc292cbc91bcb8c6c404f89ff56b1a`.
Independent proof challenges of the singleton-necklace, anti-orbit pairing,
and two-equal-orbit arguments: **PASS**. Independent final full-note/core API
audit and optimized saved-output comparison: **PASS**, reproducing all
**19,331** checks and the exact hash above. The audit separately checked the
anti-involution pairing, native block joins, split/run separation of prefix
mates, equal-odd-orbit deletion, source bindings and no-census scope. The
four owned files are frozen.
