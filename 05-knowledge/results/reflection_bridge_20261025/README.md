# Two reflections, a faithful square, and the seven-support boundary

2026-10-25. **PROVED elementary** where stated; **CITED** geometric conventions and recovered predecessors are identified; **FINITE-EXACT** receipts supplement the proofs. No literature-priority claim and no transfer of an open tournament or LRC problem.

The ten divisors really do have graph and geometric models. What cannot coexist is **all ten states + ordinary graph complement + involutive vertex/path reflection**. Requiring divisibility adds a second obstruction. We recover the existing action-level repairs, give an inclusion-faithful square model with a necessary complement twist, prove a sharp exponent-three exception under the weaker monotonicity requirement, and merge the seven-support/F1 work without identifying incompatible incidence structures.

## 1. Inheritance and portfolio

Closest proved mechanisms: [THM-4733, divisor plane versus tiling cube](../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md); the [minimal twelve-state resolution](../divisor_reflection_resolution_20261021/REFLECTION_RESOLUTION.md); the [Oct4 divisor-ideal dictionary](../collatz_four_slot_compression_20261004.md); and [Session59's explicit projective point conventions](../session59_f1_geometry_and_information_20261025/cyclotomic-support/CYCLOTOMIC_SUPPORT_DYNAMICS.md). The eight-state repair and its 32 equivariant bijections were also already in the [Klein dictionary, section 5](../black_detours_closed_blue_lines_klein_dictionary_20261020.md).

Canonical hostile: `pq,pr` are fixed by `CQ`, whereas graph complement composed with an involutive vertex reversal has no fixed graph. Corrected near miss: an equivariant bijection of finite sets need not respect divisibility, graph incidence, or directed-path typing. Least-used sidecars: ordered repeated-prime slots, opposite coordinate flags, and the distinction between rank-zero points and all support points.

**Anchor:** the user's reflection map and its strongest preserved predicates. **Niche:** the seven-point boundary quotient and its failed incidence. **Wildcard:** a genuine monoid-projective carrier and the sharp change at exponent three. LRC(14) remains an open repository anchor; this packet supplies no reduction of it.

| Live concept | What the present pull/computation changes |
|---|---|
| Klein fixed points | The old obstruction survives even for many-to-one maps; extra graph vertices do not help |
| Divisibility and gcd/lcm | Prefix ideals preserve them exactly; raw complement cannot preserve the order |
| Marked tournament predicates | The exponent-three threshold map reads closure, strength, and blue symmetry explicitly |
| Boundary quotient | Seven arises by contracting a complement orbit, losing order and some incidence |
| F1 point conventions | Seven supports and three rank-zero points coexist in one projective factor; ten comes from a product's rank-zero locus |
| Coordinate flags and graph host | These sidecars specify which ranks, complements, and incidence relations actually survive |

## 2. The original picture and the exact no-map theorem

For `N=p²qr`, write a proper nontrivial divisor as `(a,b,c)`, where `0<=a<=2` and `b,c in {0,1}`. The coordinates `f=a-b-c`, `g=b-c` give

```text
g \ f     -2     -1      0      1      2
  1              q      pq     p²q
  0       qr     pqr           p      p²
 -1              r      pr     p²r
```

The g-rows have sizes **3-4-3**; the negative/zero/positive f-sectors have sizes **4-2-4**. Divisibility ranks `a+b+c=1,2,3` separately have sizes 3-4-3, but are a different partition.

The axis reflections are `Q(a,b,c)=(a,c,b)` and `CQ(a,b,c)=(2-a,1-c,1-b)`; their product is divisor complement `C(d)=N/d`. Their fixed counts are 4,2,0, respectively.

In a graph on ordered vertices, let K be full graph complement and R be vertex reversal. Via the ordered graph/marked-path dictionary, K reverses nonpath arcs, R takes converse and reads along the reversed path, and KR reverses path arcs and then reads along that reversed path. The arc modification is essential. We repaired a historical shorthand `(T,P)->(T,P^rev)` in the Klein dictionary and recorded the mechanism in [MISTAKES.md](../../../01-canon/MISTAKES.md).

**Recovered obstruction:** R fixes an unordered vertex pair, so KR would have to complement its own edge bit. Thus KR has no fixed graph. Since CQ fixes pq and pr, there is no map intertwining `(C,Q)` with `(K,R)`, even without injectivity. Assigning either fixed-point-bearing axis reflection directly to K also fails. This is about genuine full graph/path charts, not an arbitrary Boolean subset universe.

Existing repairs remain valid: discard pq,pr for an eight-state carrier, or split each into two marked states for a minimal twelve-state carrier. Neither gives the full ten-divisor order as graph inclusion. See the predecessor links above.

## 3. A square preserving the whole divisor order

Use the square host with edges `01,12,23,03`, marking

    p1=12, p2=03, q=01, r=23.

Map `p^a q^b r^c` to the first a p-edges together with the indicated q,r edges. This is injective, and

    d|e iff I(d) is contained in I(e),
    I(gcd(d,e))=I(d) intersection I(e),
    I(lcm(d,e))=I(d) union I(e).

For gcd/lcm include the unit and N in the ambient lattice. Four edge slots are necessary and sufficient. Let J swap the two p-slots, let R swap q,r, and let K_E complement **inside the square host**. Then

    I C = K_E J I,       I Q = R I.

Both J and R are actual square reflections. The missing p-slot reversal is exactly the cost of choosing an inclusion-faithful representative of each repeated-prime occupancy. The omitted diagonals explain why this evades the full graph-complement obstruction. [Full proof, sharp all-k edge budget, and quotient loss](divisor_ideal_graph.md).

This twist is necessary: `p|C(p)`. Any monotone encoding with raw complement forces `I(p)` to be empty, which would lie below `I(q)` although p does not divide q. Thus **no order embedding** of the full proper-divisor set of `p^kqr`, k>=2, can carry C to ordinary complement.

## 4. The sharp exponent-three exception

We can weaken order embedding to an injective map preserving every divisibility implication. Such a map into full ordered-graph charts, with `C->K,Q->R`, exists for `p^kqr` **exactly when k=1 or k=3**. Even k fails the action test. For k>=4, both p and p² lie below their complements and would both have to become the empty graph. The squarefree k=1 case is the direct three-bit map.

At k=3 use four graph vertices and these six increasing edge predicates:

    (01,02,03,12,13,23)
      = ([a+2b>=3], b, [a>=2], [a+b+c>=3], c, [a+2c>=3]).

They are equivariant and recover all fourteen source divisors. The marked five-vertex tournament lift then has:

- closed marked path iff `a>=2`;
- strong tournament iff `a>=2 or b=c=1`;
- blue marked path iff `b=c`;
- transitive tournament iff the source is p.

It also reads `g=b-c` directly from edges02 and13. It creates some extra inclusion relations, so is not an order embedding. [Proof and exact controls](graph_bridge.md).

## 5. Merging 3-1-3 with seven and F1

The squarefree seven inside the ten are `p,q,r; pqr; pq,pr,qr`. Let B be the eight subsets of `{p,q,r}`. Ordinary complement acts on B, not on its seven nonempty subsets. Identify empty with full, calling the resulting class star. The descended involution fixes star and swaps each singleton with its opposite pair: **3+1+3**, cycle type `1·2³`. The same quotient identifies the empty graph and triangle in the eight graphs on three labelled vertices.

This is an exact set-with-involution quotient. It is not an ordered-lattice quotient, and the induced permutation of the seven Fano points fails four of the seven lines. The old Keller fragment's 3,1,3 inverse-fibre counts are a separate source: its first-layer coordinate reflection does not iterate to the seven second-layer leaves. [Maps, witnesses, provenance, and strongest survivors](seven.md).

A positive geometric merge uses the pointed-monoid model

    U=(P² x P¹ x P¹) minus its two opposite coordinate sections.

It has **ten rank-zero points** indexed by the divisor triples. Coordinate reversals give C; swapping the last two factors gives Q. These are actual commuting automorphisms of the model. Its P² factor simultaneously has **seven support/morphism points but only three closed/rank-zero points**. The seven-support reflection obtained by boundary contraction is not one of its coordinate permutations: their fixed counts differ. The distinction follows the established pointed-monoid convention, not an ordinary field of cardinality one. [Construction, flags, counts, and source pins](geometric_carrier.md).

For a cheap typing check, U has 61 literal A1-valued points, ten rank-zero points, and finite-field count `(1+Q+Q²)(1+Q)²-2`. The ten-label rank polynomial `3t+4t²+3t³` is different. A shared value at one cannot replace the missing map.

## 6. Verification and next boundary

Run the four standard-library checkers from this folder:

```text
python3 check_divisor_ideal_graph.py
python3 graph_bridge_check.py
python3 seven_check.py
python3 geometry_check.py
```

Repeat with `python3 -O`; gates are explicit exceptions. [VALIDATION.md](VALIDATION.md) records scopes, independent audits, and receipt comparison. The unrestricted conclusions rest on the proofs, not on the finite ranges.

The next precise question is whether a useful graph/tournament predicate survives the square's twisted complement beyond the decoded divisor order. Plain complement/order preservation is now closed by obstruction; additional numerical coincidences do not reopen it. The geometric carrier also leaves a concrete question about flag-compatible correspondences, since C exchanges standard and opposite flags rather than the chosen cells themselves.
