# Twin primes on two Berggren rays: an exact edge-to-node dictionary

**Status: PROVED elementary deductions; independently audited PASS.**
The [independent audit](../arithmetic/AUDIT_TREE.md) adds 2,892 gates through
a separate side-enumeration and factorization path. Computation is
FINITE-EXACT. Classical Euclid/Berggren structure and the existing repository
charts are inherited; no priority claim is made. Twin-prime infinitude remains
OPEN. This note identifies a subset and an actual graph map, not the whole
ternary tree with twin primes.

## 1. Inheritance, conventions, and board

Use ordered primitive triples `(a,b,c)=(odd leg,even leg,hypotenuse)` and
Euclid parameters `m>n>0`, `gcd(m,n)=1`, `m+n` odd:

    Phi(m,n)=(m²−n²,2mn,m²+n²).

Closest proved mechanisms:

- [THM-3333, Gaussian-square/Farey carrier](../../../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md)
  retains content and ordered legs. A Farey edge is not generally a Berggren edge.
- [THM-3334, parabolic spine](../../../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
  already proves the U-spine formula used below.
- [THM-3756, odd-square chart](../../../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md)
  already proves the lossless two-root chart. Keeping only the outer root
  confuses `(5,12,13)` and `(15,8,17)`.
- [THM-2422, twin-center ancestry](../../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md)
  gives centers `{4} union {6k:6k−1,6k+1 prime}`. Its additive ancestry
  statement through center `10^8` is FINITE-EXACT; the all-center claim is OPEN.
- [Unique-six theorem A](../../aliquot_power_124_20261026/aliquot/README.md)
  says `s(p²)` is perfect for prime `p` exactly when `p=5` and `s(p²)=6`.

The canonical hostile is loss of the source predicate under tree descent.
The corrected near miss is calling any countable subset of a tree a tree
isomorphism. The least-used sidecars here are the inner square root, branch
word, and total prime-factor multiplicity `Omega`, rather than distinct-factor
count `omega`.

| Live concept | Exact representation | Operation and decisive boundary |
|---|---|---|
| Twin pair | ordered primes `p<q`, `q−p=2` | retain the difference, not only product |
| Prime leg | factorization `(m−n)(m+n)` | primality forces one factor to be 1 |
| Berggren tree | branch-labelled Euclid pairs | only U preserves the prime-leg possibility |
| Semiprime fibre | coprime odd factor splits | the trivial and nontrivial splits are different nodes |
| Aliquot anchor | `s(p²)=p+1` | the perfect center is uniquely 6 |

Use the inherited branch names

    U = [ 1 −2 2;  2 −1 2;  2 −2 3],
    A = [ 1  2 2;  2  1 2;  2  2 3],
    D = [−1  2 2; −2  1 2; −2  2 3].

On our Euclid parameters these are respectively

    U(m,n)=(2m−n,m), A(m,n)=(2m+n,m), D(m,n)=(m+2n,n).

The root is `R=(3,4,5)=Phi(2,1)`. This is the `m>n` convention, not the
swapped convention in THM-3357. Branch words act successively from the root.

## 2. The two-root chart and the twin image

For positive coprime odd integers `x<y`, define

    Xi(x,y)=(xy,(y²−x²)/2,(y²+x²)/2).                 (1)

This is a bijection onto the ordered primitive triples, with inverse

    x=sqrt(c−b), y=sqrt(c+b).                         (2)

For completeness, the two odd numbers `c−b,c+b` are coprime and have product
`a²`, hence are individually squares. Their positive roots are coprime and
odd; solving recovers (1). Conversely (1) is Pythagorean; a common odd prime
would divide both roots, while 2 cannot divide its odd leg. This reproduces
the inherited THM-3756 proof and fixes the normalization.

The decisive height identity is

    c−a=(y−x)²/2.                                   (3)

Thus `y−x=2` if and only if `c−a=2`. The identity concerns the **odd leg**;
it is not the U-spine height `c−b=1`.

**Theorem T1 (twin-pair subset bijection).** For a twin-prime pair `p,p+2`,

    F(p,p+2)=(p(p+2),2(p+1),(p+1)²+1)                (4)

is a primitive triple. This is a bijection from twin-prime pairs to primitive
triples satisfying

    c−a=2 and Omega(a)=2,                            (5)

where `Omega` counts prime factors with multiplicity. The inverse is

    p=b/2−1, p+2=b/2+1,

or equivalently the two roots in (2).

*Proof.* Apply (1) to the odd coprime pair `p,p+2`. Conversely (3) makes the
roots consecutive odd integers. If `x=1`, then `y=3` and `a=3`, contradicting
`Omega(a)=2`. Therefore both factors of `a=xy` are greater than 1; total
prime-factor count two forces each to be prime. The inverse is explicit. QED.

Equivalently, all triples with `c−a=2` are `Phi(m,1)` for positive even
`m>=2`; (5) selects those with `m−1,m+1` prime. This follows directly from
`c−a=2n²`. Primality itself is not needed for the primitive construction:
every positive even `m` gives a primitive triple.

**Hostile controls.** The root `(3,4,5)` has height two but roots `1,3`, so
is not a twin image. `(63,16,65)` also has height two and `omega(63)=2`,
but `Omega(63)=3` and roots `7,9`; distinct-factor count is insufficient.
The half-sized proposal `m=(p+1)/2,n=1` fails the parity gate when `p=5`:
it gives the nonprimitive `(8,6,10)`, whose unordered primitive reduction
duplicates the root from `p=3`. Formula (4) uses `m=p+1`, not that half.

## 3. Prime legs rigidly select one whole branch

For every odd integer `t>=3`, put

    P(t)=(t,(t²−1)/2,(t²+1)/2).

The inherited U-spine formula says

    P(t)=U^((t−3)/2) R,   U P(t)=P(t+2).             (6)

**Theorem T2 (all prime-to-prime tree edges).** A primitive triple with prime
odd leg `p` is exactly `P(p)`. Every Berggren edge with prime odd legs at
both endpoints is exactly `P(p) -> P(p+2)` for a twin pair.

*Proof.* In the unique Euclid representation, primality of
`a=(m−n)(m+n)` forces `m−n=1,m+n=p`. This gives `P(p)`. The three children
of that parent have odd legs

    p+2, p(2p+1), p(2p−1),                           (7)

respectively. The last two are composite, so the first child is the only
possible prime-leg child, and it is prime exactly for a twin pair. QED.

In particular, this is an exhaustive statement about prime-to-prime edges
of the full ternary tree, not just a search on a preselected ray.

## 4. An actual graph isomorphism, with its exact loss

For `k>=0`, put

    Q_k=D^k R=((2k+2)²−1,4k+4,(2k+2)²+1).          (8)

Let `H` be the undirected subgraph induced in the Berggren tree by triples
whose odd leg is prime. Let `J` be the subgraph induced by triples satisfying
(5). For a graph `G`, its line graph `L(G)` has one vertex per edge of `G`,
with two vertices adjacent exactly when those edges share an endpoint.

**Theorem T3 (edge-to-node isomorphism).** The following map is a graph
isomorphism, including any isolated vertices of `L(H)`:

    L(H) -> J,
    {P(p),P(p+2)} -> Q_((p−1)/2)=F(p,p+2).          (9)

*Proof.* T1 and T2 give the vertex bijection. Write the U-edge index as
`j=(p−3)/2`. Its endpoints have words `U^j,U^(j+1)` and its image has word
`D^(j+1)`. Two such U edges share an endpoint exactly when their indices
differ by one; two such D nodes are adjacent exactly in the same case.
No off-ray edge joins two D-ray vertices because these vertices have their
unique tree ancestry along D. Thus adjacency is preserved and reflected.
QED.

This map extends before primality filtering: the line graph of the complete
U ray is the D ray with its root deleted, with the same index shift. It
intertwines advancing the U edge by one with applying D to the image. A
coordinate formula using the lower endpoint alone is

    F=(t²+2t,2t+2,t²+2t+2),                         (10)

and the endpoint difference is

    P(t+2)−P(t)=(2,2t+2,2t+2).

Both formulas preserve the exact twin predicate only after its explicit
primality test. The unfiltered shift does not preserve primality.

**Sharp component boundary.** Three odd primes `p,p+2,p+4` force `p=3`:
one is divisible by 3, and the two later positions cannot equal 3. The next
candidate 9 is composite. Hence `H` has exactly one component with two
edges, `P(3)—P(5)—P(7)`; all other components are single edges or isolated
vertices. Consequently `J` has exactly one edge,

    (15,8,17) — (35,12,37),

and every other vertex of `J` is isolated. This classification does not
assert whether there are finitely or infinitely many isolated vertices
of `J`; that is precisely twin-prime infinitude.

Thus (9) is a genuine graph isomorphism of typed subgraphs after the line
graph operation. It is not an isomorphism with the full ternary tree, nor
does it identify vertices of `H` with twin pairs. Isolated prime vertices
of `H` disappear under the line graph operation. The edge incidence and the
root offset are the needed sidecars.

## 5. Semiprime fibres and the aliquot center

For distinct odd primes `p<q`, the odd-leg fibre `a=pq` has exactly two
primitive triples. In the two-root chart, coprime factor splits of `pq` are
exactly `(1,pq)` and `(p,q)`. The first gives `P(pq)` on the U ray; the
second gives `Xi(p,q)`. This also specializes the classical fixed-leg count
in Tripathi, Theorem 2. Its positive nontrivial split has height

    c−a=(q−p)²/2.

Twin pairs are exactly the minimum-height-two members among these nontrivial
semiprime splits. The semiprime predicate alone does not imply twinhood:
`Xi(3,7)=(21,20,29)` has height eight. Nor does an odd-leg value specify a
node: for `a=15`, the two nodes are `(15,112,113)` and `(15,8,17)`.

Let `h=p+1` be a twin center. Formula (4) becomes

    Q=(h²−1,2h,h²+1),
    c−b=(h−1)², c+b=(h+1)².

Since the square roots are prime, the corresponding aliquot identities are

    s(c−b)=h=b/2,
    s(c+b)=h+2=b/2+2,
    s(a)=1+p+(p+2)=2h+1=b+1.                         (11)

The last identity uses distinct primes and is not valid for arbitrary
odd-leg factorizations. Formula (11) supplies exact arithmetic maps, not
only a visual analogy between a branch and an aliquot graph.

By the inherited unique-six theorem, the lower square's aliquot image in
(11) is perfect exactly for `p=5`, so the unique twin image with perfect
half-even-leg is `(35,12,37)`, with center six. Its lower square is 25.
The preceding image `(15,8,17)` has center four, and `s(9)=4`; the two
centers must not be conflated. No classification of other aliquot fibres
or proof of twin-prime infinitude follows.

There is also an exact two-sided refinement. A twin image has a perfect
aliquot image of either square `c−b` or `c+b` exactly when that square is
`25`, by the same unique-six theorem. A twin pair containing its root `5`
must be `(3,5)` or `(5,7)`. Thus the two resulting triples are precisely
the endpoints of the unique edge of `J`: `25` is the upper square of
`(15,8,17)` and the lower square of `(35,12,37)`. This identifies the shared
square and its orientation; it does not identify the two twin centers.

## 6. Connection contract and remaining question

Source: twin pairs, or equivalently prime-to-prime Berggren edges.
Target: the height-two, `Omega(a)=2` D-ray locus.
Map: (4), with inverse (2) or `b/2 +/-1`.
Preserved: positivity, primitivity, the twin predicate, branch address with
one index shift, and edge-sharing adjacency under the line graph operation.
Lost by replacing the source graph with its edge set: isolated prime-leg
nodes. Lost by retaining only `a`: the factor split and therefore the node.
Sidecars: ordered legs, two roots, multiplicities, and the branch word.
Cheapest tests: root exclusion, `63` versus `15`, `(3,7)`, and the first
failed prime-leg continuation `P(7)->P(9)`.

The remaining infinitude problem is exactly whether the D ray contains
infinitely many nodes with `Omega(a)=2`. The unfiltered ray is infinite,
but its elementary recurrence does not establish that prime predicate.
The normalized twin-center addition from THM-2422 is a separate operation:
no preservation of it by (9) is asserted.

## 7. Verification

`check_tree.py` uses explicit runtime gates, no producer imports and only the
standard library. It checks a full finite Euclid universe and a complete
ternary prefix, compares independent tree and Euclid paths, checks the two
factor-split fibres, and records the exact restricted graphs and hostile
controls. Normal and optimized runs must match `normal.json` byte for byte:

```bash
python3 -B check_tree.py
python3 -B -O check_tree.py
```

The finite bounds and counts are printed in the receipt. All-parameter
statements rest on the proofs above, not on the finite census. Source and
producer hashes, source-reading boundaries, and audit lineage are recorded
in `SOURCE_MANIFEST.json` and `SOURCES.md`.
