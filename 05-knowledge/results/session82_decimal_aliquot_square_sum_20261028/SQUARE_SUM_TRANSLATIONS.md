# Recovering square-sum meaning: an exact translated-triple / aliquot-fiber bridge

Date: 2026-10-28. Source pin: `fbc89240231927f0eb6d1039d4006fd60d8b08d6`.
Status: **PROVED elementary classification and infinite family; FINITE-EXACT full n=1000 classification.** No novelty claim. Work is outside the publisher tree; publication belongs to the parent lane.

## 1. Inheritance and correction boundaries

The prior square-sum problem is the simple graph on positive integer labels, with distinct x,y adjacent iff **x+y is an integer square**. Its interval version asks for Hamiltonian paths through every label 1,...,N. It is not asking whether entries are squares or fourth powers in a quadratic integer ring. Gerbicz's own Jan 21, 2018 comment on [OEIS A090461](https://oeis.org/A090461), checked directly today, records paths for every N>=25 and cycles for every N>=32. This remains **CITED**, not a fresh audit of the construction.

Relevant recovered, corrected routes:

* `05-knowledge/results/session32_square_monster_20261023/square-sum/README.md`: exact small path census and endpoint-aware suppression. Degree two does NOT force both incident edges when the vertex may be a path endpoint. Q23 has a path ending at degree-two vertex22. Connectivity is not Hamiltonicity.
* `05-knowledge/results/graceful_square_exchange_20261019.md`: the square-sum Q28->Q29 extension needs at least three old-edge deletions (finite exact), despite existence of both paths. A nonempty rank-gap shift destroys old square edges, unlike its graceful-label counterpart.
* Working-copy source `session70-graph-labeling/REPORT.md` (title: *Named attachments, square-sum paths, graceful fibers, and Petersen quotients*) recovers graceful upper-tail two-color connectivity, while arbitrary-fixed-color connectivity is refuted. Its sum/difference sign gauge is not a graceful-tree proof; absolute differences require the mirror information.
* `05-knowledge/results/ap_square_pythagorean_20261019.md`: every positive increasing AP square-sum triangle is uniquely
  `2g²(C²−4AB, C², C²+4AB)`, with primitive positive Pythagorean legs A>B, g>=1, and C²>4AB. Edge roots are `2g(A−B),2gC,2g(A+B)`. This is an actual reversible map, not the claim that any three Pythagorean-looking integers are a square-sum triangle.
* `05-knowledge/results/square_sum_ppt_bridges_20261019.md`: signed endpoint contacts concatenate sparse AP petals. Signs, rational square condition, and common denominator scaling are retained. This is not interval coverage.
* `05-knowledge/results/square_sum_prescribed_edge_absorbers_20261019.md`: local absorption exists, but simultaneous use needs actual host edges plus a matching condition. Numerical room alone does not supply those hosts.

Closest mechanism: square-root factorization / square reflections. Hostile: an edge plus an isolated vertex. Corrected near miss: all-order interval existence is not an open conjecture, and endpoint-free degree forcing is invalid. Least-used sidecar here: exact factor pair (s−r,s+r), which makes translation enumeration finite and complete.

Concept board: pair-sum graph; common translation; difference of squares; aliquot fibers; Pythagorean normalization; graceful sign gauge. The first four yield the new concrete result; the last two supply explicit non-implication controls.

## 2. The original triple

For n>=2 put a=n²+1, b=n+1, c=n²+n. Then

* b+c=(n+1)², so b--c is always an edge.
* a+b=n²+n+2 lies strictly between n² and (n+1)² for n>=2; hence a--b never exists.
* a+c=2n²+n+1 is a square only for special n, not identically.

At n=1000 the last sum is 2001001, strictly between 1414²=1999396 and 1415²=2002225. Thus the specified triple {1000001,1001,1001000} spans exactly one edge, not a square-sum chain. It is not an AP triangle and is not itself a Pythagorean triple of lengths.

Original aliquot sums s(x)=sigma(x)−x are respectively

    s(1000001)=10003, s(1001)=343, s(1001000)=2143960.

These are three distinct fibers. In particular the fixed b--c square edge does not force equality of aliquot sums.

## 3. Exact classification of all common translations

Fix a prospective center m and endpoints x<y chosen from the original triple. Write d=y−x>0. A common integer translation t produces the path

    (x+t)--(m+t)--(y+t)

with positive labels iff there is a factor pair u*v=d, 1<=u<=v, satisfying

    u≡v (mod2), r=(v−u)/2, s=(v+u)/2,
    r²−m−x≡0 (mod2), t=(r²−m−x)/2,
    min(a,b,c)+t>0.

Necessity follows by taking the positive edge roots r<s: their squared difference is d, so u=s−r and v=s+r. Sufficiency follows by substitution. Every possible t is obtained, with no height cutoff. This classifies common translations, not arbitrary separate changes of labels.

### Universal forbidden center, and no translated triangle

For every n>=2 there is no positive-label common translation with center b. Here endpoints a<c differ by d=n−1. The smaller edge root satisfies

    r <= (d−1)/2 = (n−2)/2,

because s−r>=1 and s+r<=d. But b+t>0 implies

    r²=a+b+2t > a−b=n²−n,

contradicting (n−2)²/4 < n²−n. Therefore **no positive-label common translation of this triple is a square-sum triangle**, since a triangle would in particular contain both edges centered at b.

At n=1000 center a is also impossible: its two edge sums have even and odd parity respectively. Their being squares forces residues 0 and1 modulo4; their difference is c−b=n²−1≡3 (mod4), a contradiction. This argument applies to every even n.

### An infinite positive family

Let n=8k, k>=1, d=n(n−1),

    r=d/4−1, s=d/4+1,
    t=(r²−(n+1)²)/2.

Both roots are odd and t is an integer. Since r>=n+1 for n>=8, t>=0. The exact identities

    (b+t)+(c+t)=r²,
    (a+t)+(c+t)=s²

give a square-sum path centered at c+t for every k>=1. Translation preserves the three pairwise label differences but not their original edge roots, arithmetic factors, or aliquot values. No universal prime or aliquot-fiber assertion accompanies this family.

## 4. Complete n=1000 result and a genuine aliquot intersection

The factor-pair classifier yields exactly **24** integer translations with all labels positive and a square-sum spanning path. All are centered at c+t. The enumerated domain includes every t>=−1000, not just positive offsets; in fact all 24 survivors have t>0. The smallest possible translation is

    t=177612,
    path 178613 -- 1178612 -- 1177613,
    edge sums 1357225=1165² and 2356225=1535².

Its factor-pair certificate is u=370, v=2700, uv=999000. Both endpoints are prime, verified by exhaustive trial division through their integer square roots. Hence

    s(178613)=s(1177613)=1,
    1178612=2²*23²*557,   s(1178612)=981406.

Of the complete set of 24 allowable translated paths, this is the **only** one whose two endpoints both lie in s⁻¹(1). This is a finite exact intersection of a graph constraint and an aliquot-fiber constraint. It is not a proof of an infinite prime-pair family, an aliquot cycle, or a universal transport of aliquot fibers.

For any two distinct primes p,q, equality s(p)=s(q)=1 is automatic, but a common square-sum neighbor h is additional data: p+h=r² and q+h=s² imply q−p=s²−r². Here the same complete difference-of-squares test constructs the neighbor and tests positivity. This is the exact preserved predicate, rather than numerical resemblance.

### Lossless recovery requires labels, not just their aliquot observations

For every c-centered translated member with endpoints p=b+t<q=a+t, their difference d=q−p=n(n−1) recovers

    n=(1+sqrt(1+4d))/2, t=p−(n+1).

The center must then be h=n²+n+t. This is a left inverse of the label-preserving family map; it supplies a precise membership test (positive integer n>=2, integral t, correct h, positive labels, and square edge sums). At our prime pair sqrt(1+4*999000)=1999, so n=1000 and t=177612 are recovered exactly. Passing only to aliquot values (1,981406,1) discards these endpoint labels; no source recovery from that tuple alone is asserted. The center aliquot value is an exact computed observation, not a demonstrated invariant.

## 5. Typed connection ledger

| Source -> target | Map and preserved predicate | Lost information / required sidecar |
|---|---|---|
| Ordered labels + chosen center -> translated square-sum path | Factor d=(s−r)(s+r), recover t | Must retain parity, positivity, center and roots; finite factor list is complete |
| This translated path -> aliquot fiber | Evaluate s at its endpoints; both prime here | Translation does not preserve factorization or s; primality separately certified |
| Pythagorean triple -> AP square-sum triangle | Inherited quadratic map in §1 | Needs primitive scale and positivity cone; unrelated arbitrary triples do not qualify |
| Bipartite sum labels -> oriented differences | Change sign on one bipartition | Distinct oriented values do not ensure distinct absolute values; no graceful-tree theorem |
| Integer square edges -> fourth powers over Gaussian/Eisenstein integers | No map supplied | Different rings and predicates; shared word “square” supplies no implication |

The current Eisenstein audit's represented 6 and54 and local criterion n not congruent5 modulo8 do not alter these integer graph tests; finite failed fourth-power searches are not nonrepresentation proofs.

## 6. Reproduction and outcome

    python3 check_square_translations.py > certificates.json
    python3 -O check_square_translations.py > certificates.optimized.json
    cmp certificates.json certificates.optimized.json

Python standard library only. Explicit runtime checks remain active under `-O`. The factor classifier is independently checked by direct root enumeration for every n=2,...,100, using the proved finite root bound. The infinite family is additionally tested for n=8,16,...,8000. All 24 n=1000 translations, roots, factorizations and aliquot values are recorded. These checks supplement, not replace, the proofs above.

Outcome: the user's triple is exactly located in the prior graph problem; complete common-translation classification provides a least positive realization whose endpoints share aliquot fiber1. No blocker. Infinite simultaneous prime-endpoint occurrence, additional nontrivial equal aliquot fibers, and prescribed-interval embedding are not resolved here. Source publisher was not modified.

## Independent audit

A separate root-enumeration method searched the proved complete root range at n=1000, reproducing all24 translations without using the factor-pair enumeration. Exact trial division checked every endpoint and confirmed the unique prime-endpoint offset177612. The public checker includes an explicit uniqueness guard. Both normal and optimized runs agree.
