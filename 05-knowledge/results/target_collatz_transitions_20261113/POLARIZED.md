# Independent modules give exact target factorization and ternary compensation

**PROVED scoped. FINITE-EXACT controls are identified below. Chan–Pak is
OPEN.** The prior symmetric norm chart excludes every7p with prime p≠7.
The construction here releases its equality constraint on the two module
responses. Its three-arc change has a measured count drift, and that drift
can be used rather than required to vanish.

## 1. The actual three-arc transport

Let M,N be arbitrary nonempty tournaments, with

    h=H(M), p=P_M(2), j=H(N), q=P_N(2).

P(2) counts ordered two-slot spanning path covers, allowing empty slots;
H(empty)=1. Thus h,j are positive odd integers and p,q are positive even
integers. Form the ordinal chain

    r -> s -> M -> t -> N -> u -> v,

with five singleton letters and a new apex a. Before the change, a beats
precisely r,t,N; afterward it beats precisely r,s,N,u. Exactly a-s,a-t,a-u
reverse. Every old internal module arc, the root r, its entire star and
its unique predecessor a are retained. The graphs have |M|+|N|+6 vertices.
They are strong: a path through the ordinal chain closes through a.

Write the root state as(H,s_r,q_r,F_r), with the same endpoint/deletion
and two-slot conventions as the inherited
[mirror proof](../sixty_three_bridge_20261113/MIRROR_MODULE.md). Put

    A=12pq+6hq+2pj−hj,
    Z=12pq+2hq+6pj−hj,
    D=pj−hq.

Then the exact states are

    before: (A, 6hq−hj, 6hq−hj, A),
    after:  (Z, (2p+h)(2q+j)−hj, (2p+h)(2q+j), Z+hj),
    Z−A=4D.                                           (1)

In particular the move is neutral exactly when p/h=q/j. Equality of
these ratios is sufficient; no isomorphism or converse relation between
M,N is required. For TT2 and C3, the responses(1,4),(3,12) give the
same count669 before and after, with defect0→3. Unit readiness is not
preserved by replacing a singleton defect with a product hj.

**Proof.** Split each HP at its unique apex. Every module has at most one
run on each side, so its nonempty-run weights are h,p−2h or j,q−2j.
Legal assignments are exactly those whose last prefix block loses to a
and whose first suffix block is beaten by a. The finite quotient-run
table is the same one already proved for the symmetric mirror, but the
two middle columns must now remain distinct:

| observer | hj | h(q−2j) | (p−2h)j | (p−2h)(q−2j) |
|---|---:|---:|---:|---:|
| before H |63|30|26|12|
| before s,q |11|6|0|0|
| after H |63|26|30|12|
| after s |24|10|10|4|
| after q |25|10|10|4|

Expanding proves(1). The root has indegree one, so F=H+q_r−s_r by the
inherited concatenation bijection. The quotient table has a complete
finite word enumeration, not a fit to samples. The additional data beyond
the symmetric formula are precisely the two ordered cross terms. ∎

## 2. The output factorization and its primitive guard

Apply the actual N then R births to the prepared root, using the
[inherited root-child maps](../two_port_return_induction_20261103/ROOTED_RESETS.md).
They add two vertices and give the balanced rooted state

    (Y,Z,Z,Y),
    Y=3Z+2hj=36pq+6hq+18pj−hj
      =(6p+h)(6q+3j)−4hj.                             (2)

The exact primitive criterion is

    gcd(h,3p(2q+j))=1,
    gcd(j,q(6p+h))=1.                                 (3)

Indeed gcd(Y,Z)=gcd(hj,Z), since Z is odd. Reduction modulo h gives
Z≡6p(2q+j), and modulo j gives Z≡2q(6p+h). Removing the factors2 from
odd h,j gives(3), in both directions. No coprimality between h and j is
assumed in this criterion.

This genuinely escapes the norm restriction, with small physical inputs:

| M | N | Z | Y | marked content |
|---|---|---:|---:|---:|
| TT1 | TT1 |63|191|1|
| TT1 | TT2 |115|347|1|
| TT2 | TT1 |123|371=7·53|1|
| Q4 | TT2 |1763|5299=7·757|1|

Here Q4=(2,4,9,3) has(h,p)=(5,32). The prime347 is impossible for the
old norm form because it is2 modulo5;371 and5299 have odd7-valuation.
These are new construction interfaces, not new global spectrum records.
The count371 was already reached by the previous companion-word repair.

The ordering of the two modules matters. If L(M)=36p²+24hp−h² denotes
the old symmetric norm output, direct expansion gives

    Y(M,N)Y(N,M)=L(M)L(N)+144(pj−hq)².                  (4)

Thus the determinant is both the exact surgery drift and the extra term
that the diagonal norm discards. A scalar norm or an unordered pair of
responses loses the coordinate needed for the asymmetric transition.

## 3. A nontransitive primitive family divisible by seven

Fix M=Q4 and let N=TT_b, b≥1. Equation(2) gives

    Y=1182·2^b+571.

Criterion(3) is equivalent to b≠1 modulo4. Divisibility by7 is equivalent
to b≡2 modulo3. Therefore b=12t+2, t≥0, supplies a uniform family

    Y_t=4728·4096^t+571=7c_t,
    c_0=757,

with primitive balanced ports and an actual retained nontransitive
module. Since v7(4096−1)=1 and7∤4728,

    v7(c_t−c_s)=v7(t−s)   for t≠s.                    (5)

The elementary binomial valuation proof is the same as the inherited
seven-adic clocks. Thus every cofactor residue modulo7^e appears exactly
once as t runs over any7^e consecutive indices. This gives actual source
families, but not all exact seven-prime targets or infinitely many prime
c_t. A prescribed target in this family is selected by testing whether
(Y−571)/4728 is a nonnegative integral power of4096, then constructing
Q4 and TT_(12t+2).

## 4. A complete finite inverse for each target in this construction

The inherited cover inequalities are

    p≥(|M|+1)h, p≥2^|M|,
    q≥(|N|+1)j, q≥2^|N|.

The first cuts each Hamilton path; the second counts every ordered
vertex partition, whose two induced tournaments each have an HP. The
second bound is explicitly used in the older
[complete inverse audit §1](../target_directed_cores_20261112/ports/COMPLETE_INVERSE_AUDIT.md).
Neither inequality is new here.

Substitution in(2) gives Y≥191hj. Hence the following is an exact inverse
description for this polarized constructor:

1. Enumerate positive odd h,j with hj≤Y/191.
2. Factor Y+4hj=xy, retaining ordered factor pairs with
   x≥13h, y≥15j, x≡h(mod12), y≡3j(mod12).
3. Put p=(x−h)/6 and q=(y−3j)/6.
4. Find actual modules with responses(h,p),(j,q). It suffices to examine
   orders at most min(floor(p/h)−1,floor(log2p)) and
   min(floor(q/j)−1,floor(log2q)), respectively.
5. Construct the prepared seven-block graph and the two births. For a
   primitive output impose(3) as well.

Every accepted pair produces the target literally. Conversely the actual
responses of every realization in this constructor occur in these finite
loops. Arithmetic candidates without modules do not count as acceptance.

There is also a simple sharp bound independent of h,j:

    Y≥90p+11,   Y≥78q+35.                              (6)

For example, replace q by2j in(2) to obtain the first bound, then use
h,j≥1; the other follows by replacing p by2h. Therefore every module
needed for a fixed target has order O(log Y). In particular, for Y<39971,
both module orders are at most8. The threshold is sharp: M=TT1,N=TT9
has q=512 and Y=39971.

The inherited complete ordinary bank through8 therefore decides this
constructor at ALL orders for Y<39971. It gives205 attained target totals
in that interval. This is a finite head of this construction, not of the
full tournament spectrum. The simpler arithmetic inverse already has
one candidate at191,347,371,5299 and none at881. For881, hj≤4 and its
residue modulo3 forces h=j=1; none of the ordered factors of885 meets
the required size and mod12 guards. Thus even the polarized construction
needs other source charts to handle that known attained target.

## 5. Exact negative ternary transport by determinant compensation

In terms of the BEFORE count A, equation(2) reads

    Y=3A+2e,  e=6(pj−hq)+hj
              =j(6p+h)−6hq.                          (7)

Consequently the exact guard e=±1 supplies the desired3A±2 transition.
This is a determinant-one condition on the two columns
(6p+h,h) and(6q,j), a precise Farey/Stern interface. It is not a claim
that arbitrary addition of such columns produces tournament responses.
Also e is always odd, so this particular mechanism cannot give3A or3A+4.

The inherited bank through8 has1922 distinct response pairs. Of their
ordered pairings, exactly nine satisfy e=±1: eight transitive diagonal
pairs and one nontrivial pair,

    M=(4,9,82,21,35,143,187,31),    (h,p)=(181,3424),
    N=(6,12,48,21,131,219,159,15),  (j,q)=(179,3416).

For these actual modules,

    D=−5400, hj=32399, e=−1,
    A=145259777,
    Z=145238177,
    Y=435779329=3A−2,
    gcd(Y,Z)=1.                                       (8)

This is a22-vertex source, three old arc reversals, and two literal
births producing a24-vertex primitive balanced output. The large path
totals follow from the universal module theorem, not an infeasible full
24-vertex census. The two eight-vertex module H values are also checked
independently by direct permutation.

The exact fixed-M and fixed-core-padding partner selectors are proved in
[FAREY_PARTNERS](history/FAREY_PARTNERS.md). For this pair, arbitrary
independent transitive padding creates no additional±1 compensation:
the displayed source is isolated in that padding language. This is a
precise stopping boundary. Renewing the determinant condition requires
changing the module responses, not just repeating this witness.

## 6. Reproduction and scope

Run `response_probe.py` to recover the pinned ordinary bank and solve the
ordered response-pair condition. Run `polarized.py` for the full symbolic
table,16 literal module pairs through output order16, the independent
eight-vertex permutation counts, and named arithmetic inverses. Run
`selector_head.py` for the complete all-order construction head below39971.

The proofs establish the universal transition and exact finite inverse.
They do not show that every remaining prime or7p enters this chart, nor
that every successful output is again an input for a desired next digit.
The point of releasing the diagonal is an actual escape from an all-order
arithmetic barrier, with target selection and its physical guards intact.
