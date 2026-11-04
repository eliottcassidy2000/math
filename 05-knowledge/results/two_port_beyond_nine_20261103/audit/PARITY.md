# Squarefree parity identities and the even-converse stopping boundary

Status: **PROVED** universal identities P1-P3 and theorem P4; **FINITE-EXACT** endpoint
census; **REFUTED** parity shortcuts; **OPEN** exact even endpoint-balance
exclusion and universal two-port reconstruction. Independent audit lane,
2026-11-03. No novelty priority is claimed.

## Inheritance and the question tested

The closest proved mechanism is Session180's
[named squarefree path matrix](../../session180_170807_returns_20261104/ports/PROOF.md),
whose entries retain every path support and both endpoints. The earlier
[odd-converse obstruction A1](../../two_port_unit_distance_20261103/algebra/PROOF.md)
uses ordered two-path covers modulo four. The hostile is a balanced odd
tournament and its converse: start/end margins can coincide although ordered
two-path information separates them. The corrected near miss here is the
stronger claim that parity alone excludes endpoint balance at even order.
The least-used sidecar is the endpoint count on every induced subset.

The live concepts are named endpoints, squarefree support, converse symmetry,
ordered two-path covers, and the distinction between parity and exact counts.
The map in this report sends the full squarefree path matrix over the integers
to its reduction modulo two. It preserves exact coefficient parity and support
incidence, but destroys the signed orientation of C. It therefore cannot be
substituted for the integral joint endpoint interface in a reconstruction proof.

## P1. A rank-one skew relation in characteristic two

Let T be any tournament on the named finite set V, A its 0/1 adjacency
matrix, and work in

    R = F_2[z_v : v in V] / (z_v^2 : v in V).

Put Z=diag(z_v), G=(I-ZA)^(-1)Z, s=G1 and e=G^t1. Each inverse is a finite
geometric series. The entry G_uv is the sum of all simple directed paths
from u to v, weighted by their vertex support, including singleton paths.
Then

    G + G^t = s s^t.                                  (P1)

Proof: for any two adjacency matrices A,B the resolvent identity gives

    G_A-G_B = G_A(A-B)G_B.

No inverse of Z is used. Taking B=A^t gives G_B=G_A^t. In characteristic
two, tournament adjacency satisfies A+A^t=J+I, including the zero diagonal.
Also GG^t=0: every product G_uk G_vk contains z_k twice. Consequently

    G+G^t = G(J+I)G^t = (G1)(G1)^t.

This proof covers all finite orders, including zero and one.

## P2. Start and end polynomials are subset-zeta transforms

Set U=product_(v in V)(1+z_v). The classical oddness of the number of
Hamilton paths in every nonempty induced tournament gives

    P:=1^t G1 = U-1.

Summing P1 over its second matrix index yields s+e=sP. Hence

    e=Us,             s=Ue.                           (P2)

The second identity follows also from U^2=1. In particular, for every u,

    E_u(T) = sum_(W subset V, u in W) S_u(T[W]) mod 2,

and the same formula with S,E interchanged holds. Here S_u and E_u count
Hamilton paths starting and ending at u; the formula includes W={u}.
This is a relation between the full family of induced endpoint counts.
It does not say that S(T)=E(T), even modulo two, for the full vertex set.

## P3. The parity of ordered two-path covers forgets direction

Let C_uv count ordered spanning covers by two nonempty paths, the first
ending at u and the second starting at v. The squarefree product gives

    C_uv mod 2 = [V] e_u s_v
               = [V] U s_u s_v
               = [V] U(G_uv+G_vu).                   (P3)

Multiplication by U before extracting [V] sums all coefficients of its
other factor. Therefore C_uv modulo two is the parity of the number of
simple paths from u to v or from v to u, over all possible vertex subsets.
In particular **C is symmetric modulo two for every tournament**.
This is an unordered all-length path observable; it is not an adjacency
reconstruction formula. It explains why the earlier converse obstruction
needs information beyond C modulo two.

## P4. Even self-converse tournaments have unequal endpoint parity

**Theorem P4.** If T has positive even order and is isomorphic to its
converse, then S(T) and E(T) differ modulo two in at least one named
coordinate. Thus the endpoint part of the observed sidecar separates
T from its labelled converse.

Proof: any tournament automorphism has odd order. Otherwise a suitable
power is a nonidentity involution, which exchanges some pair u,v while
preserving the arc between them, an impossibility. Let tau be an
anti-automorphism, an isomorphism from T to its converse. Its square is
an automorphism of odd order. Since |V|>=2, tau cannot itself have odd
order: an odd power of an anti-automorphism is still an anti-automorphism,
whereas the identity cannot reverse an existing arc. Therefore tau has
order2m with m odd, and sigma=tau^m is an involutory anti-automorphism.

There cannot be two fixed vertices of sigma, because their arc would
have to reverse while its endpoints stay fixed. At even order the number
of fixed vertices of an involution is even, so sigma has none. Reversing
the image of each Hamilton word gives S_u=E_(sigma(u)). If S=E modulo two,
the fixed-point-free pairs of sigma have equal S parity. Their sum H(T)
would be even, contradicting Hamilton-path oddness. QED.

This argument does not assume regularity or vertex transitivity. It
does require an isomorphism to the converse. A pair of different converse
isomorphism classes is outside its scope; the following eight-vertex
hostiles show why that distinction matters.

## Exact hostile controls and the remaining even-order question

The proposed shortcut

    "For even |V|, some vertex has S_u != E_u modulo two"

is **REFUTED**. Its least order in the complete inherited bank through nine
is eight. One witness has outgoing bit masks

    (4,17,42,131,45,75,31,119),

where bit v in row u means u->v. Its upper-triangle code is39887362,
with vertices0,...,7 and bits ordered (0,1),(0,2),...,(6,7). Direct counts are

    H=249,
    S=(13,22,13,29,26,29,52,65),
    E=(87,42,31,13,26,19,18,13).

Every coordinate has the same parity, but the vectors differ over the
integers. The other eight-vertex class with this property is its converse
class, code72392840. Literal all-path enumeration verifies both witnesses
independently of the bank's dynamic programming. This invalidates the first
implication of a parity proof of even endpoint imbalance. The strongest
survivor is P2; it constrains all subset endpoint counts together.

The auxiliary claim that S dot E is even at every even order is also false.
At order six, code689 has outgoing rows(34,20,9,3,13,30), H=37,
S=(5,6,5,4,7,10), E=(4,6,9,8,6,4), and odd dot product. Thus the odd total
H cannot be converted into an even-order endpoint contradiction by this
inner-product shortcut.

The exact statement "S(T)!=E(T) whenever |V| is positive and even" remains
**OPEN** in this packet. It would separate every even-order tournament
from its labelled converse already at the endpoint level. Combined with
the older odd-C theorem it would exclude every nontrivial whole-converse
collision, but would still not prove universal two-port reconstruction.

The independent complete inherited-class endpoint scout records:

| Order | Classes | Exact S=E classes | S=E modulo two classes |
|---:|---:|---:|---:|
|1|1|1|1|
|2|1|0|0|
|3|2|1|1|
|4|4|0|0|
|5|12|2|3|
|6|56|0|0|
|7|456|3|48|
|8|6880|0|2|
|9|191536|5|4628|

Completeness is inherited from the pinned representative bank, not reproved
by this scout. The separately audited complete order-ten census also finds
zero exact S=E classes; it does not count parity-balanced classes.
Exact balance need not imply regularity: the five-vertex code152 has rows
(16,5,9,3,14), S=E=(3,3,3,3,3), and degrees(1,2,2,2,3).

## Reproduction and scope

From the repository root:

```text
python3 05-knowledge/results/two_port_beyond_nine_20261103/audit/endpoint_scout.py --output 05-knowledge/results/two_port_beyond_nine_20261103/audit/endpoints.json
python3 05-knowledge/results/two_port_beyond_nine_20261103/audit/parity_check.py --output 05-knowledge/results/two_port_beyond_nine_20261103/audit/parity.json
python3 -O 05-knowledge/results/two_port_beyond_nine_20261103/audit/parity_check.py --output 05-knowledge/results/two_port_beyond_nine_20261103/audit/parity.optimized.json
python3 05-knowledge/results/two_port_beyond_nine_20261103/audit/self_converse_check.py --output 05-knowledge/results/two_port_beyond_nine_20261103/audit/self_converse.json
python3 -O 05-knowledge/results/two_port_beyond_nine_20261103/audit/self_converse_check.py --output 05-knowledge/results/two_port_beyond_nine_20261103/audit/self_converse.optimized.json
```

The independent literal checker makes99,344 checks: every coefficient of
P1-P3 for every labelled tournament at orders1-5, every inherited class at
order6, and both eight-vertex hostile witnesses, plus the dot-product hostile.
Its normal and optimized receipts agree byte-for-byte. Checks use explicit
exceptions, not assertions. The endpoint scout is separate C++ subset DP
over all inherited classes through nine; sources and the bank are hash-pinned
in its receipt. These controls corroborate the proofs and refutations; they
are not extrapolated to an all-order endpoint-balance theorem.

The separate P4 checker enumerates every tournament compatible with the
fixed adjacent-pair anti-involution at orders2,4,6 (2,16,512 tournaments),
and sixteen specified controls at order8. It counts Hamilton words by
literal permutations and verifies the endpoint transfer and parity
separation. Every self-converse even tournament can be relabelled to
one of these fixed-pair involution families, although this control does
not claim one representative per isomorphism class. Both order8 parity
hostiles are additionally tested against all8! permutations and have no
anti-automorphism, as required by P4.
