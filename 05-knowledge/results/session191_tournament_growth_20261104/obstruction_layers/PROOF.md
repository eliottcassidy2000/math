# The top layer of bounded-H obstructions and finite near-top templates

**PROVED from the cited strong-floor theorem and elementary SCC arguments;
independently reviewed.** This is a refinement of the inherited sharp
obstruction cutoff, not a new exclusion of H=7 or H=21 and not a proof of
spectrum completeness. No literature-priority claim is made.

## 1. Scope, inheritance, and notation

H(T) counts directed Hamiltonian vertex words. Put H(empty)=1. For an
integer B>=1, an obstruction means a tournament F with H(F)>B but
H(F-v)<=B for every vertex v. Hamiltonian-path insertion makes H monotone
under induced deletion, so this is precisely induced-minimality outside
the class H<=B. A tournament's strongly connected components (SCCs) have
a unique linear dominance order, and H multiplies along that order.

Let r be the positive integer with 3^(r-1)<=B<3^r. Session133 proves that
every obstruction has order at most 3r, attained by r ordered directed
triangles. Its exact-core sharpening says every other obstruction has
at most K(B)+1 vertices. That result alone does not classify all equality
cases at order 3r when K(B)+1=3r.

Inherited sources:

* [Session133 sharp obstruction cutoff](https://github.com/eliottcassidy2000/math/blob/9aa6e0836f5c92b6448af240154d639e711ec420/05-knowledge/results/session133_fractran_spectrum_refinement_20261101/spectrum/OBSTRUCTION_PROOF.md).
* [Independent 7/21 exclusion and elementary strong-deletion proof](https://github.com/eliottcassidy2000/math/blob/9aa6e0836f5c92b6448af240154d639e711ec420/05-knowledge/results/exclusion_7_21_20261103/PROOF.md).
* THM-1370 / `h-spectrum-omits-7-21-all-n`: all-order exclusions, including
  the prior-art correction. Completeness of the spectrum remains open.
* THM-4744 / `interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve`:
  the cited Busch bound H(C)^3>=5^(|C|-1) for nontrivial strong C, and the
  exact strong floor at five vertices, H>=9.

Write C3 for the directed triangle, and Q for the unique strong tournament
on four vertices. Then H(C3)=3, H(Q)=5, and max_v H(Q-v)=3. Q is rigid
and self-converse; C3 has three automorphisms and is self-converse.
These small facts are also independently enumerated by `verify.py`.

Concept board: SCC products; induced-criticality; normalized strong cost;
order deficit; finite exceptional components. The hostile to avoid is
inferring arbitrary-target existence from a bounded-threshold obstruction
theorem. The new coordinate retained here is distance below the *sharp
obstruction-order ceiling*, not distance below the strong H floor.

## 2. A multiplicative cost for every nontriangle SCC

Set kappa=125/81>1. For every strong C with s>=3 vertices,

    H(C)^3 / 3^s >= kappa^(s-3).                       (1)

For s=3 this is equality. For s>=4, Busch gives H(C)^3>=5^(s-1), and

    [5^(s-1)/3^s] / kappa^(s-3) = (27/25)^(s-4) >= 1.

Equality in (1) holds for C3 and Q; when s>=5 the displayed lower bound
is strict. For a tournament U with no singleton SCCs define

    E(U) = sum_C (|C|-3),

the sum being over all SCCs. Thus E counts excess vertices above one
triangle per SCC, not the number of nontriangle components. Multiplying
(1) gives

    H(U)^3 >= 3^|U| kappa^E(U).                        (2)

One useful stronger fact at low defect is:

    if U has a component of order >=5, then
    H(U)^3 >= 3^(|U|+1).                              (3)

For s=5 use H>=9, so H^3=at least 3^(s+1). For s>=6, Busch gives
5^(s-1)>3^(s+1): start with 3125>2187 and extend by factors 5 versus3.
Every other component contributes at least 3^s, proving (3).

## 3. Near-top stability, uniformly in the threshold

**Theorem A.** Suppose F is an obstruction of order n=3r-d, d>=0.
If it is not an ordered sum of triangles, then

    kappa^(E(F)-1) < 3^(d+1),                         (4)
    E(F) < 1 + (d+1) log(3)/log(125/81).

The logarithmic coefficient is approximately 2.5321548074. An exact
integer test for (4) is

    125^(E-1) < 81^(E-1) * 3^(d+1),

so no floating-point comparison is needed. The combined order of all
nontriangle SCCs is at most 4E(F), and every such SCC has order at most
E(F)+3. Thus a fixed d admits only finitely many exceptional SCC lists,
independently of B and r.

**Proof.** F has no singleton SCC, since deleting it would preserve H.
Choose a component C of order at least four. The strong-deletion lemma
provides v in C with C-v strong. Consequently F-v still has no singleton
SCC, has n-1 vertices, and E(F-v)=E(F)-1. By minimality and (2),

    3^(3r-d-1) kappa^(E(F)-1)
      <= H(F-v)^3 <= B^3 < 3^(3r).

This proves (4). If s>=4 then s<=4(s-3); summing proves the exceptional
order bound. Finally, s-3<=E proves the individual component bound. QED.

The pure-triangle case is completely explicit: its count and maximum
deletion count are 3^k and 3^(k-1), so minimality forces k=r and d=0.
It is not an omitted family at positive defect.

## 4. Exact classification at the sharp maximum order

**Theorem B.** Every obstruction of order exactly 3r is one of:

1. The ordered sum of r triangles, always present.
2. When r>=4 and B>=25*3^(r-3), an ordered SCC sequence containing exactly
   three copies of Q and r-4 copies of C3, in any order.

There are no others. In the second family,

    H(F) = 125*3^(r-4),
    max_v H(F-v) = 75*3^(r-4) = 25*3^(r-3).           (5)

Thus the number of ordinary isomorphism classes at order 3r is

    1 + 1_{r>=4 and B>=25*3^(r-3)} * binom(r-1,3).    (6)

**Proof.** In the nontriangle case, (4) at d=0 forces E(F)<=3:
kappa^3>3, because 125^3=1953125>3*81^3=1594323. Since
3r=3*(number of components)+E(F), the nonzero E is divisible by three,
so E=3. The possible component sizes other than triangles are therefore:

    (6), (5,4), or (4,4,4).

The (6) option fails: strong deletion in that component leaves a
five-vertex component in the singleton-free F-v of order 3r-1, so (3)
gives H(F-v)>=3^r>B. The (5,4) option fails by deleting within the
four-vertex component while preserving its strong triangle; the five-
vertex component remains, giving the same contradiction. Therefore only
(4,4,4) survives. The order constraint gives exactly r-4 triangles.

For this family, deleting from Q has maximum count 3*5^2*3^(r-4), and
deleting from a triangle, if one exists, has count 5^3*3^(r-5), which is
smaller. Its H exceeds 3^r because125>81. Thus minimality is equivalent
to the lower threshold in the statement. Distinct SCC words are distinct
isomorphism classes, because the SCC dominance order is intrinsic. There
are r-1 components in total, with three Q positions. This proves (6). QED.

First transitions:

| r | top order | B interval | new Q-family starts at | additional classes |
|---:|---:|---:|---:|---:|
| 1 | 3 | 1..2 | absent | 0 |
| 2 | 6 | 3..8 | absent | 0 |
| 3 | 9 | 9..26 | absent | 0 |
| 4 | 12 | 27..80 | 75 | 1 |
| 5 | 15 | 81..242 | 225 | 4 |
| 6 | 18 | 243..728 | 675 | 10 |
| 7 | 21 | 729..2186 | 2025 | 20 |

In particular, Session133's unique order-nine obstruction at B=21 is the
first low-threshold instance of a complete top-layer theorem. The pure
triangle obstruction is not universally unique: Q->Q->Q is an order-12
obstruction for every B from75 through80.

The labelled count follows without new counting assumptions. It is
(3r)!/3^r for the triangle class, plus, when the second family is present,
binom(r-1,3)*(3r)!/3^(r-4). Automorphisms preserve every SCC position.

## 5. Exact finite-template law at any fixed order deficit

This is an explicit consequence of Theorem A, not merely a finiteness
analogy. Fix d. An exceptional skeleton S=(C_1,...,C_k) is an ordered list
of nontriangle strong classes allowed by (4). Write

    s = sum_i |C_i|,  h = product_i H(C_i),
    q_i = max_v H(C_i-v)/H(C_i),
    t=(s+d)/3.

Discard S unless t is an integer. At order 3r-d exactly a=r-t triangles
must be inserted into the k+1 gaps of this list. If a>=1, put

    U(S,d)=h/3^t,
    L(S,d)=U(S,d)*max(1/3,q_1,...,q_k).                (7)

For theta=B/3^r, the resulting tournaments are minimal obstructions iff

    L(S,d)<=theta<U(S,d).                            (8)

When (8) holds the skeleton contributes exactly binom(a+k,k) ordinary
classes, since each weak composition of a into the k+1 gaps is distinct.
This also remains correct with repeated exceptional types: deleting all
triangle SCCs recovers the ordered skeleton and its gap lengths uniquely.
When a=0, omit 1/3 from the maximum in (7). This boundary must not be
silently included in the stabilized formula.

**Proof.** SCC multiplication gives H=3^a h=3^r U. Deletion from a triangle
multiplies H by1/3. Deletion in C_i multiplies it by at most q_i, and the
maximum is attained. Thus H>B and all deletion counts<=B are exactly (8).
The finite-template assertion follows from the bounds in Theorem A. QED.

Consequently, at each fixed d and on each open chamber between the finitely
many rational endpoints L,U in [1/3,1), the number of obstructions of order
3r-d is eventually a polynomial in r. The same holds at each endpoint,
using the indicated closed-left/open-right convention. This is a chamber
statement in the normalized threshold theta, not a single polynomial for
an arbitrarily oscillating sequence B_r. It does not require knowing the
full spectrum. Computing the finite strong-class templates at larger d
is a separate finite task, not claimed completed here.

## 6. Verification and boundaries

`verify.py` builds all tournament classes through order six
using the complete fixed-Hamilton-path cover, verifies class counts
1,1,2,4,12,56, checks every H by literal permutations, verifies all vertex
deletion maxima, and exhausts the three possible nontriangle size patterns
in the top layer. It checks finite threshold transitions and independently
constructs actual Q/C3 sums for direct subset-DP checks through order15.

The all-order proof is Sections2--5. Finite checks are controls, not a
replacement for the quantified argument. The exceptional-floor input is
Busch's cited theorem. The ordinary 7/21 nonexistence and fixed-H rational
series are inherited; no arbitrary target constructor or full spectrum
proof follows from the present theorem.
