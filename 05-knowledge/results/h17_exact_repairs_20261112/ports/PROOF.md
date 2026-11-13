# Marked slot factors and a compensated neutral block swap

**PROVED:** marked ordinal factorization, the cyclic wrapper's full rooted
formula, an exact criterion for a two-block compensation, and the literal
positive example below. **FINITE-EXACT:** the declared small block search.
The positive example changes a continuation defect at the same current
and first-child counts; its gcd remains9. It is not a primitive repair or
an arbitrary-target selector. Chan–Pak remains **OPEN**.

## 1. What transfers from the older boundary work

The relevant old mechanism is a complete, retained boundary interface:
[Session61 palette-deck, Sections1,3,5](../../session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md)
distinguishes a full section from its parity shadow and full joint
colouring relations from their pair projections. The later
[HN six-colour interface, Sections1–3](../../hn_metagraph_information_20261026/six_color/RESULT.md)
proves that an unprescribed frame can be invisible under global palette
symmetry, while a retained palette or second frame changes the actual
boundary relation. Its gluing theorem also identifies when transparent
fragments remain transparent under arbitrarily many gluings.

The transfer here is a counting-interface construction, not an assertion
that tournament arcs are unit-distance edges. The unmarked slot polynomial
is the quotient, a physical root inside one factor is the retained mark,
and the marked slot polynomials below are the lost coordinates. Ordinal
gluing multiplies these coordinates with different arguments on the two
sides of the mark. This gives a lawful compensated swap, whereas an
unmarked equality by itself would not preserve the next target.

Inheritance board: full boundary tables / a named internal root / ordinal
factor order / cyclic two-vertex wrapper / first-child target / defect
and gcd. The closest positive example is
[the neutral H17 preparation](../../eight_nine_twentyseven_20261112/constructors/SIX_RETURN.md).
The hostile is its one-direction root shift: current H stays fixed while
the next birth target changes. The least-used sidecar is the *three-slot*
marked cover evaluation, which controls the cyclic wrapper's rooted F.

## 2. Marked slot polynomials factor with a shifted argument

For a tournament T put

    P_T(t)=sum_k C_k(T) binom(t,k),        P_empty(t)=1.

For a nonempty T with physical root v define E_(T,v)(t) to count spanning
covers in t labelled path slots, where slot1 is nonempty and ends at v;
the other slots may be empty. Define S_(T,v)(t) with slot1 starting at v.
These are integer polynomials: for example

    E_(T,v)(t)=sum_(U containing v) e_v(T[U]) P_(T−U)(t−1).       (1)

Here e_v counts Hamilton paths ending at v. The definitions are for
positive integers t; the polynomial identity extends to every argument.

**Marked factorization theorem.** If T=L→K→R and v lies in K, then

    E_T(t)=P_L(t) E_K(t) P_R(t−1),
    S_T(t)=P_L(t−1) S_K(t) P_R(t),
    P_(T−v)(t)=P_L(t) P_(K−v)(t) P_R(t).                       (2)

*Proof.* Every path in an ordinal join meets its factors in their order.
Choose the restrictions to each labelled slot independently and concatenate
the nonempty pieces. For E, the distinguished slot cannot contain R,
because it ends in K. Its other t−1 slots cover R. Every slot may meet L.
The S identity is the reverse argument. Deletion preserves the ordinal
product, including an empty K−v. QED.

Thus any permutation of a fixed collection of ordinal factors preserves
the full unmarked polynomial and H. If the physical root stays in K, its
deletion polynomial also stays fixed; its E and S generally change.
Permutations entirely to one side of K remain invisible even to this
marked interface. Moving factors across the mark changes explicit
evaluation factors in (2).

There is a useful conservation law:

    sum_v E_(T,v)(t)=sum_v S_(T,v)(t)=P_T(t)−P_T(t−1).           (3)

Every cover with nonempty distinguished slot has exactly one end and one
start for that slot. The right side removes covers with an empty first
slot. Neutrality of P therefore does not make the marked values freely
assignable; their sum is retained.

## 3. Exact rooted state of a cyclic two-vertex wrapper

Adjoin X,Y with all arcs `T→X→Y→T`, and retain v inside T. Write G(T)
for this actual tournament. Its rooted state convention is `(h,s,q,F)`:
H, Hamilton paths starting at v, H after deleting v, and ordered spanning
two-slot covers whose first slot is nonempty and ends at v.

**Wrapper theorem.** For every nonempty T,

    h=P_T(1)+P_T(2),
    s=S_T(2),
    q=P_(T−v)(1)+P_(T−v)(2)−[|T|=1],
    F=2E_T(2)+2E_T(3)−E_T(1).                              (4)

The indicator is necessary: deleting the sole module vertex leaves the
single directed edge X→Y, whose H is1, not2.

*Proof.* A Hamilton path either has Y→T→X, giving P_T(1), or has the
consecutive pair X→Y, with the T vertices split between the two outside
slots, giving P_T(2). The latter includes either empty end slot. A path
starting at v has exactly that two-slot form with its first slot starting
at v, proving the s formula. The same argument after deletion gives q,
with the explicitly stated empty-module exception.

For F, classify the path ending at v by its use of X,Y.

- If it uses neither, the other path covers X,Y and the remaining T.
  This gives E_T(2)+E_T(3)−E_T(1): the subtraction corrects the empty
  remaining module, where X→Y has one path instead of two.
- If it uses Y but not X, it starts with Y; the other path ends with X.
  This gives E_T(2).
- If it uses both, X must immediately precede Y. The module portions
  before X, after Y, and on the other path give E_T(3).
- It cannot use X without Y and still end in T.

Adding the disjoint cases proves (4). QED.

For T=TT_n, n>=2, let v occupy zero-based rank j. Then

    h=1+2^n,              q=1+2^(n−1),
    s=2^(n−1−j),
    F=2^(j+1)+2*3^j*2^(n−1−j)−[j=n−1].                    (5)

At n=4 these give `(s,F)=(8,18),(4,28),(2,44),(1,69)`.
The H17 preparation is the move from rank0 to rank2. Formula(4)
extends that example to arbitrary marked ordinal factors.

## 4. A scalar criterion that really preserves the first birth target

Fix a rooted K and nonempty blocks A,B. Compare the two literal modules
`A→K→B` and `B→K→A`, inside the same cyclic wrapper. All vertex names
and internal block arcs are retained; only the directions between blocks
change. Put

    h_i=P_i(1), k_i=P_i(2), l_i=P_i(3),
    e=E_K(2), u=S_K(2), w=E_K(3),
    alpha_K=(2e+u)/(2w).

The wrappers have the same h and q, and the same full unmarked slot
polynomial. The latter follows either by compressing module runs in each
ordered path-cover slot, or by the inherited uniform module response
theorem: P_A P_K P_B is unchanged.

For a root-child birth, both R and N have new Hamilton count `h+F−s`.
The two source states have the same first-child target **if and only if**

    (2e+u)(k_A h_B−k_B h_A)
       +2w(l_A k_B−l_B k_A)=0.                            (6)

This follows directly from (2)–(4); the E_T(1) term vanishes because
the right factor is nonempty. If `h_A k_B!=h_B k_A`, it is equivalent to

    alpha_K=(l_A k_B−l_B k_A)/(h_A k_B−h_B k_A).            (7)

The first numerator/denominator is an actual secant of the two evaluator
ratios, not an independent numerical input. Equality (6) is a complete
test for this specified physical swap. Under it the state difference is
`(0,v,0,v)`, where `v=u(h_A k_B−h_B k_A)`.

Every actual rooted K satisfies `0<alpha_K<=3/2`, with strict inequality
for |K|>1. Indeed E_K(2)<=E_K(3) by appending an empty slot. Also
S_K(2)<=E_K(3): remove the starting root into the distinguished singleton
slot and put the remaining first path and the old second path into the
two free slots. If |K|>1, the first injection is not surjective, proving
strictness. This gives a cheap necessary filter for (7).

Transitive prefix/suffix preparation has an exact effect on this filter.
For K'=TT_a→K→TT_b,

    alpha_(K')=(2^(a+1) E_K(2)+2^b S_K(2))
                /(2*3^a*2^b E_K(3)).                       (8)

This is a boundary-state transformation, not a claim that arbitrary
rational slopes are realizable.

## 5. A positive compensated swap at two exact targets

Take the strong seven-vertex blocks

    A=(32,17,67,7,13,30,59),
    B=(26,44,121,112,98,65,3).

Their exact triples `(P(1),P(2),P(3))` are

    A: (99,1236,8289),       B: (105,1404,9411).

Take K=TT4 with its first physical vertex marked. Then
`(e,u,w)=(1,8,8)` and alpha_K=5/8, exactly the secant (7).
The actual twenty-vertex strong wrappers have states

    z_AB=(27775899,1111968,13893147,186463656),
    z_BA=(27775899,1038240,13893147,186389928).                 (9)

Their difference is `(0,73728,0,73728)`. Both R births therefore
hit the same exact target213127587, and their resulting states are

    R z_AB=(213127587,13893147,13893147,359034165),
    R z_BA=(213127587,13893147,13893147,358886709).            (10)

The current count, next count, deletion coordinate, root name, and child
name agree. The retained two-cover coordinate differs by147456. Before
the birth the defect `F−h−q+s` differs by the same147456. This realizes
the `(0,v,0,v)` direction in the exact one-step continuation fibre.
It complements the separate pure-deletion repair in this session;
the baselines differ, so these two examples are not a freely additive
two-control machine.

**Gcd boundary.** Both source gcds and both returned root gcds in (9)–(10)
equal9. The operation changes the defect but does not improve this gcd.
A and B each have3 dividing every ordered-cover coordinate. They remain
disjoint uniform modules under all R/N births based outside them. Every
compressed path-count term contains one cover factor from each, so9
continues to divide the scalar count. Breaking that obstruction requires
a selective contact inside one of those blocks or another source change.

For any tournament U with H(U)=m, append U as an ordinal suffix to each
whole source. All four retained rooted coordinates multiply by m, by (2).
The same swap and one R birth then give an exact-target pair at
`m*213127587` with retained F gap `m*147456`. This is an actual unbounded
family of targets (for example m=3^r); the suffixed graphs need not be
strong. It is not positive coverage of arbitrary targets.

## 6. Finite probe, controls, and unresolved selection

The probe uses the inherited complete strong bank through order7, a
singleton, and the five previously certified eight-vertex repair blocks
C,C2,C4,C24,N. Both converse gauges are retained for roots. This gives
160 distinct unmarked triples and3443 root gauges. There are23 block-pair
secants in `(0,3/2]`, and no match with a directly strong root block from
that bank. Allowing at most6 transitive vertices on each root flank and
the construction `L→root→R` with L,R empty or in the same unmarked block
bank produces two compensating block pairs.

The smallest is (9). The second uses H133/H129 blocks and a twelve-vertex
root block; its two source and first-child root gcds are27. Its complete
data are in the receipt. Neither found pair removes its uniform3 content.
The requested compensation with both outer blocks prime to3 is unclosed
within this explicitly tested grammar; no impossibility beyond it follows.

`mixed_probe.py` records every admissible secant and the smallest matching
construction for each. The block and root counts are exact subset path
and slot counts. Small literal controls check (2), (4), and the H17
boundary, while the displayed seven-block triples receive an independent
permutation/cut recount. No twenty-vertex census is used: the source
formulas give the large literal states through the proved gluing laws.

Run from the repository root:

    python3 -B 05-knowledge/results/h17_exact_repairs_20261112/ports/mixed_probe.py

The output is `mixed_receipt.json`. The exact next gate is to realize a
compensating secant and a useful defect/gcd change at the required target
and retained ports, with the same source baseline. A rational slope match
without its physical rooted block would not discharge that gate.
