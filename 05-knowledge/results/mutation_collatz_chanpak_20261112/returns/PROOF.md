# Triangle transport changes endpoints but supplies zero ternary digit

**PROVED scoped structural consequences of inherited path-cover laws.**
**FINITE-EXACT:** the declared small graph and polynomial controls.
The ordinary-mutation and Chan–Pak completeness questions remain **OPEN**
outside the stated theorems. No priority claim is made for the inherited
birth identity or ordinal-product mechanism.

## Inheritance and the observer

The closest mechanisms are the ordinary FZ
[four-return classification](../../triples_ports_sieve_20261112/ports/FOUR_RETURNS.md),
the [unique-start funnel calculus](../../polarized_path_calculus_20261103/FUNNELS.md),
and [ternary-ready roots](../../two_entry_funnels_20261103/READY.md).
The current [endpoint frontier](../../core_descent_endpoints_20261112/README.md)
warns that even the ordinary17-coordinate observer can forget a selected
endpoint-cover gate. We do not identify the small observer below with that
full selected response.

For a rooted tournament (T,r), write

    H=H(T), s=s_T(r), q=H(T-r), F=E_(T,r)(2), b=F-H.

Here P_T(m) counts ordered m-bin spanning path covers, allowing empty bins;
S_(T,r)(m) and E_(T,r)(m) require the first bin to start or end at r.
The first bin in S/E is therefore nonempty. Let H(empty)=1,
P_empty(m)=1, and P_T(0)=0 for nonempty T.

A birth N_r adds x beaten only by r: r->x and x->T-r. This is an actual
tournament birth, not ordinary FZ mutation. Marking the latest newborn,
the inherited laws are

    (H,F,s) -> (H+F-s, H+F, H-s),
    N^2: H -> H+2F = 3H+2b.                         (R0)

For an indegree-one root with sole predecessor a, the inherited defect
identity is b=q-s=s_(T-r)(a). Merely having indegree one does not force
b=1 or2.

The live board is ordinal factor / path-cover polynomial / named root /
positive readiness / arithmetic factor / actual birth. The hostile is
three available roots whose defects are all zero. The corrected near miss
is replacing the guard b=1 or2 by the weaker indegree-one condition.

## R1. The whole endpoint response of an extreme triangle transfer

Let C be a directed triangle and K any tournament, possibly empty. Set
h=H(K). For every nonnegative integer m,

    P_C(m)=m^3+2m,
    S_(C,r)(m)=E_(C,r)(m)=m^2-m+1       (m>=1).

For P, classify bin assignments by whether all three vertices occupy one,
two or three bins. This gives3m+3m(m-1)+m(m-1)(m-2).
For an endpoint root, the first bin has three vertices in one valid order,
two vertices in one valid order with m-1 choices for the last vertex,
or just r with the other two forming a transitive edge in m-1 bins.
Thus the total is1+(m-1)+(m-1)^2=m^2-m+1.

The standard ordinal-product bijection concatenates the two factors'
paths within every bin. Consequently, for T=C->K and r in C,

    P_T(m)=(m^3+2m)P_K(m),
    S_(T,r)(m)=(m^2-m+1)P_K(m),
    E_(T,r)(m)=(m^2-m+1)P_K(m-1).                  (R1)

The endpoint condition in E forbids K vertices in the first bin.
For T_before=K->C the last two formulas exchange places. Reversing
the triangle during the four-mutation transfer does not change them,
since the three labelled triangle roots have identical one-root responses.

Thus the transfer preserves the entire unrooted polynomial P, not merely
H. It exchanges entrance and exit budgets. In particular, at the three
front triangle roots,

    H=3h, s=q=h, F=3h, b=0.                       (R2)

Before the transfer their F was3P_K(2), so their defect was
3(P_K(2)-h). This is a real count-preserving source change but it lands
on the zero-digit locus. The triangle's three rotations fix K and permute
the roots, so three roots here are not three independent positive digits.
The literal names and cyclic orientation still matter for joint queries.

The first two births at any such root therefore have counts

    3h -> 5h -> 9h,                               (R3)

and rooted states

    (3h,3h,h) -> (5h,6h,2h) -> (9h,11h,3h).

The first newborn has defect h; the second has defect2h.
They give defects1 and2 only when K has h=1, equivalently is transitive.
This boundary explains the small positive sequence3->5->9 and why it
does not furnish those same digits at arbitrary triangle-factor counts.

## R2. The factor trap for the prime and7p lanes

If a birth is marked in the current prefix A of A->K, the newborn beats
every K vertex, so the result is A'->K. Inductively every sequence of
these births and arbitrary re-markings inside A preserves K as a literal
ordinal suffix and keeps

    H=H(A)H(K).                                   (R4)

When A contains the original triangle as an induced graph, insertion
monotonicity gives H(A)>=3. If H(K)>1, (R4) can never be a prime or7p
for prime p: the latter factorization would force one of its two
nontrivial realizable factors to equal7. That count is impossible by
**THM-1370**,
[h-spectrum-omits-7-21-all-n](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
This includes p=7. The statement is about this preserved ordinal suffix;
operations crossing it are outside the argument.

There is an immediate-birth dichotomy even if the first mark may lie
anywhere in C->K:

- mark in C: the count is5h;
- mark in K: rotation of the untouched C is still an automorphism after
  the birth, and acts freely on Hamilton vertex words, so3 divides the
  output count.

No such first birth reaches a prime greater than5 or any realizable7p.
In the first case7p=5h would require p=5 and h=7. In the second case
3|7p would require p=3 and the forbidden count21. The possible prime5
occurs precisely in the first case with transitive K.

Finally, either nonzero digit3H+2 or3H+4 is coprime to3. At that output
no same-count representative can have an ordinal triangle factor.
Consequently the extreme-triangle transfer cannot be the next
count-preserving return for such a branch. A future induction needs
another source family or a count-changing intermediate route with an
exact return to its original target.

## R3. Exactly two indegree-one roots: positive readiness with its guard

Use the [two-root normal form](../mutation/TWO_ROOTS.md):

    a->b->c->a,     a,b -> R,

with arbitrary incidences between c and R. Put W=R+c,
h=H(R), d=H(W), s=s_W(c), and F_W=E_(W,c)(2).
This is exactly the inherited pair of births from (W,c).
The two named roots a,b have

    root a: start=d, deletion=d, F=H(T), defect0;
    root b: start=F_W, deletion=d+F_W-s,
            F=H(T)+d-s, defect d-s.               (R5)

For root a, T-a=b->W, so its sole predecessor c cannot start a
Hamilton path in T-a. For b, a is its sole predecessor, and paths in
T-b starting at a are exactly paths of W not starting at c:
prepend a to the latter, since a beats R but loses to c.
This proves the defect d-s. The remaining coordinates follow R0.

When c has an in-neighbour in R, insertion of c into each Hamilton path
of R gives an injection into paths of W not starting at c.
More explicitly, if the path's last vertex beats
c, append c; otherwise use the last vertex along the path that beats c
and insert immediately after it. If no vertex beats c, no such insertion
is possible; that is precisely the case c is a source and defect zero.
Hence, provided c has an in-neighbour in R,

    d-s >= h.                                    (R6)

This qualification is essential. The exactly-two-root case has such an
in-neighbour; with none, all three triangle vertices have indegree one
and R1 applies.

Therefore positive defect1 or2 at b forces h=1, so R is transitive.
Write its order r1->...->rm and I={ri:ri->c}. Let I be nonempty.
The entire low-ready classification is then:

| Defect at b | Necessary and sufficient cut | H(T) |
|---|---|---|
| 1 | I={r1,...,rj}, 1<=j<=m | 1+2^(j+1) |
| 2 | I={r2}, m>=2 | 9 |

Proof. The start count s is0 if r1 is in I and1 otherwise.
Since d is odd, d-s=1 forces(d,s)=(1,0). Then W is transitive,
so I is a nonempty initial segment. There are2^j-1 nonempty subsets
of R whose path ends in I. These count the ordered two-path covers
U_I in the normal-form formula, giving
H(T)=1+2[1+(2^j-1)]=1+2^(j+1).

For d-s=2, necessarily(d,s)=(3,1). A tournament with H3 consists of
one C3 SCC and transitive singleton padding: this is the inherited
small-count SCC classification underlying THM-1370. Since R is transitive,
c belongs to that triangle. Since s_W(c)=1, no prefix precedes it.
Its other vertices must be r1,r2, with c->r1->r2->c; every later
vertex lies in the transitive suffix. Thus I={r2}, and conversely that
cut gives d=3,s=1. Here U_I=2 and H(T)=3+2(1+2)=9.

The first family includes prime source counts5,17,257, but never a7p
source: powers of2 modulo7 are1,2,4, none equal-1. The second source
count is9. These families can provide real initial branches
3(1+2^k)+2 and31, respectively; the first has k>=2.
They are not a universal ready-source bank.

There is an exact weak-source recursion obstruction at the two
indegree-one roots themselves. Apply two births to the positive-ready
root b of this normal form. The only indegree-one roots of the output
are the two newborns: every original root or exterior vertex acquires
enough new incoming arcs to have indegree at least two. The first
newborn has defect0. The latest has

    H(T)-s_T(b) = d+F_W >= 3.                     (R7)

Indeed H(T)=d+2F_W and s_T(b)=F_W. Since I is nonempty,
U_I>=1 (choose one vertex of I for the first path and a Hamilton path
of the remaining R for the second), whence F_W=h+U_I>=2 and d>=1.
Thus even with arbitrary transitive suffix padding, neither new
indegree-one root is1- or2-ready. In the first row of the table the
new positive defect is1+2^j; in the second it is6. This statement
does not on its own exclude a root of larger indegree.

At the inherited standard two-birth output of a strong source of order
at least4, no re-marking alone restores defect1 or2; see READY section4.
That theorem explains why these positive starts still require a genuine
source replacement for recursion. Three zero-ready roots and two roots
with one positive-ready root serve different roles.

The [named Collatz lift](COLLATZ_LIFT.md) gives a positive use of the
first family: delete its ready b to realize(H+1)/2 down the entire
almost-transitive power-of-two staircase. At H5, two births followed by
their reverse deletions give an actual repeatable5->17->9->5 arithmetic
cycle, with the transient graph charts retained.

## Connection audit and arithmetic boundary

| Source and map | Preserved predicate | Lost or additional data |
|---|---|---|
| K->C to reversed-C->K by four ordinary mutations | H and all P_T(m) | Endpoint budgets exchange; names and joint cyclic data remain necessary |
| Front triangle root to one/two births | Exact R0 response | H changes; the retained suffix creates R4 |
| Two-root gadget to(d,s,F_W) | Exact defects and birth counts | Does not determine arbitrary named attachment queries |
| Collatz table label to a3-divisible intermediate | Arithmetic equality | A residual tournament realizing the required quotient still needs construction |

The moving-modulus Collatz encoding can point to a3-divisible label,
but no triangle transfer creates a tournament with that label from a
different count. If the desired label is3h, installing C3->K requires
an actual K with H(K)=h. This is an explicit realization obligation,
not a consequence of divisibility. The forbidden quotient h=7 gives
the sharp hostile3h=21.

## Reproduction

Run check.py in this directory, normally and with -O. It checks every
labelled K of orders0..4 (76 inputs), the response laws through m=3,
all triangle-root birth controls, the first-birth divisibility dichotomy,
and every transitive-residual two-root cut through length7.
Selected small controls are also counted by literal permutations.
The proofs, not those finite universes, give the all-order statements.
