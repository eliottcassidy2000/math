# Missing comparisons: a defect-graded tournament grammar

**PROVED:** formulas and symmetry bounds below, with direct proofs.
**FINITE-EXACT:** full partial-orientation census through order four, literal
signed pair-orbit checks through seven, first four defect coefficients through
20. The underlying Burnside enumeration of oriented graphs is classical; no
priority claim for that enumeration or the general use of cycle indices.
The contribution here is its explicit integration with the repository's
size/attachment grammar and the sharp defect and marked-toggle guards.

## 1. Type of object and inherited hostile

A partial tournament has, on each unordered vertex pair, one of three states:
forward, backward, or unknown. Unknown is a comparison to be filled, not an
assigned winner. Let p(n,d) count vertex-isomorphism classes with exactly d
unknown comparisons, and P_n(t)=sum_d p(n,d)t^d. Thus P_n(0)=A000568(n).
At t=1 this is the classical count of oriented graphs (A001174).

The earlier [four-arc square](../../partial_four_tournament_20261011c.md)
already distinguishes four completion patterns with the same number of missing
pairs. Its lesson survives: d controls a symmetry budget but not the completion
relation. Defect in this note means *missing comparisons*, not the older
almost-automorphism disagreement count in anti_defect_burnside_s20fk.py.

## 2. Exact cycle-type polynomial

For a partition lambda=(l_1,...,l_r) of n write
z_lambda=product_l l^(m_l) m_l!. Define

    I_l(t)=t^(l/2 if l even else 0) (2+t^l)^floor((l-1)/2),
    K_lj(t)=(2+t^lcm(l,j))^gcd(l,j).

Then

    P_n(t)=sum_(lambda partitions n) (1/z_lambda)
              product_i I_li(t) product_(i<j) K_li,lj(t).       (1)

Proof: a vertex permutation acts on unordered pairs, with an orientation sign.
Inside one l-cycle, each distance below l/2 gives an orbit of length l that
returns without reversal. It has two consistent directed choices, or all l
pairs unknown, yielding 2+t^l. If l is even its antipodal orbit has length l/2
and reverses orientation on return. Only the all-unknown choice is possible,
yielding t^(l/2). Between separate cycles l,j there are gcd(l,j) pair orbits
of length lcm(l,j), each without reversal. Multiply the independent choices
and average over permutations using class size n!/z_lambda. QED.

Consequences:
- Setting t=0 kills all even vertex cycles and recovers ordinary tournament
  Burnside exactly.
- In coefficient t^d, total even-cycle vertex mass is at most 2d. Thus arbitrary
  size has only a bounded even-cycle core at fixed defect.
- Free pair orbits longer than d cannot be erased and contribute just 2.
  The partition sum can therefore be truncated before its polynomial products.
- At n=4 the coefficients, in ascending defect, are
  4,10,12,10,4,1,1. Their sum is 42; the constant coefficient is 4.

This is an exact two-parameter classification of counting strata, not a
claim that all structures with the same (n,d) have the same future behavior.

## 3. Recursive generating function

Introduce cycle variables y_l of weight l and the formal series

    F(y;t)=sum_lambda product_l y_l^m_l / z_lambda
                     product_i I_li product_(i<j) K_li,lj.

The empty partition contributes 1. Removing a distinguished l-cycle gives

    partial F / partial y_l = (I_l/l) F((y_j K_lj)_j;t).         (2)

The factorial denominator supplies the derivative factor; every remaining
j-cycle acquires exactly K_lj for its interaction with the distinguished
cycle. Equation (2), with F(0;t)=1, determines the coefficients recursively.
The desired ordinary bivariate series is

    sum_(n>=0) P_n(t)x^n = F(x,x^2,x^3,...;t).                  (3)

For explicit finite recursion, let A_n(w;t) be the degree-n specialization
F((w_j x^j)_j;t). Euler's weighted derivative gives

    A_0(w;t)=1,
    n A_n(w;t)=sum_(l=1)^n w_l I_l(t)
                      A_(n-l)((w_j K_lj(t))_j;t).             (4)

Set every w_j=1 to obtain P_n. At t=0 only odd l survive. Truncating all
polynomials modulo t^(D+1) computes defect strata d<=D without changing them.
This is an enlarged state recurrence. Fixed-order linear recurrences with
polynomial coefficients are ruled out for A000568 by its growth; see ../limits/.
Nonlinear or exponential-coefficient scalar recurrences are not excluded.
No polynomial-time improvement over partition enumeration is claimed.

## 4. Every missing comparison pays at most one binary symmetry bit

**Theorem.** If X has d unknown pairs and G=Aut(X), then

    v_2(|G|) <= d.                                           (5)

Proof: choose a Sylow 2-subgroup S of G. It acts on the 2^d labelled full
completions of X. A stabilizer in S is a 2-subgroup of a full tournament's
automorphism group. Such a group is trivial: any nontrivial finite 2-group
has an involution, and an involution of vertices swaps some pair, reversing
that pair's directed edge. Hence the action is free. Therefore |S| divides
2^d, which proves (5). QED.

Sharpness: take d ordered blocks of two vertices, leave the pair within each
block unknown, and direct every comparison between blocks forward. The
blocks have different partial outdegrees, so every automorphism preserves
each block. Aut(X)=(C2)^d, attaining equality. For fixed d any number of
additional ordered singleton blocks preserves equality.

This sharp bound is not a free bit for faithful representation. Erasure
removes a comparison and may introduce a symmetry; filling the comparison
breaks that symmetry. The completion action remembers this exact cost.

## 5. One-hole recursion and marked metagraph edges

Let e_n count tournaments with one unordered vertex pair marked, modulo vertex
isomorphism, and let v_n count tournaments with one vertex marked. Then

    2 p(n,1) = e_n + v_(n-1), n>=2.                          (6)

Proof: reversing the marked comparison is an involution on marked tournament
classes. Its orbits are precisely one-hole partial classes. A fixed marked
class is a partial class with an automorphism swapping the missing pair's
endpoints. The induced automorphism on the remaining complete tournament has
odd order. Raise to that odd order to obtain the transposition fixing all
other vertices. Thus the two endpoints have identical incidences with every
other vertex. Conversely such twins can be swapped. Collapsing the twin pair
is a bijection with vertex-marked tournaments on n-1 vertices. Burnside for
this involution gives (6). QED.

For n=2,3,4,5,6 the triples (e_n,v_(n-1),p(n,1)) are
(1,1,1), (4,2,3), (16,4,10), (88,12,50), (704,48,376).
The twin classes are the fixed points of the *marked* reversal. Unmarked
endpoint classes may coincide for additional reasons when an isomorphism
moves the marked pair, so this is not a classification of every loop in the
ordinary simple metagraph. Retain pair marks and parallel-edge incidences.

An independent hostile already occurs on four vertices. In pair order
(01,02,03,12,13,23), encode +1 as lower-to-higher and -1 as reverse.
The tournament (-1,-1,+1,-1,-1,-1) is isomorphic to its reversal of pair
(1,3), yet those endpoints are not twins. Thus an unmarked metagraph loop
need not be a fixed point of the marked reversal involution. The independent
audit checks this example directly.

The corresponding formal series satisfy

    2 P_1(x)=E(x)+x V(x),

where P_1=sum_(n>=2)p(n,1)x^n, E=sum e_n x^n (e_0=e_1=0), and
V=sum_(n>=1)v_n x^n. This relates the first defect layer, marked toggle
incidences, and pointed growth without the incorrect division by n.

## 6. Exact validation and scope

Run `python3 check_defects.py` and `python3 -O check_defects.py`; outputs must
agree. The checker independently walks actual permutation pair orbits and
orientation signs through n=7, compares full brute ternary orientations
through n=4, checks the binary budget for every resulting class, checks the
one-hole identity through n=15, and matches A000568 through n=12.
Known A000568 inputs are a regression control, not the proof of (1).
The coefficient tables through n=20 are derived by (1), not independent
complete class enumeration. A second recurrence checker validates (4).

Classical context: Harary and Palmer, *Graphical Enumeration* (1973),
[official publisher record](https://shop.elsevier.com/books/graphical-enumeration/harary/978-0-12-324245-7), chapters on orientations and cycle indices.
The formulas needed here have direct proofs, so no inaccessible theorem is
being silently imported. No assertion of global novelty is made.
