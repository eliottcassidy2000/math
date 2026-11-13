# Every nonzero half-repair return costs at least a factor five

**PROVED, scoped to identified-root contexts and the indicated repair
direction.** This is not a lower bound on all tournament constructions.

Let z(B,r)=(h,s,q,F), where h counts Hamiltonian paths, s counts those
starting at r, q=H(B-r), and F counts ordered two-path covers whose first
path is nonempty and ends at r; the second may be empty. Put E=h+F-s.
The actual birth r->v->B-r has Hamilton count E. Let

    d=(0,1,2,1),       P={F+s=h+q}.

Take a tournament A of order at least two, with entry r and output u in A.
Identify the roots of A and B and orient A-r -> B-r. The inherited port
compiler gives z' = M(A,r,u)z for every B of order at least two.
Assume Md=lambda*d with lambda nonzero. Then:

1. H(A)>=5.
2. For every actual input on P,

       H(A star B)=e_A h+(h_A-e_A)E >= h_A h >= 5h,

   where e_A counts Hamiltonian paths of A ending at its entry.
3. If the chosen returns also preserve P, every k-step iteration grows
   H by at least 5^k. This holds even when the contexts have unbounded order.

The singleton context is the identity and is excluded. A zero image of d
is not a return: it destroys the distinguished repair.

## Neutrality gives the target row

Write the seven entry statistics of A as (h_A,s_A,e_A,q_A,F_A,G_A,P_A),
with G the reversal analogue of F and P the ordered two-path-cover count
of A-r. By the proved gluing law in
[PG1](../endpoint_context_recursion_20261113/PORT_GRAMMAR.md) and
[the seven-coordinate law](../endpoint_context_recursion_20261113/OBSTRUCTIONS.md),
the first row is

    h'=q_A h+(e_A-q_A)s-(G_A-h_A-q_A+e_A)q+(G_A-q_A)F.

Its value on d is 2h_A-G_A-e_A. Thus every return requires

    G_A+e_A=2h_A.                                      (1)

On P substitute q=E-2h+2s and F=E-h+s. Equation (1) yields
h'=e_A h+(h_A-e_A)E, with nonnegative coefficients.
For every actual rooted tournament F>=q>=s: take the first path to be
the singleton root for the first inequality, and delete the initial root
of a Hamiltonian path for the second. Thus E>=h, proving the bound once
H(A)>=5 is established.

## Why contexts with H=1 or H=3 cannot return

Hamilton counts multiply over the ordered strong components. A strong
tournament of order n>=3 has a Hamiltonian cycle, hence at least n
Hamiltonian paths. With the oddness theorem this implies:

- H=1 precisely for a transitive tournament;
- H=3 precisely for TT_a -> C3 -> TT_b, with a,b>=0.

Here we use the standard Camion theorem and Hamilton-path oddness already
routed by `05-knowledge/reference/CORE-PAPERS.md`; no new classification of
general small counts is asserted. The following covers every entry and
output, including arbitrary transitive padding.

**Transitive A.** If t vertices follow r, G_A=2^t and e_A is 1 for t=0
and 0 otherwise. Condition (1) permits only t=0 or t=1.

- For t=0 the glued tournament is L->B, with nonempty transitive L.
  At an output in L, deletion has count H(L-u)h, so its difference along
  d is zero. At r, the difference is (0,0,2,P_L), where P_L=2^|L|>0.
- For t=1 it is L->birth(B,r). At an output in L, deletion is a fixed
  multiple of E, so its difference is zero. At the new birth vertex,
  deletion has count h and difference zero. At r the difference is
  (0,2,2,0) when L is empty, and (0,0,2,0) otherwise.

None is a nonzero multiple of d. Outputs with zero q difference cannot
qualify, since d has nonzero q coordinate.

**H(A)=3.** Write A=TT_a->C3->TT_b.

- If r is in the prefix, then G_A>=12, since the two-path-cover count
  of C3 is 12. This violates (1).
- If r is in the suffix, and t suffix vertices follow it, then G_A=3*2^t
  and e_A is 3 for t=0 and 0 otherwise. Again only t=0,1 qualify.
  These are L->B and L->birth(B,r), with H(L)=3. Outputs in L or at
  the new birth vertex have zero q difference. At r their differences
  are respectively (0,0,6,P_L) and (0,0,6,0), with P_L>0.
- If r is in C3, then G_A=3*2^b and e_A is 1 when b=0 and 0 otherwise.
  Thus (1) forces b=1. If a>0, every output after the prefix has s=0;
  an output inside the prefix has zero q difference, since deleting it
  leaves the same neutral core count. Neither can return d. It remains
  only the four-vertex context C3->singleton, rooted in its cycle.

Represent that last context by out-neighbour masks (10,12,9,0), entry0.
At outputs 0,1,2,3 its exact images of d are respectively

    (0,2,2,-2), (0,0,0,0), (0,0,-2,-2), (0,-1,2,7).

These follow by the finite color-word compiler (or direct enumeration of
the constant-size context). None qualifies. This exhausts H=1 and H=3.
Since tournament Hamilton counts are odd, every nonzero return has H>=5.

## Composition does not create a slow return

Port contexts are closed under composition, even if their intermediate
repair directions differ. If A1 has entry r1 and output u1, and A2 has
entry r2 and output u2, form C by gluing A2 to A1 at r2=u1, with
A2-r2 -> A1-u1. Keep r1 as C's entry and u2 as its output. Every new
vertex dominates the original B-r1, so, after the natural relabelling,

    A2 star_(u1) (A1 star_(r1) B) = C star_(r1) B.

The internal edges and every cross edge agree. Consequently any finite
word of these operations whose net matrix returns d nontrivially is
itself covered by the theorem. The factor five applies to the net
return, not to each intermediate step. Cycling through other slopes
therefore cannot produce a slow net half-return within this operation.

## Consequence and boundary

The factor five is an all-order obstruction to embedding a slow
Calkin--Wilf generator as a repeatable return in this state space. In
particular a unipotent target matrix cannot be such a return: its powers
grow polynomially, contradicting the exponential lower bound on any
actual nonzero input. Other slopes can enlarge the branching language,
but the composition lemma excludes a slow return obtained solely by
visiting those slopes. A different seam, an output chosen inside the
attachment away from its entry, or a rule that discards this repair and
supplies another one is outside the theorem.

There is a precise incompatibility with the inherited prime-lane ternary
induction in [Session 193, Theorem 2](../session193_triangulated_proof_20261105/prime_induction/PROOF.md).
For h>=3, both 3h+2 and 3h+4 are less than 5h. Thus neither ternary
child can be produced from a parent on P by a nonzero net half-return
of these contexts. This does not refute that sufficient induction
criterion: it asks only for an existential tournament witness, with no
requirement to retain this paired interface. The extra invariant is ours.

It also does not imply that every finite alphabet of returns is sparse.
Many expanding branches can have a large scalar projection. The inherited
nine-rule language has a separate weighted-branching proof of sparsity.

The positive control U has h_A=5,e_A=1 and h'=h+4E. The other eight
catalogued return matrices also have first row (e_A,h_A-e_A), checked
independently from their entry profiles. Exact controls exhaust the H=1
and H=3 templates through order10: 654 actual-output comparisons, 74
abstract-compiler comparisons, and 10 permutation comparisons. The proof
above, not the finite cutoff, supplies the all-order quantifier.

Reproduce: `python3 05-knowledge/results/strategy_reframe_20261113/growth_probe.py`.
