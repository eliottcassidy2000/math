# Five-dimensional interfaces, the odd antipode, and seven sufficient coordinates

**2026-11-11. PROVED IN THIS NOTE**, relative to the inherited odd-cycle
formula and Session251 compiler. No independent canon promotion or priority
claim. The five-core all-height covering statement remains OPEN.

## 1. A genuine five-dimensional object attached to every coloured core

For a tournament C with a distinguished apex cut, retain the residual-H-weighted
path polynomials H=d,K(z),I(z),O(z),J(z) of
[Session251](../../session251_creative_completion_20261111/ternary/PROOF.md).
Here a path ending in the apex in-set contributes to I, starting in its out-set
to O, and satisfying both to J. Define

    T_C(z) = [[d, O(z), J(z)],
              [0, d+K(z), I(z)],
              [0, 0, d]].

The [golden note](../golden/PROOF.md) proves by classifying paths that
T_(C▷D)(z)=T_C(z)T_D(z), with cuts joined in the same ordinal order.
This is not arbitrary matrix packaging: the O_C I_D term is precisely a
path starting in the first coloured block and ending in the second.

Let R be any commutative coefficient ring. The matrices

    A(R)={[[d,o,j],[0,e,i],[0,0,d]] : d,e,o,i,j in R}

form a rank-five associative R-algebra. It has an intrinsic description:

    A(R) = {X in M3(R) : X N=N X},  N=E13, N²=0.

Indeed XN=NX forces x21=x31=x32=0 and x11=x33, and imposes no other
conditions. Multiplication of its five coordinates is

    (d,e,o,i,j)·(D,E,O,I,J)
      =(dD,eE,dO+oE,eI+iD,dJ+oI+jD).                 (A1)

Over a field k, regard k³ as a module for k[ε]/(ε²), with ε acting by N.
Its Jordan blocks have lengths2 and1, so

    A(k) ≅ End_(k[ε]/ε²)(k[ε]/ε² ⊕ k).

Thus the five-dimensional carrier is an endomorphism algebra of a two-step
nilpotent module, not a numerological identification of five cores with
five solids. Its strictly upper-triangular ideal has basis E12,E23,E13,
with E12 E23=E13 and all other products of these basis elements zero.
Its semisimple quotient over a field is k×k. Equivalently, it is the
five-dimensional path algebra of two vertices and opposite arrows with
one of their two length-two composites set to zero.

**Scope.** The polynomial matrices of actual coloured tournaments occupy a
proper arithmetic subset. The ambient algebra is not a theorem that every
five-tuple comes from a tournament. Actual row masks and the cut remain the
realization sidecar. The supplied finite exact basis demonstrates that
positive evaluations of the five-core cuts span all five ambient directions.

## 2. The odd-cycle antipode kills the negative middle channel

There is a universal relation not needed in the original nine-coordinate
presentation:

    K_C(−1)=−H(C)  for every nonempty tournament C.     (A2)

**Proof.** Work in the squarefree commutative algebra
Q[x_v:v in V]/(x_v²:v in V). Put

    h_C(x)=sum_(S⊆V) H(C[S]) x_S.

The inherited odd-cycle formula, including singleton vertices, is equivalent
to

    h_C(x)=exp(sum_v x_v + 2 sum_γ x_(Vγ)),

where γ runs through directed odd cycles of length at least3, counted once
up to cyclic rotation. Expanding the exponential picks disjoint cycles;
the factorial denominator cancels the ordering of a chosen collection,
and the unused chosen vertices contribute singleton factors. Every monomial
inside the exponential has odd degree. Consequently

    h_C(x) h_C(−x)=1.

For V nonempty, its x_V coefficient is zero:

    sum_(S⊆V) (−1)^|S| H(C[S])H(C[V\S])=0.

Every directed simple path with vertex set S is a Hamiltonian path of C[S],
so K_C(z)=sum_(S nonempty) z^|S| H(C[S])H(C[V\S]). The empty-set summand
is H(C), proving (A2). QED.

This is the historical ±1/fugacity bridge made operational: odd primitive
cycles make sign reversal into a convolution inverse. The matrix T_C(−1)
therefore has middle diagonal zero. If C,D are nonempty,

    O_(C▷D)(−1)=d_C O_D(−1),
    I_(C▷D)(−1)=I_C(−1)d_D.

Signed endpoint information sees only the corresponding extreme block,
scaled by H of the other block. This explains an exact boundary erasure;
it is not an excuse to erase positive endpoint channels.

## 3. Nine coordinates compress to seven, and the seven compose

Write

    e=d+K(1), u=O(1), u'=O(−1), v=I(1), v'=I(−1),
    q=J(1)+J(−1).

Then the response state is the seven-tuple

    s(C)=(d,e,u,u',v,v',q).

For all nonnegative binary flank addresses A,B, α=A mod2 and β=B mod2,
the exact child count simplifies to

    F_C(A,B)=d(g(A)+g(B)−1)+eAB+vA−v'α+uB−u'β+q.     (A3)

The apparent αβ interaction in Session251 vanishes for every tournament,
not just these five cores. No information is lost for this query language:
(A2) determines K(−1), and only the sum J(1)+J(−1) was observed.

More strongly, the seven-tuples are closed under ordinal concatenation:

    d_CD = d_C d_D,               e_CD = e_C e_D,
    u_CD = d_C u_D+u_C e_D,       v_CD = e_C v_D+v_C d_D,
    u'_CD = d_C u'_D,            v'_CD = v'_C d_D,
    q_CD = d_C q_D+q_C d_D+u_C v_D+u'_C v'_D.          (A4)

These formulas follow directly from the positive and negative matrix
products. They are associative, and compute the same observer after any
future ordinal concatenation of nonempty cores. The empty core has a
separate identity state, since (A2) requires nonemptiness.

This seven-tuple is **not** Session250's different seven-coordinate rooted
observer. Their coordinates and legal operations differ. Equations (A3)
and (A4), rather than the equal number seven, are the connection contract.

## 4. Linear minimality for the actual response language

Seven is not merely an available upper bound for linear response summaries.
The seven functions of (A,B) occurring in (A3) are linearly independent over
Q. Set B=0: g(A), A, A mod2,1 are independent on A=0,1,2,3, using
(g(0),g(1),g(2),g(3))=(1,1,3,1). Hence any identity forces the coefficients
of d,v,v',q to vanish. Setting A=0 then separates B from B mod2; finally
A=B=1 separates AB. All seven coefficients vanish.

The finite exact artifact gives seven actual bank-core/cut signatures with
nonzero rational determinant. Thus actual responses already span seven
independent directions. No fixed linear summary of fewer than seven scalar
coordinates can recover every response in this family. This says nothing
about nonlinear encodings of integers, a fixed uncoloured core, or smaller
restricted query languages.

## 5. A concrete route forward, without reintroducing fake fives

The five core charts now live in a proved seven-dimensional response
semigroup whose positive channel is a canonical five-dimensional algebra.
The missing global problem can be stated on this smaller carrier:

- Find a target-reducing rewrite that maps every allowed odd target and
  an applicable response chart to a smaller target with a realizable chart.
- Require it to preserve (A3), not merely H or the triangle-overlap diagram.
- Bound a proper height of the full chart/address state. The recovered A4
  golden visible/conjugate construction warns that contraction of one
  displayed coordinate is insufficient.
- Use (A4) for lawful regrouping and the Catalan pentagon; do not remove the
  cross term u_C v_D+u'_C v'_D. Associativity alone would also hold for some
  incorrect contractions that omit real paths.

The result here is a smaller exact carrier and a useful algebraic language,
not the missing decreasing rewrite. The original finite coverage remains
inherited, and five-core all-height coverage remains OPEN.

## Reproduction

Run `python3 algebra/probe.py` from this packet, or pass its full path from
any directory. Only the inherited standard-library compiler is imported.
The artifact records signed-antipode controls on all labelled tournaments
through order4 and all cuts, basis determinants, exact seven-state join
checks, and direct response comparisons. These tests support the new
identities; they are not a new global coverage scan.
