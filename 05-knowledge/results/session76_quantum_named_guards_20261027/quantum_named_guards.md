# Seven named guards, three quantum bases, and an SU(3) Cartan direction

## A faithful state encoding with deliberately lossy observations

For a named tournament let `K_ij=1` when `i→j`, `K_ji=−K_ij`, and `K_ii=0`.
For `n≥2` and `0<ε<1/[n(n−1)]`,

    ρ(T)=I/n+εiK(T)

is a strictly positive trace-one density matrix. The absolute row sum of
`iK` is `n−1`, giving the eigenvalue lower bound `1/n−ε(n−1)>0`.
Relabelling acts by permutation-unitary conjugation. If two labelled
tournaments disagree on `e` arcs, then

    ||ρ(T)−ρ(U)||_HS² = 8ε²e.

This transports labelled arc disagreement, not native path-move distance.
Every real vector `v` has `vᵀKv=0`, so a real measurement basis gives the
uniform distribution for every input. In contrast, the named probes
`(e_i±ie_j)/√2` have `iK` expectations `∓K_ij` and recover the edge sign.
Grouping pairs in matchings provides `n−1` measurement bases when `n` is
even, or `n` when it is odd. These need not be mutually unbiased.

## The three specified bases

Index six coordinates as `3x+y`, with `x=0,1` and `y=0,1,2`. Let

    ω=(-1+i√3)/2,          F_yb=ω^(yb)/√3,
    Q=diag(1,ω,ω),
    X=[[1,1],[1,−1]]/√2,  Y=[[1,1],[i,−i]]/√2.

Use `B0=I6`, `B1=X⊗F`, `B2=Y⊗QF`. They are orthonormal bases. Squared
qubit overlaps are `1/2`; qutrit overlaps are `1/3`, by the elementary
three-term Gauss sum. Hence squared cross-basis overlaps are `1/6`.
This standard product construction requires no maximum-MUB theorem.

For vertical currents `K_s=[[0,S],[−S,0]]`, `S=diag(s0,s1,s2)`, a product
vector `a⊗b` with `|b_y|²=1/3` satisfies

    (a⊗b)* iK_s (a⊗b) = −2 Im(conj(a0)a1)(s0+s1+s2)/3.

Thus any collection of product bases having flat qutrit factors sees only
the sum of these three named currents. The computational basis sees none
of them. Equal-sum vertical sign patterns have fibre sizes `1,3,3,1`.

Take all arcs `i→j` for `i<j`, then reverse only `(1,4)` to obtain `T`, or
only `(0,3)` to obtain `U`. Their degree lists are respectively
`(5,3,3,2,2,0)` and `(4,4,3,3,1,0)`. The first has a source, the second does
not: they are nonisomorphic. At `ε=1/60`, both states have smallest
eigenvalue at least `1/12`; all eighteen probabilities in the three bases
coincide. Their trace distance and Hilbert–Schmidt distance both equal
`1/15`. The probe `(|0>+i|3>)/√2` instead gives `3/20` and `11/60`.
The additional basis `Y⊗I3` restores the vertical signs but is not a fourth
MUB. Zero response to a current does not mean zero outcome probability.

## An exact, small description of all measurement fibres

Write each upper-triangle sign as `k_ij∈{−1,1}`. Direct multiplication of
the displayed basis matrices shows that their complete probabilities are
equivalent to these seven integer statistics:

    p1=k01+k34,   p2=k02+k35,   p3=k12+k45,
    v=k03+k14+k25,
    sA−sB,       sB−sC,        dA−dB+dC,

where `A=(k04,k13)`, `B=(k05,k23)`, `C=(k15,k24)`, and `s,d` denote each
ordered pair's sum and difference. Both the actual real response matrix
and its rational coefficient-expanded matrix have this seven-dimensional
row space. The included script checks that identity directly from the
bases, not by assuming the statistics capture the observations.

The variables split into three disjoint two-sign blocks, a three-sign
vertical block, and a six-sign cross block. This gives a complete elementary
classification of minimal collision supports:

* Each pair sum has a swap collision; one named sign is necessary and
  sufficient. The three blocks cost three guards, with `2³` choices.
* The vertical triple has a swap between any two coordinates; two of its
  three named signs are necessary and sufficient, with three choices.
* In the cross block normalize a difference of two sign vectors by two,
  obtaining entries in `{−1,0,1}`. Equality of statistics says all three
  pair sums equal a common `t` and `dA−dB+dC=0`. If `t=±1`, the three pair
  differences are odd, so their signed sum cannot vanish. If `t=±2`, all
  six coordinates are nonzero. If `t=0`, each pair is `(r,−r)`, and
  `rA−rB+rC=0`; every nonzero solution uses exactly two groups. Thus the
  minimal cross supports are `A∪B`, `A∪C`, `B∪C`. Two guards in distinct
  groups are necessary and sufficient, with `3·2²=12` choices.

All nine minimal supports occur as actual sign-cube collisions. Because
the blocks are disjoint, the exact minimum is

    3+2+2 = 7 additional fixed named arc-sign guards,

with exactly `2³·3·12=288` optimal sets. One is

    {01,02,03,04,05,12,14}.

For that set the unqueried edge order is
`(13,15,23,24,25,34,35,45)`. The remaining real kernel is one-dimensional,
spanned by `(2,1,2,1,0,0,0,0)`. A normalized difference of two sign vectors
must belong to `{−1,0,1}^8`: a coefficient-one coordinate forces the scalar
to be integral, and a coefficient-two coordinate then forces zero. This is
discrete injectivity despite a surviving continuous kernel. A real multiple
such as one half does fit the continuous unit cube, so the integer-domain
qualification cannot be dropped.

The cross block's 64 sign assignments have 53 signatures: 44 singleton
fibres, seven double fibres and two triple fibres. Combining independent
blocks gives `3³·4·53=5724` signatures on all `2^15=32768` tournaments, with
maximum fibre `2³·3·3=72`, attained four times. Separate exhaustive checks
agreed with this factorization. A 72-element fibre also gives the binary
information lower bound `ceil(log2 72)=7`. It is not a lower bound on a
real-valued observable carrying more than one bit.

Thus the augmented signature is a lossless **fixed-order-six** carrier.
By decoding the tournament, one can reconstruct native-path legality and
outputs, converse and relabelling. No efficient-update theorem, global
minimal encoding theorem, or size-changing attachment closure is claimed.
The observations here are exact probabilities, not a finite-shot recovery
guarantee. Seven guards are not seven MUBs or seven native moves.

## The lost direction is a precisely typed SU(3) direction

Let `J0=[[0,−I3],[I3,0]]`. Realification maps an anti-Hermitian complex
matrix `A=X+iY` to

    R(A)=[[X,−Y],[Y,X]].

Its image is the `J0`-commuting part of `so(6)`, namely `u(3)`, dimension
nine. The anticommuting skew matrices are `[[P,Q],[Q,−P]]`, with `P,Q`
skew, a six-dimensional space. The resulting real SU(3) representation
splits orthogonally as

    so(6) = su(3) ⊕ R J0 ⊕ (Λ² C³)_R,      15=8+1+6.

The last summand transforms by `C↦UCUᵀ`; the determinant form identifies
`Λ²C³` with the conjugate fundamental representation. This is a statement
about representations, not a direct sum of Lie algebras.

The vertical current is `R(−iS)`. Its mean diagonal is the central `u(1)`
part; its traceless diagonal is the two-dimensional Cartan subalgebra of
`su(3)`. The collision above has `S=diag(2,−2,0)`, exactly a traceless
Cartan direction. It is not a physical colour current or a gauge-invariant
QCD observable. Tournament sign matrices are not closed under SU(3)
conjugation or Lie brackets.

The Eisenstein identification is explicit:

    a+bω ↦ (a,−a+b,−b),    ||x||²=2(a²−ab+b²).

Cyclic coordinate permutation acts as multiplication by `ω`. The difference
of its two cyclic Weyl actions acts as `1+2ω`, whose square is `−3`.
After realification its squared norm is `4(a²−ab+b²)`. The direction
`(2,−2,0)` therefore has squared skew-matrix norm 16, matching the collision.
This is a difference of Weyl actions, not a QCD propagator.

For classical background see
[Boykin–Sitharam–Tiep–Wocjan on MUBs and orthogonal Lie-algebra decompositions](https://arxiv.org/abs/quant-ph/0506089)
and the [PDG QCD review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-qcd.pdf).
The full dynamics and continuum existence of QCD are not premises here.

## Two smaller phase and modulus lessons

If `H` and its entrywise square are both complex Hadamard matrices, ratios
`z_j` from two distinct rows satisfy `Σz_j=Σz_j²=0`. If their cubes take
values `a≠b`, and `A,B` are the two fibre sums, then
`A+B=0` and `a conjugate(A)+b conjugate(B)=0`, since
`z²=z³ conjugate(z)`. Therefore `A=B=0`. Within each occupied cube fibre,
integer multiplicities of the three Eisenstein phases must be equal.
The square-Hadamard hypothesis is necessary: the Fourier-six row ratio
`(1,−1,1,−1,1,−1)` has fibre sums `3,−3` and its squared matrix is not
Hadamard. The Tao exception in the related external Fourier statement is
also essential: a separate exact calculation gives its three-plus/three-minus characters as
`(1±i√3)/4`, not zero. That calculation is not repeated by the compact scripts here. Neither fact proves the full MUB maximum theorem.

[Earlier period-fifteen/conductor-eleven work](https://github.com/eliottcassidy2000/math/blob/3858d2afca72fb97a46a3601a6ae4b1e6e99a812/00-navigation/S74.md) has a genuine but limited
arithmetic contact: the golden and Eisenstein integer rings both reduce
to `F4` modulo two, since `x²−x−1` and `x²+x+1` become the same irreducible
polynomial. Modulo three the Eisenstein quotient is ramified and nonreduced,
whereas the golden quotient is `F9`. Thus a common finite reduction does
not identify the number fields or preserve all changing-modulus dynamics.

`check_quantum_guards.py` independently constructs the three bases, checks
the collision, proves the displayed row-space identity, verifies the small
cross-block factorization and exact seven-guard recovery, and checks the
SU(3)/Eisenstein matrix identities. It is a compact reproducible core, not
a replay of every original research computation.
