# Recovering 5 and 11 without conflating their arithmetic roles

2026-10-29. **PROVED elementary deductions / CITED descent inputs /
FINITE-EXACT computations / OPEN general finiteness.**

## Inheritance and scope

Source checkpoint: the clean incoming Session96 Sha packet at
306f2d571f40, followed by audit commit 3e4fae1986 on origin/main. The worktree
was then rebased before any mathematics was filed. The shared dirty main
checkout was left untouched.

Closest proved mechanism: [Session96 marked discriminant pairing and
growth](../session96_sha_pairing_growth_20261029/tournaments/PAIRING_GROWTH_MODES.md).
Hostile: its C_3^2 self-join becomes C_9^2, so squared order need not be a
direct square of the group. Corrected near miss: global converse is not
path-native complement; no claim of native closure is inherited. The
least-used useful sidecar in this continuation is the order of the named
port, together with its alternating contraction with a second port.

Anchor: the requested Sha finiteness problem, through exact paired
recurrences and a real arithmetic calibration. Niche: the October 13
Monsky receipt and October 14 golden 5/11 split. Wildcard: whether equal
thin F1 counts mask different ambient forms. These are operationally
underconnected sources; no claim of exhaustive rediscovery counts.

| Live concept | New result's effect / required data |
|---|---|
| Square/twice-square order | Pairing and characteristic element still matter; order alone is insufficient |
| One marked port | Exact exponent lcm(d,h^2), not repeated exponent squaring |
| Joint ports | Two contractions permit new primes 5 and 11 |
| F1 coordinate skeleton | First jet distinguishes quadratic and symplectic extensions |
| Arithmetic level clock | Copy number is not a p-power descent level |
| Golden 5/11 | Ramification/splitting is a separate predicate from the Sha defect |

## R1. Sources recovered before extension

1. [Arithmetic seams, September 21](../arithmetic_seams_20260921_operations.md),
   sections 1--2, already explains doubling/squaring through equal-parent
   deletion and the transported operations under a -> a+1. It distinguishes
   graph complement from inverse operations. We inherit that distinction.
2. [Four-state tournament recursion, October 4](../tournament_recursive_four_state_20261004.md)
   and its sources THM-291, `mode-b-multilinear-recursion`, and THM-447,
   `skew-sylvester-doubling`, distinguish endpoint padding, order join,
   substitution, and skew doubling. Their vertex-size formulas do not
   determine what happens to a paired discriminant group.
3. [Finite Selmer receipts, October 13](../bsd_selmer_growth_20261013a.md),
   section 2, supplies the exact Monsky matrices used below. Its applications
   of an owner-supplied full-BSD theorem are conditional and are **not**
   dependencies here. The [Tunnell terminals](../bsd_tunnell_terminals_20261013c.md)
   also retain 5 as a zero-coefficient control and 627=3*11*19 as a distinct
   nonzero terminal. No conclusion about 11 is borrowed from that product.
4. [Golden eleven, October 14](../strong_balance_pisano_golden_eleven_20261014.md),
   section 2 and MISTAKE-614, distinguishes the ramified prime 5, split
   prime 11, and inert primes in Q(sqrt(5)). Its icosahedral/eleven-square
   analogies do not identify the corresponding objects. We reuse only
   elementary period facts checked again here.
5. [Session59 signed flags](../session59_f1_geometry_and_information_20261025/signed-flags/THEOREM.md)
   retains whole counting polynomials before q=1. The present first-jet
   comparison is a new application of that existing research move.

Targeted searches included Shafarevich, Cassels--Tate, F1, doubling/squaring,
the growth modes, their correction entries, and BKLPR/random alternating
matrices. The incoming Session96 results were read as current truth rather
than rediscovered from older logs. No claim of literature-wide priority.

## R2. An actual arithmetic calibration at 5 and 11

For E_n: y^2=x^3-n^2x with n a positive odd squarefree integer, use the
classical Monsky interface from [Watkins et al., section 4.1](https://dms.umontreal.ca/~andrew/PDF/RankQuadTwist.pdf).
For a single prime n, the pure full 2-Selmer dimension is the nullity over
F_2 of

    [[ [2/n], [2/n] ], [ [2/n], [-2/n] ]],

where [a/n] is zero or one according as a is a quadratic residue or
nonresidue. This **cited arithmetic theorem** incorporates local solubility;
our computation of a binary kernel is not a replacement proof of it.

The resulting matrices are

    M_5  = [[1,1],[1,1]],  nullity 1,
    M_11 = [[1,1],[1,0]],  nullity 0.

The point (25/4,75/8) on E_5 is nontorsion by Nagell--Lutz, since its
coordinates are nonintegral on an integral short Weierstrass equation;
see [Sutherland, MIT 18.782, Theorem 24.21](https://www.ocw.mit.edu/courses/18-782-introduction-to-arithmetic-geometry-fall-2013/f37bbbd7fb8c7ee4e077cc5aa1d2a220_MIT18_782F13_lec24.pdf).
Its curve equation is checked with exact fractions.

Finite Kummer exactness ([Poonen's descent notes](https://math.mit.edu/~poonen/f01/weakmw.pdf)) gives

    dim Sel_2(E_n) = rank E_n(Q) + 2 + dim Sha(E_n)[2].

Thus rank E_5(Q)=1, rank E_11(Q)=0, and both Sha[2] groups vanish.
Every nonzero element of a 2-primary torsion group has a nonzero multiple
of order two. Therefore

    Sha(E_5/Q)[2^infinity] = Sha(E_11/Q)[2^infinity] = 0.  (R1)

These are classical examples reverified as calibration. No general Sha
finiteness conjecture or supplied low-corank theorem is used in (R1).
This computation alone says nothing about their odd-primary parts.

The number 5 here is a curve parameter; it is not the coefficient prime
of the descent, which is 2. Likewise C_5^2 in J3 is a different object:
it is the discriminant of a specified tournament lattice, not Sha(E_5).

## R3. Golden clocks: the old 5/11 distinction survives with its own type

Let F=[[0,1],[1,1]]. Its characteristic polynomial is t^2-t-1, of
discriminant 5. Over F_5 it is (t-3)^2, and F=3I+N with N^2=0,N!=0.
The scalar has order four; the nontrivial unipotent part has order five.
Hence F has order 20. Over F_11 the distinct roots are 4 and 8, of orders
5 and 10, so F has order 10. These recover the repository's Pisano clocks.

This is a ramified-versus-split distinction of a particular quadratic
algebra. Both primes are odd; neither plays the exceptional role of 2 in
an antisymmetric Q/Z pairing's self-pairing. For x of odd order,
2b(x,x)=0 and odd-order torsion force b(x,x)=0. A nonsingular odd-primary
part therefore has hyperbolic paired factors and square order.

## Connection contracts and stopping boundary

| Source -> target | Map / preserved fact | Lost information / restoration |
|---|---|---|
| Marked lattice -> repeated join | J1 block formula / inverse and port order | Determinant drops exponents; retain A and q |
| Two marked lattices -> attached lattice | B=UV^T / J5 Pfaffian | Marginal ports drop contractions alpha,beta |
| Quadratic or symplectic geometry -> q=1 count | Polynomial evaluation / 2^n | Ambient form lost; first jet separates counts, Q or omega restores predicate |
| E_5,E_11 -> binary kernel | Cited Monsky theorem / full 2-Selmer dimension | Kernel without source equation/local theorem is not a Selmer group |
| Golden matrix -> modular clock | Exact reduction / matrix order | Splitting type does not encode Sha or tournament composition |

## Incoming comparison: a pairing repair retains the old clock

Incoming commit f41b1d1bfd supplies a decisive control from our previous
[Lucas cokernel tower](../stabilizer_necklace_kernels_20261028/COKERNELS.md).
Its depth-two group G has Smith factors 1, seven 2s, and 38. The order
4864=2^8*19 has an odd 19-valuation, so G cannot carry a perfect
antisymmetric Q/Z pairing. A clocked cokernel alone is not a Sha carrier.

There is a **classical**, explicit paired construction which preserves G
as a named half rather than pretending the old pairing exists. Define

    H(G)=G + G^dual,
    b((x,chi),(y,psi))=psi(x)-chi(y).

This is perfect alternating, has order |G|^2, and contains G as a specified
Lagrangian. An automorphism t of G extends canonically to
(t,(t^-1)^dual). Thus the old cyclic clock survives along with its dual
phase. For any nonsingular integral presentation C, the block matrix

    [[0,C],[-C^T,0]]

presents this hyperbolic paired group, up to the harmless overall sign of
the pairing. In the nine-coordinate circulant example, simultaneous
cyclic translation on both halves preserves the matrix. The resulting
Smith factors are two 1s, fourteen 2s, and two 38s, with order 4864^2.

This gives a precise doubling/squaring operation: double the presentation
dimension, square its discriminant order, and **retain the original half
and dual clock**. Unlike self order join it is a direct hyperbolization.
It repairs the finite pairing obstruction but supplies no arithmetic
realization or Sha finiteness theorem; it also produces no nonalternating
characteristic defect. This comparison is classical algebra applied to
the specific incoming object, not a new theorem about hyperbolization.

The other incoming recurrence packets concern charged observation arity
and source-conditioned affine carries. Their lesson here is consistent:
arity, copy count, Fourier phase and p-power torsion depth are different
clocks. None of those new results is used as an arithmetic premise.

The used protocol card is `Canonicalize mechanics; retain the predicate`.
Its counterindication is visible here: matching thin counts do not preserve
quadratic isotropy. The niche recovery card motivated checking the older
arithmetic receipt before proposing another abstract finite-group model.

The general finiteness conjecture remains OPEN. Poonen--Stoll's [square or
twice-square theorem](https://math.mit.edu/~poonen/papers/sha.pdf), Corollary 9,
concerns a finite Sha (or finite nondivisible quotient) under principal
polarization. A nonalternating pairing can have square order. This finite
shape theorem cannot bound a divisible radical or the prime/exponent
support of the actual arithmetic group. The stopping criterion and hostile
finite-prefix examples in [Session96](../session96_sha_pairing_growth_20261029/selmer/SELMER_TOWER_FINITE_CERTIFICATES.md)
remain the operative boundary. This session constructs no arbitrary
native-complement carrier and makes no LRC/Collatz claim.
