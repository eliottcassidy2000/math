# Eight-square signs as an affine code; a five-rectangle obstruction at sixteen

**PROVED:** the rectangle criterion, four-sign recovery in the fixed eight-dimensional XOR basis, and the explicit five-equation obstruction below. **FINITE-EXACT:** all coefficient checks in `verify_norms.py`. The old repository already established the 16 valid Fano orientations and their affine Hamming-code structure. We do not claim those counts or the classical dimension restriction as new. This note gives a common local equation, an explicit small obstruction, and an exact connection to the new paper's signed norm identity. It does not prove the paper's crossing-cost estimate or transfer its factoring algorithm to tournaments.

## Inheritance and the connection being made

The earlier [Paley/Fano/octonion note](../../collatz_mod6_20260922_paley_fano_octonion_design.md), S10–S16, constructs octonion multiplication from a chosen Fano plane and tournament signs. The [Fano code note](../../arithmetic_braids2_20260917_fano_code.md), FC1–FC4, proves that its 16 coherent gauges form an affine `[7,4,3]` Hamming coset. The hostile consists of the other 112 line orientations: the same unsigned Fano support does not imply composition. The new [factoring paper](https://hexagonmath.org/pdf/2610.00600v1), §2 and §4, cites an eight-square identity in its counting engine.

The map here takes a signed XOR multiplication table to a binary rectangle constraint system. It preserves the entire Euclidean norm-composition identity, not just support or an example. Forgetting signs discards the solution to that system. Four multiplication probes restore it at dimension eight, **assuming** it belongs to the valid family. A failed rectangle yields a sparse explicit norm defect. A five-equation contradiction blocks the same signed-XOR construction at dimension sixteen and every higher power of two.

This is not the two-bit unsigned quotient from Session111: its four bits are coordinates of a different affine solution space. The mathematical connection is the need to retain signed constraints, not equality of dimensions or group orders.

## 1. A necessary and sufficient rectangle condition

Let `V=F_2^d`, `q=2^d`, and let a real algebra have basis `(e_a)_(a in V)` with

    e_a e_b = (-1)^F(a,b) e_(a+b),
    F(0,b)=F(a,0)=0.

Addition of indices is XOR. Let `N(x)=sum_a x_a^2`. There is no assumption of associativity. Then `N(xy)=N(x)N(y)` for all real vectors **iff**, for every nonzero `t` and all `a,b`,

    F(a,b)+F(a+t,b+t)+F(a,b+t)+F(a+t,b) = 1   (mod 2).       (R)

**Proof.** The coefficient of `x_a^2 y_b^2` in the norm is one. A mixed monomial `x_a x_c y_b y_e`, with `a!=c`, `b!=e`, can arise only if `a+b=c+e`, equivalently `t=a+c=b+e!=0`. Its coefficient is

    2[(-1)^(F(a,b)+F(c,e)) + (-1)^(F(a,e)+F(c,b))].

This vanishes exactly when (R) holds. Monomials with only one unequal pair do not occur, because XOR is cancellative. These are all polynomial coefficients, proving both implications.

Consequently the valid sign tables form an affine binary code. An incoherent table has a violated rectangle; taking `x=e_a+e_(a+t)` and `y=e_b+e_(b+t)` gives `N(x)N(y)=4` but `N(xy)` equal to 0 or 8. Thus every failure has a two-term witness on each side, not merely a distant numerical defect.

## 2. Reduce to one sign on each projective line

Equation (R) with `a=b=0` gives `F(t,t)=1`. With `t=a+b`, for distinct nonzero `a,b`, it gives `F(a,b)+F(b,a)=1`. With a zero row, it gives `F(a,b)+F(a,a+b)=1`. Therefore for every line `{u,v,w}` with `u<v<w` and `u+v=w`, there is one bit `z_uvw` such that

    F(u,v)=F(v,w)=F(w,u)=z_uvw,

and the three reverse products have exponent `z_uvw+1`. These facts follow from norm composition; they are not extra imposed assumptions. Conversely any line-bit assignment defines a normalized multiplication table, and (R) tests whether it composes the norm.

## 3. Four named signs recover the eight-square table

For `q=8`, name the seven line signs in order

    A=z_123, B=z_145, C=z_167, D=z_246,
    E=z_257, F=z_347, G=z_356.

Here `F` in this display is a single line bit, not the table function. Three instances of (R), with `(t,a,b)` equal to `(1,2,4)`, `(2,1,4)`, `(4,1,2)`, force

    D+E+F+G=1,   B+C+D+E=1,   A+C+D+F=1.

Their solution is

    E=1+B+C+D,
    F=1+A+C+D,
    G=1+A+B+D.                                             (C)

All remaining rectangles follow. One fully finite verification is to substitute (C) into (R) for the 112 rectangles after identifying swaps of each row/column pair; the checker also independently expands every polynomial coefficient of each of the 128 line tables. Alternatively, the table `(A,B,C,D,E,F,G)=(0,0,0,0,1,1,1)` composes the norm by those same coefficient checks, and the following gauge construction produces all 16 possibilities.

For a function `s:V->F_2` with `s(0)=0`, change basis by `e_a'=(-1)^s(a)e_a`. Then

    F_s(a,b)=F_0(a,b)+s(a)+s(b)+s(a+b).

Its added part has zero rectangle sum, or equivalently the diagonal sign change is an isometry transporting multiplication. Thus it preserves composition. The kernel consists exactly of additive functions `s`; these have dimension three. The seven gauge signs therefore have a four-dimensional image. Fixing `s(1)=s(2)=s(4)=0` removes the kernel; choosing `s(3),s(5),s(6),s(7)` produces all 16 tables. This is the inherited Hamming-coset family in the present coordinate convention.

Four probes `(1,2),(1,4),(1,6),(2,4)` reveal `A,B,C,D`, hence all signs by (C). Four binary probes are information-theoretically necessary to distinguish 16 promised-valid tables, including adaptive probes. These probes do **not** certify validity of an arbitrary table: an unobserved line may have been changed. Certification uses the rectangle equations or equivalent syndrome checks.

In dimensions two and four there are respectively one and two normalized norm-composing XOR tables. Dimension one is the real scalar case. These are the real/complex/quaternion/octonion sizes, but the assertion here concerns only this explicit XOR-basis model.

## 4. A five-equation contradiction in dimension sixteen

Use the forced line-bit convention from §2 for the 35 lines on nonzero `F_2^4`. The following five necessary rectangle equations are enough to contradict one another. All additions below are modulo two; subscripts denote an unordered line, with commas preventing decimal ambiguity.

    z_(1,4,5)+z_(1,6,7)+z_(3,4,7)+z_(3,5,6) = 0,
    z_(1,6,7)+z_(1,14,15)+z_(6,9,15)+z_(7,9,14) = 0,
    z_(1,4,5)+z_(1,14,15)+z_(4,10,14)+z_(5,10,15) = 1,
    z_(3,5,6)+z_(3,9,10)+z_(5,10,15)+z_(6,9,15) = 0,
    z_(3,4,7)+z_(3,9,10)+z_(4,10,14)+z_(7,9,14) = 0.

They are (R) at `(t,a,b)=(2,5,4),(8,1,7),(11,4,1),(12,6,5),(13,4,3)`, respectively, after replacing each ordered product by its forced line sign. Every variable on the left occurs exactly twice. Summing gives `0=1`, a contradiction.

Hence **no normalized real signed-XOR algebra in dimension sixteen composes the Euclidean norm**. Every `F_2^d` for `d>=4` contains a four-dimensional subspace; restricting a proposed multiplication table to it yields the same impossible table. This proves the obstruction at every higher power of two. We do not claim this five-equation certificate is minimal, nor that it replaces Hurwitz's theorem for arbitrary bilinear products.

The certificate has ten line variables. The checker independently tries all 1,024 assignments to these ten bits and confirms at least one of the five exact rectangles fails. More importantly, the displayed XOR sum itself is a short proof valid without trusting a search or an exhaustive 35-variable computation.

## 5. What this does and does not transfer to the paper

Eight-square composition is a coherent **signed** tensor identity. It can be encoded by linear parity constraints, recovered under a valid-family promise, and locally falsified. The old tournament/Fano construction therefore provides an actual coefficient-level source for the identity cited by the factoring paper.

It does not supply the paper's network regrouping, term-count bound, shared-bit cost, machine compilation, or `11/12` contraction. Nor does failure of a sixteen-square identity rule out other higher-dimensional crossing decompositions or faster factoring methods. Only the direct normalized signed-XOR norm-composition extension is excluded.

For a particularly small sign-erasure hostile, take `i^2=-1`, `x=1+i`, `y=1-i`: correct multiplication gives `xy=2`, but replacing the structure constant of `i^2` by `+1` gives zero. Coherent gauge changes preserve the norm; arbitrary removal of minus signs does not.

## Reproduction

Run `python3 verify_norms.py` and `python3 -O verify_norms.py`. Both use explicit exception gates. The universe includes all 131 line tables at dimensions 2,4,8, with exact norm-polynomial coefficient expansion, four-probe reconstruction, all 16 gauge choices, and the ten-variable contradiction certificate. No random test or floating arithmetic is a proof dependency.

## 6. Exact transport from the earlier divisor-sum Hamming bank

Session108's [sigma-symmetry bank](../../session108_sigma_symmetry_20261030/README.md) uses primes `(2,5,7,19,29,59,89)` whose `sigma(p)=p+1` squareclasses are the seven nonzero vectors on the prime basis `(2,3,5)`. Here is an explicit coordinate identification, rather than an analogy between two appearances of a Hamming code:

| Line bit | A | B | C | D | E | F | G |
|---|---:|---:|---:|---:|---:|---:|---:|
| Selected arithmetic prime | 19 | 2 | 59 | 29 | 5 | 89 | 7 |
| Squareclass of p+1 | 5 | 3 | 15 | 30 | 6 | 10 | 2 |

For a line orientation `z`, let `m(z)` be the product of the primes in its selected columns. Distinct primes ensure multiplicativity, so the three exponent-parity bits of `sigma(m(z))` are exactly the left sides of the three equations in §3, in the order `(2,3,5)`. Therefore

    this XOR line algebra composes the Euclidean norm
       iff sigma(m(z)) has squareclass 30
       iff sigma(m(z))/30 is an integer square.             (A)

For example, choosing only `D=1` gives `m=29`, `sigma(m)=30`, and a valid eight-square sign table. The gauge difference between two valid tables lies in the zero-syndrome Hamming kernel, precisely the old square-sigma selection condition. For this difference one uses symmetric difference of prime selections, **not ordinary multiplication with repeated prime factors**.

This gives a literal bijection between 16 arithmetic selections in one fixed bank and the 16 signed tables. It preserves the three-bit defect syndrome and its one-line correction rule. It does not identify the integers with octonions, preserve numerical divisor sums or Hamiltonian-path values, or extend (A) to arbitrary prime sets. The numerical squareclass 30 depends on the stated orientation convention; a coherent change of base table shifts the affine code description but not the validity predicate.
