# A Weil-shaped exact formula for the determinant-seven Farey address filter

Let L=[[1,1],[0,1]], R=[[1,0],[1,1]], and evaluate words at (1,1).
Distinct words are distinct primitive positive vectors, by subtractive
Euclid. For the H=23 projective transfer M=[[11,4],[23,9]], a primitive
address vector (x,y) loses content 7 under M iff x+y≡0 mod 7.
Set N_n to be the number of such addresses of exact length n.

## Exact entrance series and sharp square-root error

Reduce positive words to the projective line over F₇, with initial
state 1 and kernel state 6=−1. Reversing the order in which addresses
are read is a bijection on all length-n words. Thus
N_n=(A^n)_(6,1) for the adjacency A=P_L+P_R, where each P is a
permutation matrix. Direct eight-state calculation gives

    det(I−zA)=(1−2z)(1−z)(1+2z²)(1+z+2z³+4z⁴),
    ∑ N_n z^n = (4z³−4z⁴)/((1−2z)(1+z+2z³+4z⁴)).

The checker verifies the determinant from Newton identities for traces,
and the entrance series by recurrence/finite coefficients. Put
s=√17 and t_±=(−1±s)/2. The quartic denominator factors as
(1−t_+z+2z²)(1−t_-z+2z²).
Each reciprocal pole from the quartic has absolute value √2,
since |t_±|<2√2. More explicitly, for n≥0 one has

    |N_n − 2^n/8| ≤ 4·2^(n/2).

To see the actual constant, subtract 1/(8(1−2z)) from the series.
The remaining two numerators over the quadratic factors are
(−17−23s+(153−17s)z)/272 and
(−17+23s+(153+17s)z)/272, respectively. Coefficients of
1/(1−t z+2z²) equal (√2)^n sin((n+1)θ)/sin θ, where
cos θ=t/(2√2). For either t, |t|<31/12 and sin θ>2/5.
The sum of the resulting error constants is less than
(5/2)(46s+306/√2)/272<4. No quartic pole cancels: the original
numerator's roots are 0 and 1. Thus there is not an o(2^(n/2))
error. For instance, the normalized error is a sum of four
nonzero geometric oscillations with distinct unit bases; its
mean square tends to the sum of squared amplitudes.

There is also a literal finite trace/closed-orbit explicit formula.
If B_d is the number of primitive directed edge necklaces of period
d in the labelled two-edge transition graph, then
tr(A^n)=∑_(d|n)d B_d, because a primitive necklace of period d has
exactly d distinct marked starts after repetition. Hence

    det(I−zA)⁻¹=exp(∑_(n≥1) tr(A^n) z^n/n)
              =∏_(d≥1)(1−z^d)^(-B_d).

These identities are formal power series statements. The 2^n main
term and square-root fluctuations for N_n are for address length,
not height or a count of tournament modules.

## Structural criterion and exact failure at a prime

On any finite projective quotient, write S:x↦−1/x and T=S L.
Projectively S²=T³=1 and L=S T, R=S T⁻¹. With P the orthogonal
projection onto the T-fixed subspace, A=S(3P−I). Let U=P S P
act on the range of P. For a unit vector u in that range with Uu=tu
and |t|<1, set v=Su−tu≠0. Then Pv=0 and Sv=(1−t²)u−tv.
The plane spanned by u,v carries the characteristic polynomial
z²−3tz+2. For t=±1 the line gives the eigenvalue ±2. The remaining
orthogonal complement is in ker(P), where A=−S has eigenvalues
±1. These orthogonal components account for the full space by
taking an orthonormal eigenbasis of U (the associated v are mutually
orthogonal). Thus a square-root bound requires the nontrivial
compression eigenvalues t to satisfy |t|≤2√2/3. The Euler product
alone gives no such inequality.

Indeed at p=433 the same projective address graph is *not*
Ramanujan in this sense. The T action on P¹(F_433) gives 146 orbits.
Let D be the diagonal of their sizes, and C_(ij) count elements of
orbit i mapped by S into orbit j. Then U has matrix D^(-1/2) C
D^(-1/2) in the orthonormal orbit basis. The exact symmetric
integer matrix 20000 C−18857 D has two positive eigenvalues
(certified by fraction-free elimination; the signs of the pivots
are the signs of the successive leading-minor ratios, by congruence).
Thus U has at least two eigenvalues greater than 18857/20000,
which is strictly larger than 2√2/3 by squaring integers. One is
the constant eigenvalue 1, which is simple by transitivity of S,T:
equality in norm forces a common fixed vector. The other supplies
a nonconstant A-eigenvalue greater than √2 by the displayed
quadratic. This is a prime-field counterexample to extrapolating
the p=7 result. The failure is visible in the closed-orbit
explicit formula itself: for even lengths 2m the real eigenvalues
contribute nonnegatively, and every nonreal eigenvalue has modulus
√2 by the quadratic blocks. If λ>√2 is the certified eigenvalue,
then tr(A^(2m))−2^(2m) ≥ λ^(2m)−434·2^m. Hence even the global
closed-orbit trace error fails a square-root bound on this prime
quotient.

Uniform expansion is a weaker assertion and does not contradict
this. [Bourgain–Gamburd, Annals 167 (2008)](https://annals.math.princeton.edu/2008/167-2/p07)
give a uniform spectral gap for the fixed symmetric generators
S,T,T⁻¹ modulo primes. Applied to a T-fixed unit vector orthogonal
to constants, the Rayleigh test yields a uniform upper bound
t≤1−δ. For a lower bound, tensor the projective action at odd p
with the sign action modulo 2 (S has sign −1, T has sign +1).
This is a representation of the group modulo 2p and has no fixed
vector by CRT. The arbitrary-modulus expansion theorem of
[Bourgain–Varjú](https://arxiv.org/abs/1006.3365), for the same fixed
generators generating SL₂(Z) projectively, applies to this tensor
representation and gives −t≤1−δ. Thus one can bound nontrivial
eigenvalue moduli strictly away from 2, but not by √2 uniformly.

## The arithmetic-geometric analogy and its boundary

The source discussion is Alpöge's [announcement](https://alpo.ge/remark.pdf),
the [independent Sutherland Lean audit](https://math.mit.edu/~drew/grh/),
and the [Arakelov perspective](https://hexagonmath.org/pdf/2611.00833v1).
The latter gives a slope-method interpretation: a theta-jet evaluation
through special CM cycles, with quantitative equidistribution/zero
estimate and arithmetic Hilbert–Samuel inputs. It explicitly does not
verify the whole argument. It does not identify an F₁ geometry or a
Frobenius functor with those inputs. No such realization is asserted
here either. The strongest exact resemblance in the current theorem is
a closed-orbit determinant plus a genuine square-root bound where a
specific finite-dimensional spectral estimate is provable. The prime
433 example shows why such analytic or geometric estimates cannot be
substituted by a formal F₁-like Euler product.

Replay: run farey-spectrum-check.py with Python 3.
