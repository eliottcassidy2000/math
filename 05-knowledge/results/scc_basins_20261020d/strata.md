# Every fixed black-deficit stratum is rational

Status: **PROVED**, conditional only on the explicitly proved inputs linked
below, with finite formal-language checks in `strata.py`.
This classifies the form of all fixed-deficit counting series; it does not
supply every finite strong-centre coefficient or determine black wiring.

The source objects are ordinary tournament SCC words and self-converse
palindromic SCC words. The maps retain vertex size, Hamiltonian-path orbit
mass, centre blue mass, and the number of nonrigid half-components. They
forget individual Hamiltonian words and complement destinations. The
canonical hostile remains `R4 R4`: its black deficit is20, not4.

Dependencies are the all-order rigid alphabet `{1,C3}`, the literal SCC
factorization in [centres, Section2](centers.md), the half-supported
permutation obstruction and regular rotation-orbit connectivity bound in
[the imported deficit-two proof](../deficit_two_centers_20261020/session15-deficit-centers/PREFIX_CONNECTIVITY_CUTOFF.md),
and the general black-deficit cutoff proved in [the general lane](general.md).
The ordinary-mass cutoff below uses only the even-order automorphism
obstruction; it does not assume self-converseness.

## 1. There are finitely many strong tournaments of each fixed orbit mass

Let `h(T)=H(T)/|Aut(T)|`, an odd positive integer. Let G_h(z) count ordinary
isomorphism classes of strong tournaments with this mass.

**Theorem 1.** `G_1=z+z^3`. For every odd h>=3, let L be the least power
of two strictly greater than2h. Then G_h is a polynomial with

    deg G_h <= h(2L-1) < 8h^2.                              (1)

The statement is a computable finite bound, not a claim that all its
coefficients have been enumerated.

**Proof.** The h1 assertion is the inherited rigid-strong classification.
Suppose `n>=2L` and vertex connectivity is at least L+1. A strong tournament
has a strong induced L-vertex subtournament S, by the inherited
vertex-pancyclicity input. Choose any vertex v outside S. At least L/2
vertices of S beat v, or at least L/2 are beaten by v. In the second case
replace the whole tournament by its converse; h and connectivity are
unchanged, and the first case holds.

The complement of S is strong, so it has a Hamiltonian path beginning
at v. Choose a Hamiltonian cycle of S and cut it immediately after each
of the at least L/2 vertices beating v, then append the fixed complement
path. The resulting Hamiltonian words have relative permutations that
are nonidentity powers of an L-cycle on S, hence of even order. Two such
words cannot be in the same Aut(T)-orbit: the unique vertex permutation
taking one full word coordinatewise to the other would be that even-order
permutation, impossible for a tournament automorphism. Thus `h>=L/2>h`,
a contradiction. We have proved `kappa(T)<=L`.

Unlike the black-deficit clique lemma, this argument does not require
disjoint reflected supports or a global anti-automorphism; using `n>=2L`
merely gives a uniform convenient cutoff.

For any Hamiltonian cycle let m be the order of its automorphic rotation
subgroup and q=n/m. The n cyclic cuts give exactly q distinct path orbits,
so `q<=h`; m is odd. If q>=2, the imported rotation-orbit connectivity
bound gives `ceil(m/2)<=kappa<=L`, hence `m<=2L-1` and
`n<=h(2L-1)`. If q=1, the induced regular tournament gives
`ceil(n/3)<=kappa<=L`, so `n<=3L<=h(2L-1)` for h>=3. Orders below2L
already satisfy (1). Finally `L<4h` for odd h>=3, proving the quadratic
estimate. No assumption about the full automorphism group being cyclic
is made.

For example, the bounds for h3 and h5 are45 and155 respectively. They
do **not** prove `G_3=z^5`: that would still require a classification or
finite exclusion of every additional strong mass-three class within the
bound. The known C5 term nevertheless gives `G_3(z)>=z^5` coefficientwise.

## 2. Exact rational SCC words at fixed multiplicative mass

Put `A(z)=1/(1-z-z^3)` and let W_q(z) count all ordinary tournament
classes of mass q, including the empty word when q=1. Unique SCC
decomposition gives

    W_q(z)=sum_(j>=0) A(z)^(j+1)
                 sum_(m_1...m_j=q; each m_i>=3 odd)
                      G_(m_1)(z)...G_(m_j)(z).             (2)

For q=1 the empty product gives W_1=A. For q>1 the sum is finite.
Writing Omega(q) for the number of prime factors of q with multiplicity,
there can be at most Omega(q) nonrigid SCCs. Hence W_q is rational with
denominator dividing

    (1-z-z^3)^(Omega(q)+1).                                 (3)

All numerator ingredients are finite polynomials by Theorem1. The
factorizations are ordered; for example the two orders of masses3 and5
are distinct terms. This is a statement about ordinary classes, since
the left half of a self-converse SCC palindrome is an ordinary class.

The first multiplicative layers are

    W_3=A^2 G_3,
    W_5=A^2 G_5,
    W_9=A^2 G_9+A^3 G_3^2,
    W_15=A^2 G_15+2A^3 G_3 G_5,
    W_27=A^2 G_27+2A^3 G_3 G_9+A^4 G_3^3.                (4)

The square, cube, and fourth powers of A record independent rigid gaps
around one, two, or three nonrigid SCCs. They do not denote additive
black-deficit powers.

## 3. Every fixed black-deficit layer is rational

For a strong self-converse centre Z write `c=h(Z),d=b(Z)`. Include the
empty centre with c=d=1. Let `C_(c,d)(z)` count these centres. The general
cutoff in `general.md` proves that each positive-deficit strong-centre
series D_beta is a polynomial. If c>d, `C_(c,d)` is a coefficientwise
subseries of the finite polynomial D_(c-d); if c=d, the pure-blue
classification leaves only empty, singleton, C3, and C5. Thus every
`C_(c,d)` is a polynomial.

Let B_beta(z) count self-converse ordinary classes with black half-edge
deficit beta. The centre law gives

    B_beta(z)=sum_(q,c,d: q^2c-qd=beta) W_q(z^2) C_(c,d)(z).

For positive beta, the q1 term is `A(z^2)D_beta(z)`. The q>1 corrections
are finite because q divides beta and

    q(q-1)<=beta,
    c<=beta/[q(q-1)],
    1<=d<=c.                                                (5)

Combining this with (2) proves the following all-order result.

**Theorem 2.** For each even beta>=2, B_beta(z) is rational. Its
denominator divides a power of `1-z^2-z^6`, with exponent at most

    1+floor(log_3 Q_beta),
    Q_beta=(1+sqrt(1+4beta))/2.                              (6)

An integer-only formulation chooses the largest j such that
`3^j(3^j-1)<=beta`, then uses exponent j+1. Empty strata have the zero
series. Odd black deficits are impossible. The beta0 series is the
already proved pure-blue rational series.

Equation(6) gives a constant-coefficient eventual recurrence of order at
most six times that exponent. It does not assert a fully reduced
denominator, a minimal recurrence order, or noncancellation at every
complex denominator root. Individual centre coefficients still need
their own finite computations. For beta2, the imported complete theorem
does give the sharper exact result `B_2=z^4/(1-z^2-z^6)`.

This is stronger than a finite catalogue of low layers: every fixed
black deficit has an eventual recurrence, even before all its finite
centre coefficients have been determined.

## 4. Exact thresholds for higher pole powers

Let rho be the positive root of `rho+rho^3=1`, and let lambda=1/rho.
Set

    beta_j=3^j(3^j-1),       j>=1.

**Theorem 3.** B_(beta_j) has a pole of order exactly j+1 at
`z=sqrt(rho)`. No smaller positive deficit can have a pole of order
j+1 there. Its top-order contribution is exactly

    (1+z+z^3) G_3(z^2)^j / (1-z^2-z^6)^(j+1),             (7)

with all other terms having strictly smaller pole order at that point.

**Proof.** A term of pole order j+1 in (2) requires at least j nonrigid
half-components, each of odd mass at least3. Hence q>=3^j, and (5)
requires beta>=q(q-1)>=beta_j. At equality, q=3^j and every one of the
j nonrigid components has mass3. The centre equation becomes
`qc-d=q-1`; since `d<=c`, it forces c=d=1. The allowed centres are
exactly empty, singleton, and C3. This gives (7).

The known strong C5 realizes mass3, so `G_3(rho)>0`. All expansions before
putting them over a common denominator have nonnegative coefficients;
their top poles at a positive point cannot cancel. Formula(6) supplies
the matching upper bound. This proves both exactness and the minimal
deficit threshold.

The first thresholds are:

| deficit | pole order | polynomial prefactor on either parity subsequence |
|---:|---:|---|
| 6 | 2 | linear |
| 72 | 3 | quadratic |
| 702 | 4 | cubic |
| 6480 | 5 | quartic |

These give a precise answer to how squares, cubes, and fourth powers
enter the black-deficit problem. They arise when the half contains one,
two, or three nonrigid mass-three SCCs. They are not `D_2^2,D_2^3,D_2^4`,
whose atom-count interpretation is different.

More precisely, writing `b_(beta,n)=[z^n]B_beta(z)`, singularity expansion
of(7) gives positive constants

    b_(beta_j,2n)   ~ K_j n^j lambda^n,
    b_(beta_j,2n+1) ~ (1+rho) K_j n^j lambda^n,

    K_j=G_3(rho)^j /
        [j! (rho(1+3rho^2))^(j+1)].                        (8)

The cubic `1-t-t^3` has a unique root of smallest modulus, the positive
rho; this follows from the triangle inequality and the coprime exponent
set `{1,3}`. Its zero is simple. On splitting the series into even and
odd parts in `t=z^2`, the respective leading numerators are `G_3(t)^j`
and `(1+t)G_3(t)^j`. Every other term has lower pole order, proving (8).
Thus even the limiting odd/even ratio `1+rho` is determined without
assuming that C5 is the only strong mass-three atom. The finite polynomial
G_3 controls the common amplitude, not this ratio or the pole order.

These asymptotics concern a fixed black deficit. They do not determine
the user’s distinct matching sequence M_k/k!, nor the unrestricted
tournament count, nor the connectivity of nodes at that deficit.

## 5. Exact formal controls and scope

`strata.py` compares two different enumerations on the genuine restricted
strong alphabet `{1,C3,R4,C5,Paley7}`: direct ordered-word dynamic
programming by vertex size and multiplicative mass, and formula(2) with
strong-mass polynomials `G_3=z^5,G_5=z^4,G_9=z^7` for this restricted
alphabet only. It checks all sizes through30 and masses through81.
It also checks the first pole thresholds arithmetically without floating
point and verifies the symbolic first multiplicative layers.

The restricted choices for G_m are finite controls, not an assertion
that the unrestricted strong-mass alphabets have been classified. The
proofs of rationality and (8) retain the unknown finite G_m polynomials.
