# Marked joins: exponent saturation and the two-port prime threshold

2026-10-29. **PROVED in integral alternating lattices / FINITE-EXACT controls.**
Same-agent derivation and audit; no external novelty claim. These groups are
not asserted to be realized as Tate–Shafarevich groups.

## J1. An all-copy recurrence retaining the marked port

Let M be a nonsingular integral alternating matrix of positive even size,
v an integral column, A=M^-1, q=Av, d the exponent of coker M, and h the
order of [v] in coker M. Equivalently d is the least denominator of A and
h the least denominator of q. Let S_r have zero diagonal, +1 above it,
and -1 below. For r>=1 put

    K_r = I_r tensor M + S_r tensor vv^T.

For a tournament dominance matrix and v=1 this is the order join of r
copies. For general v it is an integral skew lattice, not necessarily a
tournament. Then

    K_r^-1 = I_r tensor A + S_r tensor qq^T,                 (J1)
    det K_r = (det M)^r,
    K_r^-1(v,...,v) = (q,...,q).

For every r>=2 the exact, r-independent exponent is

    exp(coker K_r) = lcm(d,h^2).                            (J2)

The repeated marked port still has order h. In particular, iterative self
order join, r=2^k, saturates its exponent after the first doubling. Squaring
the exponent at every subsequent doubling is false.

**Proof.** Alternation gives v^TAv=0 and v^TA=-q^T. Multiplying the two
displayed matrices makes the cross terms cancel; the remaining quadratic
term vanishes. The determinant lemma, applied with E=I_r tensor v, gives
det K_r=(det M)^r det(I_r+S_r E^T(I_r tensor A)E)=(det M)^r.
The repeated-port identity follows from q^Tv=0.

The least positive integer making an inverse matrix integral equals its
cokernel exponent. Diagonal blocks of the inverse have denominator d.
Write q=a/h in lowest simultaneous terms. For each prime ell dividing h,
some a_i is not divisible by ell; the corresponding diagonal entry of
qq^T has denominator with ell-valuation 2v_ell(h). All its entries have
denominator dividing h^2. Thus qq^T has least denominator exactly h^2.
It occurs as an off-diagonal block whenever r>=2. This proves (J2).

The mechanism is a zero self-contraction of one alternating port. It is
not a generic theorem that low-rank perturbations cannot accumulate debt.

## J2. Full elementary-primary group recurrence

Assume coker M = (C_p)^(2g), g>=1, and [v] is nonzero, with p any prime.
Write delta_r=r mod 2. Then, as finite alternating paired groups,

    coker K_r ~= (C_p)^[2r(g-1)+2delta_r]
                 directsum (C_(p^2))^[r-delta_r].           (J3)

The isometry is an unmarked classification, not a canonical choice of
generators or a promised transport of every port. Formula (J1), rather
than just (J3), retains the named port and pairing.

**Proof.** The base exponent is p. Over F_p, ker M has dimension 2g and
v is not in im M: otherwise v=Mx+py would imply Av=x+pAy integral.
Conversely Av integral would put v in im M modulo p. Since im M is the
orthogonal complement of ker M under the ordinary dot product,
x -> v^Tx is a nonzero functional on ker M.

The kernel equations for K_r read

    Mx_i + v sum_j (S_r)_ij t_j = 0,   t_i=v^Tx_i.

Reducing each equation modulo im M forces S_r t=0, then Mx_i=0. Each
prescribed t_i leaves 2g-1 free coordinates. Consecutive equations for
S_r imply t_(i+1)=-t_i. The first equation then leaves one free coordinate
for odd r and none for even r, including in characteristic two. Hence

    dim ker(K_r mod p)=r(2g-1)+delta_r.

By J1 the order is p^(2gr), and by J2 every cyclic factor has order p or
p^2 (for r=1 only p). If a,b are their multiplicities, a+b is that nullity
and a+2b=2gr; solving gives (J3). The perfect alternating pairing is the
one induced by K_r^-1. Hyperbolic peeling classifies such pairings by
their cyclic factors, as in the inherited Session96 pairing proof.

For g=1, denote a hyperbolic pair C_m^2 by H_m. The recurrence becomes

    G_0=0, G_1=H_p, G_(r+2) ~= G_r orthogonal-sum H_(p^2), (J4)
    sum_(r>=0) type(G_r) z^r = (1+u z)/(1-w z^2),

where u marks H_p and w marks H_(p^2). Thus

| copies r | group |
|---:|---|
| 1 | C_p^2 |
| 2 | C_(p^2)^2 |
| 3 | C_(p^2)^2 + C_p^2 |
| 4 | C_(p^2)^4 |
| 8 | C_(p^2)^8 |

The +2 here is **two copies**, not endpoint padding by two vertices.
Endpoint padding itself preserves the discriminant pairing, by
[Session96, section 3](../session96_sha_pairing_growth_20261029/tournaments/PAIRING_GROWTH_MODES.md).
Keeping those two meanings separate is essential.

## J3. Two interacting ports are the sharp threshold for new primes

For nonsingular alternating M,N, a rank-one cross block B=uv^T gives

    det [M,B; -B^T,N] = det M det N,

since u^TM^-1u=0 makes the Schur correction zero. This includes all
one-port order joins; their prime support cannot grow.

For a rank-at-most-two factorization B=UV^T with two columns in U,V, put

    alpha = u_1^T M^-1 u_2,  beta = v_1^T N^-1 v_2.

Schur elimination and the two-column determinant lemma give

    det [M,B; -B^T,N] = det M det N (1-alpha beta)^2,
    Pf  [M,B; -B^T,N] = Pf M Pf N (1-alpha beta).          (J5)

For the Pfaffian sign, elimination is determinant-one congruence and the
rank-two update is linear in alpha; its value at alpha=0 fixes the sign.
Equivalently expand the Pfaffian with its two updated columns. The
intermediate factor is rational; the final Pfaffian remains integral.
The required data are **joint pairings of two named ports** on each side.

Here are genuine tournaments proving rank two is sufficient. Use the same
four-vertex block on each side:

    M = [ 0  1  1  1 ]
        [-1  0 -1  1 ]
        [-1  1  0 -1 ]
        [-1 -1  1  0 ].

It has Pfaffian 3 and coker M=C_3^2. Put e=(1,1,1,1),
w=(-2,0,-2,-2), and B=e e^T+w v^T. Every entry of B below is +/-1:

| v | rows of B | alpha=e^TM^-1w | beta=e^TM^-1v | new group |
|---|---|---:|---:|---|
| (0,1,1,0) | (1,-1,-1,1); (1,1,1,1); (1,-1,-1,1); (1,-1,-1,1) | -2/3 | -2/3 | C_5^2 |
| (1,1,1,0) | (-1,-1,-1,1); (1,1,1,1); (-1,-1,-1,1); (-1,-1,-1,1) | -2/3 | 1/3 | C_11^2 |

Their Pfaffians are 9(1-4/9)=5 and 9(1+2/9)=11. Each cross block has
rank exactly two: its explicit factorization bounds rank above, and a
nonzero 2-by-2 minor bounds it below. The paired group of order p^2 is
C_p^2. Thus the appearance of 5 and 11 here has an explicit operation
and preserved pairing, rather than a comparison of unrelated constants.

## Scope and decisive boundary

An arbitrary alternating cokernel is already a classical Sha model; see
[BKLPR, section 3](https://math.mit.edu/~poonen/papers/rst_distribution.pdf).
Our result concerns a structured marked family within that category.
Its bounded exponent and fixed prime support under one-port joins exclude
it from representing every finite alternating group. In particular a
base with exponent p cannot reach H_(p^3) by these joins alone.

The copy parameter r is not the exponent level in an arithmetic Selmer
tower. These constructions provide neither localization kernels nor a
bound on actual Sha. A transition from combinatorics to arithmetic still
needs a source curve/abelian variety, local conditions and certified maps.
