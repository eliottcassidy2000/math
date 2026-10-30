# Scaled counts, lost Fourier phases, and the signed crossing carrier

**PROVED elementary deductions / FINITE-EXACT named finite universes /
CITED external identities. Same-agent review; no global novelty claim.**
Notation and each consumer are fixed below. A residue, a zero decision,
a witness, an isomorphism class and a weighted contribution are distinct.

Concurrent precedence: Session112 published the same general effective-
modulus theorem, its Burnside application and a dyadic cutoff at
`fea3b07542` while this packet was in progress. C1-C4 are credited
rederivations/applications, with a slightly smaller elementary constant
in C2, not distinct discoveries. Its later norm-rectangle and defect
results are integrated in S4. See RECOVERY for pinned paths and scope.

## C1. The exact information in scaled residue observations

Let D>0 be known, let m_1,...,m_s be positive moduli, and observe

    r_i = D c mod m_i,      0 <= r_i < m_i.

Put L=lcm(m_i) and

    Q=lcm_i(m_i/gcd(m_i,D))=L/gcd(L,D).                 (1)

**Theorem.** Two integers c,c' give the same observations iff
c=c' mod Q. Thus every fibre is one residue class modulo Q.
The observations determine every c in [0,B] iff Q>B. They distinguish
zero from every positive c<=B under the same condition. On[-B,B],
zero detection needs Q>B, while full reconstruction needs Q>2B.

**Proof.** Equality of the i-th residues is m_i | D(c-c'), equivalent
to m_i/gcd(m_i,D) | c-c'. Intersect the divisibility conditions.
The equality in(1) follows prime by prime:

    max_i max(v_l(m_i)-v_l(D),0)
      =max(max_i v_l(m_i)-v_l(D),0).

The interval statements follow from the smallest nonzero period Q.
At the failure boundary the pair 0,Q is an explicit collision. QED.

**Constructive decoder using only small residues.** Put d_i=D mod m_i
and g_i=gcd(d_i,m_i). The necessary image condition is g_i | r_i.
When q_i=m_i/g_i>1, the factor d_i/g_i is invertible modulo q_i, so

    c = (r_i/g_i)(d_i/g_i)^(-1) mod q_i.                (2)

Generalized CRT merges these compatible congruences. Neither a
factorization of D nor its expanded integer value is needed if its
residues can be evaluated. A multiplicative expression for D is
therefore a sufficient sidecar. Without D, arbitrary division of r_i
is invalid.

For a vector of counts, the same statement applies coordinatewise.
Residues are closed under addition and multiplication, but a bounded
integer reconstruction must update its height: B1+B2 or B1*B2.
Closed modular arithmetic does not preserve an obsolete height budget.

## C2. A scale-free zero certificate for a dyadic multiplier

Suppose D=2^E, E>=0, and 0<=c<2^n. The odd moduli

    3,5,7,...,2n+5                                       (3)

suffice to decide c=0 from the residues of D c, independently of E.
There are n+2 such moduli. They also determine c uniquely.

**Proof.** Set M=n+3 and Q=lcm(1,3,...,2M-1). The odd part of
binom(2M,M) divides Q: every prime power occurring in that binomial
coefficient is at most2M. For example, the usual factorial-valuation
formula bounds each valuation by floor(log_l(2M)).

The factorial identity v2(j!)=j-s_2(j) gives
v2(binom(2M,M))=s_2(M). Consequently

    Q >= binom(2M,M)/2^s_2(M)
      >= 4^M/[(2M+1)2M] > 2^n.                         (4)

For the last inequality use 2^(n+6)>(2n+7)(2n+6):
it holds at n=0 and the right-hand side grows by a factor less than 2
at each increment. Each modulus in(3) is coprime to D, so C1 applies.
This proof needs no prime-number theorem or factorization oracle. QED.

**Exact limit.** For a general known D, odd moduli need not suffice
at this height. If D=lcm(1,...,H), every observation modulo m<=H
is zero for both c=0 and c=1. No theorem can remove arbitrary content
without changing the modulus bank or retaining more information.

The Hexagon prose describes a power-of-two scaling. Its actual generic
Lean grid multiplier is

    t=ceil(lambda/3),
    R=binom(G,2)^2,   A=G^4-R,
    mult(G,lambda)=(t+1)^A 2^(t R).                     (5)

Deep levels multiply such factors. For G=2,lambda=4, this is
3^15*2^2=57395628. Thus C2 is valid for the dyadic subkernel, but is
not an unconditional optimization of the full submitted implementation.
C1 is the correct general replacement criterion. See SOURCE_AUDIT.md
for the exact files and the distinction between this generic parameter
test and the large schedules chosen by the top-level machine.

## C3. Burnside division is the same content problem

For a finite group G acting on a finite set X,

    S=sum_(g in G) |Fix(g)|=|G| a,   a=|X/G|.

Computing S modulo the m_i determines a modulo
L/gcd(L,|G|), not automatically modulo L. This is C1 with D=|G|.
In a single prime-power channel, a mod l^b can always be recovered
from S mod l^(b+v_l(|G|)); the extra valuation is necessary for
arbitrary integer a. If m is coprime to |G|, ordinary modular inversion
is legal.

**Literal tournament application.** Fix n unnamed old vertices and r
individually named attachment ports, whose internal arcs are prescribed.
The Q=2^r colors record each old vertex's complete incidence vector
to those ports. Permutations act on old vertices only. If lambda is
an odd-parts partition of n, define

    z_lambda=product_j j^(m_j) m_j!,
    e(lambda)=sum_i(lambda_i-1)/2 + sum_(i<j)gcd(lambda_i,lambda_j).

The classical cycle-index formula is

    A_n(Q)=sum_(lambda odd partition of n)
                   2^e(lambda) Q^(length lambda)/z_lambda.          (6)

An even permutation cycle would reverse some pair and fixes no
tournament. With odd cycles, internal arc orbits number(lambda_i-1)/2
and cross-cycle arc orbits number gcd(lambda_i,lambda_j). A fixed
coloring is constant on each cycle. This proves(6).

Our modular path computes the integer numerator n! A_n(Q), applies
(2), and reconstructs under the labelled height
2^(binom(n,2)+rn). Independent rational summation of(6) agrees for
n=0..20 and r=0,1,2. Uncolored n<=12 also matches the inherited
A000568 prefix. This is a reproducibility check of an existing
formula, not a new sequence term, compact family classifier or
complexity improvement.

**81/27 hostile.** Native groups of orders81 and 27 have1620 and 1728
regular binary masks, hence20 and 64 free orbits. Both raw counts are
zero mod3. Keeping mod243 for the first and mod81 for the second
recovers respectively2 and 1, their orbit counts mod3. The same
cyclic-triple set does not supply the missing group-order division.

## C4. Prefix search preserves witnesses only with the exact fibre

For any finite nonempty set of fixed-length binary witnesses, counts
of prefix fibres determine its lexicographically least member: test
the zero-child fibre, choose it if nonempty, otherwise choose the
one-child. The current prefix always has a completion.

The paper applies this classical self-reduction to divisors d>=2 of N.
Each such divisor has exactly one partner e with de=N, and each
deterministic multiplication computation has one internal filling.
This uniqueness is why the network count has the required multiplicity.

The verifier independently implements the reduction for 2<=N<=500,
using C2 on synthetically scaled exact divisor counts. It does not
implement the paper's fast counting engine or benchmark factoring.
For tournament orbits, a prefix of labelled arcs need not be invariant;
one cannot divide each prefix count by a fixed group order without
handling its actual stabilizers or orbit-saturated membership.

## F1. Exact Fourier failure on odd-base plus-one rings

Let b>=3 be ANY odd integer, a>=1, m=2^a, T=2m, and

    N=b^m+1=2u,   u odd.

Consider F:(Z/N)^T -> (Z/N)^T,

    (Fx)_j=sum_(i=0)^(T-1) b^(ij) x_i.                 (7)

**Theorem.** F has kernel and cokernel isomorphic to
(Z/2)^(T-1). Every nonempty fibre has exactly 2^(T-1) elements.
The odd CRT component is fully invertible; the modulo-two component
retains only the sum of the input parities.

**Proof.** Since b^m=-1 mod N and m is a power of two, the
geometric product

    sum_(i=0)^(T-1) X^i = product_(j=0)^(a) (1+X^(2^j))

shows that sum_i b^(ij)=0 mod N for 0<j<T: choose the factor that
raises b^j to an odd multiple of m. Thus F(b)F(b^-1)=T I.
Modulo u the integer T is a unit, proving invertibility there.

Modulo2 every entry of F is 1. This rank-one matrix has kernel the
even-weight binary words, dimension T-1, and cokernel dimension T-1.
CRT proves the assertion. Explicit kernel generators are
u(e_i+e_(T-1)) for 0<=i<T-1. QED.

**Optimal reconstruction.** F supplies x modulo u and sum_i x_i mod2.
The parities of any T-1 named input coordinates determine the last
parity and hence x modulo N by CRT. No auxiliary datum with fewer
than 2^(T-1) values can separate a full fibre; T-1 bits are optimal
in this universal reconstruction sense.

This also applies to composite odd b. Primality is not a hypothesis.
The even exponent is essential to the stated two-primary structure:
we used N=2 mod4. The original paper's base2 modulus2^m+1 is odd,
so T is already a unit and this lost-parity component is absent.

| b | m | N | T | missing bits |
|---|---:|---:|---:|---:|
| 3 | 2 | 10 | 4 | 3 |
| 5 | 2 | 26 | 4 | 3 |
| 7 | 2 | 50 | 4 | 3 |
| 3 | 4 | 82 | 8 | 7 |
| 3 | 8 | 6562 | 16 | 15 |

## F2. Retaining the clock gives a nilpotent algebra

Represent the parity word as a polynomial in X modulo X^T-1.
In characteristic two,

    F2[X]/(X^T-1) = F2[epsilon]/(epsilon^T), X=1+epsilon.             (8)

The Fourier shadow keeps only the constant epsilon coefficient,
the augmentation sum_i x_i. Its kernel is the nilpotent ideal
(epsilon), with dimension T-1. Cyclic convolution is multiplication
in(8), so the discarded information is an algebraic ideal, not random
rounding error.

A one-step cyclic shift multiplies by 1+epsilon, which has EXACT
order T: its T-th power is 1, and its T/2-th power is
1+epsilon^(T/2), not 1. The nonzero hidden parity word 1+X has a full T-cycle
for T>=4. Its lift u*(1+X) modulo N has zero odd component and
zero Fourier output at every shift. Repeating read-only
Fourier observations can therefore never recover its phase.

One named coefficient queried after every shift distinguishes every
parity word. In the presence of this future consumer, all T-1 missing
bits are necessary. The quotient is nevertheless closed for consumers
that ask only for total parity under convolution. This is a precise
consumer-dependent boundary.

Frobenius sends epsilon^i to epsilon^(2i). Its j-th iterate on the
whole parity algebra has image dimension max(T/2^j,1), with the
quotient understood integrally for j<=log2(T); thereafter the image
is the constants. On(epsilon) its image dimension is
max(T/2^j-1,0). It halves the supplied shift clock at each nonterminal
step while destroying coefficients; a clock alone does not recover
those coefficients.

Doubling T gives the exact memory recurrences

    K(2T)=2K(T)+1,      H(2T)=2 H(T)^2,                 (9)

where K=T-1 is the missing-bit count and H=2^K the number of states
in a transform fibre. This supplies a literal doubling/squaring
connection and a3->7->15 chain. It is not an identification with the
different seven-prime Hamming code.

## S1. The paper's crossing identity selects the old Hamming code

The external Cross.lean gives an explicit octonion sign
sigma(x,y) on x,y in F2^3. It uses the XOR product index and defines

    R=[x+y=v+u] sigma(x,y)sigma(v,u),
    C=[x+u=v+y] sigma(x,u)sigma(v,y),
    R+C=2[x=v][y=u].                                   (10)

All signs, indices and four named ports are retained. The verifier
replays all 4096 quadruples directly from the polynomial sign formula.

The seven nonzero points carry the seven XOR Fano lines. Flip the
orientation of any subset e of these lines; this gives128 multiplication
sign tables and 128 regular tournaments. This family, and its16 good
versus112 defective members, are INHERITED from the September Fano
work. The new comparison uses exactly the external crossing formula.

Label each line by its nonzero normal vector a in F2^3 and put

    syndrome(e)=sum_a e_a a.                            (11)

**FINITE-EXACT equivalence on all 128 line tables:** the complete
crossing identity holds iff syndrome(e)=0. There are 16 such masks.
Every other mask fails on exactly 96 of the4096 four-port inputs.
These are the inherited[7,4,3] Hamming-code masks.

**Mechanism for sufficiency.** Replacing basis vector e_x by
(-1)^f(x)e_x changes sigma(x,y) by
(-1)^(f(x)+f(y)+f(x+y)). In each supported R or C term the internal
label occurs twice, so both terms gain the same external factor.
On the nonzero right-hand side that factor is 1. Thus every basis-sign
gauge preserves(10).

On line orientations this gauge is the Fano point-line incidence map.
Its binary rank is 4, and its image is ker(syndrome): the three normals
through any point sum to zero. Rank and all remaining masks are checked
independently in the finite enumeration. The converse classification
here is explicitly the complete finite check, not an unproved
classification of arbitrary larger multiplication algebras.

## S2. The smallest carrier depends on the requested output

For the operation of flipping any future set of NAMED lines, and the
question whether(10) still holds, the smallest deterministic carrier
has exactly EIGHT states: the syndrome in F2^3.

Equal syndromes give the same answer after every future flip.
If s!=t, flip any mask of syndrome s; the first state becomes valid
and the other remains invalid. Hence no two syndromes can be merged.
Literal converse flips all seven lines and has syndrome zero, so it
is closed on this carrier. Keeping only today's valid/invalid flag
is insufficient: masks1 and 2 are both invalid, but flipping line1
makes just the first valid.

**Re-marking closure.** Let A be any element of GL(3,2), acting on the
seven nonzero point labels. A line normal transforms by A^(-transpose).
Relative to the fixed reference sign table, the action on flip masks is

    e -> P_A e + b_A,   syndrome(b_A)=0.

The translation b_A is necessary: re-marking a reference orientation
does not act by a plain permutation of its orientation bits. The
transported reference table still satisfies (10), since A preserves
XOR and equality; S1 therefore puts b_A in the good code. Consequently

    syndrome(e) -> A^(-transpose) syndrome(e).

All 168 linear re-markings act faithfully on the eight syndromes.
Together with arbitrary named line flips they induce AGL(3,2), of
order 8*168=1344. The eight-state lower bound already follows from
line flips, so adding these re-markings does not enlarge the minimum.
Exact enumeration gives a reference-table stabilizer of order 21 and
an eight-element GL orbit; adjoining converse reaches all 16 good
tables. This concerns the supplied XOR Fano geometry. It does not
assert closure under arbitrary tournament substitutions or arbitrary
permutations that forget that geometry.

If the consumer requests the ACTUAL R branch coefficient at named
ports, all 128 line masks must be retained. For a line containing
x,y,x+y, choose u=0,v=x+y. Then R=sigma(x,y), recovering that line's
orientation. Seven such queries recover every mask. The eight-state
legality carrier is therefore not a complete weighted crossing carrier.
Under a promise that only the 16 good tables can occur, four bits suffice,
as shown in concurrent Session112 NORM_SIGNS. Our 128-state statement
allows arbitrary future named flips, which can leave that valid locus.

## S3. Cancellation is essential, even on a zero target

At x=y=0 and u=v=1, the external formula has R=-1,C=1 and target0.
Replacing contributions by absolute values gives2, not 0. The four-port
identity is exact over signed rings; it is not a positive probability
decomposition. An externally nonnegative integer witness count can
still have signed internal terms. Applying the compression to LRC or
Collatz needs the actual witness count and its height, not an average
of absolute pieces. No uniform target-entry or convergence theorem
follows from(10).

## S4. Crossing error and norm error have the same squared energy

This comparison integrates Session112's concurrent NORM_SIGNS rectangle
criterion and NORM_DEFECT energy formula with the literal Cross.lean
port convention. Their 24-rectangle and 384-energy conclusions are
credited inputs/recoveries, not separately claimed discoveries here.

Let V=F2^d and let ANY sign table sigma:V*V->{-1,1} define
e_x e_y=sigma(x,y)e_(x+y). Normalization, alternativity and a Fano-line
promise are not needed. Write N for the Euclidean square norm, define R
and C as in (10), and define the crossing error

    E(x,y,u,v)=R+C-2[x=v][y=u].

For random vectors X,Y with all coordinates mutually independent,
centered, and of variance one,

    sum_(x,y,u,v) E(x,y,u,v)^2
       = expectation [(N(XY)-N(X)N(Y))^2].              (12)

Only second moments are required. **Proof.** The diagonal coefficients
cancel because sigma(x,y)^2=1. Every remaining monomial has the form
X_x X_v Y_y Y_u with x!=v, y!=u and x+v=y+u=t!=0. Its coefficient is

    2[sigma(x,y)sigma(v,u)+sigma(x,u)sigma(v,y)] = 2E.

Each unordered pair of row indices and unordered pair of column indices
gives four ordered port assignments, all with the same E. Their total
contribution 4E^2 is exactly the coefficient squared. Distinct monomials
are orthogonal in expectation: any distinct index pair leaves at least
one centered coordinate to the first power. Each monomial's squared
expectation is one. Summing proves (12). QED.

The bijection from rectangles to ordered ports is explicit:

    (t,{a,a+t},{b,b+t}) -> (x,y,u,v)=(a,b,b+t,a+t),

with the four choices of representatives a,b. Thus the concurrent 24
failed rectangles give 96 failed ordered port assignments. At each
failure E=+2 or -2, giving energy 96*4=384. In the supplied eight-
dimensional Fano family, validity is equivalent to zero energy; every
invalid table has energy 384. Individual signed coefficients still
require more than the syndrome, as S2 shows.

The verifier independently expands norm polynomials for every Fano mask
and compares their coefficient energy with the four-port sum. Additional
arbitrary sign tables at dimensions 1,2,4,8,16 check the unrestricted
identity. Equation (12) is a proof, not extrapolation from those tests.
It describes error of an identity; it does not make a defective crossing
replacement legal inside an exact counting engine.

## R1. Recovered 12/13 clocks and the distinction from ramification

This is an inherited mechanism, replayed independently here, from Q1-Q3
of [TRACE_CARRIER](../cursor_cm_phase_20261028/TRACE_CARRIER.md).
For O=Z[u], u^2-a*u-1=0, Tr(r+s*u)=2r+a*s, the primitive pair
(Tr x,Tr(u*x)) has determinant a^2+4 and Smith factors
(gcd(2,a),(a^2+4)/gcd(2,a)).

For a=5 modulo 13, u^12=1: the entire sampled trace history has 13
values on 169 ring states, while the primitive pair recovers all 169.
For a=2 modulo 11, u^12=-1: sampling alternates and has 11 values on
121 states; the primitive pair recovers all 121. For a=5 modulo 29,
the discriminant vanishes, and even the primitive pair has only 29
values on 841 states. Another power clock cannot fix that trace radical.

The earlier [Ramanujan/Zagier recovery](../cursor_cm_phase_20261028/RAMANUJAN.md)
separates the reciprocal two-sheet quotient of the Fricke function from
this finite-ring trace loss. We inherit that distinction; no new modular
form or singular-modulus identity is claimed here. The current Fourier
example extends the same operation-sensitive question to a nilpotent
ideal of dimension T-1, rather than identifying its rings or its clock
with the quadratic example.
