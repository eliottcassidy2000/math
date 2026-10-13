# Finite Selmer growth certificates and a tournament obstruction

2026-10-13. **PROVED elementary interfaces + CITED arithmetic inputs +
FINITE-EXACT computations.** The full BSD consequences below are conditional
applications of the owner's supplied low-corank BSD theorem, accepted as an
input for this session. No independent proof of general BSD is claimed.

## Inheritance and the exact target

The closest proved repository mechanism is
[THM-3646, Berggren fixed-107 Mordell rank exactly three](../../01-canon/theorems/THM-3646-berggren-fixed107-mordell-rank-exactly-three.md):
an authenticated Selmer upper bound and independent rational points together
give an exact rank. Its corrected near miss is essential: a two-point
subgroup can be 3-saturated and still miss a third independent direction.
The rank-three curve there lies outside the low-corank input used here.
[THM-3643, real quadratic class number one](../../01-canon/theorems/THM-3643-berggren-fixed107-real-quadratic-class-number-one.md)
is a useful sidecar, not itself a Selmer or rank computation.

The new target is a **finite receipt for the low-corank hypothesis**, rather
than replacing it by the much weaker assertion that few rational points
have been found. A second lane connects actual quadratic-residue relations
to the repository's tournament switching operation. The live concepts are
finite growth, divisible tails, full local conditions, switching gauges,
and retained arithmetic labels.

The supplied *Exact Birch–Swinnerton-Dyer Formula from Low Selmer Corank*,
Theorem 1.1 (`C:/Users/Eliott/Downloads/exact-bsd-low-selmer-corank.pdf`),
states: if the full q-infinity Selmer corank at some prime is at most one,
then algebraic and analytic ranks equal that corank, Sha is finite, and
the complete BSD formula holds. This is **CITED, supplied accepted input**.
The full formula uses the total real Neron period, the full Mordell–Weil
regulator, torsion squared, and all finite Tamagawa factors. No numerical
leading term is evaluated in this packet.

Its proof mechanism distinguishes an integral comparison from a unit:
horizontal divisors must also be removed. The related supplied two-primary
paper retains marked local conditions and full lattice indices. Here the
analogous missing coordinates are the rational torsion correction and the
finite part of Sha: a rank statistic alone discards both.

## 1. Low corank has a finite growth certificate

Let E/Q be an elliptic curve, q a prime, r its rational rank, and s its full
q-infinity Selmer corank. For m >= 0 define

    u_m = log_q #Sel_(q^m)(E/Q) - log_q #E(Q)[q^m],   u_0=0.

These are full multiplication-by-q^m Selmer groups, not isogeny Selmer
groups. Standard finite Kummer exactness gives

    u_m = mr + log_q #Sha(E/Q)[q^m].                        (1)

No finiteness conjecture for Sha is needed. Cofinite generation and the
structure theorem for finitely generated Z_q modules give, noncanonically,

    Sha[q^infinity] = (Q_q/Z_q)^delta + F,
    F = direct_sum_i Z/q^(e_i),       s=r+delta.

Consequently

    u_m = ms + sum_i min(m,e_i),                           (2)
    u_(m+1)-u_m = s + #{i:e_i >= m+1}.                     (3)

**PROVED.** The following are equivalent:

1. s <= 1;
2. for some m >= 1, u_m < 2m;
3. for some m >= 0, u_(m+1)-u_m <= 1.

Indeed, each quotient u_m/m and each difference is at least s. Conversely,
the finite sum in (2) is eventually constant, and the count in (3)
eventually vanishes. The differences form a nonincreasing sequence of
nonnegative integers with eventual value s. This is an exact
semidecision interface: a successful finite receipt proves low corank;
failure at any prescribed finite depth does not prove high corank.

Exact orders at two adjacent levels suffice; compatible transition maps
are not an extra hypothesis for this numerical implication, because the
groups' full Selmer identities supply (1). Authentication of those groups
and their local conditions remains essential. One may also use a rigorous
upper bound at level m+1 and a rigorous lower bound at level m. Subtracting
two upper bounds is not a valid upper bound for the difference.

**Hostile model.** Take s=1 and F=(Z/q^k)^2. Then u_m=3m through level k,
exactly like a corank-three group with no finite part. The two-level
certificate first succeeds at (k,k+1), whereas the one-level average first
succeeds at m=2k+1. This is an abstract cofinitely generated group control,
not a claim that every such model is realized by an elliptic curve.
Arbitrarily deep agreement of initial statistics can hide the decisive
finite-versus-divisible distinction.

**BSD application, conditional on the supplied theorem.** Each authenticated
receipt in items 2 or 3 triggers Theorem 1.1. It is unnecessary to know Sha
is finite before constructing the receipt; finiteness is a consequence of
the cited theorem after low corank is established.

## 2. An exact curve family using Monsky's matrix

For positive odd squarefree n with prime factors p_1,...,p_k, let

    E_n: y^2=x^3-n^2 x,
    [a/p]=(1-(a/p))/2 in F_2.

Define A_ij=[p_j/p_i] for i != j, choose A_ii so every row sums to zero,
and put D_a=diag([a/p_i]). The cited Monsky formula is

             [ A+D_2     D_2     ]
    M_n =    [                  ],                        (4)
             [ D_2       A+D_-2 ]

    dim Sel_2(E_n/Q) - dim E_n(Q)[2] = nullity(M_n).

Source: Watkins, Donnelly, Elkies, Fisher, Granville and Rogers,
[*Ranks of quadratic twists of elliptic curves*, section 4.1](https://dms.umontreal.ca/~andrew/PDF/RankQuadTwist.pdf).
The displayed model is rationally isomorphic to their twist model; n is
odd and squarefree throughout. Rational 2-torsion has dimension two,
including n=1 with an empty Monsky matrix.

Thus a binary matrix of nullity at most one provides the m=1 receipt in
section 1. For a single prime p, the four residues modulo eight give

| p mod 8 | D_2 | D_-2 | nullity M |
|---|---:|---:|---:|
| 1 | 0 | 0 | 2 |
| 3 | 1 | 0 | 0 |
| 5 | 1 | 1 | 1 |
| 7 | 0 | 1 | 1 |

The p=1mod8 row is inconclusive for low corank, not a rank-two conclusion.
These classical small-prime consequences are not claimed as new BSD
theorems. The contribution is the explicit finite certificate interface
and its integration with the accepted input.

**FINITE-EXACT universe.** All 4,056 positive odd squarefree n <= 10,000
were included, with no rank, prime-residue or point-search filter.

| Pure 2-Selmer dimension | Number | Least n |
|---:|---:|---:|
| 0 | 1,018 | 1 |
| 1 | 1,868 | 5 |
| 2 | 988 | 17 |
| 3 | 164 | 205 |
| 4 | 18 | 1513 |

Exactly **2,886 curves** have a low-corank certificate at this first level.
The other 1,170 are unresolved by this test, not BSD counterexamples.
Each matrix rank was checked independently by enumerating its entire
binary kernel. The script authenticates squarefreeness and complete prime
factorization by exact trial division. It does not compute arbitrary
Selmer groups outside this cited matrix family.

## 3. A natural tournament detects the obstruction

Assume every p_i is 3 modulo 8. Then D_2=I and D_-2=0. The equations
M_n(x,y)=0 eliminate y and give the explicit isomorphism

    ker(A^2+A+I) -> ker(M_n),
    x |-> (x,(A+I)x).                                    (5)

Because X^2+X+1 is irreducible over F_2, its kernel is naturally a vector
space over F_4 with the primitive cube root acting by A. In particular the
pure Selmer dimension is even. In this subfamily, a first-level low-corank
receipt is therefore exactly invertibility of A^2+A+I.

Quadratic reciprocity says A_ij+A_ji=1 for distinct i,j. This supplies an
intrinsic tournament: vertices are the prime factors; i -> j means
[p_j/p_i]=1. There are no ties. Its row-parity diagonal is retained in A.
The target is the kernel in (5), not chromatic number or a geometric
embedding. Erasing the prime labels loses n; erasing their residue class
loses the fixed diagonal matrices in (4).

For three vertices, a transitive tournament makes A triangular with
diagonal entries zero or one, so A^2+A+I is invertible. A cyclic tournament
has A=I+C with C the 3-cycle permutation matrix. Its all-one vector has
A-eigenvalue zero, and its complementary two-dimensional space is exactly
ker(A^2+A+I). Hence

    transitive triple -> pure Selmer dimension 0,
    cyclic triple     -> pure Selmer dimension 2.          (6)

Actual arithmetic witnesses are

    1947 = 3*11*59: transitive, nullity 0;
     627 = 3*11*19: cyclic, nullity 2.

The latter is an explicit next higher-descent obligation. The cyclic
kernel does not by itself distinguish finite Sha classes from surviving
rank or divisible Sha. A complete Sel_4 order would make the adjacent
growth test available.

**2026-10-13 continuation:** [Tunnell terminals](bsd_tunnell_terminals_20261013c.md)
now ground E17 and E627 without computing Sel_4: exact central-value
nonvanishing and the classical CM theorem give rank zero and finite Sha.
The same second gate settles724 of this packet's1170 first-gate residuals;
the unchanged4056-source union covers3610. The first-level counts above
remain correct and retain their narrower scope.

Every prescribed finite tournament occurs for primes all 3 modulo 8:
choose a new prime in specified nonzero square/nonsquare classes modulo
the previous primes and in class 3 modulo 8, then use CRT and Dirichlet's
prime-in-progressions theorem. Reciprocity fixes the reverse arcs. This
is an existence statement, not a claimed frequency distribution of prime
tuples.

## 4. Even-order switching preserves this Selmer statistic

This uses the same cut reversal as
[THM-474, tilings are switching classes](../../01-canon/theorems/THM-474-tilings-are-switching-classes.md),
but the preserved arithmetic predicate and parity restriction are new to
this packet. It does not identify tournament isomorphism with switching;
[THM-1415, canonical star quotient](../../01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md)
already records the severe loss in that quotient.

Let k be even, let 1 be the all-one vector, and let u be the indicator of
the switched vertex set. For the row-zero binary tournament matrix,

    A+A^t = 11^t+I,       A1=0,       1^t A=1^t,
    A' = A+u1^t+1u^t+(sum u)I.                            (7)

Set W=1-perp and L=span(1). Since k is even, L is contained in W. Both
are A-stable. On V/W, A acts as 1; on L, it acts as 0. The polynomial
Phi_3(X)=X^2+X+1 equals one at both 0 and 1, so

    ker Phi_3(A) = ker Phi_3(A acting on W/L)

up to a canonical isomorphism. Explicitly, a kernel vector lies in W;
the quotient is injective on it. If Phi_3(A)w=a1, replace w by w+a1 to
lift a quotient kernel vector uniquely. On W/L, (7) changes A only by
(sum u)I, and Phi_3(X+1)=Phi_3(X). Thus:

**PROVED.** For even k, cut switching preserves nullity(A^2+A+I).

Normalize the last marked vertex to a sink, using the unique switching
mask with last bit zero, and delete it. If B is the row-zero matrix on
the remaining k-1 vertices, the normalized matrix is

             [ B+I   1 ]
    A_norm = [         ].                                (8)
             [  0    0 ]

Equation (8) preserves the Phi_3 nullity once more. The smaller tournament
and the retained k-1 switching bits recover the original exactly. Therefore
the labelled distribution satisfies, for every even k and every d,

    #{T on k vertices: dim ker Phi_3(A_T)=d}
      = 2^(k-1) #{T on k-1 vertices: dim ker Phi_3(A_T)=d}. (9)

This is a lossless representation when the mask is kept, and a sufficient
quotient for this particular statistic when it is discarded. At odd k it
fails: switching one vertex of a cyclic triangle makes it transitive and
changes the nullity from two to zero. The parity restriction is necessary.

The complete labelled atlas gives:

| k | nullity 0 | nullity 2 | nullity 4 |
|---:|---:|---:|---:|
| 1 | 1 | 0 | 0 |
| 2 | 2 | 0 | 0 |
| 3 | 6 | 2 | 0 |
| 4 | 48 | 16 | 0 |
| 5 | 624 | 360 | 40 |
| 6 | 19,968 | 11,520 | 1,280 |

Appending a universal source or sink preserves the Phi_3 nullity at either
parity by the same block computation. This gives infinite arithmetic
extensions preserving a low-corank certificate, with primes chosen as
above. It does not preserve arbitrary higher Selmer groups, Sha, regulators,
or analytic leading terms of the changing curves.

## 5. Reproduction and the next decisive obligations

Run from the repository root:

    python 04-computation/experiments/bsd_selmer_growth_20261013a.py
    python -O 04-computation/experiments/bsd_selmer_growth_20261013a.py

[Exact output](bsd_selmer_growth_20261013a.json) records 247,561 always-active
checks: 1,705 abstract group models; all 4,056 arithmetic matrices with
independent exhaustive kernel counts; all labelled tournaments through
six vertices; reversible marked-sink receipts; and named positive and
hostile cases. Arithmetic and abstract-tournament universes are kept
separate. Matrix input authentication is arithmetic, not a random test.

The next BSD task is **not** another nullity census: compute a complete
4-Selmer certificate or the relevant higher pairing for n=627 (and n=17),
then apply the difference test with exact torsion correction. For the
arithmetic tournament recursion, derive what additional marked pairing
data must accompany switching to preserve a higher-descent predicate.
For the rank-three THM-3646 curve, the supplied low-corank input does not
apply; a genuinely higher-rank theorem or separate leading-term argument
is required.

The connection to the other current lanes is specific. A successful
finite certificate must retain the data needed for its terminal predicate:
original-source charge for Collatz, complete geometric relations for
plane graphs, characteristic-zero derivative closure for Jacobian, and
full Selmer growth here. Similarity of their scalar summaries does not
transfer any conjecture's conclusion.
