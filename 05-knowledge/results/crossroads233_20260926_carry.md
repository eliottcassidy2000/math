# Depth 233: collision decoding and a pointwise height certificate for rank

**Status: PROVED, independently audited, with FINITE-EXACT controls.**
The inherited depth-233 collision is established. The new statements below are
self-contained elementary proofs, not a proof of Collatz or of G2. No novelty
claim is made for eventual multiplicative independence of distinct affine forms.
The root coordinator independently audited sections 2--4; the geometry lane
independently audited sections 3--4 and checked the positive-integer type
boundary in section 2. The latter is explicit below.

## 1. Inheritance, board, and exact scope

The closest proved mechanism is **THM-4490**,
[`THM-4490-collatz-affine-word-sunit-specialization.md`](../../01-canon/theorems/THM-4490-collatz-affine-word-sunit-specialization.md):
subword carries control incident prime banks, giving full rank for almost all
sources through a growing prefix, using a cited S-unit bound. The canonical
hostile here is the two-source, same-slope collision at depth 233 in
[`crossroads_20260926_flow.md`, section 8](crossroads_20260926_flow.md).
The corrected near miss is fixed-slope injectivity: it holds at all heights
through depth 31, fails at depth 233, and is not an input to the surviving
polynomial carry capacity. The least-used sidecar is the *magnitude*, rather
than merely the prime support, of each incident carry.

The live board is: ordered carries; rational-base block decoding; valuation
relations; gcd overlaps; and source height under a no-descent condition.
The anchor is a pointwise rank statement beyond almost-all counting. The niche
is what actually happens when a collision gadget is composed. The wildcard is
the user's triangle-inequality/unique-relation suggestion, tested on the exact
valuation matrix rather than on a graph with arbitrarily oriented edges.

The current G2 question concerns independence on every injective chronological
orbit; it remains OPEN. Full rank is not a certificate of descent, and an
additive source-to-hub collision is not a multiplicative relation among the
nodes on one path.

## 2. The collision does not amplify by naive block concatenation

Let `T(n)=n/2` for even n and `(3n+1)/2` for odd n. Put

    A=3^153, B=2^233,
    u=A^(-1) mod 2^80 = 186937257649965781719819,
    N=2^153 u-1, M=N-4, v=(A u-1)/2^80.

The inherited actual words are `w_0=1^153 0^80` and

    w_1=(110)^51 01110111100111110111011110110011111111110010011110000100101111010011001101010101.

They take N and M to v, respectively. Both have 153 odd letters, every
nonempty prefix has slope above one, and their carries are

    C_0=A-2^153, C_1=5A-2^153=C_0+4A.

Replacing the two sources by `N+tB,M+tB` preserves the words and gives the
common endpoint `v+tA` for every integer t>=0. Therefore no large-height
condition can restore source-to-hub injectivity at this fixed depth.

**Block decoding lemma.** More generally, let A,B be positive integers with
`gcd(A,B)=1`, and let a finite
integer digit set D have distinct residues modulo A. Consider fixed-length
affine blocks

    F_d(x)=(Ax+C_0+Ad)/B, d in D.

Among all length-t block strings, at most `#D` integer sources can end at one
fixed integer endpoint. In fact their digits after the first block must agree.

To prove this, the difference of two total normalized carries is

    sum_(j=0)^(t-1) (d_j-d'_j)(B/A)^j.

An integer source difference at a common endpoint requires this expression
to be an integer. If the latest differing index is s>=1, multiplication by
`A^s` and reduction modulo A gives
`(d_s-d'_s)B^s=0 mod A`, contrary to the residue assumption. Thus only the
first digit can differ. No distributional model or parity approximation is
used; imposing actual integer guards can only remove members of the family.

For the 233 blocks, D={0,4}. Exactly `2^(t-1)` distinct total-carry residues
modulo `A^t` occur, each with two symbolic words. All block strings retain
the strict prefix slope property. The script checks all 2046 strings through
ten blocks, and obtains and replays a positive integer representative for all
126 strings through six blocks. One collision therefore supplies a two-way
fibre, not an exponentially growing fibre under this particular composition.

The source-to-target map is an ordered block word to its normalized carry in
`Q/Z`. It preserves the possibility of an integer source difference at a
common endpoint. It discards the endpoint's integrality and the actual parity
guards, which are retained separately for the finite positive controls. The
bound concerns this two-block alphabet, not all words of length 233t.

## 3. A valuation triangle inequality gives an effective rank certificate

**Gcd dominance lemma.** If positive integers `a_0,...,a_(r-1)>1` are
multiplicatively dependent, then there exists an index i with

    a_i divides product_(j!=i) gcd(a_i,a_j).             (1)

Consequently the strict inequalities

    a_i > product_(j!=i) gcd(a_i,a_j) for EVERY i        (2)

are sufficient for multiplicative independence.

Choose a nonzero integral relation `product a_i^c_i=1`, and choose a supported
i of maximal `|c_i|`. For every prime p, its valuation equation implies

    |c_i| v_p(a_i) <= sum_(j!=i)|c_j|v_p(a_j)
                   <= |c_i| sum_(j!=i)v_p(a_j).

Thus `v_p(a_i)<=sum_(j!=i)v_p(a_j)`. Truncating each term on the right at
`v_p(a_i)` preserves that inequality: either one term already reaches that
value, or no term is changed. This proves (1), prime by prime. The proof
identifies one maximal-coefficient supported node; it does **not** show that
every node violating (1) has zero coefficient in every relation.

For an actual positive odd Collatz trajectory let

    m_0=n, m_(j+1)=(3m_j+1)/2^k_(j+1),
    S_0=0, S_j=k_1+...+k_j, k_j>=1.

For i<j define the positive integer subword carry

    C_(i,j)=sum_(ell=i)^(j-1) 3^(j-1-ell)2^(S_ell-S_i).

Then

    2^(S_j-S_i)m_j=3^(j-i)m_i+C_(i,j).                 (3)

Since all nodes are odd, `gcd(m_i,m_j)` divides `C_(i,j)`. Set

    Q_i=product_(j!=i) C_(min(i,j),max(i,j)).

For r=N nodes, applying (2) to `a_i=3m_i` shows that

    m_i > 3^(N-2) Q_i for EVERY i                      (4)

is sufficient for full first-slot rank N, and hence for independence of the
actual pairs `(-3m_i,3m_i+1)`. A relation of pairs would give, on taking the
absolute value of the first coordinate, a forbidden relation among `3m_i`.

For any fixed finite valuation word this is already an explicit height
cutoff: because `m_i >= 3^i n/2^S_i`, it is sufficient that

    n > max_i [3^(N-2) Q_i 2^S_i/3^i].                (5)

All quantities on the right are exact computable rationals determined by the
word, so this proves fixed-word eventual independence without an S-unit
theorem. It does not provide the strong global count of exceptional sources
that the S-unit method supplies.

## 4. Pointwise no-dip theorem, with an explicit threshold

**Theorem.** Let N>=2. Suppose the first N odd nodes, indexed
`0,...,N-1`, satisfy `m_i>=n=m_0`, with n>=N. Define

    H_N = 2^(N-1) (N-1)! 3^(N(N-2))
          / 2^((N-1)(N-2)/2).                          (6)

If `n>H_N`, then these N first slots and their N actual pairs are
multiplicatively independent. The bound is rational, not necessarily integer.
For example `H_2=2`, `H_3=108`, `H_4=39366`, `H_5=86093442`.

For ell<N, the exact product identity and the no-dip hypothesis give

    2^S_ell <= 3^ell product_(h<ell)(1+1/(3m_h))
             <= 3^ell(1+1/(3n))^ell < 2*3^ell.        (7)

The final strict inequality follows from `ell<N<=n` and `exp(1/3)<2`.
Using (7) term by term in (3), and `S_i>=i`, gives

    C_(i,j) < 2(j-i)3^(j-1)/2^S_i
             <= 2(j-i)3^(j-1)/2^i.                    (8)

Consequently `Q_i<R_(N,i)`, where

    R_(N,i) = 2^(N-1) i! (N-1-i)! 3^E3 / 2^E2,
    E3=(N-1)(N-2)/2+i(i-1)/2,
    E2=i(N-1)-i(i+1)/2.

The exact successive ratio is

    R_(N,i+1)/R_(N,i)
      =[(i+1)/(N-1-i)] 3^i 2^(i-N+2).                (9)

These ratios strictly increase with i, so log R is strictly convex and its
maximum is at an endpoint. The late endpoint is at least the early endpoint,
since their ratio is `(3/2)^((N-1)(N-2)/2)`. Thus

    Q_i < R_(N,N-1)
        = 2^(N-1)(N-1)!3^((N-1)(N-2))
          /2^((N-1)(N-2)/2).

Now `n>H_N=3^(N-2)R_(N,N-1)` implies (4). This proves the statement.
For N=2 the two endpoint bounds agree; the proof retains strictness in (8).

Writing `alpha=log_2 3`, Stirling's formula in (6) gives

    log_2 H_N=(alpha-1/2)N^2+N log_2 N+O(N).           (10)

Therefore for every fixed

    0<c<1/sqrt(log_2 3-1/2)=0.960047311978290...,

every sufficiently large odd n whose first
`N=floor(c sqrt(log_2 n))` odd nodes stay at least n has full rank N.
The endpoint constant c equal to the displayed value is not proved: the
positive `N log N` correction remains.

**Finite no-descent-horizon sidecar.** Put `ell=ceil(N log_2 3)+1` and assume
`n>=max(N,2^ell)`. If the half-Collatz path does not descend below n through
`floor(log_2 n)` half-steps, then `S_N<=ell`. Otherwise the first ell
half-steps have at most N odd letters, and all intermediate values are at
least n. Their exact product expansion gives

    T^ell(n)/n <= (3^N/2^ell)(1+1/(3n))^N
                <= exp(1/3)/2 < 1,

a contradiction. In particular all first N nodes lie within the known
no-descent horizon. Thus, with the additional condition `n>H_N`, the finite
theorem applies to that customary hard-source set. For N of order
`sqrt(log n)`, the sidecar's extra lower bounds hold for all sufficiently
large n.

This is pointwise for every sufficiently high no-dip source, whereas
THM-4490 gives a longer `sqrt(log H log log H)/8` prefix for almost all
sources in every interval of length H. The statements have different
quantifiers; neither subsumes all of the other. Neither eliminates a
particular putative divergent orbit, which could satisfy all these rank
conditions indefinitely.

**Descent-or-rank corollary.** For every odd n>max(N,H_N), either one of the
first N odd nodes is below n, or their N first slots are independent.
An infinite injective positive orbit consequently has arbitrarily long
consecutive full-rank odd blocks and infinite rank of its entire first-slot
family. Its odd nodes tend to infinity, since injectivity eventually removes
every finite set of positive integers. Future tail minima also tend to
infinity and are attained. Choose such a minimum above max(N,H_N) and apply
the theorem to the next N nodes. This proves unbounded rank, not G2's
independence of every chronological prefix, and supplies no contradiction
to divergence.

## 5. Controls, non-equivalences, and the graph boundary

The depth-233 collision has a particularly useful rank control. Each of the
153 odd source nodes on **each** of its two paths has a private prime relative
to that path. The script proves this without factoring any large integer:
repeatedly strip gcds with every other node, and check the remaining cofactor
is greater than one. Both initial sources and all later odd nodes are
coprime to 3, so every surviving private prime also certifies the corresponding
first slot `3m_i`. Hence both 153-column first-slot valuation matrices have
full rank, despite the common endpoint and same total affine slope.

The raw gcd criterion on the unscaled odd nodes succeeds on all 153 nodes of
the M path, and on 152 of 153 nodes of the N path (failure at index 17).
Private primes certify all nodes of both. These are sufficient certificates,
not claimed equivalent characterizations.

Cheap hostile controls are essential:

- `(6,10,15)` is multiplicatively independent but has no private prime in any
  coordinate and fails strict gcd dominance in every coordinate. Neither
  certificate is necessary.
- `(2,3,6)` has a one-dimensional relation space, generated by `(1,1,-1)`.
  Uniqueness of a relation does not make it vanish.
- For `(2,4)`, the matrix `log gcd(a_i,a_j)` is
  `log(2)*[[1,1],[1,2]]`, positive definite, but the multiplicative rank is one.
  Positivity of this natural gcd kernel cannot replace the valuation rank.
- Sources `n=(4^a-1)/3` can be arbitrarily large but go directly to 1, then
  repeat 1. Height without a no-dip condition does not give a uniform growing
  chronological rank theorem when valuation words are allowed to vary.
- Height alone does not repair the additive collision: its lifts above give
  actual collisions at arbitrarily high sources.

The valid source-to-target map for the user's triangle suggestion sends an
integer multiplicative relation to its coordinatewise prime valuation
equations; a maximal coefficient and the triangle inequality then imply (1).
It preserves the entire integral relation condition. Passing further to the
symmetric gcd overlaps loses signed valuation dependence, which is why only
a sufficient certificate is claimed. The necessary sidecar is the diagonal
size of **every** node, not just a selected large one.

There is no intrinsic antisymmetric binary relation here. Orienting pairs by
their numerical sizes produces a transitive tournament with a Hamiltonian
path for the dependent tuple `(2,3,6)` as well as for independent tuples.
Such a Hamiltonian path therefore retains the order of values and loses the
target predicate. No tournament conclusion is imported.

## 6. Literature boundary and reproducibility

Eventual multiplicative independence under large translations is already a
literature subject. Dubickas--Sha, *Multiplicative dependence of the translations
of algebraic numbers*, prove eventual independence of pairwise distinct
algebraic numbers translated by sufficiently large integers; see their
[primary preprint](https://arxiv.org/abs/1608.05458) and
[published paper](https://ems.press/content/serial-article-files/38742).
Ostafe--Sha--Shparlinski--Zannier study multiplicative dependence of rational
function values and dynamical orbits in their
[primary preprint](https://arxiv.org/abs/1706.05874).
These are background, not dependencies for the elementary proof above.
The present claim is the explicit Collatz incident-carry cutoff and its
pointwise no-dip consequence; no priority claim is made for the general idea.

Run [the exact script](../../04-computation/experiments/crossroads233_20260926_carry.py):

    python 04-computation/experiments/crossroads233_20260926_carry.py
    python -O 04-computation/experiments/crossroads233_20260926_carry.py

The explicit universes are: all 2046 two-block words through ten blocks;
126 actual positive concatenation replays; all 8835 multisets of lengths
2,3,4 from integers 2..20; 4950 exact bound ratios through N=100; all 1364
valuation words with N=2..6 and k_i=1..4 specialized above their exact
fixed-word cutoffs; all odd
sources 3..100000 with N=2..5, of which 66405 prefixes obey the no-dip guard;
12007 guarded finite-horizon sidecar checks with N=2..4; and 29 manufactured
high no-dip controls through N=30. There are 51339
above-threshold rank certificates including the high controls. All checks
use exact integers or rational numbers and survive optimization; the final
asymptotic decimal alone is floating-point display. The
[output](crossroads233_20260926_carry.out) is the exact normal stdout; the
normal and optimized stdout were separately compared and are byte-identical.

The decisive stopping point for Collatz is clear: these bounds constrain
possible multiplicative relations, while a putative divergent path may be
fully independent. The rank coordinate must enter an additional theorem
that forces descent or excludes a positive orbit. No such implication is
proved here.
