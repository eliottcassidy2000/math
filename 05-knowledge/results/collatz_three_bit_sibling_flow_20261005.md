# Three bits, three colours, and a finite kernel for sibling flows

2026-10-05. **PROVED, scoped:** source disintegration, entropy and sibling
scaling, inverse-fibre summation, the equivalent base-flow criterion, and
the finite-kernel compiler. **FINITE-EXACT:** the supplied kernel and checks.
**OPEN:** a positive summable base weight on every base, the uniform refuel
bound, and universal Collatz. No claim of literature priority or canon
promotion is made.

The useful 3+1 connection is operational: the two source measures differ
exactly where sibling insertion crosses a binary height boundary. The
adjusted measure removes that discrepancy, allowing every infinite inverse
fibre to be summed and replaced by one checked base edge. A finite grounded
kernel consequently certifies infinitely many inputs. Covering all bases
remains the precise outstanding requirement.

Artifacts: [program](../../04-computation/experiments/collatz_three_bit_sibling_flow_20261005.py),
[exact output](collatz_three_bit_sibling_flow_20261005.out), and
[finite kernel](collatz_three_bit_sibling_flow_20261005.json).
Two independently checked supplements develop the
[three-colour lift](collatz_three_colour_lift_20261005.md) and
[2,10,42 neighbour macros](collatz_2_10_42_neighbour_macros_20261005.md).

## 1. Inheritance and the concept board

Anchor: source-labelled paid coverage. Niche: the exact integer measures.
Wildcard: marked Zeckendorf colours and the Berggren sequence2,10,42.
The board is **source entropy / climb and cofactor / marked colours /
sibling quotient / affine carry / refuel budget**.

The closest mechanism is the summable incoming-flow equivalence in
[effective prefix mass](collatz_effective_prefix_mass_20261005.md), P6,
and the gamma source measure in
[algorithmic Collatz measure](algorithmic_collatz_measure_20261005.md), §4.
The canonical hostile is9→7→11→17→13: its endpoints have the same binary
height and climb register. The corrected near miss is identifying a
contracting coefficient with actual descent while dropping its carry.
The least-used coordinate is the borrow in sibling insertion on the source
index, together with the minimal inverse-fibre representative.

The sequence is inherited from
[ternary Berggren](ternary_berggren_20260925.md), §4:
H(j)=2(4^j−1)/3 conjugates phase increment to height h→4h+2 and is a
3-adic isometry. Its ordinary nonnegative image remains sparse. Residue
completeness of that phase map is not coverage of all ordinary heights.
The applied meta-patterns are **Separate local support from bounded-height
coverage** and **Retain local incidence through counting**. The source
label, sibling depth and finite ROOT graph remain attached to every mass
calculation; no new repository-wide meta-pattern is promoted here.

## 2. Three bits split into the climb register and its cofactor

For a positive odd source n, let m=(n+1)/2, and let ell be binary bit length.
The inherited weights are

    w(m)=2/4^ell(m),
    mu(n)=w((n+1)/2),       nu(n)=8/(3*4^ell(n)).              (1)

Both sum to1. The gamma-index blocks have mass2^(-k-1), k>=0.
For nu, n=1 has mass2/3, and the odd numbers of bit length l>=2
have total mass(2/3)2^-l.

**C1. Exact product law.** Write uniquely

    n=2^h*t-1,   h=v2(n+1)>=1,   t odd>=1.

Then

    mu(n)=q(h)*nu(t),       q(h)=3/4^h.                    (2)

Indeed ell((n+1)/2)=h-1+ell(t), so both sides are8/4^(h+ell(t)).
The map is a bijection between positive odd n and these pairs(h,t).
Consequently h and t are independent under mu. The entropy calculation is

    E_q h=4/3,
    H(q)=8/3-log2(3),
    H(nu)=log2(3)+1/3,
    H(mu)=H(q)+H(nu)=3 bits.                              (3)

For the second identity involving nu, E_nu ell(n)=5/3 and
-log2 nu(n)=2ell(n)+log2(3)-3. The entropy addition follows from the
bijection and independence, not from equating three colours with three bits.

The q in(2) is exactly the distribution used as the squared-Haar tilt
q(a)=3/4^a in the previous likelihood construction. Its coordinate here is
h=v2(n+1), not the next Collatz valuation a=v2(3n+1). At n=1 these are1
and2; at n=3 they are2 and1. Equality of laws does not equate these observables.

This product law is compatible with an actual operation. While h>=2,

    U(n)=(3n+1)/2=2^(h-1)*(3t)-1,
    (h,t) -> (h-1,3t).                                   (4)

The climb register decreases while the cofactor grows. At its exit a new
register can be large; dropping t would lose the information controlling
that replenishment. The9→13 hostile rules out every positive correction
mu(n)F(h(n)) or nu(n)F(h(n)) satisfying strict discounted edge payment,
because both baseline weights and h agree at the two endpoints.

There is also a literal three-head-plus-continuation parser for gamma indices:
heads1,2,3 have masses1/2,1/8,1/8, and the remaining mass1/4 obeys

    w(4q+r)=w(q)/16,   q>=1, 0<=r<4.

Conditional on continuation, q has the original gamma law and r is uniform.
Its entropy equation is H=7/4+(2+H)/4, again giving3. Regrouping binary
digits this way does not turn the three heads into equiprobable colours.

## 3. The comparison factors repair genuine sibling insertion

Let D be the source indices m that are powers of2. The inherited comparison is

    nu(2m-1)/w(m)=4/3 on D, and1/3 off D.                 (5)

Thus mu(D)=2/3 and nu(D)=8/9. Equivalently, on odd source labels,
nu=(1/3)mu+mu restricted to the Mersenne sources n=2^j−1.
The shared climb marginal is q, but its smallest ordinary representative
has conditional probability2/3 under mu and8/9 under nu. A climb-only
observation loses exactly this height distinction.

The genuine sibling map is

    S(n)=4n+1,       U(S(n))=U(n),
    S^k(n)=4^k*n+(4^k-1)/3.                              (6)

In gamma indices its action is

    M_k=(S^k(n)+1)/2
       =4^k*(m-1)+(a_k+1),  a_k=2(4^k-1)/3.             (7)

The remainder a_k+1 has base-four digits22...23. Retaining the borrow
m−1 is essential; plain concatenation onto m gives a different integer.
For k>=1 the exact weight laws are

    nu(S^k n)/nu(n)=16^-k,
    mu(S^k n)/mu(n)=16^-k * (4 if m in D, else1).         (8)

Proof: ell(S^k n)=ell(n)+2k. For gamma indices, M_k has bit length
ell(m)+2k−1 at a dyadic m, and ell(m)+2k otherwise. Also M_k is odd and
at least3, hence never dyadic. Substituting(5) cancels exactly the factor4.
This is the arithmetic role of the measure comparison: it makes the
same-successor operation have a constant mass scaling.

Equation(8) alone does not make a supplied source a sibling of a certified
source. The base and depth must actually decode and the base must be grounded.

## 4. The neighbours of2,10,42 expose opposite payment behaviour

The [neighbour supplement](collatz_2_10_42_neighbour_macros_20261005.md)
proves all of the following, with k>=1 and a_k=2(4^k−1)/3:

* a_k/2=(4^k−1)/3 has3(a_k/2)+1=4^k, hence is explicitly rooted.
* a_k+1 has the exact odd-step word
  (1,2 repeated k−1 times,2+epsilon_k), where epsilon_k=1 for even k
  and2 for odd k. It ends at(3^k+1)/2^epsilon_k and shares that endpoint
  with the smaller child3^(k−1). The macro pays; that child is not thereby
  supplied with a ROOT proof.
* For k>=2, a_k−1 has the exact word(2,1 repeated2k−2 times,2), ending
  at(3^(2k−1)−1)/2. It resets h to1 while expanding the source.

The rooted ray has exact mu-mass19/30 and nu-mass32/45. This proves a
large atomic mass of a sparse certified family, not a natural-density claim.
The adjacency is useful precisely because the same short recurrence encodes
both a ROOT mechanism and the large refuel bill. The offset must be retained.

## 5. Summing every inverse fibre leaves one base edge

Let B be the positive odd bases satisfying v2(3b+1) in{1,2}. Every positive
odd n has a unique decomposition n=S^k(b), b in B, with

    k=floor((v2(3n+1)-1)/2).                              (9)

An independent decoder repeatedly replaces n by(n−1)/4 when n=5 mod8.
Every step decreases n, and the stopping condition is exactly membership
in B. Thus this decoder is total without assuming Collatz convergence.

For an odd target y, the inverse fibre is empty if3 divides y. Otherwise
its unique base is

    b0(y)=(2y-1)/3 if y=5 mod6,
    b0(y)=(4y-1)/3 if y=1 mod6,                         (10)

and all its predecessors are S^j(b0(y)), j>=0. This follows by solving
3n+1=2^a*y: a has the required parity modulo2, and adding2 to a applies S.

**C2. Exact incoming sums.** For y>1, with the incoming operator killed at1,

    K nu(y)=(16/15)nu(b0(y))                            (11)

when b0 exists, and0 otherwise. The analogous multiplier for mu is19/15
if(b0+1)/2 is dyadic, and16/15 otherwise. These are convergent geometric
sums, not truncated orbit simulations. In particular

    K nu(y)/nu(y) in {0,4/15,16/15,64/15}.               (12)

The target5 has ratio64/15 and target7 has4/15. The baseline nu is
therefore not a discounted flow. The root target1 is excluded from K;
if it were retained while deleting the root source, its ratio would be1/15.

Under nu the sibling depth k is independent of its base, with
Pr(k=j)=(15/16)16^-j and normalized base weight(16/15)nu(b). Consequently
any finite certified set of distinct bases C has infinite sibling closure
of exact mass(16/15)sum_(b in C)nu(b). Keeping the unique base avoids
double-counting multiple certificates of the same family.

## 6. An exact base criterion and the bill for changing patterns

Fix rational0<r<1 and0<rho<1. Let g be positive on every b in B and summable.
Extend it by

    f(S^j b)=r^j*g(b), including f(1)=g(1)>0.            (13)

This is summable with total(sum_B g)/(1−r). The root value is bookkeeping;
no contraction of the root self-loop is imposed.
For each b in B\{1}, decode its actual endpoint as

    U(b)=S^{k(b)}(G(b)),

where k(b)>=0 and G(b) belongs to B. In this display k(b) is a function
of b; S is iterated k(b) times on G(b). Every nonroot base has U(b)>1.

**C3. Fibre elimination.** The full killed inequality Kf<=rho*f is
equivalent to the individual base inequalities

    g(b)<=rho*(1-r)*r^{k(b)}*g(G(b)),    b in B\{1}.     (14)

Indeed the row at y=U(b)=S^k(c) has incoming sum g(b)/(1−r) and value
f(y)=r^k*g(c). Every target with predecessors has exactly one such minimal
base; zero-incoming rows impose no condition. Multiple base edges can
share c, but their target rows S^k(c) are distinct. This is why(14) has
no additional sum over those edges.

Take any finite bound W>=sum_B g, and set R(b)=log(W/g(b)). Then(14) is
the explicit budget condition

    R(b)-R(G(b)) >= log(1/[rho*(1-r)]) + k(b)*log(1/r).  (15)

The first charge pays the full incoming fibre and strict progress; the
second pays the sibling depth used to change the base. Because R>=0 and
the first charge is positive, an infinite sequence of these paid changes
is impossible. This gives a concrete refuel obligation rather than an
expected gain under an unspecified distribution.

**C4. Exact equivalence, not a solution.** For each fixed rational r and rho in(0,1),
such a positive summable g exists if and only if every positive integer
reaches1. Sufficiency follows either from(15) and common-future transport,
or from the inherited summable-flow argument on actual U-orbits.

For necessity, assume convergence. Every G-path terminates at1: taking U
strictly decreases actual root rank, and subsequently removing sibling
depth preserves that rank or reduces it further when the base is1.
Enumerate bases b_j, j>=1, including1. For a path
b_0=b_j,...,b_T=1 set d_b=rho*(1−r)*r^{k(b)}, choose
epsilon_j=2^-j/(T+1), and contribute at its i-th vertex

    epsilon_j * product_(l=i to T-1) d_(b_l).

Each path contributes total at most2^-j. The sum g is positive at every
base and summable. Each contribution through a nonroot b also occurs at
G(b) divided by d_b; therefore g(b)<=d_b*g(G(b)). Under the assumption,
this construction is computable with a tail bound sum_(j>J)2^-j. It uses
the terminating paths and cannot be invoked to establish their existence.
There is no independent no-go theorem against fixed r.

## 7. A finite certificate compiles infinitely many grounded inputs

The finite form of the preceding construction is useful without an
assumption about other inputs. Supply a finite set C of bases, containing1,
with checked equations U(b)=S^k G(b), parents in C, and every parent chain
ending at1. The verifier checks those equations, integer types, and the
finite dependency graph. It does not trust discovery search or an assumed
rooted seed. Every base and every sibling S^j(b) is then rooted by repeated
actual common-future transport.

Finite path contributions give g>0 on C and g=0 off C. Extending by(13)
gives a summable discounted flow supported on the infinite union of these
sibling families. Its positivity is intentionally restricted to this
certified union; promoting it to all positive inputs would be invalid.

The supplied experiment starts with all48 bases<=127 and closes their
explicitly checked G-paths. Its retained kernel has88 bases, largest2051,
and maximum parent-chain length42. It certifies all odd inputs<=127 and
infinitely many additional siblings, with exact source masses

    mu(C_sibling)=1946497/1966080,
    nu(C_sibling)=5878417/5898240.                       (16)

These remain strictly below1. They are probabilities under the stated
atomic priors, not percentages of integers in natural density. For
r=1/16,rho=1/2, every supported base inequality passes and sum_C g<1.
The kernel contains literal edge equations; the construction script can
reconstruct the weights from its finite paths.

The callable `compile_root_word(entries,n)` also generates a literal strict
first-hit ROOT word for any supplied n in the certified sibling union,
without an orbit search. Inserting depth j into a nonroot base's first
edge adds2j to that edge's valuation and retains the suffix. For base1,
j=0 has the empty ROOT word and j>0 has the single valuation2+2j. These
rules also reconstruct each base from its checked parent equation. The
compiler rejects unlisted bases, rather than treating them as rooted.

## 8. What the Zeckendorf connection contributes

The [colour supplement](collatz_three_colour_lift_20261005.md) recovers the
ordered extra-unit alphabet1_K,1_B,1_R,2,3,5,... and the separately defined
F2^2 charge. The owner stripe, charge, and independently decorated source
words are kept distinct. On the independent four-state lift, zero projects
to0 and three nonzero colours to1, so a binary word of cost A with r ones
has3^r lifts and likelihood3^r/2^A. This realizes the previous slope
identity as a finite-fibre count.

The aggregate charge has matrix J4−I4 and eigenvalues3,−1. It compresses
counts exactly but loses the word order and carry. One further unrestricted
symbol makes this charge uniform independently of every earlier finite
guard. Thus a colour restriction alone supplies no new paid guard.
At intermediate shortcut cuts the exact boundary factor is3^(p_end−p_start).
For common-future receipts the relevant coefficient is the ratio of the
source and child likelihoods; the paid J rule supplies an explicit example
whose source likelihood exceeds1 while the relative coefficient is243/256.

| source -> target | preserved | lost unless retained |
|---|---|---|
| n -> (h,t) | exact mu product law and valuation-one update | nothing; dropping t loses refuel data |
| marked Zeckendorf representation -> value | integer | ordered unit marker and carry |
| decorated word -> binary word | likelihood through fibre count | colour history |
| word -> aggregate charge | four-state convolution count | order, source guard and affine carry |
| n -> sibling base | next U-target | depth; restore it for weight and actual source |
| finite base kernel -> infinite sibling flow | ROOT proofs and exact discounted inequalities | positivity at unlisted bases |

The next productive target is an explicit summable g satisfying(14) for
every base, or a sound expanding family of finite kernels with a proved
atomic residual tending to0. The source entropy, colour lift and sibling
normalization make the coordinates and charges exact. Neither currently
proves that residual decay or the global refuel inequality.

## 9. Reproduction and hostile controls

    python -X utf8 -B 04-computation/experiments/collatz_three_bit_sibling_flow_20261005.py
    python -O -X utf8 -B 04-computation/experiments/collatz_three_bit_sibling_flow_20261005.py

Both runs agree:195,544 exact checks, with source disintegration on8192
positive odds, six sibling depths per source, direct valuation enumeration
versus sibling inverse fibres through odd target1023, exact geometric-tail
completion, and the finite kernel. Malformed edge, root, duplicate,
boolean-depth and floating-valuation controls are rejected. The kernel
JSON is regenerated deterministically. Infinite statements above have
algebraic proofs; these finite tests do not establish universal coverage.
The strict ROOT compiler additionally replays all88 bases at sibling
depths0 through6, with no earlier ROOT visit; its longest word has43
accelerated edges. Each package received an independent proof and code audit.
