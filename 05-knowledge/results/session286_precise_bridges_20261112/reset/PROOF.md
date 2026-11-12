# Rerooting a central spin fiber: an exactly solvable long-range Ising system

**Status: PROVED scoped identities and extrema; independently reviewed in the
session audit. Chan–Pak remains OPEN.**

## Inheritance and the precise bridge

The physical rooted fusion, h13 core and central AA/BC trade are inherited
from [Session282](../../session282_chanpak_20261112/wildcard/PROOF.md).
Its [three-strand compiler and arithmetic reset](../../session282_chanpak_20261112/wildcard/ROOT_RESET.md)
are inherited as well. The sibling [automata lane](../automata/PROOF.md)
classifies the whole four-letter spin fiber. This note combines the operations
rather than presenting the old reset as new.

The source is a physical fixed-central-interface family. The map reroots at a
specified vertex and prepends a rooted directed triangle. Its target is an
explicit binary quadratic energy, with all pair interactions negative. The
preserved predicate is the exact Hamiltonian-path count of the resulting
literal tournament. The information **not** preserved is the old root's
continuation interface and its simple normalized numerical-semigroup action.
Keeping the new marked vertex and three-strand state repairs that loss.
A graph cut computes the minimum exactly, and the coefficients in fact give
both extrema in closed form. This is a real transfer of an optimization
method, not a claim that an Ising thermodynamic limit proves count coverage.

## 1. The physical spin cube

Let D have outgoing masks (18,8,3,5,14). Attach a common root z=0 using cut
masks A=9, B=3, E=18, C=24. The four six-vertex profiles in coordinates
(H,h,Q,U,V,d,e) are

    A=(37,13,104,46,71,8,5), B=(33,13,104,45,63,7,5),
    E=(37,13,104,53,64,7,6), C=(43,13,104,54,72,8,6).

For k>=3 take binary spins s_1,...,s_(k-2), set s_0=s_(k-1)=s_k=0, and let
letter i correspond to (s_(i-1),s_i): 00=A,01=B,11=E,10=C. Fuse the k letters
in increasing order; all arcs between different deletion blocks point from
earlier to later. Call this tournament T(s). Its final letter is A.

Let r be local vertex1 in that final block, global label 1+5(k-1). Reroot at
r and fuse a rooted directed triangle on the left. The result R(s) has5k+3
vertices and, by the inherited exact fusion law,

    F_k(s)=H(R(s))=H(T(s))+2V_r(T(s)).                 (1)

This is an actual construction at every k. It has no fractional or signed
physical multiplicities. The original H(T(s)) is a known central spin energy
with nearest-neighbor term -18*13^(k-2); (1) supplies the new long-range terms.

## 2. Exact quadratic response

Use zero-based variable indices i=0,...,k-3. Put c_k=F_k(0), and

    S(a,b)=[-34104960*471^a*104^b +3945600*104^a*471^b
             -140160*104^(a+b)]/367,
    a_i=2*13^(k-2)+2*S(k-i-3,i).

Then

    F_k(s)=c_k+sum_i a_i s_i+sum_(i<j) b_(j-i) s_i s_j,            (2)

where

    b_1=-1658880*104^(k-4)-18*13^(k-2),
    b_d=-55296000*104^(k-d-3)*471^(d-2),   2<=d<=k-3.             (3)

The b1 formula is used only when k>=4. Every pair coefficient is strictly
negative. The baseline c_k is explicitly evaluated by the18-state integer
compiler below; no enumeration of the2^(k-2) spin cube is needed.

### Why degree two is exact

In each deletion block the binary cut vector is

    x_i=x_A+s_i(x_B-x_A)+s_(i-1)(x_C-x_A).

The two difference vectors have disjoint supports, and x_E=x_B+x_C-x_A
entrywise. Thus every orientation bit incident with z is affine in the spins.
After the triangle insertion all other arcs are fixed, because the marked
final block is A. A Hamiltonian path uses at most two arcs incident with z.
Its count is therefore a multilinear polynomial of degree at most two in
those cut variables, hence degree at most two on the spin cube. Singleton
and pair evaluations determine (2) completely. This argument rules out hidden
higher interactions; it does not infer degree from a numerical fit.

## 3. Finite integer certificates for the all-order formulas

Retain the old root z. V_r counts an ending-r path and a disjoint complementary
path, possibly empty. Removing z splits these into three ordered strands
P,Q,R. The first strand may end at an inneighbor of z; the second must begin
at an outneighbor; the third is unrestricted. The state is(P status,Q status,
R occupied), of sizes3,3,2. Within a block distribute the five vertices among
these strands and count their literal directed paths. This is precisely the
inherited compiler, restated in `transfer.py` and independently checked by
literal18-vertex Hamilton-path counts.

Let F_X be its18-by18 integer transition for letter X, acting on column
states; F=F_A. Let u be the empty state and L the final-A ending-r functional.
The following exact certificates are checked on all matrix entries:

    F(F-13I)^2(F-104I)^2(F-471I)=0,
    Z=F_C F_B-F^2,  K=Z^2,
    FK=KF=104K,  ZK=KZ=0,  LKu=-27648000,
    Z F^m Z=471^m K  for m=0,...,5.                         (4)

The degree-six annihilator extends the last identities to every m>=0.
For adjacent spins put

    R=F_C F_E F_B-FZ-ZF-F^3.

Here **FR=RF=104R is false**. The valid boundary identities are

    LRu=-829440,  LRF=104LR,
    LF^aR=104^a LR for a=0,...,5.                            (5)

Again the annihilator extends the left equality to every a>=0. These are
observed row identities, not an unsupported full-matrix eigenclaim.
Finally, the36 exact checks

    L F^a Z F^b u=S(a,b),   0<=a,b<=5,                       (6)

prove the formula for all a,b>=0: each side obeys the same annihilator in
each index. The proposed expression consists only of13/104/471-compatible
exponentials (indeed only104 and471 here).

A single spin replaces AA by BC; (6) gives a_i. Two nonadjacent spins insert
two Z's, and (4) gives the distant term in (3). Adjacent spins replace AAA
by BEC; (5), together with the original central nearest-neighbor action,
gives b1. This proves (2)–(3). The number471 is the ordered three-path-cover
total of D, not an arbitrary fitted growth constant.

## 4. Exact extrema over an exponentially large family

**Theorem.** For every k>=3, the unique minimizing spin vector has

    s_i=1 for 0<=i<floor(k/2), and s_i=0 thereafter.

The unique maximizing vector is its complement on the k-2 free variables.
Thus both extremizing words are explicit and have length k.

Proof. Write a=k-i-3,b=i, P=104^a471^b, N=471^a104^b, M=104^(a+b).
Always M<=N and13^(k-2)<=13M. The ratio rho=471/104 satisfies rho<5 and
rho²>20.

If i<floor(k/2), then b-a<=1, so P<5N. The linear coefficient satisfies

    a_i <= [(7891200*5-68209920)/367+26]N <0.

Every pair coefficient is negative, so turning on this spin strictly lowers
F_k in every setting of the other spins.

If i>=floor(k/2), then b-a>=2, so P>20N. Sum the magnitudes of all incident
pair terms. The distance>=2 geometric sums are bounded above by

    55296000(P+N)/(367*471).

There are at most two adjacent terms, bounded by
3317760M/104+36*13^(k-2). Therefore the smallest possible marginal, with
all other spins on, is strictly greater than

  {[(7891200-55296000/471)*20-(68209920+55296000/471)-280320]/367
      -3317760/104-442} N >0.

Both displayed signs are certified rational inequalities. Turning on any
such spin strictly raises F_k in every setting of the others. The assertions
follow by changing coordinates one at a time. No asymptotic exception or
unproved ferromagnetic phase argument is used. ∎

More generally, (2) is exactly reducible to an integer minimum-cut problem.
Write w_ij=-b_ij>0 and double the energy:

    2F=2c+sum_i (2a_i-sum_j w_ij)s_i
                     +sum_(i<j) w_ij|s_i-s_j|.

Positive unary coefficients give arcs to the sink; negative ones give arcs
from the source and their constant offset. Bidirectional edges of capacity
w_ij encode the final sum. The cut identity is proved by the four binary
cases, and `flow.py` implements it. This supplies a polynomial-state exact
method even for variants where the explicit marginal-sign separation fails.
The verifier compares cuts with exhaustive spin evaluation at small k and
with the proved explicit extrema through k100. This is not a claim that
arbitrary prescribed energy membership reduces to a minimum cut.

## 5. A thirteen-adic inverse: the reset separates every spin state

**Theorem.** For fixed k>=3 the map s->F_k(s) is injective. A prescribed integer
can be tested for membership in this physical fiber, and its unique spin word
recovered if present, using polynomially many exact arithmetic operations in k.
In particular the reset fiber has exactly2^(k-2) distinct scalar counts,
although its unreset central shadow has only linearly many count values.

Proof. The k3 one-bit case is direct, with nonzero linear coefficient−165094.
For k>=4 put n=k-2. Recover variables from the two ends. Suppose those at
depths less than t are already known; subtract from F_k(s)-c_k all linear and
pair terms involving only known variables. The residual is divisible by13^t.
Multiply its quotient modulo13 by the inverse of

    lambda_t=8^t*471^(k-t-5) mod13.

Negative exponents at the small boundary are interpreted in the nonzero
multiplicative group modulo13. At t0 the two unknown outermost variables x,y
satisfy the residual equation

    digit=3x+2y-6xy mod13.

The four possible values are0,3,2,12, all distinct. At later two-variable
layers the equation is

    digit=(3-6*s_(n-1))*x+(2-6*s_0)*y mod13.

The coefficient pairs, according to the already known outermost bits, are
(3,2),(10,2),(3,9),(10,9). Their four subset sums are respectively
{0,3,2,5}, {0,10,2,12}, {0,3,9,12}, {0,10,9,6}: all distinct. If a single
central variable remains (t>=1), its coefficient is
5-6*(s_0+s_(n-1)), which is one of5,12,6 and is nonzero modulo13.
Thus every layer has at most one binary solution.

For completeness, these digit equations follow from (2)–(3), not a fitted
residue observation. A distant pair(i,j) has13-adic valuation
n-1-(j-i)=i+(n-1-j). The adjacent coefficient has valuation n-2, with the
same leading residue needed at the k4 and k5 boundaries. At depth t>=1 all
unknown-unknown pair terms have valuation at least2t>t. A boundary-to-unknown
term contributes modulo13 after division only when its known vertex is the
opposite **outermost** variable. The linear term at the current left boundary
has leading coefficient3 and the right2 in the displayed normalization;
the middle term104^(k-3) has higher valuation. For a lone center both linear
contributions survive and sum to5. The constants are explicitly

    (-68209920/367)*471²=3,
      (7891200/367)*471²=2,
                 55296000=6                       modulo13.

Every further layer vanishes in the current digit. These facts prove the
induction. At the end, evaluate the exact integer polynomial and accept only
if it equals the prescribed target; compatible low digits alone are not
acceptance. Formula construction, at mostceil(n/2) digit steps, and a final
quadratic evaluation give a polynomial arithmetic algorithm. ∎

This is a precise passage from a collapsing numerical-semigroup shadow to a
faithful p-adic code under a **physical** reroot operation. It does not use
p-adic completion points as physical tournaments. The inverse decides one
specified family; it does not assert that arbitrary odd targets are accepted.

## 6. Why arithmetic escape is not a renewable selector

For the all-A word the old reset theorem gives
F_k(0)=4*3^(k-4) mod13 for k>=4. The original count is divisible by13^(k-2).
Thus the construction genuinely removes arbitrarily much13-adic content.
But the semigroup coordinate from the old central fiber becomes a long-range
quadratic, depending on the **positions** of excursions, not just their
number and total length. At k3 the words AAA and BCA give reset counts
16,801,195 and16,636,101, rather than a unit step. The inverse above now
proves that no two spin states collide at fixed k, but that is separation,
not a small spacing guarantee. Their new-root deletion
counts are66,689 and66,695, so they no longer share a receiver.

The whole family is still power-sparse. For each k it has at most2^(k-2)
outputs, and every output is at least13^k: the original tournament has an
induced ordinal sum of k copies of D and its path count cannot decrease
under vertex insertion. One elementary injection inserts the new vertex
at the first valid Hamilton-path location in each old path; deletion of
that vertex recovers the input path. The left triangle only increases count
by2V_r. Hence all these outputs below X number O(X^(log(2)/log(13))).

This is a sharp conceptual separation: rerooting breaks an exponentially
large common divisor but leaves an exponentially large **height** and loses
the renewable interface. Exact optimization of a physical spectrum does not
imply interpolation of every integer between its extrema. Chan–Pak remains
open; a target-bearing return map after reset is still missing.

## 7. Reproduction and exploration boundary

Run `python3 verify.py` and `python3 -O verify.py` from this directory (or use
the repository-relative script path). `core.py` and `transfer.py` restate the
inherited literal definitions locally; no mutable project imports or external
packages are needed. The normal/optimized receipts must agree.
Initial probes found gcd2 among the reset trade coefficients at bounded k;
that observation is **not** promoted to an all-order gcd theorem or coverage.
The exact coefficient and extremum theorems replace that weak numerical
signal with an auditable mechanism.
