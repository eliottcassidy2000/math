# The complete row-sixteen repair quotient and the repeated weight-24 payer

**Status: PROVED IN THIS NOTE.** The row-sixteen bracket continuation is
conditional, not an actual row-sixteen lift of THM-4438's source family.
The later [all-row depth theorem](ALL_ROW_DEPTH.md) proves the exact
homogeneous depth kernels and source-weight staircase at every row in
its stated range; finite exact matrices provide independent controls. JC(2), DC(2), chart entry,
polynomial termination, and Chan–Pak remain open. No independent canon
promotion or literature-priority claim is made.

## 1. What advances

The previous session specified an uncomputed row-sixteen quotient. We compute
that **entire response map**, without reconstructing the inherited nonlinear
source debts. There are seven independent bracket debts, only one of which
can be paid by any prefix-preserving valuation-sixteen Hamiltonian. Its
primitive response observable is

    4495 s0 + 870 s2 + 540 s4 + 600 s6 + 1008 s8.             (1)

The other six debts have canonical representatives: the coefficients of
x^11,...,x^16 in the bracket discrepancy. In particular the monomial y^8 is a
visible payer of weight 24. Row fifteen used p^3 y^6, also of weight 24:
one need not raise the weight to purchase the *next scalar* compatibility.
But no combination of the nine new source coefficients can touch the six
other debts. This is a precise separation of purchased response from
remaining memory, not a dimension-count prediction.

We also compute the homogeneous depth selector at row sixteen. The
[all-row continuation](ALL_ROW_DEPTH.md) now proves its general form:
for every n>=3 the homogeneous terminal tangent kernel has degree bound
floor(n/2)+2, hence dimension floor(n/2)+3. Whenever
its inhomogeneous compatibility equations are satisfied, its terminal
fiber is **A^11**, not an extrapolated A^10. Its existence is not proved.

The next actual-source task has therefore been reduced to six explicit
high-degree residuals and the inhomogeneous depth obstruction after an
explicit scalar payment. It is not necessary to repeat a symbolic
left-nullspace search for the row-sixteen bracket.

## 2. Inheritance, concept board, and source boundaries

The closest mechanisms are the source-normal Student operator from
THM-4308/THM-4315, the full associated-graded source diagonal of THM-4415,
and the restored-memory conic and rational boundary of THM-4426/THM-4438.
The hostile is THM-4415's two active debts with only one same-row response.
The corrected near miss is treating an unpaid same-row class as a genuine
geometric impossibility before restoring earlier memory. The underused
sidecar is the **homogeneous depth kernel**, independent of the very large
inhomogeneous formulas.

The live concept board is: Student moments; top-row depth kernels; complete
source filtrations; memory versus repair cokernels; and positive operation
cones on the Chan–Pak interface. The anchor is JC row16, the niche is the
exact truncated Stein quotient, and the wildcard is positive-cone reachability
rather than linear-span reachability.

Source pins are in SOURCES.md. The particular inherited data used are:

* p=t(1+x^2 t), y=x t p, q=-(x^2+6)/2;
* A0=1+x^2/4 and C0=-3x/4-x^3/8;
* changing the last tangent by theta changes the next bracket by
  L_m(theta)=q' theta-q theta'/m;
* the projected depth spaces P_d are generated as described in §4;
* THM-4438 supplies an actual row15 solution family over its named G_m,
  with seven free valuation15 response parameters and ten terminal
  tangent parameters. We inherit that theorem in its stated scope.

Session253, moments/PROOF.md §8, specifically requests the complete
valuation16 packet and its full quotient map. Its §4 centralizer bridge
also motivates retaining a lawful operation carrier rather than mistaking
abstract linear reachability for positive tournament reachability.

## 3. A general truncated Student quotient theorem

Let K be any field of characteristic zero, let m>=1 and 0<=d<m be integers,
and write K[x]_{<=e} for polynomials of degree at most e. Define

    L_m(f) = -x f + (x^2+6) f'/(2m),
    mu_0=1, mu_(2j+1)=0,
    mu_(2j)= product_(i=1)^j 6(2i-1)/(2m-2i+1).

For N>=d+1, define Q_(m,d,N):K[x]_{<=N} -> K^(N-d) by

    Q(F)=(sum_(j=0)^(d+1) mu_j [x^j]F,
          [x^(d+2)]F,...,[x^N]F).                         (2)

Then the sequence

    0 -> K[x]_{<=d} --L_m--> K[x]_{<=N} --Q--> K^(N-d) -> 0

is exact. In particular (2) is a complete, not merely necessary, obstruction
ledger for solving L_m(theta)=F with degree(theta)<=d.

**Proof.** On a monomial,

    L_m(x^j)=(j/(2m)-1)x^(j+1)+(3j/m)x^(j-1).             (3)

Its leading coefficient is nonzero for 0<=j<=d<m, so L_m is
injective. The defining recurrence is

    mu_(j+1)=6j mu_(j-1)/(2m-j).

It makes the first coordinate of Q vanish on (3); all higher coordinates
vanish by degree. The coordinates of Q are independent: the first is
nonzero on 1 and the others select distinct higher monomials. Its kernel
has dimension d+1, equal to the dimension of the image of L_m. QED.

No integration or ordered field is needed for this argument. Over R the
same mu_j are moments of the normalized density (x^2+6)^(-(m+1)) wherever
these moments converge; the proof above only uses rational recurrences.

**General source response.** If a source packet has leading coefficients
1,x,...,x^k with k<=d+1, its image under (2) is exactly one-dimensional,
with observable sum mu_j s_j. It is not legitimate to infer that the full
cokernel is one-dimensional: its dimension is N-d.

### All-row weight staircase, after the depth continuation

[ALL_ROW_DEPTH.md](ALL_ROW_DEPTH.md) proves the missing homogeneous
hypothesis: for every row m>=5 the preceding terminal tangent kernel is
K[x]_{<=d_m}, where d_m=floor((m-1)/2)+2. The complete valuation-m source
packet has leading monomials x^b, 0<=b<=floor(m/2), all within the low
part of this quotient. Its visible b are exactly the even ones. Thus

    w_min(m)=2m-2 floor(m/4).                              (3a)

For k>=2, w_min(4k-1)=w_min(4k)=6k. The corresponding payers are
p^3 y^(2k-2) and y^(2k). This explains the repeated weight24 payment at
rows15/16 as part of an all-row parity staircase, not an accidental
coincidence. The full bracket quotient dimension is floor(m/2)-1,
so floor(m/2)-2 other debts are untouched by the entire same-row packet.
These homogeneous statements do not assert existence of compatible
prefixes or vanishing of any later inhomogeneous debt.

## 4. Exact terminal depth kernels at rows fifteen and sixteen

This finite rational computation is independent of the inhomogeneous
source boundary. Here is the complete universe, so it can be reproduced
without importing an inherited audit script.

For d in {2,3} and n in {15,16}, use coordinates

    (i,r), 0<=i<=n, 0<=r<=i+d.

The projected depth module is the column span of the polynomials

    x^a u^b p^c y^e mod t^(n+1),
    u=x^2 t,  a,b,c,e>=0, a+b<=d, b+c+2e<=n.

A column (a,b,c,e) has entries binomial(c+e,j) at

    (i,r)=(b+c+2e+j, a+2b+e+2j), i<=n.                    (4)

The top-row tangent increment is t^n theta A0' in the d=2 module and
t^n theta C0' in the d=3 module. For theta of degree <=n, apply every
left-nullspace vector of the two column matrices to this increment,
and stack the outputs as D_n. This tests membership in both projected
depth modules, preserving every earlier row.

The exact matrices in probe.py give:

    n=15: D_15 has shape 91 x 16, rank 6,
          kernel = span(1,x,...,x^9).
    n=16: D_16 has shape 105 x 17, rank 6,
          kernel = span(1,x,...,x^10).                     (5)

The frozen output gives dimensions of both generating matrices, the full
stacked shape for n=16, and one nonzero 6x6 minor for each selector. The
script verifies that the first ten, respectively eleven, columns vanish
identically and that the last six columns are independent. These facts
prove (5) over Q and hence over every characteristic-zero field. This is
an exhaustive finite-dimensional rational linear-algebra certificate,
not a numerical rank estimate or a scan over boundary parameters.

Since an affine solution space of a fixed linear map is either empty or
a translate of its kernel, a nonempty row15 depth fiber has the first
kernel in (5), and a nonempty row16 depth fiber is A^11. None of this
proves the inhomogeneous row16 fiber is nonempty.

## 5. The row-sixteen theorem, with all hypotheses

Fix any of THM-4438's admissible row15 source parameters over K, including
a point of its G_m and its seven free relative response coefficients.
Choose any origin in the resulting nonempty row15 terminal depth fiber.
Let F(x), of degree <=16, be the row16 discrepancy

    inherited source coefficient - predicted bracket coefficient

at this origin, before adding valuation16 source terms. The inherited
source-normal bracket formula guarantees this degree bound. The same
following statement applies to any source-normal prefix with these
properties; no special identity of the G_m debt is used.

Add the **complete** prefix-preserving valuation16 source packet

    R16=sum_(b=0)^8 s_b p^(16-2b) y^b.                     (6)

A monomial p^a y^b has t-valuation a+2b and leading coefficient x^b.
Distinct b have distinct leading monomials, so no cancellation can hide
an earlier valuation. Thus (6) is exhaustive modulo valuation>16.

By (5), every terminal15 variation is theta of degree <=9. Its effect on
the predicted bracket is L_16(theta). Therefore the bracket problem is
exactly

    L_16(theta)=F+sum_(b=0)^8 s_b x^b.                    (7)

Apply §3 with m=N=16 and d=9. Write F_j=[x^j]F. Equation (7) has a
solution if and only if

    F_11=F_12=...=F_16=0,                                 (8)
    ell(F)+s0+(6/31)s2+(108/899)s4
          +(120/899)s6+(1008/4495)s8=0,                   (9)

where ell(F)=sum_(j=0)^10 mu_j F_j, and

    (mu0,mu2,mu4,mu6,mu8,mu10)
      =(1,6/31,108/899,120/899,1008/4495,54432/103385).

This proves (1). It also shows that changing the chosen origin does not
change any of the seven ledger coordinates: the discrepancy changes by
an element of im L_16.

If (8) holds, choose all s_b=0 except

    s8=-(4495/1008) ell(F).                               (10)

There is then a unique theta solving (7). More generally (9) cuts the
nine-dimensional source packet to an affine space A^8, with theta
uniquely determined for each source choice. These are bracket solution
fibers, not fibers after imposing row16 depth.

**Sharp weight.** On (6), the source weight is 2a+3b=32-b. The smallest
weight available is 24 at b=8, and mu8 is nonzero. Thus y^8 is the least
weight visible payer inside the frozen-prefix packet. It also proves a
successive scalar-control statement: after the row15 weight24 payment
from THM-4438, a second independent weight24 coefficient can adjust the
row16 scalar without disturbing any row<=15. It does not repair (8).

**Hostile and positive controls.** Each F=x^k for 11<=k<=16 fails (8)
and is unpayable by the entire source packet, regardless of scalar
payment. Conversely every F of degree<=10 is paid by (10), with a unique
tangent. The script checks a target with all coefficients 1,...,11
nonzero, as well as all six individual high-degree hostiles. Odd source
coefficients, including s7, are invisible to (9), but may still affect
later rows or depth. Invisible here does not mean globally irrelevant.

## 6. What remains on the actual boundary, and how to compute it

The result is deliberately a response theorem. We have **not** computed
F_11,...,F_16 on the actual rational G_m family. They might vanish,
reduce to inherited equations, or impose genuinely new constraints.
No nonvanishing or vanishing claim about them is made.

A minimal next computation is now completely specified:

1. Retain all seven free row15 response coefficients, the rational G_m
   parameter, and the consumed lower-memory intercepts from THM-4438.
2. Take the terminal15 origin with its ten free low theta coefficients
   zero and its six selected high coefficients as solved by depth.
3. Form F=G16-predicted_G16. Record the six polynomials F_11,...,F_16,
   with denominators and the exact parameter domain. This is the entire
   unpayable bracket ledger, not a subset of raw nullspace tests.
4. If these vanish on a specified locus, use (9) with eight free source
   parameters, solve the unique ten-coordinate tangent, append the
   particular determinant row16 and all seventeen new theta16 variables.
5. Test all inhomogeneous projected depth residuals. If compatible, the
   terminal depth fiber is automatically A^11 by (5). Retain every new
   source equation; do not replace failure by a claim about JC(2).

This is a smaller and more transparent stopping criterion than treating
all seventeen bracket equations as a black box. If a debt survives,
restore earlier source memory before claiming a geometric obstruction.

## 7. Cone-sensitive transfer to Chan–Pak: a useful boundary, not a lift

The concrete transferable mechanism is **quotient first, then the legal
operation cone**. In the JC packet the coefficients s_b are unrestricted
in K. Over R, if instead they were constrained to be nonnegative, the
response observable in (9) has image exactly [0,infinity): every even
mu_j is positive, odd coordinates are invisible, and s0 realizes every
nonnegative value. Thus a residual scalar ell(F)>0 cannot be paid by
nonnegative source additions, even though the unrestricted response map
is surjective. This assertion follows directly from (9), not analogy.

The relevant Chan–Pak algebra from Session252/253 has coordinates
(d,e,o,i,j) represented by [[d,o,j],[0,e,i],[0,0,d]]. For two such matrices,
the central coordinate of a product is

    d*j' + o*i' + j*d'.                                   (11)

On the nonnegative interface cone every term in (11) is nonnegative.
Therefore replacing a lawful positive operation family by its vector
span can lose a one-sided obstruction just as it does in (9). A proposed
Chan–Pak repair should retain its exact cut-mask/source realization and
its attainable central increments, rather than relying on the existence
of the central E13 direction in the abstract algebra.

Source: the JC finite-row map (7) and the actual Chan–Pak interface
multiplication. Target: a test on a proposed positive repair family.
Map: send each legal repair to its compatibility increment (JC) or
central interface increment (Chan–Pak). Preserved: the increment and its
sign when the operation itself preserves it. Lost: source realizability,
noncentral coordinates, nonlinear cross effects, and absolute Hamilton
multiplicity. Required sidecars: actual tournament/cut-mask operators,
unchanged target coordinates, and reachable increments. Cheapest test:
find a legal positive repair with the needed sign before using a
linear-span certificate. There is no claimed map from a JC Hamiltonian
to a tournament and no new Chan–Pak coverage theorem here.

## 8. Reproduction and validity boundaries

From the repository root run:

    python3 05-knowledge/results/session254_four_frontiers_20261111/jacobian/probe.py

SymPy exact rational arithmetic is the only nonstandard dependency.
probe.out freezes the results. The finite depth checks are exhaustive
for their explicitly defined matrices. The general quotient theorem
has a proof in §3, supplemented by exact controls for m=2,...,20 and
0<=d<m; those finite tests do not supply its quantifiers.

The checks concern new response/depth-kernel claims, not an audit of
inherited theorems. We use the repository meta-pattern “Compute the repair
quotient before testing the residual defect.” The result shows why
same-row source abundance and fixed-weight scalar payments must retain
an independent unpayable-debt ledger.
