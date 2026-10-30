# Random selectors: coupon collection plus a disappearing orientation phase

**Status:** PROVED, conditional only on the explicitly inherited universal-group
classification U1/U2/U7; independent finite checks below. No claim about the
full automorphism group of an arbitrary single tournament. No priority claim.

## Inheritance and the new question

The source is `universal_ramified_closure_20261030/PROOFS.md`, U1/U2/U4/U6/U7,
at main commit `b51ef5cd739375e6475ceafd42e36e86858c70ca`. Fix prime `p>=5`,
`k>=2`, and `v_p(D)=1`, and use

    R=(Z/p^k)[delta]/(delta^2-D), q=p^k, N=q^2,
    rho(x+y delta)=x mod p.

U1 identifies the full directed-color group H, of order `p^(3k+p-2)`.
U7 identifies the full unsigned-color group `W=H semidirect C2^2`.
U2 says a named pin set meeting `j>=2` residue fibres has H-stabilizer
`C_p^(p-j)`. Its base condition is exactly meeting all p fibres.
U7 also exhibits a tilted p-point W-base, while the integer transversal
retains conjugation. These are inherited results, not proved again here.

The new question is quantitative and algorithmic: how frequently does an
ordinary random transversal accidentally retain the unsigned phase, and how
much longer does random point selection take? The answer distinguishes a
coupon-collection cost from a depth-sensitive algebraic exceptional set.

Closest mechanism: the U2 occupancy stabilizer. Hostile: U7's integer
transversal. Corrected near miss: meeting every fibre kills H but need not
kill W. Least-used sidecar: the actual delta coordinate, not just occupancy.
The concept board is occupancy / shared orientation / affine equations /
random stopping time / exact counting.

## S1. The reflection axes are bent affine lines

Write `Q=p^(k-1)`. The fixed sets of elements in the coset Hc, where c is
ring conjugation, are either empty or exactly a set

    L(s,b,h) = {(x, s*x+b+Q*h(x mod p)) : x in Z/q},       (S1)

where `s,b` run over residues modulo Q and `h:F_p -> F_p` is arbitrary.
Every set (S1) is the fixed set of exactly one element of Hc. Consequently
there are exactly

    A = Q^2 p^p = p^(2k+p-2)                              (S2)

such reflections, and each fixes q vertices, q/p in each residue fibre.

**Proof.** An Hc element has the form

    f(z) = (u+v delta) conjugate(z) + beta + eta*h(rho(z)),

where `(u+v delta)` is a principal norm-one unit and eta spans the top
socle `Q delta F_p`. In particular u=1 mod p, so u+1 is invertible.
The two fixed-point equations, writing beta=(B_x,B_y), give

    (u-1)x-Dv*y+B_x=0,
    v*x-(u+1)y+B_y+socle(rho(x))=0.

Solve the second for y. The coefficient of x in the first then vanishes
because `u^2-Dv^2=1`; a socle term is annihilated by D. Thus either
`B_x=Dv*B_y/(u+1)` and every x has one fixed y, or no point is fixed.
The resulting graph is (S1) with slope `v/(u+1)` and intercept
`B_y/(u+1)`, both initially modulo q. Their highest digits can be absorbed
into h, leaving s,b modulo Q.

Conversely each slope s is realized by the norm-one unit
`lambda=(1+s delta)/(1-s delta)`: its denominator is a unit, and lambda
has residue1, hence is a clock power by U1. The reflection
`z -> lambda conjugate(z)` fixes `y=s*x`. Translating the axis by b delta
and adding an appropriate top-socle fibre shift realizes every h. Its
square is the identity, as it fixes every point on this graph and acts by
reflection in the transverse coordinate (alternatively, its square lies
in H and fixes a transversal, hence is identity by U2).

The parametrization of axes is injective: equality at all x=0 mod p and
x=p forces the slopes equal modulo Q, then the intercepts equal modulo Q,
then each h value equal. More explicitly, subtract the two graph functions
at x and x+p to get `(s-s')p=0 mod q`; the chosen residues force s=s'.
Their equality at x=0 then gives b=b' modulo Q, and equality at a point of
each residue gives h=h'. Two Hc elements fixing the same axis have quotient
in H fixing all fibres, so they are identical. This proves S1-S2. QED.

These are genuinely coupled axes, not p unrelated choices of local
reflection. The parameters s and b are shared across every residue fibre.
The nonlinear top digit h is the remaining fibrewise freedom.

## S2. Exact orientation-phase survival under random pins

Draw an ordered r-tuple of vertices independently and uniformly from R,
allowing repeats. Let J be the number of residue fibres visited. Conditional
on `J=j>=2`, the probability that the W-stabilizer has an extra reflection is

    P(extra phase | J=j) = p^(2k+j-2-kr).                 (S3)

The stabilizer itself is `C_p^(p-j)` if no phase survives and
`C_p^(p-j) semidirect C2` otherwise, where the involution acts by inversion.
Its orders are respectively `p^(p-j)` and `2 p^(p-j)`. Thus its distribution
is exact, not merely a union bound. In particular, conditioned on a minimal transversal (r=j=p),

    P(not a W-base) = p^(-(k-1)(p-2)).                    (S4)

For p=5,k=2, exactly 1/125 of random transversals retain the phase.
The integer transversal is in this exceptional set; the tilted source
transversal is outside it.

**Proof.** An element in either of the two W/H cosets involving global
reversal acts on residues as `a -> -a+t`; it cannot fix two distinct
residues. Only H and Hc can fix such a tuple. The H-stabilizer has size
`p^(p-j)` by U2. If any Hc element fixes the tuple, the entire set of such
elements is one coset of that stabilizer and has the same size.
Each such element is a reflection by S1. Its residue action is the identity,
and conjugation negates the top socle, so it inverts every remaining
unoccupied-fibre shift. This gives the stated semidirect-product type.

The total number of ordered tuples with J=j is

    M_j = (p)_j S(r,j) (N/p)^r,

with Stirling number S(r,j). Each of the A reflection axes has q/p points
in each fibre and contains `(p)_j S(r,j)(q/p)^r` such tuples. Double count
pairs (reflection, fixed tuple), dividing by `p^(p-j)` because every bad
tuple has exactly that many fixing reflections. Dividing again by M_j
leaves `A*p^(-(p-j))*q^(-r)`, which is (S3). QED.

For example, after two different-residue pins every pair retains a phase:
S3 equals1 for r=j=2. A third pin in a new fibre retains it with probability
`p^(1-k)`; a third pin in an already visited fibre retains it with probability
`p^(-k)`. An additional point in an old fibre is a stronger phase test,
although it does not advance coupon collection. This is a precise tradeoff
between two purposes of a selector, not a claim that fewer than p pins form
a base.

## S3. Exact base probability and expected stopping time

Let `C_r=p! S(r,p)/p^r` be the probability that r independent residue draws
cover all p fibres (C_r=0 if r<p). Then, for r>=p,

    P(H-base after r draws) = C_r,
    P(W-base after r draws) = C_r [1-p^(2k+p-2-kr)].      (S5)

The second formula is zero for r<p by definition, not by multiplying an
undefined or inadmissible probability factor.

Let T_H and T_W be the first times the collected pins are respective bases.
Then

    E[T_H] = p H_p,
    E[T_W] = p H_p + Delta_(p,k),                        (S6)

where

    Delta_(p,k) =
      p! p^(-k(p-2)-2) / product_(j=1..p)(1-j/p^(k+1)). (S7)

**Proof.** The H-base condition is U2's coupon condition. Apply S3 at j=p
for S5. For a positive integer stopping time,
`E[T]=sum_(r>=0)P(T>r)`. Subtracting the tails in S5 gives

    E[T_W]-E[T_H] = p^(2k+p-2) sum_(r>=0) C_r p^(-kr).

Use the ordinary generating function
`sum_r S(r,p)z^r = z^p/product_(j=1..p)(1-jz)` at
`z=p^(-(k+1))` to obtain S7. The familiar directed coupon mean follows
by splitting the wait for each new fibre into geometric waiting times
of means p/(p-j). QED.

Depth k has no effect on the directed selection time, but the extra unsigned
phase cost decreases exponentially with k. This is a structural difference
between point selection and the growing bank-coordinate count: neither
should be substituted for the other.

The random procedure is a Las Vegas selector construction only after a
valid membership/base test is supplied. S5-S7 count its mathematical stopping
time; they do not charge bit operations or assert an improvement in the
factoring paper's algorithm. One explicit base test with all fibres occupied
is: choose two pins of different residues, solve the possible bent-axis
slope and intercept modulo Q. First require y-s*x-b = 0 modulo Q
for every pin, then require (y-s*x-b)/Q modulo p to be constant within
each residue fibre. If both conditions hold the tuple is not a W-base;
otherwise it is. This is exact
and uses the literal coordinates, not an unlabeled tournament alone.

## S3b. A one-move repair of every bad transversal

Given any transversal with one pin in each fibre, apply the two-pin phase
test. If it is bad, replace any one pin z by z+delta. The resulting set is
always a W-base and still has exactly p pins.

Indeed at least two other fibres remain represented. Their difference has
unit real coordinate, so they determine the possible slope and intercept
of every containing axis modulo Q. The original pin satisfied that equation
modulo Q; adding delta changes its second coordinate by1, so it no longer
does. Since Q is divisible by p, this is not a top-socle adjustment. No
reflection can fix the repaired transversal, and occupancy still kills H.

Thus a uniform random transversal can be repaired with at most one pin move,
and the expected number of moves is exactly `p^(-(k-1)(p-2))`. This extends
the source's specific tilted integer transversal to every transversal. It
is a change in chosen selectors, not a claim that existing measurements can
be altered or that the changed pin has the same color observations.

## S4. Sharp depth-one hostile

At k=1 U7 gives the imprimitive unsigned group
`W=(D_(2p)) wr D_(2p)` on p fibres of size p. A named set is a base iff each
fibre contains at least two DISTINCT vertices: an empty or singleton fibre
retains an internal dihedral symmetry, and two distinct points kill its
internal dihedral action. Once all fibres are individually pinned the outer
action is also trivial.

For r iid vertex draws the exact shallow base probability is

    r!/(p^2)^r [z^r] (exp(pz)-p exp(z)+(p-1))^p,          (S8)

or equivalently

    sum_(a+b+c=p) p!/(a!b!c!) (-p)^b (p-1)^c
                      ((ap+b)/p^2)^r.                  (S9)

The EGF for assignments to a single fibre with at least two distinct
vertices is `exp(pz)-1-p(exp(z)-1)`. Taking its pth power and extracting
labelled samples proves S8 and its multinomial expansion S9.
For r=0 take 0^0=1 in S9. In particular p pins never suffice, whereas
S4 gives a positive base probability with p pins at every k>=2.
This is the precise boundary preventing a careless k=1 substitution.

## Verification, source contract and boundaries

Run `python3 check.py`; no third-party dependency is needed. It independently
implements literal conjugation-plus-fibre-shift maps fixing 0 and1 at p=5,
k=2, tests their fixed sets against the graph equation, exhausts all625
third pins and all625^2 ordered fourth-pin choices, and compares both the
number of surviving reflections and the conditional probability with S3.
It exhausts all3125 low-coordinate minimal-transversal patterns, of which25
are bad. It also compares the two-pin phase test with direct clock/conjugation
parameter enumeration in four literal-ring families, including D=10 and
k=3, and tests the one-move repair on bad transversals. It checks exact
rational occupancy, bounds and expectation identities
in nine families p=5,7,11 and k=2,3,4. The shallow wreath counting identity
is exhaustively checked in smaller abstract fibre actions p=2,3,r<=6; those
are tests of S8 only, not assertions that the ramified clock works at p=3.

The source object is the colored ramified ring; the target is a named pin
sample. The preserved predicate is triviality of its pointwise stabilizer.
Occupancy alone loses a reflection phase; the bent-axis residual restores
it. The cheapest hostile is the integer transversal, contrasted with the
single tilted pin in U7. No result here identifies full automorphism groups
of arbitrary individual tournament banks, enumerates all tournaments, or
imports the factoring preprint's complexity theorem as an independently
verified fact. The proof uses only elementary finite counting after the
named inherited group classification.
