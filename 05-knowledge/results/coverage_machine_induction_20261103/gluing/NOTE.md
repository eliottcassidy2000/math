# Cyclic gluing with two old ports: differential closure and exact divisor escape

**PROVED scoped operator and arithmetic theorems; FINITE-EXACT controls.**
These supply a reusable gluing operation and a complete criterion for its
prime-divisor escape. They do not prove prime/7p spectrum coverage.

The closest proved mechanisms are the
[eighteen-polynomial birth carrier](../../source_selection_closure_20261103/GRADED_OBSERVER.md)
and [Session159's unmarked wrapper operator](../../session159_constructions_20261103/profile_operator/NOTE.md).
The hostile is C3: unrestricted iteration of this wrapper never erases its
factor3. The corrected near miss is treating coefficient primitiveness as
target coverage. The least-used coordinate is the number of anonymous path
components joined on both sides of the two new wrapper vertices.

Anchor = a legal source-changing gluing step; niche = exact divisor escape;
wildcard = the modular cellular automaton of cover evaluations. The concept
board is old physical names / anonymous components / ordered cover content /
wrapper histories / modular states / actual source guards.

## G1. An exact wrapper retaining both old names

For a nonempty tournament T define

    G(T) = T -> a -> b -> T.

All old arcs and names are retained. The two new vertices are fresh, a->b,
every old vertex beats a, and b beats every old vertex. The output is strong.
Fix two distinct old physical ports u,v; the output ports are still u,v.

For a spanning unordered directed-path cover, its skeleton s is the multiset
of component endpoint pairs involving u or v; write * for every other
endpoint and remove the (*,*) pairs. If m such anonymous components were
removed, let c_(s,m) count covers with that signature, and set

    P_s(z)=sum_m c_(s,m) z^m.

There are eighteen skeletons, as in the inherited carrier. The two-port
cyclic wrapper has an exact transfer

    P_new = P A(z) + P' B(z) + P'' C(z),                 (G1)

with nonnegative integer polynomial matrices, each of degree at most two.
The receipt gives all312 nonzero entries:156 in A,120 in B,36 in C.
No dimension bound on polynomial degree is imposed.

**Proof and explicit compiler.** Treat each path of an old cover as an
atomic directed component. Any new cover is determined by choosing:

1. zero or one old component to enter a;
2. zero or one old component to follow b;
3. whether the edge a->b is used.

If the same old component is selected in both1 and2 and a->b is used,
the result is a directed cycle and is rejected. Otherwise the selected
arcs join atomic paths into paths, and all unselected components survive.
Deleting a,b recovers the unique old cover and all three choices, proving
the bijection.

At most two distinct anonymous old components are selected. Selecting j
distinct such components contributes (m)_j choices, hence P_s^(j), with
j<=2. All resulting new anonymous components lie in the at most two paths
containing a or b; their number is at most two. This proves both coefficient
positivity and the degree bound.

For a directly reproducible extraction, let K_m(s) be the resulting output
polynomial for one input skeleton with exactly m anonymous components.
Then its row matrices are

    A=K_0,
    B=K_1-z K_0,
    C=(K_2-z^2 K_0-2z B)/2.

The extraction uses subtraction only to identify coefficients; the
selection proof proves that the final coefficients count actual choices.
The script also checks K_m=z^m A+m z^(m-1) B+m(m-1)z^(m-2)C for m=0..6.

Second derivative order is necessary for this free signature kernel:
from an empty skeleton with m anonymous paths, the coefficient with m-1
output components is m(m-1). An order-one operator with nonnegative
polynomial powers cannot produce that quadratic coefficient at z^(m-1).
This is a statement about this kernel representation, not a minimality
claim over arbitrary nonlinear encodings of actual tournaments.

A chronological birth compiler can temporarily retain a third name a:
append the sink a; then append b losing only to a and beating T; then
forget a,b as names. There are88 three-port skeletons,64 independent
endpoint statuses and24 possibilities pairing two distinct ports.
The direct second-order operator avoids maintaining that temporary carrier.
For |T|>=3 a two-birth implementation which keeps both u,v at every stage
and forgets the first newborn cannot make the second newborn distinguish
that forgotten vertex from a forgotten old vertex. This is only a failure
of that immediate two-port compiler; longer native programs are not ruled
out, and |T|=2 does not have the same obstruction.

The old seventeen positive-birth queries are insufficient for G: TT2 and
TT3 with their first two vertices marked have the same inherited observer,
but their wrapped counts are5 and9. The equality is the inherited
ordinal-sink extension law, not a new finite inference.

## G2. The unmarked projection and the right content invariant

Set

    Q_T(z)=sum_s z^|s| P_s(z)=sum_k p_k(T) z^k,
    C_k(T)=k! p_k(T),
    F_T(r)=sum_k p_k(T)(r)_k=sum_k C_k(T) binom(r,k).

Here C_k counts ordered nonempty k-path covers, while F_T counts r ordered
slots allowing empties. Projection of (G1) gives exactly Session159's law

    Q_G=(z+z^2)Q''+(3z+2z^2)Q'+(z+z^2)Q,
    F_G(r)=r^2 F_T(r)+r[F_T(r+1)-F_T(r-1)].           (G2)

Hamilton count is F_T(1), and H(G(T))=F_T(1)+F_T(2).
Projection forgets the endpoint names and paired endpoint data. The
eighteen polynomials and literal source must accompany it when later
operations use those names.

The gcd of all coefficients of the eighteen polynomials is always1:
the cover of all singleton vertices contributes exactly1 in skeleton
{(u,u),(v,v)} at anonymous degree |T|-2. Thus full coefficient content
cannot replace the old seventeen-query continuation divisor.

Define the **fixed slot divisor**

    d_slot(T)=gcd_(r>=1) F_T(r)=gcd_(1<=k<=n) C_k(T),  n=|T|.  (G3)

The equality follows from the unitriangular binomial evaluation matrix at
r=1,...,n, with F_T(0)=0. In particular d_slot divides H(T) and n!.

**Uniform-context theorem.** If T is substituted as one homogeneous block
in any fixed tournament context with fixed other blocks, then

    d_slot(T) divides d_slot(output).                 (G4)

To prove this for any r-slot output cover, cut each slot at every transition
between T and the exterior. Read the maximal T-runs in slot order and in
path order within each slot. They form an ordered k-path cover of T.
For each fixed exterior pattern the attachment directions are uniform, so
either every such ordered cover is admitted or none is. Thus

    F_output(r)=sum_k a_(r,k) C_k(T),  a_(r,k)>=0 integers.

Every output cover has a unique pattern, which proves (G4). This includes
cyclic gluing, the wrapper G, and ordinal extension by arbitrary fixed
tournaments. An isomorphic repositioning of homogeneous module vertices
also preserves d_slot. A general native source change need not.

For example C3 has d_slot=3, F2=12 and wrapped H=15. No sequence of
uniform-context substitutions retaining this block can remove that3.
Yet an indegree-one birth at a chosen C3 root gives H5: a pointed birth
breaks the homogeneous-block condition. Actual H9 and H15 sources in the
receipt also have d_slot3, despite full polynomial coefficient gcd1.
This is why adding more observation coordinates alone is not a reset.

## G3. Exact wrapper content and a bounded prime escape test

Let h_j=H(G^j(T)), and put B_r=(r-1)!F_T(r), r>=1. Equation(G2) becomes

    B'_r=B_(r+1)+r^2 B_r-r(r-1)B_(r-1),              (G5)

where the r=1 last term vanishes. Therefore

    h_j=B_(j+1)+an integer linear combination of B_1,...,B_j.

The first n histories and B_1,...,B_n are related by a unitriangular
integer matrix. Equivalently h_j=sum_(r<=j+1) a_(j,r)F_T(r), with
a_(j,j+1)=j!; some lower coefficients may be negative.

**Exact finite content theorem.**

    d_wrap(T):=gcd_(j>=0) h_j
       =gcd_(0<=j<n) h_j
       =gcd_(1<=r<=n) (r-1)!F_T(r).                  (G6)

Indeed F_T is a degree-n polynomial with F_T(0)=0. The recurrence
Delta^(n+1) F_T=0 expresses every F_T(m), m>n, in terms of the preceding
n+1 samples. Multiplying by (m-1)! makes every coefficient on an earlier
B_r integral, since (m-1)!/(r-1)! is an integer. The first step m=n+1
uses F_T(0)=0. Induction shows that every B_m is an integer combination
of B_1,...,B_n. Combine this with the unitriangular history relation.

**Prime escape criterion.** For every prime ell,

    ell divides every h_j  iff  ell divides d_slot(T).       (G7)

If ell divides d_slot, (G4) proves the forward persistence. Conversely,
the expression F_T(r)=sum_k p_k(T)(r)_k makes F_T an integer-coefficient
polynomial, so F_T(r+ell)=F_T(r) modulo
ell. If ell does not divide d_slot, some r in1,...,ell-1 has F_T(r) nonzero
modulo ell. Choose the first such r. The triangular history formula then
gives h_(r-1)=(r-1)!F_T(r) nonzero modulo ell. Hence an escape is found
within ell-2 wrappers. Equation(G6) also bounds the first escape by n-1.
Thus the bound is min(n-1,ell-2).

In particular, if ell>n, it cannot divide d_slot because C_n=n!.
Every such prime divisor of H(T) is removed by at most n-1 wrappers.
Small primes are the only possible permanent wrapper obstruction.
This is a constructive factor-escape result, not a target-value selector:
the new count need not be prime or equal to any requested integer.

The prime supports of d_wrap and d_slot agree, but their valuations should
not be identified without proof. The formal positive odd polynomial
r^7+4r^5+4r^3+18r has d_slot9 and weighted content27. No tournament
realization of that polynomial is asserted; it is a hostile to a purely
polynomial argument claiming equality of those contents.

## G4. A literal modular cellular automaton

Fix M>=2 and record v_i=F_T(i) modulo M, i in Z/MZ. The same integer
coefficient property makes this a lossless carrier for all slot evaluations
modulo M. The wrapper acts by the radius-one update

    v'_i=i^2 v_i+i(v_(i+1)-v_(i-1)) modulo M.         (G8)

The array has cyclic boundary conditions. It is a spatially periodic linear
cellular automaton. If a translation-uniform rule is wanted, retain the
fixed label i modulo M in each cell and update (i,v_i) by the same local
formula while leaving i unchanged. Legal initial label arrays have adjacent
labels differing by1. Ordinal sum of tournaments becomes pointwise
multiplication of these arrays, by F_(A->B)(r)=F_A(r)F_B(r).

This map preserves all modular slot observations and every iteration of G.
It destroys integer magnitudes, polynomial degree, source isomorphism,
ordered port data, and legal source-changing guards. Retain the actual
tournament and its eighteen polynomials to lift back to exact constructions.
In particular an orbit modulo49 does not choose an exact 7p output. The
fixed-divisor theorem and the counter-carrier collision in the companion
packet show two different reasons that computational richness is not
spectrum coverage.

### G4b. Mod49 reduces exactly to value/derivative registers

The49-cell state admits a smaller exact representation. Put

    V_j=F_T(j) mod49, j=0,...,6,
    D_j=F'_T(j) mod7, j=0,...,6.

V_0=0 for every nonempty source. The remaining six values and seven
derivatives are13 scalar registers, with different moduli. For every
integer a and j in0,...,6, the integer-coefficient Taylor expansion gives

    F_T(7a+j)=V_j+7a D_j mod49.                       (G9)

All higher terms contain49. This reconstructs the entire49-cell state
and applies also to negative a. The exact wrapper update is

    V'_j=j^2 V_j+j[F_T(j+1)-F_T(j-1)] mod49,
    D'_j=2j V_j+j^2 D_j+V_(j+1)-V_(j-1)
          +j[D_(j+1)-D_(j-1)] mod7.                 (G10)

In the value update use(G9) at the boundaries: F_T(7)=7D_0 and
F_T(-1)=V_6-7D_6 modulo49. In the derivative update all indices are
cyclic modulo7 and values are reduced modulo7. This is simply the
derivative of(G2), so it retains the information needed by every future
wrapper, not merely the next Hamilton residue. Ordinal gluing has the
exact product update

    V_j(A->B)=V_j(A)V_j(B) mod49,
    D_j(A->B)=D_j(A)V_j(B)+V_j(A)D_j(B) mod7.         (G11)

The [inherited parity law](../../session159_constructions_20261103/cycle_profiles/NOTE.md)
F_T(-r)=epsilon F_T(r), epsilon=(-1)^|T|, compresses further:

    V_(7-j)=epsilon[V_j-7D_j] mod49,
    D_(7-j)=-epsilon D_j mod7,      j=1,2,3.

For even order D_0=0. Thus, retaining the parity sign, the three V_j and
three D_j at j=1,2,3, together with D_0 in odd order, suffice: at most
seven scalar registers. The wrapper preserves parity, while ordinal sum
multiplies the two parity signs. No minimality claim is made.

These registers resolve the entire mod49 wrapper calculation, including
whether a displayed count is divisible by7 but not49. They do not encode
exact magnitude, primality of H/7, old port names, or availability of a
native reset. Their dimension is unrelated to the seventeen exact
two-port continuation queries. A successful modular orbit still needs
the actual integer source and a target-directed construction.

## G5. Reproduction and boundary

Run `python3 -I -B check.py` and `python3 -I -B -O check.py`; stdout is
`check.json`, identical in both runs. No repository code is imported.
All checks use explicit runtime conditions.

The universe is all74 labelled tournaments of orders2,3,4 with named
ports0,1. For each, literal permutation/cut enumeration gives the complete
marked path-cover vector before and after G; a separate subset-endpoint DP
recounts H. Additional checks cover126 formal signatures with m=0..6,
all unmarked projections,4514 cellular-automaton cell updates at
M=2,3,7,49,370 prime-escape controls, twelve wrapper histories per source,
ten small nontrivial uniform cyclic contexts, and the three explicit
fixed-divisor hostiles. There are296 mod49 jet wrapper/49-cell/parity
reconstructions and74 literal ordinal-product jet controls. The finite
checks support the all-order bijections
and algebra; they are not a scan for a larger spectrum interval.

The constructive induction now has an actual source-changing G step that
retains both physical ports, a precise interface cost, a finite criterion
for every prime-factor obstruction under its iterations, and a modular
machine representation. It still needs a guarded choice among such steps
that lands on each requested prime/7p target and leaves a reusable source.
