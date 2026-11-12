# Central count intervals, symmetric deletion cores, and literal integer inputs

**Session276. PROVED from the inherited rooted fusion law, with FINITE-EXACT
literal seeds. Chan–Pak remains OPEN.** General central-matrix transport is
inherited, not a novelty claim: see Session252 and
`chanpak_five_core_bridges_20261110/algebra/FIVE_MODULES.md` (its B95 switch).
The regular-C5 balanced count set `{41,43}` is already recorded in
THM-4123, `THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice.md`.
Here it is made into a small complete rooted central packet, a symmetric-core
supplier, and an exact positive target consumer with unrestricted integer inputs.

## 1. A systematic supplier of central packets

Let U be a vertex-transitive tournament of odd order n, let h=H(U), and let
Q=sum_{S subset V(U)} H(U[S])H(U[V\S]), with H(empty)=1. Add root r and let
A be its outgoing set, of fixed size a. For the rooted fusion coordinates of
Session264, every choice of A has

    root-deletion count h, two-cover coordinate Q, k=0,
    U_r = h+a(Q-h)/n,             V_r = h+(n-a)(Q-h)/n,
    d_r = h-ah/n,                 e_r = ah/n.                  (1)

Only the Hamiltonian count H(r+A) depends on the arrangement of A. Thus two
such cuts with counts differing by delta have matrix difference delta E14.
Any consecutive odd run in this Johnson-layer count image is an actual
central interval bank, not merely a scalar interval with unrelated roots.

Proof. Write s_i(S) and t_i(S) for paths in U[S] starting and ending at i.
A path starting at r first takes an outgoing neighbor i, then a path in U[S];
its complement contributes H(U\S). Its contribution is therefore linear in
the incidence vector of A, with coefficient

    alpha_i=sum_{nonempty S} s_i(S) H(U\S).

Vertex transitivity makes all alpha_i equal; their sum is Q-h. The analogous
ending coefficients have the same sum and are also constant. For alternating
path-length weights the coefficient is

    gamma_i=sum_{nonempty S} (-1)^|S| s_i(S) H(U\S).

Their total is -h, since the full alternating two-cover sum is zero by
complement cancellation at odd n. This proves (1), including the root-alone
term h. All coordinates are integers because they count actual objects.
The singleton deletion n=1 obeys the same identities; it supplies no
nonconstant count bank. No symmetry converse is asserted.

## 2. Small literal central banks and a content-free one

Root is 0 in outgoing-mask notation. Two six-vertex inputs are

    C41=(56,13,25,48,34,6), C43=(52,13,24,49,34,6).
    (h,Q,U,V,d,e)=(15,120,78,57,6,9).

Their full matrices differ by exactly 2E14. Reversing the chosen cut sizes
interchanges U,V and d,e and retains the same two counts.

For the cyclic seven-vertex deletion with steps {1,2,3}, h=175 and Q=1988;
its outgoing-size-two menu is `{559,561,575}`. The first two are again a
central +2 pair. The deletion's common coefficient (Q-h)/7 is 259.

The order-nine circulant with steps {1,7,6,5} (code1 in `central.py`) has

    h=3357, Q=50004, (Q-h)/9=5183, gcd(3357,5183)=1.

Its outgoing-size-four menu contains `{15381,15383}`. This central pair has

    (h,Q,U,V,d,e)=(3357,50004,24089,29272,1865,1492).

It has no common content in h and the symmetric endpoint coefficient. That
removes the simple content obstruction of the C5/C7 examples; it does not
prove any interval overlap or congruence universality.

The exact cyclic-core experiment covers every binary orientation choice for
steps {±1,...,±(n-1)/2}, for n=3,5,7,9,11, and every root cut. These are
circulant presentations, NOT all tournament classes. The longest consecutive
odd run in any fixed-size cut menu has length2 for n=5,7,9 and length1 for
n=3,11. Thus more vertices/symmetry do not automatically thicken this bank.
No global nonexistence of a longer central interval is claimed.

## 3. Reusable interval transport and the exact amplification invoice

If a central bank has H-values b+2{0,...,L}, and fixed rooted fusion contexts
P,R have root-deletion counts h_P,h_R, its image is exactly an arithmetic run
with spacing 2h_P h_R and length L+1. Every other rooted coordinate agrees
across the entire image. This is the inherited identity

    M(P) E14 M(R)=h_P h_R E14.                              (2)

In particular pointed contexts whose deletion is transitive have h_P=h_R=1,
so they translate the central interval without changing its step2 or width.
There is no finite fuel restriction: any finite word of such contexts is
allowed, of arbitrary length.

The exact catch is important. Fusing m identical C41/C43 banks gives only

    B_m+2*15^(m-1){0,...,m},

not a longer step2 interval. The m controls commute in the central direction,
but each pays the deletion count of every other core. Arbitrarily many unit
contexts cannot remove that invoice. This is an instance of the already known
central-copy spacing law, not a new no-go for other operations.

## 4. An exact positive target consumer with unrestricted integer inputs

Define f on the nonnegative integers by

    f(0)=1,
    f(2a+1)=f(a),
    f(2a)=f(a)+a+(a mod2)       for a>=1.                  (3)

It is the Hamiltonian count of a rooted tournament L(A) whose deletion is
transitive and whose outgoing incidence word is the binary expansion of A.
Leading zero bits add neutral initial sources and do not change f. Its rooted
coordinates include h=1, U=1+A, and d=1-(A mod2). Formula (3) follows by
appending one old vertex to the transitive deletion: an outgoing final arc
adds no Hamiltonian paths, while an incoming final arc adds U-d. Thus (3)
is an exact inherited one-apex digital count, not a fit or heuristic.

For an arbitrary pointed nontrivial core C, place L(A) on the left and the
converse of L(B) on the right using root fusion. Write

    a=U_C-h, b=V_C-h, c=d_C-h, e=e_C-h.

Direct multiplication of the inherited carrier gives, for EVERY A,B>=0,

    F_C(A,B)=H_C+h(f(A)+f(B)-2)+aB+bA+QAB
                    -c(B mod2)-e(A mod2).                (4)

The empty binary word means the rooted singleton identity. At A=0 or B=0
one may use that identity rather than a padded nontrivial zero word; (4)
still holds because its associated first-row or last-column entries vanish.
It is crucial that A,B are arbitrary nonnegative integers: they are not
members of the sparse one-way digit language of Session275. On the other
hand f(A),f(B) are paid digital corrections, not removable counts.

For a symmetric deletion core with 0<a_out<n, all correction terms of (4)
are nonnegative: h>0, a,b,Q>0 and c,e<0. Consequently a target N in a central
bank starting at H0 has a finite, exact inverse search

    0<=A<=(N-H0)/b,  0<=B<=(N-H0)/a,  QAB<=N-H0.          (5)

For every pair in this explicit hyperbolic domain compute the translation in
(4), then test whether N minus that translation belongs to the central bank.
A positive test constructs literal vertices and arcs. A negative test rejects
only this core-bank sandwich family. It is not a test of the full spectrum.
The integer bounds supply absolute anchoring absent from an argument about
residues or count differences alone. They do not establish successful inversion.

For C41/C43, (4) is the explicit union of pairs

    {41,43}+15(f(A)+f(B)-2)+63B+42A+120AB
                                      +9(B mod2)+6(A mod2).

Every translation here is divisible6. Hence this entire unrestricted A,B
family still omits multiples of3. It cannot by itself close Chan–Pak.

For the C9 pair, (4) becomes

    {15381,15383}+3357(f(A)+f(B)-2)+20732B+25915A+50004AB
                                       +1492(B mod2)+1865(A mod2).

The earlier content3 obstruction no longer follows. Exact target inversion,
not the gcd, is the remaining question.

## 5. Exact lattice content and a stronger local sieve

The gcd of ALL translations in (4), over arbitrary A,B>=0, is exactly

    G=gcd(2h,2a,2b,a-c,b-e,Q).                              (7)

Proof. Rewrite each one-variable contribution using A=2q+epsilon. Since
f(A) is odd, every term is divisible by the displayed gcd. Conversely the
translation values at (1,0),(2,0),(3,0),(0,1),(0,2),(0,3),(1,1) generate
b-e,2b,2h,a-c,2a and Q by integer differences. Thus no smaller lattice
contains every translation. This use of integer differences establishes
the gcd, NOT physical subtraction or exact interval filling.

For C41, G=6. For the C7 outgoing-size-two core, G=2, but its two-point
central bank STILL has a global obstruction modulo7. Since
h=175,Q=1988,a=518,b=1295,c=-50,e=-125, (4) reduces modulo7 to

    F_C(A,B)=H_C+(B mod2)-(A mod2).

For H_C in {559,561}, the only residues are {0,1,2,5,6}; residues3,4 never
occur at any word length. This is a useful exact hostile example against
the implication “primitive spacing plus a central +2 pair gives coverage.”

The C9 outgoing-size-six pair14207/14209 also retains content3, despite
gcd(h,(Q-h)/9)=1. The cut size itself is a required sidecar. In contrast
the outgoing-size-four pair15381/15383 has G=2 and gcd(a,h)=1. Taking A=0
and B=2t realizes EVERY odd residue modulo2h=6714: the f term vanishes
modulo2h, leaving H_C+2at. This is an all-length exact modular statement,
not exact coverage of any prescribed integer. It explains why the
size-four source is strictly better than the superficially similar
size-six source for this consumer.

More generally, modulo any odd divisor m of h, (4) reduces to an explicit
quadratic polynomial in A,B plus their parities. Its exact infinite-input
image is obtained from A,B modulo2m, or by taking unions of linear cosets
with each parity fixed. This gives a terminating local obstruction sieve
without enumerating ever larger tournament words. Passing that sieve does
not establish integer realization.

## 6. Carries are a resource, not just noise

The one-way Session275 bank erased mixed carries and consequently produced
superincreasing holes. Formula (4) retains a lawful digital correction while
making its unrestricted integer arguments explicit. In particular

    f(2^k-2)=2^(k-1)+1,   f(2^k-1)=1       (k>=2).

Increasing A across that binary carry resets f by 2^(k-1). In the sandwich,
for fixed B, the two neighboring arguments therefore differ by

    Q B+b - h*2^(k-1) - e.                               (6)

The last term uses that A changes from even to odd. This difference can have
either sign and arbitrarily large magnitude. Thus scalar order is not the
binary argument order; monotone envelope claims would discard the useful
reset. A quantitative overlap proof must compare the real ranges on both
sides of this reset with their absolute baselines. Equation (6) identifies
one concrete carry surface on which to attempt that proof.

The finite menu experiment alone does not supply overlap. Section8 obtains
an explicit finite-width overlap; unbounded interval width is still open. Future work
should combine the primitive C9 central packet with carries, or physically
change the deletion core while retaining enough endpoint coordinates to join
translated intervals. Neither more bit strings nor a finer modulus proves
that obligation.

## 7. A positive radix return map and an absolute-target decision tree

Put g=f-1, so g(0)=0 and g(2x+r)=g(x)+(1-r)(x+(x mod2)). The consumer (4)
has the form

    Phi(A,B)=c0+h(g(A)+g(B))+qAB+uA+vB+e0(A mod2)+f0(B mod2),

with c0=H_C, q=Q, u=V_C-h, v=U_C-h, e0=h-e_C and f0=h-d_C. These seven
coefficients are nonnegative, with h,q,u,v strictly positive for our banks.
Splitting A=2a+r and B=2b+s, r,s in{0,1}, returns EXACTLY the same family:

    c0'=c0+ur+vs+qrs+e0*r+f0*s,
    h'=h, q'=4q,
    u'=2u+2qs+h(1-r), v'=2v+2qr+h(1-s),
    e0'=h(1-r), f0'=h(1-s).                               (8)

This is algebraic substitution in the displayed recurrence, with no discarded
floor, parity, offset or interface. The physical interpretation is retained:
the branch bits reconstruct actual integers A,B, which reconstruct actual
transitive-deletion flanks and hence an actual tournament. The intermediate
coefficient state need not itself be a tournament and is not claimed to be.

For an absolute target N, discard a state once c0>N. If u>N-c0, its only
possible remaining A is0; likewise for v and B. A state with c0=N accepts
only remaining A=B=0. Along every branch, u and v at least double, so after
at most

    max(0,1+floor(log2((N-c0_initial)/min(u_initial,v_initial))))

splits (when N>c0_initial), both are greater than the initial budget; no
nonzero tail can survive. Thus the complete target tree always terminates.
There is no hypothesis that N is below a fixed computational bound. Running
both central initial constants H0 and H0+2 tests the complete two-point bank.
A failed search is an exact rejection for this sandwich family, not a missing
value in the full tournament spectrum.

This is an explicit finite program for every target, not a finite-state
coverage theorem: q grows as4^depth, the absolute constant varies, and the
number of surviving states need not be uniformly bounded. Quotienting these
states by residue or projective scale would lose exactly the height/anchor
needed for the decision. The remaining proof ambition is to find an invariant
collection of **accepting** target states or a compatible family change that
makes every allowed target admit a branch. Termination alone does not do so.

The implementation `radix.py` independently compares this bitwise decision
with bounded direct integer inversion for both central seeds. One-apex
recursion, symmetry supply, central transport and carry-aware inverse now
form one exact construction pipeline rather than an analogy between scalars.

## 8. A genuinely gapless overlap at every height

The radix consumer permits more than a negative diagnosis. Fix B and write
A=2a. Its even-to-odd increment is exactly

    F_C(2a+1,B)-F_C(2a,B)=QB+(V_C-e_C)-h(a+(a mod2)).       (9)

Every even nonnegative number occurs as a+(a mod2): value2u occurs at
a=2u and, for u>=1, at a=2u-1. Hence a prescribed increment delta has
arbitrarily large nonnegative solutions A even, B>=0 if and only if

    gcd(Q,2h) divides V_C-e_C-delta.                       (10)

Necessity is immediate. For sufficiency solve QB-2hu=delta-(V_C-e_C).
Every integer solution yields another by B+=2h/g and u+=Q/g, where
g=gcd(Q,2h). Both variables therefore become positive together; setting
A=4u or4u-2 gives actual binary flanks. This is a **physical** overlap
criterion, not just a gcd inference: the explicit solution is substituted
into the exact positive constructor (4).

Take C7=circulant(7,0), with steps {4,5,6}, the converse orientation of
the core used in Section2. Its outgoing-size-three central cuts are

    C641=apex(C7,{0,2,4}), C643=apex(C7,{0,1,3}),
    (h,Q,U,V,d,e)=(175,1988,952,1211,100,75).

Their matrices differ by2E14. For every integer t>=0 set

    B=7+25t,  A=170+568t       (alternatively A=172+568t).

Then (9) equals2. The four literal tournaments formed by independently
choosing core641/643 and left flank L(A)/L(A+1), with common right flank
converse(L(B)), have counts

    N, N+2, N+2, N+4.                                      (11)

Thus two physical two-point digit ranges overlap at their middle value and
form a gapless three-point odd interval at unbounded absolute counts. For
the first choice of A, its starts are

    2,568,295; 47,837,731; 149,542,367; 307,740,703; ...

These examples are not the proof: equations(9)-(10) give every t. Literal
nineteen-vertex DP independently checks all four t=0 constructions. Because
A is positive and even, A and A+1 have the same binary length and differ
only in their last bit. The corresponding left-flank graphs differ in one
root-incident arc; all other vertices and arcs can be retained. The first two such single-edge moves also have explicit common Hamiltonian
frames, reconstructed independently from the intersection digraphs in
`OVERLAP.json`. Common frames at arbitrary t are not asserted here.

Earlier sessions already constructed consecutive odd triples. This theorem
claims the explicit central-interface carry-overlap mechanism and its exact
all-height source parameters, not the first existence of high triples.
The concurrent carry lane extends the approach using a compensating diagonal
move to four-point intervals; its separate proof retains that operation.

There are two sharp renewal invoices:

* A rectangular attempt to install independent A and B increment switches has
  mixed difference exactly Q. If both lower edges have count jump2, the two
  opposite edges jump2+Q. Independent digit controls cannot be assumed here.
* After an even-to-odd A move of jump2, the next same-coordinate unit move
  jumps at least2+h*a for symmetric cores, where A=2a. Indeed it is exactly
  2-2E+h(2a+2+f(a+1)-f(a)), where E=h-e_C is between0 andh;
  f(a+1)>=1 and f(a)<=a+1 give the bound. The latter upper bound follows
  directly from(3) by induction (for odd a use f(a)=f((a-1)/2)).

So (11) is not an initialized corridor continuing forever by immediate unit
moves. It provides an unbounded family of **real overlaps**, but neither their
width nor their attachment to the small-count corridor grows. Compensating
diagonal moves, re-rooting or a wider compatible central bank remain genuine
new tasks rather than consequences of the first overlap.

## 9. Distant-source connections and provenance

Seeded archive sampling (seed276411 among tracked reflection filenames
containing carry/coin/knapsack/exchange/radix/Zeckendorf) selected
`summand-graph-fermat-zeckendorf.md` and
`lrc14-endpoint636-exchange-zero-response-and-recursive-ideals-20260831.md`.
The first's explicit correction warns that repeated summation needs actual
independent operands; it suggested retaining the complete central-root
interface before adding digits. The second keeps a whole consumer-response
class rather than one successful scalar; here the central fiber serves that
role. Neither its LRC theorem nor rejected historical arithmetic assertions
are mathematical dependencies here.

The transformation ledger is: symmetric deletion core -> fixed-cardinality
apex cuts -> full fusion-interface fiber -> translated exact count bank.
It preserves rooted fusion continuations, but loses full path-cover data and
therefore does NOT preserve arbitrary vertex-substitution contexts. The
coordinator's independent cover lane records that distinction explicitly.

## Reproduction

`python3 -B .../digits/central.py` rebuilds the complete stated cyclic menus
from endpoint subset DP. `python3 -B .../digits/replay.py` tests the new
formulas against literal path counting, the binary boundaries, and integer
inversion. `radix.py` checks the positive return map and complete target trees. `overlap.py` constructs the all-height overlap family and independently counts its first four literal nineteen-vertex tournaments using `literal_dp.cpp`. The all-order claims are proved above, not inferred from the tests.
