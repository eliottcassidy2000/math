# Six-state arithmetic carriers and actual constant-target tournament returns

**Status: PROVED elementary statements below; FINITE-EXACT receipts.** This is
an audited research packet, not a canon promotion. Chan–Pak remains **OPEN**.
Date: 2026-11-12. The new scope is a complete integral monomial family, an
autonomous marked-cut controller, and its actual strong tournament lift. It is
not a classification of all tournament metagraphs or a factorization into the
directed cycle-opening primitive.

## Inheritance and concept board

[HISTORY.md](HISTORY.md) records the recovery and correction lineage. Closest
mechanisms: THM-4146's signed order-six lift, THM-4849's actual circulant/Collatz
map, THM-3166's all-order ordinal product, and the independently audited
Session36 module-context classification. The underused sidecar is THM-3023's
Newton-ratio fixed-point recurrence; its frontmatter still requires an
independent hostile audit, so the recurrence used here is proved afresh.
The canonical hostile is the sign lost under projectivization; the corrected
near miss is confusing the quadratic three-cycle with a rational scalar six.
MISTAKE-696's base-point qualification for THM-4849 is retained.

The live board has five objects: a rooted tournament and full H; a marked
ordinal cut; a positive monomial return; a signed order-six carrier; and
clocked Collatz words. The sibling compiler's six-dimensional minimal state
is a rank statement, different from any order-six statement here.

## S1. Complete rational period classification

On (Q*)^2 put L(a,b)=(b,b/a). Its coordinate sequence is

    a, b, b/a, 1/a, 1/b, a/b, a, b.

Consequently L^3(a,b)=(1/a,1/b) and L^6=id. Exactly:

* (1,1) is fixed;
* the other three pairs in {−1,1}^2 form one exact three-cycle;
* every other rational pair has exact period six.

A fixed point forces a=b=1. A point fixed by L^2 satisfies a=b/a and b=1/a,
so a^3=1 and hence a=b=1 over Q. A point fixed by L^3 has a^2=b^2=1. These
exhaust the proper divisors of six. Every nontrivial positive rational pair
therefore has exact six, and no point has exact four or five. This is a
**two-dimensional birational map**, not the one-variable polynomial x^2+c.

## S2. Complete integral classification by prime-wise hexagons

Fix positive rational C and a positive rational sequence (h[i]), and impose,
for every integer i,

    h[i−1] h[i+1] = C h[i].                           (1)

With a=h[0]/C and b=h[1]/C, the recurrence uniquely forces

    h/C = (a,b,b/a,1/a,1/b,a/b), repeated.             (2)

In particular h[i]h[i+3]=C^2. If all six entries are integers, C^2 is an
integer; rationality of C then forces C to be a positive integer.
For p^c||C let x=v_p(a), y=v_p(b). The six valuations are

    c+x, c+y, c+y−x, c−x, c−y, c+x−y.

They are nonnegative if and only if

    max(|x|,|y|,|x−y|) <= c.                          (3)

For p not dividing C this forces x=y=0. Conversely independent integer
choices in these finitely many hexagons produce every positive integral
state exactly once. Thus (3) is a necessary and sufficient classification.

The radius-c hexagon has 1+3c(c+1) points: its jth shell has 6j points. This
lattice count is inherited from the 2026-11-02 Collatz geometry packet; the
new application is the prime-wise integrality criterion. The number of
ordered states (h[0],h[1]) is exactly

    M(C) = product over p^c||C of [1+3c(c+1)].         (4)

There is one constant state, h[i]=C. Every other state has exact six by S1,
so there are (M(C)−1)/6 nonconstant unmarked cycles. Each factor in (4) is
1 modulo six. All six entries can be odd iff C is odd, and every integral
state at odd C is odd: at p=2 all six valuations can vanish only if c=x=y=0.

C=85=5*17 has49 ordered states and eight nonconstant six-cycles. One is

    425,1445,289,17,5,25,425.

For C=5^u17^v the state count is [1+3u(u+1)][1+3v(v+1)]. Every entry has
only prime factors5 and17. Because tournaments of counts5 and17 exist,
each entry has an ordinal-product tournament realization. This classifies
arithmetic states, not every tournament having those counts, and creates
no new prime factors in the individual entries.

## S3. Autonomous marked-cut return at constant full target

Let nonempty tournaments A,B have a=H(A), b=H(B), and choose m>=1. Retain
2m named copies of each, A1,...,A(2m), B1,...,B(2m). Set

    X_m = {(r,s) in Z^2: 0<=r,s<=2m, m<=r+s<=3m},
    F_m(r,s) = (2m−s,r+s−m).                          (5)

A state (r,s) means the actual ordinal sum

    A1 ... Ar B1 ... Bs | A(r+1) ... A(2m) B(s+1) ... B(2m), (6)

omitting empty lists. The bar is a marked cut. Internal labeled tournaments
are retained. The operation applies (5) and sorts each side by (6). This is
an autonomous rule; no externally supplied six-phase clock is needed.

Center the coordinates by x=r−m,y=s−m. Then the action is

    Q(x,y)=(-y,x+y),  Q^3=−I,  Q^6=I.

The domain becomes |x|,|y|,|x+y|<=m, which Q preserves. Only the center is
fixed; every other state has exact six because Q−I and Q^2−I are invertible,
and Q^3v=v forces v=0. Hence X_m has1+3m(m+1) states and m(m+1)/2
nontrivial six-cycles. After six operations the same labeled tournament and
marked cut return. Exact period six belongs to this **named-resource and
cut state**. Forgetting the cut can identify phases: at m=1, (2,1) and
(2,0) both have block order A1,A2,B1,B2 and hence the same labeled
unmarked tournament, but their cuts differ. No exact-six orbit is asserted
in the unmarked tournament quotient.

An ordinal-sum Hamiltonian path traverses the blocks in order and consists
of one internal Hamiltonian path per block. This bijection gives

    H_full=a^(2m)b^(2m),     h(r,s)=a^r b^s.           (7)

Full H is constant at every phase, although interblock arcs change. Since
Q+Q^(-1)=I, successive prefix counts obey (1) with C=(ab)^m. If a,b are not
both1, every noncentral prefix sequence has least period six. A shorter
period d=1,2,3 would make the nonzero functional (log a,log b) vanish on
the Q-orbit of (Q^d−I)(x,y); that orbit spans R^2, a contradiction.
Repeated entries within the six values are allowed.

For m=1 the boundary cycle is

    A1 A2 B1 | B2       a^2 b
    A1 B1 B2 | A2       a b^2
    B1 B2 | A1 A2       b^2
    B1 | A1 A2 B2       b
    A1 | A2 B1 B2       a
    A1 A2 | B1 B2       a^2.

The center A1 B1 | A2 B2 is fixed. Taking A=C3[1,1,2], H=5, and
B=C3[1,1,4], H=17, gives prefixes425,1445,289,17,5,25 and full target7225.
All six20-vertex full counts are independently checked by literal subset DP.

A persistent leading rooted anchor R preserves its physical root, all its
internal edges, and its root indegree. Full H and root-start count become
H(R)a^(2m)b^(2m) and H_r(R)a^(2m)b^(2m); deletions confined to R retain the
same factor. These changes do not alter anchored defect during the phase
cycle. The raw carriers are weak. Moreover a prime full target in (7), even
with an anchor, forces a=b=1 because a^2b^2 divides the target. Increasing m
changes the target unless a=b=1. Thus this finite return alone is neither a
strong ready-source reset nor the missing prime-target bridge.

## S3a. Strongification and the complete module observer

The all-order product law below is inherited from
[THM-3166, falling-factorial order-join transform](../../../../01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md).
The contextual law is inherited from
[THM-3121, substitution kernel](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md)
and [Session36, Theorem A](../../session36_information_curl_contexts_20261023/context/RESULT.md),
with its independent audit. The new application is the autonomous hexagon
return with a persistent physical root. Short proofs make the application
reviewable; the older observer classification is not claimed as new.

Let P_T(z) count ordered z-tuples of possibly empty directed paths partitioning
V(T), empty paths counting1. Splitting each path at an ordinal cut is a
bijection, so

    P_(A→B)(z)=P_A(z)P_B(z).

Thus every resource phase M has the same full polynomial

    P_M(z)=P_A(z)^(2m)P_B(z)^(2m).                    (8)

Equality P_M=P_N preserves the full Hamiltonian count in every fixed uniform
module context. It also preserves all endpoint and deletion counts for roots
outside the module. Indeed write

    P_M(z)=sum_j c_j(M) binom(z,j),

where c_j counts ordered partitions into j nonempty directed paths. For a
fixed outside arrangement a global Hamiltonian path meets the module in j
runs, and its internal choices are exactly c_j(M). Uniform outside adjacency
makes that number independent of which internal paths were chosen. Rooting
or deleting outside only restricts the outside arrangements.

The inherited observer theorem is an iff. The contexts
C_t[M]=C3[M,TT_t,TT_t], C_0[M]=M, give triangular observations in c_j,
with maximal j=t+1 and leading coefficient(t!)^2. Every gap between module
runs needs one vertex from each outside block. At the maximal value the word
is M A B M A B ... M, with t! independent assignments of A and B vertices.
The nonzero diagonal recovers all c_j successively. This is precisely the
previously audited Session36 classification. It concerns an unmarked module;
it does not promise neutrality for a root inside M or individually addressable
quotient slots.

A six-vertex hostile is M1=C3→C3 versus M2=C3[1,1,3]→TT_1. Both have H=9,
but P_M1(2)=144 and P_M2(2)=168. Their one-vertex-per-side cyclic wrappers
have H=153 and177, independently checked in the receipt. This illustrates
the same scalar-state loss as the older rooted9-to107/73/17 example, without
claiming to be that different compiler construction.

Now retain two fixed outside vertices X,Y and form W(M)=C3[1,1,M]:

    X→Y, Y→M, M→X.

W is strong even when M is weak: any two module vertices u,v are joined by
u→X→Y→v. Y has indegree exactly one from the same physical X at every phase.
A direct path decomposition gives

    H(W(M))=H(M)+P_M(2),
    start_Y(W(M))=H(W(M)−Y)=H(M).                     (9)

For the first identity, three one-M-run arrangements contribute3H(M); the
two-M-run arrangement contributes P_M(2)−2H(M). Starting at Y forces all of
M and then X. Deleting Y leaves M→X. Thus the rooted defect
H(W−Y)−start_Y(W) is zero.

The controller F_m therefore supplies an **actual strong, target-preserving,
repeatable six-step return with the same physical root and unique port**
on the marked resource-and-cut carrier. The same cut-forgetting collapse
persists after strongification; exact period six is not claimed for the
unmarked rooted graph alone.
Writing alpha=P_A(2), beta=P_B(2), the full target is

    N_m=(ab)^(2m)+(alpha beta)^(2m).                  (10)

It is a module permutation, not a proved sequence of directed cycle-opening
moves. Its complete contextual neutrality also explains the obstruction:
external rooted readouts and defect are unchanged, so it cannot repair a
missing cover coordinate. A defect-changing repair must break the uniform
module interface, root inside it, or introduce another nonneutral operation.

For the5 and17 cores, direct two-bin enumeration gives

    (a,alpha)=(5,32),  (b,beta)=(17,224).

At m=1 this22-vertex strong carrier has

    N_1=85^2+7168^2=51,387,449=29*1,771,981,
    start_Y=H(W−Y)=7225.

Its marked internal prefixes remain425,1445,289,17,5,25. The factor29 is a
hostile to interpreting the sum-of-squares form as a primality certificate.
The22-vertex number uses proved identity(9), with independent local P2
counts; it is not reported as full22-vertex subset DP. Independent full
subset DP checks all six16-vertex carriers with A of count5 and B=C3 of
count3: each full count is147681 and both root-start and deletion counts
are225, with the same indegree-one root and port.

There is an exact near-square/Fermat link, with a sharp limit. a,b are odd;
alpha,beta are even because swapping two bins is a fixed-point-free
involution. N_m is therefore an odd square plus an even square. If prime,
it is1 modulo4. A common divisor of ab and alpha beta forces compositeness.
If m has an odd divisor d>1, the expression factors as U^d+V^d into two
nontrivial positive integer factors. Thus a prime in this family requires m
to be a power of two and the two bases to be coprime; neither is sufficient.

To obtain the structural unit-square branch1+n^2 by setting the first
displayed square equal to1 requires ab=1. Other accidental representations
of N_m as1+n^2 are not excluded by this argument. Then A,B are
transitive by the H=1 criterion in the current correction to
[THM-209, independence polynomial](../../../../01-canon/theorems/THM-209-independence-polynomial.md),
so P_A(2)=2^|A| and P_B(2)=2^|B|. The branch reduces to

    N_m=1+2^[2m(|A|+|B|)].                            (11)

For singleton A,B and m=1,2,4,8 it gives17,257,65537,4294967297. Primality
requires the displayed exponent to be a power of two by the same odd-factor
argument. This is an actual strong rooted realization of Fermat-type
near-squares, not general n^2+1, new Fermat primes, or Landau's conjecture.
The structural restriction is H=1 forcing transitivity and hence a
power-of-two two-bin count.

## S4. Positive storage for a signed order-six lift

Incoming-surface credit: commit c022b73c3c, read via `git show origin/main`,
contains [Session293, sixth-mode carrier](../../session293_cycle_returns_20261113/six_history/PROOF.md).
Its conserved-positive-mass theorem gives a stronger six-hidden-state lower
bound without assuming a globally finite-order integer update, and it credits
[THM-2619, projective kernel and retained-target obstruction](../../../../01-canon/theorems/THM-2619-psl2f13-seven-edge-norm-minimal-projective-kernel-and-retained-target-obstruction.md)
for the older nonnegative-inverse mechanism. The elementary lemma
below is a narrower recovered control, not a novelty claim. Session293 also
constructs a different rational Möbius map of projective exact order six;
that source should be used for the one-dimensional rational-map variant.

Let a rational2x2 matrix A satisfy A^2−A+I=0. Then A^3=−I and A^6=I.
A pointed cone K cannot contain a whole nonzero A-orbit, because it contains
both v and−v. Nor can a linear readout have the same positive value on that
orbit: its values at v and−v are opposite. Hence a positive full H cannot
be read directly as a constant from the unshifted signed two-state carrier.

A sharp narrower lower bound also holds. If N is a nonnegative integral
dxd matrix, N^6=I, and a rank-two rational P satisfies PN=AP, then d>=6.
Indeed N^{-1}=N^5 is nonnegative. A matrix and its inverse both nonnegative
must be monomial: the off-diagonal zero entries of their product force
disjoint positive supports, and invertibility forces one per row and column.
Integral inverses force these entries to be1. Hence N is a permutation
matrix with cycle lengths dividing6. Its quotient A has primitive sixth
roots as eigenvalues, so N must contain a six-cycle. Equality d=6 is attained
by the six-cycle permutation with quotient columns v,Av,...,A^5v.

This does not bound nonlinear positive encodings, arbitrary affine storage,
or the minimal six-state Hamiltonian compiler. The latter is an independent
rank theorem. S3 translates signed coordinates into positive resources and
retains their complement; that reservoir preserves H, not the number six
alone.

## S5. Three Collatz clocks and the actual tournament map

On positive integers distinguish raw C_-(n)=n/2 if even,3n−1 if odd; Terras
T_-(n)=n/2 if even,(3n−1)/2 if odd; and odd-only U_-(n)=oddpart(3n−1).
The three familiar cycles have these lengths:

| Minimum | Raw | Terras | Odd-only |
|---:|---:|---:|---:|
|1|2|1|1|
|5|5|3|2|
|17|18|11|7|

The third odd cycle is17,25,37,55,41,61,91 with valuations
1,1,1,2,1,1,4, whose sum is11. Seven odd raw steps yield total length18.
All full words are in the receipt. Exhaustiveness of these familiar positive
cycles is not proved here.

An equivariant one-step map from a directed m-cycle onto an n-cycle exists
iff n divides m. Applying the mth iterate proves necessity; reduction of
phase modulo n proves sufficiency. After an initial phase choice it is
unique, with fibres of size m/n. Thus the raw18-cycle can map six-to-one
to the quadratic three-cycle, while raw2 and5 cannot. It preserves neither
numerical value nor parity: the18-cycle has seven odd positions, which
cannot be a union of six-element fibres. Under the Terras clock only the
cycle at minimum5 has suitable length3; under the odd-only clock none has
length divisible by3. Three cycle components themselves are fixed by
iteration, not a three-cycle of points.

THM-4849 provides a genuine different map. A circulant tournament on odd
q=2k+1>=3 has connection set S with S disjoint from−S and
S union−S=(Z/qZ)\{0}. The indicator word has q letters and k ones. Its
Terras3x+1 rational periodic point has positive numerator and denominator

    D=2^(2k+1)−3^k>0.

For prime q the word is primitive and has exact period q. Reflection x→−x
conjugates T_+ to T_−, so the corresponding3x−1 rational cycle is negative.
More generally every3x−1 parity word of this same shape has fixed point
−c_w/D<0. Hence the familiar positive3x−1 cycles cannot be these half-density
circulant cycles. The dictionary retains the marked parity word and its
base-point condition; it does not identify a cycle point with H(T).
Independent enumeration at q=3,5,7 gives H={3},{15},{175,189}, respectively.
The q7 values share(q,k,D)=(7,3,101), proving that the arithmetic shape and
unreduced denominator do not determine the full Hamiltonian count.

## S6. Recovered quadratic orbit and literature boundary

The recovered orbit is

    f(x)=x^2−29/16,
    −7/4 → 5/4 → −1/4 → −7/4.

The repeated triple(−7/4,−1/4,−7/4) is not an exact three-cycle and omits5/4
from this recovered orbit. At this specific parameter, the full affine
rational preperiodic set is {±1,±3,±5,±7}/4. An odd-prime denominator forces
negative valuations to double. At2, preperiodicity requires v2(x)=−2:
smaller valuations double and larger ones jump to−4. The real escape bound
(2+sqrt33)/4 then forces |4x|<=7. All eight candidates are checked directly.
Thus this parameter has no rational exact six-cycle.

THM-4146's interpolation in y=4x is M(y)=(y−13)/(y+3). The normalized
matrix B=(1/4)[[1,−13],[1,3]] satisfies B^3=−I,B^6=I, but M has projective
order3. Its central six is not an extra scalar orbit of f; S4 adds the
positive-readout obstruction.

For every rational c, rational periods4 and5 of x^2+c are excluded by Morton
and Flynn–Poonen–Schaefer. Stoll's exclusion of period6 assumes analytic
continuation and functional equation of the relevant L-series and BSD rank
equality for a specific genus-four Jacobian. These are quadratic-polynomial
statements, not assertions about every rational dynamical system. Primary
sources checked2026-11-12:

* [Morton, Arithmetic properties of periodic points of quadratic maps II](https://matwbn.icm.edu.pl/ksiazki/aa/aa87/aa8721.pdf).
* [Flynn–Poonen–Schaefer, Cycles of quadratic polynomials and rational points on a genus2 curve](https://people.maths.ox.ac.uk/flynn/arts/art11.pdf).
* [Stoll, Rational6-cycles under iteration of quadratic polynomials](https://www.mathe2.uni-bayreuth.de/stoll/papers/Xdyn06-revised.pdf).

S1 gives an explicit different rational dynamical system with six but no
four or five, without contradicting any of these polynomial results.

## S7. Five-period exchange, golden ratio, and the odd-count obstruction

For positive rationals set K(a,b)=(b,(1+b)/a). Direct substitution gives

    a,b,(1+b)/a,(1+a+b)/(ab),(1+a)/b,a,b.

Thus K^5=id. A fixed state satisfies a=b and a^2=a+1, whose positive solution
is the irrational golden ratio. Since5 is prime, every positive rational
state has exact period5. This is the elementary five-period exchange
recurrence; no cluster-algebra theorem is needed for the identity. Deleting
the additive1 changes K into L, with positive rational period6. This is a
change of operation on the same space, not a conjugacy.

If a,b are odd integer full tournament counts, (1+b)/a has even numerator
and odd denominator in lowest terms. If integral, it is even. Every tournament
has odd full H, so this exchange cannot directly act on successive full
counts. The simplest example is1,1,2,3,2. A proposed five-core/exchange
construction needs another readout or an auxiliary coordinate. The golden
ratio likewise cannot itself be an integer full count. The monomial map
avoids this parity obstruction at odd C, but its prime-target and module
neutrality boundaries remain.

## Verification and stopping boundary

Run from the repository root (all checks remain active under optimization):

    python -B -O 05-knowledge/results/cycle_opening_six_return_20261112/dynamics/verify.py
    python -B -O 05-knowledge/results/cycle_opening_six_return_20261112/dynamics/strong_wrapper_control.py
    python -B -O 05-knowledge/results/cycle_opening_six_return_20261112/dynamics/exchange_control.py
    python -B -O 05-knowledge/results/cycle_opening_six_return_20261112/dynamics/audit_compiler.py

The main receipt covers7396 rational pairs, all reservoirs1..200 by independent
divisor search, resource radii1..12, all circulant words at q3/5/7, all three
specified Collatz clocks, and six full20-vertex ordinal carriers. Strong
controls count all six16-vertex wrappers, full and rooted, plus the equal-H9
hostile. Exchange controls cover1849 positive rational pairs and256 odd-input
pairs. The sibling compiler audit uses independent rational elimination and
permutation/weak-cut enumeration rather than importing its counting code.

The infinite results are proved above; no finite bound is extrapolated.
The sibling [module return](../compiler/MODULE_RETURN.md) now proves one
uniform successor replacement that changes the rooted state at the required
full target. What remains missing is an indefinitely iterable construction
with the coverage needed for Chan–Pak. The contextually neutral controller
here cannot supply that state change by itself. No claim closes
Chan–Pak, Landau, or the positive3x−1 cycle classification.
