# Session50: A000568 from fibre arithmetic and symmetry-corrected growth

24 October 2026. **Exact aggregate reconstruction, proved marked-parent transfer,
and finite-exact corrected predictions at orders eight and nine.** This is not
a closed formula in the order alone, a recurrence on the naked blue/black line
graph, or a claim to outperform known Burnside formulas.

## 1. What information the prediction actually uses

A000568 counts ordinary tournament isomorphism classes. The converse-merged
metagraph counts a different object. The session separates four data levels:

1. Word-cube size, reflection counts, edge totals and other aggregate constraints.
2. Authenticated ordinary fibre masses, obtainable from the correct merged word
   masses and their parity. A loop contributes two literal port incidences.
3. Actual parent representatives, attachment subsets and parent automorphism
   actions, retained as native marked incidence data.
4. Owner-permutation labels on reobservation arrows and their closed walks.

The aggregate decoder uses level2. The corrected growth predictors use level3.
A native parent representative can be read from its fixed-path word, but simply
forgetting the word/owner data and keeping the blue/black graph does not supply
these inputs. No proof below says that arbitrary quotient summaries do.

Known degree-local Burnside engines remain exact arithmetic baselines. The new
question is which native metagraph information can replace child classification,
and which symmetry data cannot be discarded. [A000568](https://oeis.org/A000568)
is already a computable classical sequence; it is the proposed compressed
metagraph explanation, not existence of an exact counting formula, that remains
open here.

## 2. Odd fibres give unusually strong modular count information

The [moment theorem](mass-predictions/THEOREM.md), with
[independent audit](audit-moments/AUDIT.md), inherits the fibre identity

    h(T)=H(T)/|Aut(T)|,

where h is a positive odd integer. For n>=2, fixed-path words number
N=2^m, m=(n-1)(n-2)/2. Define the collision moments

    K_r = sum_words h(word)^r = sum_ordinary_classes h^(r+1).

On a merged vertex of word mass mu, h=mu for odd mu and h=mu/2 otherwise.
This uses the authenticated fibre rule, not an arbitrary graph degree or a
loop-deleting Laplacian diagonal.

For k>=1 let Q_k(x)=product_(j=0)^(k-1)(x-(2j+1)). At every odd h,
Q_k(h) is divisible by2^k k!, while Q_k(0) is odd. Consequently the first k
moments K_0,...,K_(k-1) determine the unknown class count modulo

    2^(k+v2(k!)) = 2^(2k-popcount(k)).

Newton interpolation proves this exponent is sharp within the stated family
of degree-at-most-k integer polynomials with odd constant coefficient. This is
not an optimality theorem for every conceivable reconstruction method.

If the exponent is at least m, the residue and1<=A<=2^m already determine A.
Thus a polynomial-length vector of modular moments is a universal exact count
certificate. Generating that vector without the underlying census remains a
separate problem.

### A short order-eight certificate with no count input

The isolated decoder accepts only K_0,...,K_7, eight integers. K_0=2^21 is known
from n, so seven positive-degree moments are additional data. They give

    A000568(8) = 6880 modulo32768.

No finite mass range, atlas maximum, support list or target count is needed.
For the short explicit upper certificate put

    q(x)=(x^2-1250x+465000)/50000000,
    P(x)=2q(x)-x q(x)^2+(1-x q(x))^2.

For all real x>=1,

    P(x)-1/x = (x-1)(1-x q(x))^2/x >=0.

The aggregate sum of this degree-six polynomial uses only K_0,...,K_6 and is
exactly

    496641366973194539/15625000000000 <32768.

Since A=sum_words1/h is positive, the residue and this bound force6880.
The seventh-degree word moment enters only the modular residue. An independent
decoder, complete small-order word censuses, and a two-file-only run verify the
boundary: the moments were extracted from an existing atlas, but neither that
atlas nor its class count is an input to the decoding calculation. This is
held-out aggregate recovery, not prediction from n alone.

## 3. Marked-parent growth has an exact unpointing projector

The [growth theorem](growth-counts/THEOREM.md), with
[independent audit](audit-growth/AUDIT.md), starts from every ordinary parent T
of order n and its full automorphism action on attachment subsets S. A state
p=(T,[S]) represents one vertex-pointed child of order N=n+1. Its stabilizer is
h_p=|Aut(T)_S|=|Aut(U)_marked-root|.

Re-rooting the constructed child U at each of its vertices requires only
order-n deletion identifications. Let D[p,q] count the resulting roots of type q.
Then

    D^2=N D,                 A000568(N)=trace(D)/N=rank(D),
    |Aut(U)|=h_p D[p,p].

The proof uses child classes to explain the matrix blocks, not as algorithmic
input. Each block repeats the vector of actual vertex-orbit sizes and has trace
N. Generation never canonicalizes an order-N unpointed child.

Converse gives a commuting permutation C, with self-converse count trace(CD)/N.
The weighted even/odd projectors recover merged and converse-pair counts. These
are spectra of a marked-incidence operator, not of the original complement-line
metagraph. Endpoint-only growth has a different state set and needs the proved
anti-endpoint kernel; ordinary converse does not preserve end-root eligibility.

The full-D verification recovers at order eight:54,256 vertex-pointed states,
trace55,040,6880 ordinary classes,176 self-converse classes and3528 merged classes.
Marked totals must not be divided by N without the actual orbit correction.
Diagonal-only computation and singleton pruning are proved opportunities; the
full-D engine itself was not advertised as a benchmarked faster implementation.

## 4. A different correction can eliminate most overlap tests

For a tournament U and vertex v, let b(U,v) be the size of its stable directed
refinement cell and r(U,v) its actual automorphism-orbit size. Define

    e(U)=sum_vertex_orbits O (b(U,O)-|O|),
    E_N=sum_ordinary_U e(U),
    U_N=sum_pointed_states p b(U_p,root).

Refinement cells contain entire automorphism orbits, so e>=0, and exactly

    A000568(N)=(U_N-E_N)/N.

This E_N is not the earlier root-moving symmetry defect
Delta_N=sum_p(D[p,p]-1). The three sums satisfy

    U_N=|P_N|+Delta_N+E_N.

For example order eight has |P_8|=54,256, Delta_8=784, E8=28, and U8=55,068.
The next sections determine E structurally rather than fitting it to a known A.

### The inherited cell machinery and its new counting extraction

The [cell-error theorem](cell-holonomy/THEOREM.md), with
[independent audit](audit-cell-holonomy/AUDIT.md), explicitly credits Session39's
affine connections and Session40's odd-cell, coprime-cell and bounded-core results.
The signed-kernel automorphism formula and216-to16 matching-class calculation
were already proved there; they are not discoveries of this session.

Stable cells induce regular tournaments and therefore have odd size. Coprime
cell sizes force a uniform crossing. In the singleton/triangle carrier, signs
of matching permutations control the rotational automorphism kernel. The new
error extraction is e(T)=6u, where u is the number of triangle cells in unbalanced
signed components. Translation phases can distinguish classes without changing
this kernel. Thus the data needed to count the orbit correction can be strictly
smaller than the data needed to reconstruct complete class wiring.

Strong components have disjoint total-outdegree intervals. Refinement and vertex
orbits restrict to them, making e additive over their unique ordered sequence.
For ordinary class series A(z) and strong-component error series E_s(z),

    E(z)=A(z)^2 E_s(z).

This is an exact weighted SCC identity, not a new scalar recurrence for A.

## 5. Order eight without child canonicalization or reroot tests

The odd stable-cell partitions of a strong eight-vertex tournament permit at
most two interacting triangle cells, or a regular five-cell with uniform
coprime attachments. Their cells are actual automorphism orbits. A seven-cell
plus singleton would have uniform dominance and cannot be strong. Hence

    E_s,8=0.

The independently complete regular-seven census has three classes, with
Aut orders21,7,3. Only the last is not vertex-transitive and contributes
E_s,7=14; lower strong errors vanish. The SCC identity gives E8=28.

The [isolated refinement predictor](refinement-predictor/THEOREM.md) now needs
only the complete parent atlas and automorphism actions through seven. It
computes U8=55,068 from attachment orbits and stable root-cell sizes, computes
the parent error14, and returns

    (55,068-28)/8=6880.

It uses no child-order canonicalization, child automorphism search, deletion or
reroot-isomorphism test. Its independent audit reconstructs the parent atlases
from a different exhaustive smaller-order construction, verifies every supplied
full automorphism action, and runs the blind predictor with only its declared
input file. The target order-nine request is explicitly rejected by this
order-eight implementation rather than extrapolating the vanishing theorem.

## 6. Order nine: a complete structural error correction

The [order-nine theorem and predictor](order-nine/THEOREM.md), with
[independent audit](audit-order-nine/AUDIT.md), take the next genuine growth step.
Their input is the complete order-eight parent atlas, recovered from authenticated
native metagraph representatives. Parent distinctness, full automorphism actions
and labelled orbit mass2^28 are verified independently.

The complete list of odd stable-cell partitions of nine yields the following
error contributions. These are not residuals obtained by subtracting a supplied
A000568(9).

1. **Regular nine:630.** Each almost-regular eight-parent (four scores3, four
   scores4) has a unique regular-nine extension. The85 such parents are exactly
   its vertex-pointed states. Only765 order-eight deletion identifications give
   regular reroot trace135, hence15 regular-nine classes and
   E_regular9=9(85-15)=630. No order-nine canonicalization is used.
2. **A seven-cell inside a strong cyclic three-block quotient:14.** The bad
   regular-seven cell with two singleton cells retains the earlier error14.
3. **Three nonregular triangle cells:72.** Unbalanced matching cycles give the
   four inherited rigid types, each with error18. A cyclic majority-orientation pattern
   would make the whole tournament regular and was already counted in item1;
   the nonregular pattern is transitive and has no such overlap.
4. **Nonstrong padding:42.** The bad seven-component with two singleton strong
   components has three positions, each contributing14.

Every other stable-cell case has zero error. Therefore E9=758.

The regular-sector reduction itself works at every odd order N=2m+1. An
almost-regular2m-parent has exactly m scores m-1 and m scores m; its unique
regular extension directs the new root toward the score-m vertices. Thus
P_regular,N is the number of these ordinary parent classes,
A_regular,N=trace(D_regular)/N, and
E_regular,N=N(P_regular,N-A_regular,N). The same root-transfer proof applies.
This is a conditional exceptional-sector algorithm requiring the complete
parent data, not a closed scalar recurrence for all regular tournaments.

The bulk calculation then uses refinement only on1,716,608 pointed attachment
states. Its root-cell histogram is

    size1:1,712,964; size3:3,491; size5:48; size7:20; size9:85.

It gives U9=1,724,582, and hence

    A000568(9)=(1,724,582-758)/9=191,536.

The only reroot-identification subproblem is the85-state regular part above.
The bulk scan does not classify or canonicalize order-nine tournaments. A
separately written implementation begins from the original parent source,
reconstructs automorphisms by a different search and computes refinement by
incoming-neighbor splitting; it obtains the same complete result. Python and
compiled controls agree. The known target value is used only for validation.

This proves the count from full order-eight parent data without first building
the unpointed order-nine atlas or metagraph. It does not produce complete
order-nine wiring, and it is not a count obtained from order eight's node total
alone. The recorded compiled scan time excludes input preparation and other
work and is not a universal speed comparison with known counting formulas.

## 7. Why the remaining information barriers are real

The [relaxed realization theorem](identifiability/THEOREM.md) and
[information-boundary audit](audit-boundaries/AUDIT.md) show that cube symmetry,
local signature support and several aggregate constraints are insufficient.
At n=6, explicit C,R partitions satisfy native V=34,N=1024,B=64,M2=48724,
allowed local signatures, loop totals and stated connectivity requirements,
but the exact relaxed ambiguity is A in {48,56}. They contain explicit score
witnesses proving that they are not tournament-isomorphism partitions. They
are not competing actual tournament counts and do not refute the native56.

The [owner-labelled reobservation theorem](inertia-counts/THEOREM.md), also
independently audited, retains permutation labels on internal observer changes.
For a component with b normalized words and closed-walk group G,

    b|G|<=H.

Equality certifies the whole class fibre and all its automorphisms. Partial
components give sum_B b_B|G_B|/H_B <= A <= number of components. Finite unrestricted-
span relocation saturates through order seven without canonicalization in
generation. No all-size connectivity theorem is inferred from those checks.

Fixed-width contiguous-window reobservation provably misses information at
arbitrarily large orders. Cyclic three-block inflations have locally unique
observer windows while retaining a genuine order-three global symmetry.
Bounded span is not the same restriction as moving only one vertex: an end-to-end
single-vertex relocation can span the whole observer.

## 8. What this session establishes, and the next open step

The two counting routes have an exact interface, at fixed N>=1 and integer s>=0.
For a pointed state p, let
L_p count Hamiltonian paths ending at its marked vertex, h_p be the parent
attachment stabilizer, r_p=D[p,p], and H_p the child's total path count. Then

    K_s = sum_p (L_p/h_p) (H_p/(h_p r_p))^s.

To verify it, fix a child U. Its rooted vertex orbits have sizes r_p and
h_p=|Aut(U)|/r_p, while sum_p r_p L_p=H(U). Its total contribution is therefore
(H(U)/|Aut(U)|)^(s+1), as required. This identity connects the proved endpoint
fibres and isotropy extraction; it is not a new implemented fast moment engine.
In particular a global correction E_N does not supply every individual r_p.

The advances have three different strengths:

- A small exact arithmetic certificate reconstructs a count from authenticated
  mass aggregates.
- Full marked parent incidence gives an all-order unpointing identity.
- Structural refinement-error classification turns that identity into concrete
  order-eight and order-nine predictions without general child classification.

No claim is made that the new counts themselves were previously unknown, that
classical Burnside, factorial divisibility or affine gauge theory is new, or that
this is the fastest numerical evaluator of A000568. Inherited and newly applied
mechanisms are separated in each packet.

The remaining research target is a compressed native source of moments or of
marked attachment/refinement data, together with a structural strong-error law
beyond the proved orders. In particular, further interacting five- and larger
cells and new regular cores cannot be assigned the triangle rule by analogy.
Recovering the order-nine count does not automatically provide the order-nine
parent atlas needed for an order-ten run.
Neither a naked margin table nor a guessed continuation of zero strong error
supplies that law. The complete graph, observer and port-realizability questions
remain distinct from the successful count recovery.
