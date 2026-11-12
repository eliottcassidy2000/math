# A complete central fiber is an Ising chain with a numerical-semigroup edge

**Session286. PROVED relative to the inherited rooted-fusion law.** The
finite attachment table is independently reproduced from literal tournaments.
The all-length statements follow from the displayed finite automaton and
algebra, not a length census. Chan–Pak remains OPEN.

## 1. Inheritance and the precise bridge

[Session282's rational receiver carries](../../session282_chanpak_20261112/wildcard/PROOF.md)
classified the three-letter sublanguage `(A|BC)*` of a physical central fiber.
Here we use **all32** root cuts of its fixed deletion core, and classify the
entire fiber, not just that sublanguage. The additional letter opens an
unbounded negative-gap direction and turns the exact spectrum into a truncated
numerical semigroup. The fusion law and carry localization are inherited.

The distant inspirations are the exact spin representation of a *different*
tournament statistic in [THM-559, c3 line-graph energy](../../../../01-canon/theorems/THM-559-c3-is-a-line-graph-ising-energy.md),
and the [five-core binary-energy localization](../../five_core_web_20261111/root/ENERGY.md).
Neither statement is imported to prove the result here. Our map is explicit:
root-cut words in one full-central fiber -> binary spins on a chain; the
preserved observable is the **normalized difference** of Hamilton path
counts, not the full path count or all tournament structure. Required sidecars
are the baseline, module length, factor13^(k−2), and distinguished root.
Unlike the earlier energy-localization example, these spins have only two
states; unlike the older c3 energy, this one directly controls an H difference.

The bridge makes finite-state/Ising and coin-semigroup tools genuinely
applicable, while displaying exactly why their completeness does not prove
Chan–Pak.

## 2. Physical alphabet and finite certificate

Let D have outgoing bit masks `(18,8,3,5,14)` on0..4. It is strong and has
h=H(D)=13 and Q=sum_S H(D[S])H(D−S)=104. Add a root z, with cut c indicating
which vertices z beats. Let T_c denote this actual six-vertex tournament.
Use the seven-coordinate packet `(H,h,Q,U,V,d,e)` and inherited matrix

    M(T) = [h U−h d−h H−h; 0 Q 0 V−h; 0 0 0 e−h; 0 0 0 h].

Rooted fusion identifies roots and directs every left nonroot vertex to every
right nonroot vertex. Its matrix is the product. For nonempty words the root
is the same distinguished vertex throughout. Write

    A=T_9, B=T_3, E=T_18, C=T_24.

Their varying coordinates `(H,U,V,d,e)` are

    A=(37,46,71,8,5), B=(33,45,63,7,5),
    E=(37,53,64,7,6), C=(43,54,72,8,6).

The companion `attachment_table.json` gives all32 exact packets. It is a
finite proof certificate obtained both by subset dynamic programming and by
independent permutation enumeration; every attachment has at most6 vertices.
No numerical isomorphism or sample filter determines the alphabet.

For a word R of length k, compare R with A^k. S282 gives exact ordinary
receiver carries, with p=1,q=8,g=13:

    x'= (x+V(R_i)−71)/8,    y'=8y−(U(R_i)−46),
    |x|,|y|<13;             x_0=y_0=x_k=y_k=0.

The full central condition also requires first e=5 and final d=8. Applying
these explicit integer rules to the32-row table yields exactly the following
reachable states (a preceding letter is retained for the action):

|state|outgoing letter and next state|
|---|---|
|start|(3,b), (5,f), (9,a)|
|a=(0,0,9)|(3,b), (5,f), (9,a)|
|c=(0,0,24)|(3,b), (5,f), (9,a)|
|b=(−1,1,3)|(18,e), (20,g), (24,c)|
|e=(−1,1,18)|(18,e), (20,g), (24,c)|
|f=(0,−3,5)|none|
|g=(0,−3,20)|none|

The accepting states are a,c. The dead f,g do not accept. Thus trimming and
merging equivalent a,c and b,e gives the complete language

    L = (A | B E* C)*.                                      (1)

At length0 the identity/singleton is treated separately. For every k>=1,
(1) is **exactly** the words over all32 attachments whose full packet agrees
with A^k except possibly H. This is an all-length classification in a specified
physical universe, not a universal tournament grammar.

Each transition in the finite table is a finite arithmetic deduction from
literal packets; the bound |x|,|y|<13 is the all-length carry theorem that
makes this finite table exhaustive.

### Complete classification across all32 constant baselines

The same finite carry test classifies the fiber of T_c^k for EVERY cut c.
For28 of the32 choices it trims to just the repeated letter T_c: the full
central fiber is rigid at every length. The four exceptions are below; each
has precisely the grammar `(A|BE*C)*`, with the letters reinterpreted by row:

|baseline A|B|E|C|excursion action J for BE^rC|
|---|---|---|---|---|
|9|3|18|24|2−16r|
|13|7|22|28|34+16r|
|18|24|9|3|34+16r|
|22|28|13|7|2−16r|

For each baseline the verifier recomputes the entire finite reachable and
co-reachable carry graph, with its own endpoint values, rather than assuming
that symmetry transfers the first graph. All discarded states fail acceptance.
This finite graph certificate proves the statement for all word lengths.
The respective B/E/C edge weights are(−52,−16,54), (−52,16,86),
(78,16,−44), and(78,−16,−76); the three weights in a row telescope to
the displayed excursion action. Baseline/return transitions have zero action.

Thus the two other fibers have normalized energy `j=17m+8r`, or
`j=17 sum s_i−9 sum s_i s_(i+1)`. Their exact spectrum is

    {0} union {17m+8r: m>=1,r>=0,2m+r<=k}.

This is `0 union (17+<8,17>)` truncated at
`17 floor(k/2)+8(k mod2)`, by the same resource proof used below. The signed
and positive-energy sectors therefore share the same coin-semigroup boundary
but have different baselines, endpoint contributions and spin fields.
No equality of their absolute H spectra is asserted. The rest of this note
uses baseline9 unless explicitly stated otherwise.

## 3. Exact central action and binary spins

S282's local integer action J satisfies `H(R)−H(A^k)=13^(k−2) J` for k>=2.
On the trimmed automaton its weights are

    A:0, B:−52, E:−16, C:54.

Consequently one excursion BE^rC contributes J=2−16r. If a word contains
m excursions and r total E letters, set j=m−8r. Then

    H(R) = H(A^k)+2*13^(k−2) j.                         (2)

The remaining k−2m−r letters are monomers A. For k1 only A occurs and its
count is37; the negative exponent in(2) is not used.

Equivalently pin binary spins s_0=s_k=0. At step i use

    00 -> A, 01 -> B, 11 -> E, 10 -> C.

Every internal binary string `(s_1,...,s_(k−1))` corresponds to exactly one
word in the fiber. If N is its number of ones and R its number of one-runs,
then m=R, r=N−R. Hence the exact short-range lattice-gas energy is

    j=9R−8N = sum_(i=1)^(k−1) s_i
                    −9 sum_(i=1)^(k−2) s_i s_(i+1).     (3)

This is a physical **Ising-chain Hamiltonian difference** with field1 and
nearest-neighbor coupling−9 in0/1 convention. It is not a claim that the
ordinary Hamiltonian path-count spectrum is the spectrum of all Ising models.

There are exactly2^(k−1) rooted isomorphism classes in this fiber. Indeed,
deleting the root gives the ordered SCC sequence D->...->D. A rooted
isomorphism preserves each SCC position, and the four letter packets are
pairwise different. Thus it preserves the entire word. This also explains why
the spin model counts actual rooted classes, not only constructions.

The baseline is inherited from S282:

    H(A^k)=13^(k−2)[169+312k+1914(8^k−7k−1)/49+40(k−1)]  (k>=2).

The quotient by49 is integral. Keeping this baseline prevents a normalized
energy theorem from silently becoming an absolute-target theorem.

## 4. Complete spectrum, Apéry residues, and exactly56 persistent holes

For k>=2 the normalized spectrum is exactly

    J_k={0} union {m−8r: m>=1, r>=0, 2m+r<=k}.            (4)

For k>=3 its smallest and largest elements are

    L_k=17−8k,             U_k=floor(k/2).

At k2 the spectrum is{0,1}; L_2=1 is only the least nonempty-excursion
value, while the isolated monomer value0 is smaller. We retain this L_2
in the union formula below, which is valid at k2 as well.

With u=m−1 and a=k−2m−r, equation(4) becomes

    j−L_k=17u+8a,         u,a>=0, 2u+a<=k−2.             (5)

In fact the resource bound in(5) follows from the upper endpoint:

    J_k = {0} union [L_k + (<8,17> intersect [0,U_k−L_k])]. (6)

To prove the nontrivial inclusion, write t=17u+8a with a>=0; replacing
(u,a) by(u+8,a−17) gives a representative with0<=a<=16.
For k even, `t<=U_k−L_k=17(k/2−1)` yields
`2u+a=(2t+a)/17<=k−2+16/17`, hence2u+a<=k−2.
For k odd the upper bound is17((k−3)/2)+8, giving
`2u+a<=k−3+32/17<k−1`, again the required integer bound.
The converse is immediate. This is an exact finite-cost inverse, not a
heuristic Frobenius approximation.

The Apéry representatives of <8,17> modulo8 are

    0,17,34,51,68,85,102,119.

Thus its positive holes are precisely the nonnegative t whose residue i
modulo8 satisfies t<17i. There are sum_(i=0)^7 2i=56 holes, largest111.
For k>=17 the whole gap set lies below U_k−L_k and the extra j=0 already
belongs to the semigroup part (write8k−17=119+8(k−17)). Therefore:

**Theorem.** For every k>=17, J_k is the integer interval[L_k,U_k] with
exactly56 translated semigroup holes removed. In particular it contains
EVERY integer from129−8k throughfloor(k/2), and

    |J_k| = 8k+floor(k/2)−72.                            (7)

This is a full all-scale spectrum theorem in the fiber. It proves both the
persistent holes and the growing complete interior, rather than extrapolating
from a census. In ordinary H units this interval has step2*13^(k−2), not2.

## 5. Exact target inverse and multiplicities

Given k>=2 and desired ordinary count N, subtract the actual baseline and
reject unless divisible by2*13^(k−2). Put j=(N−H(A^k))/(2*13^(k−2)).
Accept j0 by A^k. Otherwise test t=j−L_k, requiring0<=t<=U_k−L_k and
`t>=17(t mod8)`. Start with u=t mod8 and a=(t−17u)/8, then repeatedly
replace(u,a) by(u+8,a−17) until0<=a<=16. Equation(5) gives a physical word,
for example `B E^r C (BC)^u A^a`, r=k−2−2u−a. This inverse is necessary
and sufficient for the ENTIRE specified fiber, and constructs the graph.

For a fixed normalized integer j, a shortest nonempty realization length is

    k_min(j)=2j                  if j>=1,
    k_min(j)=2j+17 ceil((1−j)/8)  if j<0,
    k_min(0)=1.

This follows by solving m=j+8r>=1 and minimizing2m+r=2j+17r.
This length formula is for the normalized energy across different lengths;
(2) changes its scale and baseline with length.

The bivariate program/rooted-class generating function is rational:

    F(z,t)=1/(1−z−z²t/(1−z t^(−8))).                    (8)

It counts length and normalized energy, including the empty identity.
Equivalently it is the(0,0) entry of `(I−zK(t))^(-1)` where
`K(t)=[[1,t],[1,t^(-8)]]`. Formula(8) has finite Laurent coefficients at
every z-degree. From run counting, the multiplicity of j at length k is

    sum binom(k−r−m,m) binom(r+m−1,m−1),               (9)

over m>=1,r>=0, m−8r=j,2m+r<=k, plus1 for j0. The first binomial places
m excursions among monomers; the second distributes the r E letters.

Hence exponentially many actual rooted classes share only linearly many
full packets. For k>=17, some full-packet fiber contains at least

    ceil(2^(k−1)/(8k+floor(k/2)−72))

distinct rooted isomorphism classes. This is a maximum-fiber statement,
not a lower multiplicity at each target.

Under uniform choice among these2^(k−1) classes the internal spins are
independent fair bits. For k>=3, direct covariance counting in(3) gives

    E[j]=4−7k/4,            Var(j)=(337k−832)/16.        (10)

Thus the thermodynamic perspective gives exact finite-size moments as well
as exact support. No central-limit, large-deviation or mixing hypothesis is
needed for these formulas; none is being used to assert nonempty fibers.

## 6. General carry languages have an exact integer-flow selector

There is a second precise imported-tool interface, valid beyond this example.
For any finite physical alphabet sharing h,Q, S282's carry theorem produces
a finite directed graph recognizing equal U,V and the appropriate endpoint
flags. Store the previous letter pair to make the central action an integer
edge weight. Fix a source and one acceptable terminal state (or take their
finite union). Each edge is an actual pair of physical letters.

**Finite-flow selector.** A length-k accepted walk of action J exists iff
there is an edge support F whose underlying undirected graph on incident
vertices is connected and contains source and terminal, and integers n_e>=1
for e in F, zero otherwise, satisfying

    out_n(v)−in_n(v) = [v=source]−[v=terminal],
    sum_e n_e=k,                  sum_e weight(e)n_e=J.       (11)

For source=terminal the balance is zero, with the zero-length walk treated
separately. Thus a finite union over supports and terminal states gives an
EXACT integer linear feasibility problem, not just a necessary gcd test.
If a flag implementation uses an initial edge, it consumes the corresponding
one physical letter and is included in k, as in our start state above.

Proof: a walk supplies(11) and connected support. Conversely replace each
edge by n_e parallel copies; add terminal->source if different. The resulting
balanced weakly connected directed multigraph has an Euler circuit. To see
this without an extra hypothesis, its condensation DAG could have no edge
leaving a source component: balance would force its outgoing count to equal
its zero incoming count. Weak connectedness therefore leaves one component.
The usual unused-edge splicing yields a directed Euler circuit. Remove the
added edge to get the required source-terminal trail. Reading its edge labels
constructs the physical module word/pair. ∎

Connectivity cannot be dropped: take an edge source->terminal of weight0
and an isolated loop of weight1. The balance equations with each edge used
once give length2/action1, but there is no such source-terminal walk. This is
a minimal missing-sidecar example for naïve circulation or coin-gcd transfers.

The finite support list can be enormous; no polynomial-time claim follows.
For a **fixed baseline word at each length**, first constrain the left tape
to that exact word: take the product of the pair-labelled carry graph with
the position automaton for the supplied baseline, retaining only edges whose
left letter matches the next baseline letter. Apply(11) to this product graph,
not the unrestricted pair graph. Constant-baseline graphs used above already
enforce this constraint. Only then add the supplied baseline's actual H and
factor h^(k−2). A regular left language alone does not fix the absolute H;
if it admits several left words, retain their actual counts separately or
further constrain the left tape. If h>1 and the alphabet's deletion count is
h, deleting the fused root leaves h^k Hamilton paths, and insertion into any
path gives H(full)>=h^k. Hence a prescribed count N only requires
k<=floor(log_h N), plus its length1 boundary. This yields a finite exact
all-length decision procedure in that specified baseline/fiber universe.
It neither makes the dependence on N efficient nor bounds the deletion
interfaces needed by an arbitrary tournament. General varying-left comparison
languages do not have a supplied absolute baseline and must retain it.

## 7. Boundary and next constructive question

At every k>=3 the absolute mesh is at least26. No amount of symbolic
mixing or semigroup saturation IN THIS SAME fixed-core fiber can fill all
odd H targets. The normalized semigroup is saturated away from its boundary,
but the physical normalization has unbounded arithmetic content. S282's
reroot reset demonstrates that this content can be removed outside the fixed
fusion chart. The substantive next question is whether a reset retains a
controlled, sufficiently wide part of the explicitly classified spin energy
while returning a usable receiver. A generic Ising tool does not answer that.

The supplied exact classification is useful precisely because it separates:
finite-state compatibility, full gap support, multiplicity/typical behavior,
physical scaling, and the absolute baseline. Only the first three are closed
here; Chan–Pak target anchoring remains open.
