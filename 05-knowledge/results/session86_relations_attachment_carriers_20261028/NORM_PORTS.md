# A closed norm-one port carrier, and exact partner-query costs

2026-10-28. **PROVED normalization and separation; FINITE-EXACT query tests.**
This is a scoped positive complement to `FIELD_RESULTS.md`: arbitrary port
alignment is not recovered, but an explicitly normalized port model does
have a canonical joint response.

## 1. Marked K-lines have canonical norm-one representatives

Keep E=F_(q^2), K=F_q and the actual field structure, q=2^m. For x!=0 set

    N(x)=x^(q+1),    nu(x)=x/sqrt(N(x)),
    U={u in E^*:N(u)=1}.

Squaring is a bijection of K, so the square root is unique. If k in K^*,
then N(kx)=k^2N(x), hence nu(kx)=nu(x). Also N(nu(x))=1, and the only
norm-one point on xK is nu(x): any other candidate differs by k with
k^2=1, which in characteristic two forces k=1. Therefore nu induces a
canonical bijection

    E^*/K^*  -> U.

Since nu(xy)=nu(x)nu(y), this is a group isomorphism. The carrier U is
cyclic of order q+1, closed under multiplication, inversion and Frobenius.
It is **not** closed under general vector addition. If an attachment
language allows arbitrary addition of ports, its closure is all E because
U spans E (proved below): 64 or 256 states, not 9 or 17. It has 9 ports in F64
and 17 in F256; these counts are not the numbers of nonzero field vectors.

The normalization is extra *declared semantics*, using an already retained
norm map. It deliberately removes the K^*-scale. It does not discover an
unknown physical or arithmetic scalar on an arbitrary port. If a marked
line is already supplied, its normalized port is uniquely computable and
requires no partner queries at all.

## 2. A genuinely determined joint table

For u,v in U, v^q=v^(-1), so the reconstructed trace form gives

    B_1(u,v)=Tr_(E/F2)(u/v).                                (1)

Thus line labels, the field/norm structure, and the declared normalization
suffice to calculate this joint table; independently estimated marginal
histograms are not being promoted to joint data.

The complete rows of (1) distinguish all normalized ports. To prove it,
U spans E over F2. Indeed if z generates U, its degree over F2 is the
multiplicative order of 2 modulo q+1, namely 2m. The order divides 2m. No
positive d<m can give 2^d=1 modulo 2^m+1 by size; for m<=d<2m it would
force 2^m+1 to divide 2^(d-m)+1, again impossible by size. Hence the first
2m powers of z form a binary basis and lie in U. For u!=v, nondegeneracy
of B_1 then gives some w in U with B_1(u,w)!=B_1(v,w).

Consequently no quotient identifying two members of U preserves **all**
partner responses. In this precise sense U is a minimal faithful carrier
for normalized-port identity under all partner queries. This is not a
universal minimal attachment carrier, nor a claim about arbitrary external
port states. Its group multiplication is not an encoding of every
possible composition rule from another application.

## 3. How many allowed partners identify an unknown normalized port?

Now change the information model: u in U is unknown, and the allowed
binary observations are B_1(u,w), where each fixed partner w must itself
lie in U. We ask for a fixed, nonadaptive set of distinct partners whose
response vector identifies every u in U. This is different from computing
nu(x) when the original line xK is already known.

| Field | Ports | Information lower bound | Fixed U-partner minimum | Adaptive U-partner worst case | Number of minimum fixed sets | One minimum fixed set |
|---|---:|---:|---:|---:|---:|---|
| F64 | 9 | 4 | **6** | **5** | 27 | `{1,6,11,20,26,58}` |
| F256 | 17 | 5 | **5** | **5** | 714 | `{1,15,26,36,85}` |

The counts are among subsets of the explicitly marked U, not up to group
equivalence. `norm-ports.json` records all minimal subsets, complete binary
pairing rows, example separating response vectors, and the exhaustive
subset-search counts at every smaller cardinality. The search tests all
subsets of size 0 through the attained minimum; no heuristic optimizer is
used. Normal and `python -O` outputs are identical.

### Why F64 needs six instead of the naive four

Here U is cyclic of order 9 and splits into three cosets of its subgroup
mu_3. In the chosen polynomial marking the cosets are

    {1,58,59}, {6,26,28}, {11,20,31}.

Each is the three nonzero vectors of a binary plane: if omega has order 3,
then 1+omega+omega^2=0. Distinct members of the same coset have ratio of
order 3, whose absolute trace in F64 is 1. Ratios from different cosets
have order 9, hence degree 6 and minimal polynomial

    Phi_9(t)=t^6+t^3+1,

whose t^5 coefficient is zero. Their trace is therefore zero. Equal
arguments also give zero. Thus (1) is the disjoint union of three 3x3
blocks having zero diagonal and ones off diagonal. The planes are mutually
orthogonal and together span E.

A partner in one block gives zero on both other blocks. Fewer than two
partners from any one block cannot distinguish its three ports: with one
partner, two ports both respond 1; with zero, all three respond 0. At least
two partners from each block are therefore necessary, for a total of six.
Taking exactly two from each is sufficient: all nine responses are
distinct, with nonzero support in their own probe block. There are
binom(3,2)^3=27 such choices. This proves the exact F64 minimum without
relying only on enumeration.

The six-versus-four gap belongs to this restricted partner-observation
model. It does not say that nine labels contain six unconstrained bits,
nor that arbitrary binary encodings require six bits. The six-query
minimum is explicitly **nonadaptive**; adaptivity changes it as follows.

### Adaptivity saves one U-partner query in F64

The exact adaptive worst-case cost is **five**. For an upper bound, query
one partner from each of the three blocks. If a 1 appears, its block and
two remaining candidates are identified; one further partner query
separates them, using at most four queries in total. If all three answers
are 0, exactly three candidates remain, namely the queried partner from
each block. At most two further queries distinguish those three, giving
five in the worst case.

Four cannot suffice. The 0 branch of any first query leaves seven ports:
one from its block and all three from each other block. Every possible
second query has a 0 branch containing at least five of those ports.
Two remaining binary queries can distinguish at most four, so every
four-level decision tree fails on that branch.

The checker independently solves the exact minimax recurrence on all 512
F64 candidate subsets S:

    d(S)=0 for |S|<=1,
    d(S)=1+min_w max(d(S intersect O_w),d(S minus O_w)),
    O_w={u:B_1(u,w)=1},

where queries having an empty branch are omitted. It records every subset
cost and an explicit optimal decision tree, and verifies its answer on
all nine ports. For F16 the same computation uses all 32 subsets and gives
four; its K5 interpretation below is repeated equality testing. For F256
the counting lower bound five and the fixed five-probe witness already
prove adaptive optimality five, without a larger DP.

Thus the F64 cost sequence is six for fixed U-partners, five for adaptive
U-partners, and four when arbitrary linear partners are allowed. The last
model reaches the counting bound even nonadaptively; adaptivity cannot
reduce it further.

### F256 reaches its counting bound

There are 17 ports, so four binary probes can give at most 16 response
vectors and cannot suffice. The displayed five-partner set has 17 distinct
responses, verified in the receipt, establishing the sharp minimum five.
Exhausting all binom(17,5)=6188 five-subsets gives exactly 714 successful
sets. This is a finite marked-field result, not an all-m query formula.

### The F256 joint table is a Paley graph, not a tournament

Choose zeta=36 in the specified F256 marking; it generates U of order 17.
Index the ports by i in Z/17 via i->zeta^i. Direct exact evaluation gives

    B_1(zeta^i,zeta^j)=1 iff i-j is a nonzero square modulo 17.

The residue set is `{1,2,4,8,9,13,15,16}`. This is precisely the Paley
**undirected graph** on F17, with its exponent marking explicitly supplied.
Because -1 is a square modulo 17, the adjacency is symmetric, as an
alternating binary pairing must be. It is not the directed Paley-tournament
construction at 7. The receipt checks every one of the 17^2 entries.
This is a structural interpretation of the computed table, not an
identification of F256 addition with F17 addition: the latter labels the
multiplicative cyclic port group by exponents.

## 4. Enlarging the allowed probe objects saves two queries in F64

If partners may be arbitrary w in E, rather than only normalized w in U,
the exact minima are **4 for F64 and 5 for F256**. F256 already reaches
its counting bound with U-partners. For F64 the four probes

    {3,12,17,33}

have nine distinct response vectors on U; their common B_1-annihilator
is the plane `{0,2,40,42}`. Four is minimal because three binary answers
can encode only eight outcomes.

Here is an exhaustive kernel certificate and its structural explanation.
Let D={u+v:u,v in U,u!=v}. A linear response map is injective on U exactly
when its kernel avoids D. In F64, |D|=36. The complete list of 651 binary
two-dimensional subspaces contains exactly 36 avoiding kernels, all saved
in the receipt. With the orthogonal three-plane decomposition above, D
consists exactly of vectors supported on one or two of those planes.
Therefore every nonzero vector of an avoiding kernel has all three
components nonzero. Each projection from the kernel onto a component plane
is an isomorphism. Using the first projection to identify the kernel with
that plane, the other two projections are independently in GL(2,2), giving
6^2=36 avoiding kernels. Their B_1-annihilators have dimension four and
supply four-probe separating maps.

Thus the six-to-four improvement is an exact change of allowed query
objects, not the creation of information. Closing the port vocabulary
under addition makes those new probe vectors available, while increasing
the state carrier from U to all E. The small multiplicative carrier and
the additive carrier serve different attachment languages. Neither is
silently substituted for the other.

## 5. The F16 control connects the five-port model without conflating fifteen

Using alpha^4+alpha+1, the norm-one set is
`{1,8,10,12,15}`. It is a canonical transversal of the five F4-lines in
the fifteen nonzero binary vectors, not the full fifteen-vector space.
Its pairing matrix has ones exactly off the diagonal: the complete graph
K5. Therefore four fixed U-partners are necessary and sufficient to distinguish
its five ports (five minimum sets). Adaptivity still needs four in the
worst case: each response 1 merely excludes the queried port, while a 0
identifies it; an all-1 transcript must exclude four of five candidates. Two unqueried ports would otherwise
have identical all-one responses. Allowing arbitrary field partners drops
the exact minimum to three, achieved by `{1,2,4}` and forced by the counting
bound. The common annihilator is `{0,1}`. All 15 one-dimensional kernels
were tested; exactly five avoid the ten nonzero port differences.

The control field model is certified independently in `norm_ports.py`:
powers of alpha run through all fifteen nonzero residues and return to 1.
The same certification is applied to F64/F256 in that script. The F16
control is recorded separately under `controls`, leaving the requested
F64/F256 cases under `cases`.

## 6. A residual phase can be invisible to the chosen normalized observable

The two short-orbit F64 lines in `FIELD_RESULTS.md` normalize to 58 and 59.
These elements lie in F4 and are fixed by F^2. Hence the residual order-three
Frobenius stabilizer of either exceptional line does not move its canonical
norm-one port. It can still act on other unnormalized points of that line:
F^2(3)=17 was the earlier one-port phase hostile. Thus form/frame ambiguity
does not destroy every possible joint observable; the answer depends on
which port state is retained and which normalization is authorized.

## Verification

Run:

    python norm_ports.py > norm-ports.json
    python -O norm_ports.py > norm-ports-optimized.json
    cmp norm-ports.json norm-ports-optimized.json

The script imports only the elementary field arithmetic and rank/exception
helpers from the independently audited `check_field_forms.py`. It verifies normalization
for all nonzero field elements, scale invariance, multiplicativity for all
ordered nonzero pairs, closure of U, the ratio formula, full row separation,
every candidate normalized-partner subset up to the optimum, all 651
possible two-plane kernels for the unrestricted F64 probe problem, and the
F16 control. The F256 Paley identification is checked entry by entry. The adaptive
F64/F16 claims use complete candidate-subset DP plus explicit decision-tree
verification. Its assertions are
explicit exceptions and survive optimized Python.
