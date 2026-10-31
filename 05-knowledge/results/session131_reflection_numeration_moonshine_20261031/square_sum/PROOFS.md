# Square-sum reflection words and sharp residuewise seam bounds

**Session131. PROVED elementary theorems; FINITE-EXACT witness verification.**
No new global Hamiltonian existence theorem is asserted. Gerbicz's inherited
all-order existence for N>=32 remains credited to its original source.
The new main result is the sharp *residuewise uniform lower bound* below,
not attainment at every order. It is attained at the eight orders80..87.

## 0. Inheritance and the typed connection

Write G_N for the simple graph on1..N with edge xy iff x!=y and x+y is a
square. Session125's `SUCCESSOR_PHASE_AND_CARRY.md` section2 already identifies
the four odd-square channels modulo8 and uses them to prove a six-defect
obstruction for a permanently position-colored family. We retain those
channels, but now count the endpoint incidences of **all** odd-square path
components of a single Hamiltonian cycle. Vertex-label parity, position
parity, edge-root parity and digit-position parity are different types.

Session125's rooted-forest bridge is a bijection from fixed-length CG words
to rooted forests of an ordered path, preserving nonzero digit count. Its
forests are not Hamiltonian boundary forests in G_N. The CG/Zeckendorf
numeration proofs are likewise inherited; section5 below supplies an explicit
odd/even-digit interleaving map, not an inference from equal counts.

The overnight packet l correction overrides a graph-complete reading of
packet j: cutting prescribed ports of selected local Hamiltonian cycles
need not enumerate every global cycle on a cyclic assembly. Our cut below
starts from the **actual global cycle** and retains every resulting path,
including isolated vertices with two endpoint incidences. It does not assume
that each local boundary state recloses into a single source cycle.

## 1. Main theorem: even-square seams have a sharp eight-residue floor

For a Hamiltonian cycle C of G_N, let t(C) be the number of edges xy whose
positive square root sqrt(x+y) is even. If N>=8, then

| N modulo8 |0|1|2|3|4|5|6|7|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Lower bound on t(C) |4|7|8|7|8|7|8|7|

All eight bounds are sharp as uniform bounds depending only on N modulo8:
`witnesses.json` supplies actual Hamiltonian cycles for N=80,81,...,87
attaining the corresponding entries. This does **not** say that every N
in a residue class attains its floor. For example the inherited unique
G32 cycle has sixteen even-root edges, not four.

### Proof by endpoint flows

Odd squares are1 modulo8. Removing every even-root edge from C therefore
leaves a spanning linear forest, all of whose paths lie in one of

    A={0,1}, B={2,7}, C={3,6}, D={4,5} modulo8.

There is at least one removed edge: otherwise a single cycle could not
visit all four nonempty blocks. Removing t>0 distinct edges from a cycle
produces exactly t paths. An isolated vertex is a path with two distinct
*incidences* at the same physical vertex. Let k_A,...,k_D count these paths;
then each k_i>=1 and t=k_A+k_B+k_C+k_D.

Let n_r count the labels congruent to r modulo8, and put

    δ_A=n_1−n_0, δ_B=n_7−n_2, δ_C=n_3−n_6, δ_D=n_5−n_4.

Each odd-root edge in a block uses one vertex of each of its residues.
If a block has m such edges, its endpoint-incidence count in residue r is
2n_r−m, while its number of paths is k=n_r+n_q−m. Consequently

    E_1=k_A+δ_A, E_0=k_A−δ_A;
    E_7=k_B+δ_B, E_2=k_B−δ_B;
    E_3=k_C+δ_C, E_6=k_C−δ_C;
    E_5=k_D+δ_D, E_4=k_D−δ_D.             (1)

The removed even-root edges pair these incidences. Their square sums are
0 or4 modulo8. Among odd residues they pair {1,5} with {3,7}. Among even
residues, {0,4} pair only internally, and {2,6} pair only internally. Set

    X=k_A+k_D, Y=k_B+k_C,
    Δ_L=δ_A+δ_D, Δ_R=δ_B+δ_C.

It follows that

    X−Y=Δ_R−Δ_L,
    X≡Δ_L (mod2), Y≡Δ_R (mod2), X,Y>=2.   (2)

The imbalance pairs for N modulo8 are

    (Δ_L,Δ_R) = (0,0),(1,0),(1,−1),(1,0),
                 (0,0),(1,0),(1,−1),(1,0).

For residue0 or4, (2) gives X=Y even and t>=4. For odd N, it gives
X odd, Y even, X=Y−1, hence X>=3,Y>=4 and t>=7. For residue2 or6,
X,Y are odd and X=Y−2, hence X>=3,Y>=5 and t>=8.

For residue4, t is a multiple of4. If t=4, then all k_i=1 and
δ=(1,−1,1,−1). Block B's single path has both endpoints in residue2.
An even-root seam from residue2 can go only to residue2 or6. No other
block has residue2 endpoints, and block C's two endpoints are residue3,
not6. Thus B's endpoints would have to be sewn to each other, producing
an isolated cycle, contrary to connectedness of C. Therefore t!=4 and
its next possible value is8. This completes the universal proof. QED.

### What the theorem preserves and what it loses

The flow vector records all four path counts and all eight endpoint
residue counts. It gives necessary arithmetic and connectivity obstructions.
It does not tell which actual labels are endpoints, which endpoints belong
to the same path, whether the proposed seam sums are squares, or whether
all paths become one cycle. The complete sidecar is the named endpoint
matching plus actual labels and the union-connectedness test. This is
exactly where the packet l correction matters.

## 2. Marked path endpoints: a reusable boundary version

Let P be a Hamiltonian **path** in G_N, N>=8, and let b_r be the number
of its two original endpoints in residue r. Cut its t even-root edges.
The resulting forest has t+1 paths and equation(1) still holds. The two
original endpoint incidences are left unpaired, so equation(2) becomes

    X−Y=Δ_R−Δ_L+b_1+b_5−b_3−b_7,
    X≡Δ_L+b_0+b_4 (mod2),
    Y≡Δ_R+b_2+b_6 (mod2),
    t+1=X+Y, all k_i>=1.                 (3)

If both original endpoints lie in one residue block, that block has at
least two cut-forest paths: otherwise its one path already joins the two
original endpoints and cannot visit the other blocks. This is an additional
necessary connectedness constraint, not permission to forget label pairings.

Equation(3) is proved by subtracting b_r from each endpoint pool before
pairing. For original nice-pair ports1 and8, the closing edge is odd-root;
for ports1 and3 it is even-root. Thus the endpoint-aware version aligns
with, but does not replace, the inherited position-color certificate.

## 3. The odd-square graph has exactly four components from order44

Let H_N keep only odd-root edges of G_N. For N>=44, its components are
exactly A,B,C,D intersected with[1,N].

An explicit base proof uses H40. For each a=1,2,3,4 its two components are
these paths (their consecutive sums are9,25 or49):

    [24+a,25−a,a,9−a,16+a,33−a],
    [32+a,17−a,8+a,41−a].

These eight paths partition1..40 and exhaust the odd-square edges there.
Arrival40+a joins9−a via49 and41−a via81, merging exactly the two paths
of its residue block. Hence H44 has exactly the stated four components.

Every arrival n>=13 has an earlier odd-square neighbor. Indeed for odd
s>=5 the integer intervals

    [(s²+1)/2, s²−1]

cover every integer from13 onward: the next left endpoint is at most
s² because s²−4s−5>=0. An n in this interval has neighbor s²−n in[1,n−1].
Every later vertex therefore joins an already present component in its
own residue block, and different blocks can never merge. Induction proves
the claim. This component theorem is separate from Hamiltonicity of H_N's
components; connectedness alone is not a Hamiltonian-path certificate.

## 4. Reflection-root words: exact return and collision conditions

For a positive integer root s, define the partial reflection

    R_s(x)=s²−x,
    D_(s,N)={x∈[1,N]:s²−x∈[1,N], 2x!=s²}.

On D_(s,N), R_s is a fixed-point-free involution and its two-element
orbits are exactly the edges of root s. Odd s flips label parity; even s
preserves it. It is not generally an automorphism of the entire square-sum
graph: R_5 sends the edge1--3 to24--22, whose sum46 is nonsquare.

For a root word s_1,...,s_m, define c_0=0 and c_i=s_i²−c_(i−1).
Starting at x, the successive labels are

    x_i=(-1)^i x+c_i.                     (4)

Thus closure at m has two different types:

* m even: c_m=0, with x initially free;
* m odd: 2x=c_m, with at most one possible initial integer.

The labels for i=0,...,m−1 lie in[1,N] exactly when the integer x lies
in the intersection of

    [1−c_i,N−c_i] for even i,
    [c_i−N,c_i−1] for odd i.              (5)

Injectivity is equally explicit. If two equal-parity indices have equal
c_i, every x is invalid. Otherwise equal-parity indices never collide.
For every even i and odd j, exclude x=(c_j−c_i)/2 when it is integral.
Equations(4)--(5), closure, and these exclusions are necessary and
sufficient for a rooted simple cycle of length m>=3. For m=N they are
necessary and sufficient for a Hamiltonian cycle of G_N. Rotation and
reversal still identify multiple rooted descriptions of the same cycle.

This is a complete word validator, not a method proving that a word exists
for every N. Closure alone is too weak: the word(3,3,3,3) has c_4=0 but
repeats its start after two steps, so it never gives a simple four-cycle.

Two valid realizations of the same even root word differ by an alternating
shift x_i -> x_i+(-1)^i h. An odd word is rigid. More generally square
scaling plus alternating shift gives x_i -> r²x_i+(-1)^i h, preserving
all consecutive square sums; an odd closing edge forces h=0. Bounds and
injectivity must still be rechecked. This is the phase-aware mechanism
behind the inherited Gerbicz block formula, not a uniform label translation.

For comparison, an increasing integer-affine map T(x)=ax+b preserving
every square-sum edge on the positive integers must be T(x)=r²x. Apply it
to1 and s²−1 for all s>=2: a s²+2b=t_s². Then the integer differences
t_(s+1)−t_s converge to sqrt(a), forcing sqrt(a)=r integer. Now the integer
t_s−rs=2b/(t_s+rs) tends to0, so is eventually0 and b=0. The converse
is immediate. Square scaling conjugates reflections:

    R_(rs)(r²x)=r²R_s(x).

It embeds a chosen carrier, but its cycle is not Hamiltonian in the larger
full initial-interval graph. CG growth is not such an embedding: the
canonical shift sends1->3 and3->8, destroying the edge1--3.

## 5. Explicit nonconsecutive odd/even Fibonacci rows

Use F_0=0,F_1=1 and Zeckendorf coefficients z_k for k>=2, so z_k∈{0,1}
and z_k z_(k+1)=0. Split into interwoven rows

    e_j=z_(2j+2), o_j=z_(2j+3), j>=0.

The exact legality conditions are

    e_j o_j=0 and o_j e_(j+1)=0 for every j. (6)

Conversely any finitely supported pair of binary rows satisfying(6)
interleaves into a Zeckendorf word. There is no additional within-row
nonconsecutivity rule: neighboring entries of one row are separated by
one Fibonacci position. The omitted sidecar in a row-count comparison is
precisely the alignment and the cross-row rung/diagonal exclusions(6).

With A_j=F_(2j+2), its value is exactly

    V=Σ_j e_j A_j + Σ_j o_j(A_(j+1)−A_j),

so its signed even-weight coefficients are

    c_j=e_j+o_(j−1)−o_j, o_(-1)=0.         (7)

This is the explicitly indexed version of Session125's signed conversion.
The ordinary integer parity is

    V mod2 = Σ_(k not divisible by3) z_k mod2,

not parity of j, digit count, or just a terminal digit. On canonical
numerals the actual square reflection is the conjugate

    z -> Zeckendorf(s²−V(z))

on the canonical encodings of D_(s,N). It is an involution because decoding and
canonical re-encoding are inverse. It preserves nonconsecutivity by
normalization, while its integer parity follows root parity as above.
This is not digit-string reversal or raw bit complement: for s=3,
1 (digit1) reflects to8 (digit10000).

The source/target map preserves actual integer values and square-sum
adjacency only when the square-root label s and the exact value equation
are retained. Constraint(6) alone does not give a square-sum edge or a
Hamiltonian continuation. The CG split into *two even-index-only binary
summands* is another inherited construction and is not these two interwoven
index-parity rows.

## 6. Verification and scope

`check.py` uses only Python's standard library, explicit exceptions (not
assert), and data relative to its own location. Normal and optimized runs
must agree. It independently validates every cycle and seam, checks the
path-forest endpoint counts and all incidence pairing data, and proves
finite optimality at80..87 by comparison with the elementary universal
bound. It tests the reflection compiler against direct walking over a
bounded complete word universe, the H44 base and all arrivals through
2000, and the two-row bijection over all binary words of lengths0..18.
It also checks conjugated digit reflections and hostile cases.

The numerical optimizer was only a witness finder. Neither its floating
arithmetic nor its optimality message is used in the proof or portable
checker. Source paths, Git blobs and SHA256 pins are in `source_pins.json`.
No all-order attainment, prescribed-port completeness, induced sparse
carrier equivalence, or new proof of Gerbicz's existence theorem is claimed.
