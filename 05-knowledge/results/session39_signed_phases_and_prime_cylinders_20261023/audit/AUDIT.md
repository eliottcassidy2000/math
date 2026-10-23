# Session 39 independent audit: first checkpoint

Status: PROVED algebraic checks; FINITE-EXACT independent enumeration.
No Git mutation was performed by the auditor.

## 1. Point sections are not synchronized frames

Let connected graph edges carry bijections g_ij of a common fiber F, with
g_ji=g_ij^{-1}. Choose a root and spanning-tree transport T_i from the root
to i. Every section has x_i=T_i(x_root). A non-tree edge imposes
x_root=H_ij(x_root), where H_ij=T_j^{-1}g_ij T_i. Thus sections are in
bijection with the intersection of the based cycle fixed sets. A synchronized
full frame requires all H_ij to be identity, a strictly stronger condition.
Here full frames permit the entire affine/permutation transport group. If frames must preserve each cell's fixed cyclic orientation, every reflection edge is already an obstruction, even on a tree. The rotation-only gauge classification instead retains these edge signs.
For F_3 an affine reflection has exactly one fixed point; a nonzero
translation has none. Two individually soluble reflection constraints may
still conflict if their fixed points differ.

Independent enumeration used all 6^5=7,776 affine assignments on edges
01,02,03,12,23, comparing all 3^4 point assignments with tree-chart fixed
sets. Counts: no section 5,616; one section 1,944; three sections 216.
This was executed independently in JavaScript, without importing worker code.

## 2. The girth-five phase-blindness theorem survives audit

For named cyclic triple cells with matching/complement blocks on a marked
support graph of girth at least five, varying affine phases does not change
the exact incidence row multiset E in either orientation, while retaining
fixed signs, polarities, and uniform background blocks.

Proof mechanism: c_v is constant within each cell and phase-independent,
since triangle counts involve at most three cells, whose support is a forest.
For e_uv, if the two cells are adjacent, no third cell is marked-adjacent to
both, so the only variable is their direct relative coordinate. If the two
cells are nonadjacent, at most one third cell is a marked common neighbor,
so there is at most one length-two transported relative coordinate. Other
third-cell contributions are constant. In the adjacent case the outgoing
restriction uses the same direct coordinate; in the nonadjacent case it is
uniform. Summing the multiset over the target cell therefore eliminates phase.
The same proof works after reversing all arcs.

Independent JavaScript reconstruction of the five-cell cycle witness gave
vertexwise E equality for phases 0,1,2, Hamilton-path counts 45,132,777;
45,150,825;45,102,765, cell scores (9,10,7,4,5), and cyclic-triangle vertex
incidences (21,17,25,17,21). The distinct H values certify nonisomorphism
without needing cell markings. Every motif using fewer than the support
girth distinct cells is phase-blind by forest gauge; fewer than girth
vertices is a sufficient, weaker condition.

## 3. Signed valuation words and prime cylinders

For a word a_1,...,a_k>=1, let A_j=sum_{i<=j}a_i and
B_k=sum_{i=0}^{k-1}3^{k-1-i}2^{A_i}. The final oddness congruence
3^k n+epsilon B_k=2^{A_k} mod 2^{A_k+1} specifies exactly one odd
starting residue. To prove sufficiency, decompose
B_k=3^{k-j}B_j+2^{A_j}C_j. The endpoint congruence forces every prefix
numerator divisible by 2^{A_j}. All intermediate values are therefore
integral. Backward recurrence from the odd endpoint forces all preceding
values odd, hence all prescribed valuations exact.

For a closed word D=2^A-3^k is odd and nonzero. Divisibility D|B and
epsilon D>0 give the positive integral fixed point n=epsilon B/D, which
is odd; the same argument proves exact admissibility. These statements
do not imply global convergence or classify all cycles.

Independent exhaustive residue checks covered 240 signed words of lengths
one through four with letters one through three, including 12 integral
positive closed-word controls. Every case passed. Prime infinitude and
the asymptotic proportion 2^{-A} use Dirichlet and PNT in arithmetic
progressions as cited dependencies, not as experimental conclusions.
A primality algorithm using only finite valuation-stream observations
cannot terminate correctly: each finite prefix is shared by primes and
composites. This excludes access to the original integer or other side data.

## 4. Critical linking and characteristic refinements

For nonsingular symmetric integral A, D=coker A has
b([x],[y])=x^T A^{-1}y mod 1. If c_i=A_ii mod 2, then
q_c([x])=(x^T A^{-1}x-x^T A^{-1}c)/2 mod 1 is well-defined:
the shift x->x+Az changes it by x^Tz+(z^TAz-z^Tc)/2, an integer.
Its polar form is b. This affine quadratic refinement is not the naive
homogeneous bosonic twist for an odd lattice.

The map ker(A mod2)->D[2], u->[Au/2], is an isomorphism without assuming
the 2-primary group has exponent two. Its pairing is u^TAv/4=q1(u,v)/2.
The exponent-two hypothesis is needed only to replace D[2] by the whole
2-primary group. Thus equal q1 cannot be separated by 2-linking in an
exponent-two carrier, though a characteristic refinement can add information.

For K_2=(Z/2)^3358+(Z/4)^4+Z/256, let g have order256 and odd self-linking
unit u. Every order256 element is ag+h with a odd and ord(h)<=4. Multiplying
self-linking by256 gives a^2u plus terms divisible by8. Since a^2=1 mod8,
u mod8 is an unmarked invariant. The supplied energy arithmetic checks:
14,960,072,960/256=58,437,785=153 mod256=1 mod8.
Scaling A by3 preserves every dyadic kernel and mod-two form but multiplies
this unit by3, separating the 2-primary pairings. The full group can acquire
odd torsion; this is a scaled-graph control, not a same-capacity counterfeit.
The underlying large graph energy was independently recomputed from all 460,699 frozen edge-ledger rows on 3,528 vertices. The integer current matched the supplied witness entrywise, and the edge and witness SHA-256 hashes are recorded in terminal_receipt.json.

## 5. Flux algebra and scope

K=[[3,4],[4,3]]=[[0,1],[1,0]]+3tt^T=-I+4tt^T with t=(1,1),
det K=-7 and t^TK^{-1}t=2/7. Three-flux attachment changes fermion
exchange parity, but the bosonic integer-Hall precursor is an additional
physical assumption. It is not a derivation of an actual phase or a claim
that every cited source uses this convention.

For even-order tournament skew S, S mod2=I+J is invertible and det S is
an odd Pfaffian square. The even symmetric off-diagonal BF matrix
[[0,S],[S^T,0]] has signature zero, discriminant group
coker S plus coker S^T, and odd fourth-power order. It cannot encode a
nontrivial dyadic critical-group sector. This is an exact algebraic mismatch,
not a physical no-go theorem. The opposite homogeneous forms
q_+(x)=x^2/512 and q_-(x)=-x^2/512 on Z/256 are not isometric, because
that would require an odd square equal to -1 modulo512, already impossible
modulo8.

## Limits

The independent checks establish the listed finite universes only. They do
not prove LRC, Collatz convergence, a universal complete graph invariant,
physical realization of K, or completeness of any dyadic sidecar. External
physics source attribution remains the flux worker's separately documented
primary-source task.

## Reproduction and additional structural audit

Run from the workspace root:

- python3 session39-audit/audit_exact.py
- python3 -O session39-audit/audit_exact.py
- python3 session39-audit/audit_graph.py
- python3 session39-audit/audit_terminal.py

Matching receipt.json, receipt.optimized.json, graph_receipt.json, and
terminal_receipt.json are saved. Normal and optimized exact-check receipts
are identical; the graph DP runs the backwards endpoint recurrence,
independently of the worker's forward recurrence.

Rotation-gauge enumeration independently gave 16 triangle classes
(stabilizer one:4; stabilizer three:12) and 144 diamond classes
(stabilizer one:72; stabilizer three:72). A symbolic count also reproduces
these: flat sign systems are 2^(v-1) choices; flat phase orbits have
root-translation stabilizer3, whereas nonflat systems have stabilizer1.

In the cycle family, the scores are s_0=3g-6, s_i=3g-2-3i for
1<=i<=g-2, and s_last=5. They are distinct, so cells are intrinsic.
Any automorphism rotates each cell and connected matching support forces
the rotations equal. Thus the automorphism group is exactly C3, and
cycle translation distinguishes all phases even without Hamilton-path counts.
With intrinsic cell colors supplied, individualizing one vertex resolves
its cyclic triple in one directed-WL round; singleton colors then resolve
each adjacent matching-linked cell in one further round. Therefore
1+eccentricity(root cell) rounds suffice. Without preloaded cell colors,
include degree-refinement preprocessing. Strong connectivity follows from
strong cells and bidirectional support links.

For the doubled-edge switch ab,cd -> ac,bd, diagonal entries stay fixed and
Delta(u^TAu)=4(u_a-u_d)(u_b-u_c). Set alpha=u_a+u_d and beta=u_b+u_c
over F2. Then Delta q_c=alpha beta/2, whose polarization is
alpha tensor beta + beta tensor alpha. If alpha=beta nonzero, the bilinear
form is unchanged but the refinement changes by the nonzero linear
character alpha/2. The supplied quadratic-value count vectors give
Gauss sums -2-2i and 2-2i, ruling out isometry of those refinements.

The flux note was independently inspected after its final edit: it
explicitly distinguishes the cited four-vortex convention from the
alternative conditional three-flux EFT and distinguishes q_t(e1)=-2/7
from the representative exchange phase exp(-3 pi i/7).

## Native partition and characteristic addenda

The signed triple-cell carrier is forced, not merely an ansatz, for any
stable equitable tournament partition whose cells have sizes one and three.
Each triple has constant internal outdegree and three arcs, hence is a
cyclic triple. A triple-triple block has constant row and column sums,
equal by double counting; degrees zero/three are uniform, and one/two are
permutation/complement blocks. Triple-singleton blocks are uniform. An
exact encoding also retains the singleton tournament and uniform block
orientations. Named cells admit rotation gauge; unnamed cell permutations
would require a further quotient.

For the native merged metagraph, the inherited odd fibre size h yields
full degree h on self-converse nodes and2h otherwise. Removing loops
changes full degree by an even number, so the grounded diagonal modulo2
is exactly the retained self-converse indicator c. Therefore
z^TAz=c^Tz mod2 and c supplies a characteristic refinement. Under the
explicit formal choice K=A, the local lattice excitation Az has exchange
sign (-1)^(sum_SC z_i). This says nothing about microscopic electrons,
charges, or realization of a physical phase. The flux addendum retains
these qualifications.

The concurrent incoming packet bc3c88ec033 already contains the native
cyclic256 squareclass and characteristic convention. Its generator has
105/256 rather than153/256:9^2*153=105 mod256, while27^2*105=1 mod256.
This is an explicit abstract cyclic-form isometry, not verification that
saved potentials in different graph labelings are scalar multiples.
The linking worker correctly credits the earlier mechanism and identifies
the four newly detected n5 quadratic switches as its finite increment.

## First-checkpoint publication boundary

The marked-cell producer is forthcoming and is not among the approved producer pins in this first checkpoint. The mathematical statements above and independent audit_graph.py controls are included as independent audit results, without requiring any unpublished producer artifact. This checkpoint approves only the hash-pinned flux, characteristic, signed-cylinder, and linking files in PRODUCER_PINS.json. Later closed-target and marked-cell producer additions require separate audit pins.

For a relocated packet, run audit_terminal.py --edges PATH_TO_SESSION38_WIRING/edges.tsv.gz --witness PATH_TO_LINKING/n8_terminal_lift.json. Other audit scripts are standalone.
