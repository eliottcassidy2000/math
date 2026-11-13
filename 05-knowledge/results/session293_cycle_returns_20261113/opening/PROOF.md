# Neutral cycle-opening is renewable, but its sleeping seam has a price

**PROVED IN NOTE; independent review requested.** This is a classification of a
restricted but actual physical return, not a proof of Chan–Pak. The inherited
mechanisms are THM-1900 (neutral insertion iff an SCC cut), THM-002 (odd-cycle
packing formula), THM-4876 (cycle-weight deletion increment), and
`triples_ports_sieve_20261112/squares/PROOF.md` §9 (directed cycle opening).
The one-arc join and near-join laws in THM-4865 / HYP-9457 are close antecedents
of §3 below; the positive separator derivation here is a two-sided application,
not a claim to discover near-joins.

## 1. Complete neutral-opening classification

Let T be a tournament, a,b,v,t four distinct old vertices, and a→v an old
arc. Adjoin z with a→z→b, retaining all old arcs. Let O be the number of
Hamilton paths of the extension starting v, ending t, and using a→z→b.
Suppose the full count is unchanged: H(T+z)=H(T). Then O>0 **if and only if**

* T=A→B, where A is its first SCC and is nontrivial;
* z is placed exactly at the cut A→z→B;
* a,v belong to A, b,t to B;
* A has a directed Hamilton cycle using a→v, and B has a Hamilton path b→t.

Moreover, with c_av(A) denoting the cycles using that oriented edge,

    O = c_av(A) h_bt(B),              H(T+z)=H(A)H(B).

Proof. THM-1900 places z at an SCC cut P→z→Q. The inherited inverse of
cycle opening partitions T into a directed cycle through a→v and a path
b→t. Since b lies in Q, that directed path cannot visit P. Since a lies in
P, its directed cycle cannot visit Q. Consequently the cycle spans P and
the path spans Q. A tournament admitting a spanning directed cycle is
strong, so P must be exactly the first SCC A. The converse follows by
opening each such cycle and joining it through z to each b→t path. The
partition and inverse are unique, giving the factorization. This proof does
not assume every edge in a strong tournament belongs to a Hamilton cycle.

## 2. A genuine indefinitely repeatable marked return

Fix such A,B,a,v,b,t. For exact conditioned renewal assume additionally
that every Hamilton path of B ending t starts b, that is, e_t(B)=h_bt(B).
A sufficient physical condition is that b is a source of B. Put

    T_m = A → z_m → z_(m-1) → ... → z_1 → B,  m≥0.

At m=0 use b as the source-path start; at m≥1 use b_m=z_m. Each birth
of z_(m+1) is a directed cycle opening using the SAME a→v cycle edge,
fixed final endpoint t, and the current path-start port b_m. Its conditioned
count is always c_av(A)h_bt(B); its full count is always H(A)H(B).
For m≥1 the conditioned suffix path is z_m,...,z_1 followed by a path
of B ending t; the additional hypothesis makes its start b automatic.
Without that hypothesis, the first opening count is c_av(A)h_bt(B),
whereas every later opening count is c_av(A)e_t(B). Thus the endpoint-column
support condition is also NECESSARY for equality of these nonzero sector
counts throughout this conveyor. It cannot be silently discarded when
handing the port to a newborn.

All old-old Hamilton endpoint counts are unchanged: every Hamilton path
traverses A, then the unique seam chain, then B. All newly inserted seam
vertices have zero global start and end counts. Thus this return preserves
more than the scalar H. It renews its operational b-port, but it does NOT
retain that port as the same physical vertex: it deliberately hands the
port to the newest seam vertex. This distinction is essential.

**Fixed-name obstruction.** If all four old names a,v,b,t must remain the
same, a second positive neutral opening birth is impossible after the first,
provided all earlier vertices and arcs are retained. The first newborn z_1
is a singleton SCC after A and before b. A later positive neutral opening
would again have to occur immediately after the first strong SCC A, but
its source path starting at the old b cannot visit z_1, which precedes b.
The source cycle cannot visit z_1 either. This contradicts spanning.
Thus moving the operational port is not cosmetic: it is exactly what makes
the all-length conveyor possible.

Conversely every positive full-H-neutral opening birth is of this form at
its current cut. Therefore a sequence consisting solely of induced neutral
births cannot change its strong core through cycle opening. This conclusion
is stronger than merely saying that one proposed gadget fails, but much
narrower than ruling out deletion/replacement or arc-changing returns.

## 3. The seam is invisible now, exponentially visible after feedback

For arbitrary nonempty tournaments A,B, choose a∈A,b∈B and m≥1. Begin with
A→TT_m→B, and reverse only the cross arc a→b to b→a. Call the result J_m.
Define the nonnegative integer responses

    L_a(A) = Σ_{S⊆A,a∈S} s_a(A[S]) H(A\S),
    R_b(B) = Σ_{U⊆B,b∈U} e_b(B[U]) H(B\U),

where s_a counts paths starting a, e_b paths ending b, and H(empty)=1.
Then the exact full count is

    H(J_m)=H(A)H(B)+2^m L_a(A)R_b(B).                 (1)

Proof using THM-002. Old odd cycles and their packings all lie inside A or
B and contribute H(A)H(B). A new cycle must use b→a, and therefore consists
of a path starting a on S⊆A, a possibly empty increasing subsequence of
TT_m, and a path ending b on U⊆B. At most one such new cycle is present in
a disjoint packing. For each S,U, exactly 2^(m−1) subsets of the middle chain
make its total order odd. Its OCF weight is 2; the complement's packing
weight sums to H(A\S)H(B\U), since unused middle vertices are singleton SCCs.
Summing proves (1). The assumption m≥1 is essential: there is no equal
split of the empty chain's subsets by parity.

There is also a direct Hamilton-path proof, independent of the OCF identity.
Paths avoiding b→a traverse A, the whole middle chain, then B. For m≥1 none
uses the removed a→b arc, so these give H(A)H(B). A path using b→a has the
unique form

    A_left, M_left, B_left ending b, a starting A_right, M_right, B_right.

The two A-pieces partition A; the two B-pieces partition B; and each middle
vertex independently chooses left or right, in its inherited increasing
order. All empty pieces are allowed except A_right and B_left. Summing the
outer path choices gives L_a R_b; assigning the middle vertices gives 2^m.
This is an actual bijection, not just a coefficient identity.

At m=0 the same backedge sector remains L_aR_b, but the avoiding sector
loses e_a(A)s_b(B), the old paths using a→b. Therefore the exact boundary is

    H(J_0)=H(A)H(B)+L_a(A)R_b(B)−e_a(A)s_b(B),

recovering the inherited two-base one-arc-join law in its orientation gauge.
The jump from m=0 to m=1 includes this restored endpoint product; replacing
the empty seam by a formal nonempty identity would miss it.

Both responses are positive (use the singleton subsets), so every added
neutral seam vertex doubles the feedback excess above the baseline. The
old-old endpoint matrix alone therefore does not certify neutrality for
future arbitrary attachments. The newborn vertices are visible to future
cycles even though no current global Hamilton path starts or ends there.

If the outer product h=H(A)H(B) is fixed while cores/ports may be replaced,
a fixed activated target c requires

    L_a(A)R_b(B) = (c−h)/2^m.

Consequently m≤v_2(c−h), and advancing m by one requires exact halving of
the positive response product. This is a sharp integral budget for this
specific reactivation language, not an unrestricted count-coverage bound.
It identifies what a successful growing-defect return must change: the
baseline, the activated interface, or the closure operation, rather than
merely insert more H-neutral length.

## 4. Smallest transparent example and scope

Take A=C3 and B=TT2, a→v a cycle edge, b the first and t the last vertex of
B. Then H=3 and the conditioned opening count is 1 at every seam height.
Here L_a(A)=3 and R_t(B)=2, so feeding the last suffix vertex t back to a
produces H=3+6·2^m for m≥1 (15,27,51,...). The target-neutral chain therefore
exists and renews, but closing its return is not neutral.

This separates three notions often conflated: a support-level cycle is
opened; the current full H and endpoint matrix return; and a future feedback
operator also has the same output. Only the first two are automatic in this
classified family. Directed cycle-opening remains useful, but it must be
combined with exact response compensation to make new Chan–Pak targets.
