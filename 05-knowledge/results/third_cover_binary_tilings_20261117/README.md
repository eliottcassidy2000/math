# A third cover match and the two-tiling Hadamard chart

**Proved.** There are literal tournaments B and B' with
P_B'(q)=7P_B(q) at q=1,2,3 but not 4. This gives a ×7 switch
that survives two arbitrary common cyclic wraps. This witness
is not prime-useful: its lower H and P(2) share 3. Separately,
the two-tiling idea gives an exact chart for framed sign matrices.
Neither the literal metagraph complement pairing nor its
literal mirror pairing gives a Hadamard matrix at order ≥4
in this chart. Ordinary Hadamards still occur in the chart
for every power of two, via Gray-coded Sylvester matrices.

## 1. From a reflexive binary matrix to two framed tournaments

Take the entrywise map M=2B−J from 0/1 to ±1. It is a bijection
of entries, not linear over the usual rings. Orthogonality in
sign coordinates says that every two distinct binary rows have
Hamming distance n/2. Now require B_ii=1 and, in a fixed order,
B_i,i+1=1, B_i+1,i=0. The upper triangle gives a tournament T,
and the *complement* of the transposed lower triangle gives
a second tournament U. Both are decoded along 1→2→…→n.
There are two independent copies of the ordinary marked-path
tiling cube: every nonconsecutive pair is free in each copy.
The diagonal 1's are not steps in the Hamiltonian path. They
set the matrix normalization; the adjacent off-diagonal
entries prescribe the path. Rédei gives a path for an
individual tournament but not a common labelled path for
an arbitrary pair of tournaments.

Let t_ij,u_ij be their ±1 orientation signs for i<j.
Then M_ij=t_ij and M_ji=−u_ij. For a given frame the exact
Hadamard test is, for all i<j,

    0 = t_ij−u_ij
        + sum_(k<i) u_ki u_kj
        − sum_(i<k<j) t_ik u_kj
        + sum_(k>j) t_ik t_jk.

This is the row inner product split into the two diagonal
contributions and the three intervals. The pair-products
eta_ij=M_ij M_ji are −1 exactly where T,U agree. Simultaneous
sign switching DMD preserves those products and the diagonal.
After jointly reordering rows and columns, the (+,−) chain
can be imposed by recursive choices of d_(i+1) iff the
undirected eta=−1 graph has a Hamiltonian path. The alternative
(+,+) chain uses the eta=+1 graph; with this decoding convention
then the two directed paths run opposite ways. This is a
criterion on actual framed matrices, not a claim about
unframed metagraph nodes.

**All-order positive frames.** On F_2^m use the sign matrix with
entry (-1)^(x·y+y·y), obtained from Sylvester by column signs.
Its diagonal is +1 and eta_xy=(-1)^wt(x+y). Put the vertices
in reflected binary Gray order; consecutive differences
have weight one. The recursive simultaneous switch above
makes the upper chain +1 and lower chain −1 without losing
orthogonality. This works for every m≥1, but it is not a
circulant construction over Z/(2^m).

## 2. Literal complement and mirror pairs are obstructed

Use the [Klein operations from the p²qr divisor fold](../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md):
kappa flips every free tile (the complement move along the
fixed Hamiltonian path), and sigma mirrors the tile addresses
(converse in reversed path labels). Colors or endpoints in the
merged metagraph are statements after quotienting presentations;
here the input is the actual pair of words in one common frame.

For any M with diagonal +1, write M=I+L+K with L symmetric
zero-diagonal and K skew. Hadamard is equivalent to

    [L,K]=0,     (I+L)^2−K^2=nI,

by comparing MM^T and M^TM. L records disagreements and K
agreements. This immediately forbids a single disagreement
pair in orders n>2: the symmetric 2×2 L block on that pair
is invertible, hence commutation kills all K entries from
it to outside vertices, which would have to be ±1.

**Complement pairing U=kappa(T).** K is the fixed skew
tridiagonal path matrix (upper +1, lower −1). Its first
coordinate vector is cyclic: K^j e_1 reaches exactly j steps
with nonzero last coefficient for 0≤j<n. Any matrix commuting
with K is therefore a polynomial of degree <n in K (choose
a polynomial matching its action on e_1, then commute to
match the entire cyclic basis). Since L is symmetric and
K^T=−K, that polynomial has only even powers. Even powers
do not connect opposite path parities. But for n≥4 the
nonconsecutive pair (1,4) in the complement gluing is a
disagreement, so L_14=±1. Contradiction. This rules out a
literal complement-move witness to a Hadamard gluing, not
all gluings of other representatives of the same endpoints.

**Mirror pairing U=sigma(T).** Let R reverse the vertex order.
Then RMR=2I−M. If M is Hadamard, comparing the two row Gram
matrices forces M+M^T=2I; hence U=T and K=M−I is skew and
anticommutes with R. For n=2m use the normalized plus/minus
basis e_i±e_(n+1−i), i≤m, for R. In this basis K has only
an off-diagonal block A_ij=K_ij−K_i,n+1−j and its negative
transpose. In A each diagonal entry is ±1 and off-diagonal
entries are even. Thus each row of A has squared norm 1 mod 4.
Skew Hadamard would require AA^T=(n−1)I, i.e. n≡2 mod 4,
impossible for a real Hadamard of order n≥4. This also
excludes the fixed (blue) mirror case T=U=sigma(T), not
arbitrary pairs of two mirror-fixed tilings.

**Neither relation is necessary for compatibility.** At n=4
there are three (+,−)-framed gluings, checked by replay (two
literal diagonal pairs and one mixed pair). The mixed pair
has free signs (1,−1,1) and (−1,−1,−1) in order (1,3),(1,4),(2,4).
They are two presentations of the same strong 5-path class.
So even a class glued to itself need not be literal T=U.
If two Hadamard matrices differ in k rows, a changed row
differs in d≥n/k entries: its difference vector has squared
length 4d and is orthogonal to all unchanged old rows, and
its inner product with each of the k remaining old rows has
magnitude ≤2d. Parseval in the old row basis gives
4nd≤k(2d)^2. Hence the compatible pairs form a code with
minimum entry Hamming distance at least n; a single tile
edit, or changing two paired cells, cannot stay inside it
when n≥4. This is a framed-code bound, not a metric bound
on the quotient metagraph.

**All four quotient possibilities are physically visible.** At order
eight the replay gives literal Hadamard gluings whose decoded
endpoints have the indicated path counts:

| Endpoint counts | Relation in the merged metagraph |
|---|---|
| 405, 405 | same node |
| 565, 637 | blue line only |
| 527, 441 | black line only |
| 103, 253 | distinct nonadjacent nodes |

The relation is checked through *all* Hamiltonian frames of the
first endpoint and its converse, with the inherited line rule:
complement its free tiles, discard labels and identify converse;
a mirror-fixed starting tiling gives a blue line. The black/blue
line is witnessed in a **different** frame from the actual
Hadamard pair, as it must be by the complement obstruction.
This finite example answers the possibility question; it does not
classify all compatible pairs or all metagraph edges.

## 3. The third cover repair is physical

Write the profile of a physical block as (P(1),P(2),P(3)).
For ordinal composition, optional cover counts multiply
in **each** color: each color concatenates its local paths
in block order. Take the prime valuations of the three
components, as one finite integer column indexed by (q,prime).
Then an exact scaling by 7 is the target column with entry
1 at (1,7),(2,7),(3,7) and zero elsewhere. A **signed**
integer combination is physically implementable: put the
positive blocks on one side and negative blocks on the other.
No formal inverse of a tournament is needed.

The [literal certificate](witness.json) contains the two
lists of ordinal factors. Subset path DP and optional
cover convolution verify all factors in [replay.py](replay.py).
Both assembled arms have order 137, and their exact q=1,2,3
product ratios are 7,7,7. The q=4 ratio is not 7. A
search used small physical profiles, restricted to a smooth
prime palette, and integer lattice reduction on their
valuation columns. Common singleton and strong factors were
cancelled from both sides afterwards. The certificate
requires no optimizer or trust in a numerical search.

By the [one-color-per-wrap theorem](../packing_hadamard_flatness_20261117/README.md#5-a-concrete-cost-for-the-7-proof-route),
after one common C3[A,·,1] wrap the pair still matches at
q=1,2 and first fails at 3. After two common wraps it still
matches at q=1 and first fails at 2. After a third it fails
at q=1. The common arm A at each wrap may be any nonempty
tournament and may change between wraps. Thus the third
match is genuinely possible and buys exactly one additional
nesting, not a perpetual return.

**Remaining obstruction for the prime route.** In this pair
3 divides both lower H and P(2), so all lower Hamilton
outputs of the first cyclic reader are divisible by 3.
More strongly, the [universal module divisor theorem](../atomic_symmetry_descent_20261117/README.md)
already makes divisibility by 3 hereditary from H alone
for uniform wrapping. Removing common graph factors does
not remove the unwanted common prime. The useful next
version must reach the same valuation target with a lower
profile primitive at (H,P(2)); we do not have that here.
The previous two-color prime-useful switch is still a
separate valid result.

## 4. Fano, Berggren, and the opposite meanings of seven

The [Fano triangle-move archive](../fano_tournament_move_quotients_20261021/README.md)
already distinguishes genuine triangle moves from the marked
complement line, and identifies two Fano designs in Paley7.
Here is an extra exact packing read: use Paley steps
{1,2,4} on Z/7. Translates of {0,1,3} and of {0,2,3} are
two disjoint sheets of seven directed triangles each. In
each sheet a pair of points lies on one line (the six
nonzero differences of each base block occur once), and
any two lines meet. Thus **within a sheet no pair contributes
an x² packing term**. Between sheets each line has exactly
one disjoint partner. Indeed each of its three points lies
on three lines of the other sheet and each pair of points
on one; no line is common, so exactly 9−3=6 other lines
meet it. This gives precisely seven disjoint pairs and
explains the term 7x²y in G_Paley7. Its H is 189, not 7 or
21. Seven here is an orbit size and a matching multiplicity,
not a possible scalar number of Hamiltonian paths.

For N=p²qr with distinct primes, the prior F=U+S work has
F=10 proper divisors, U=3 prime proper divisors and S=7
squarefree proper divisors. Since the primes are *already
inside* the squarefree set, the identity actually says that
the nonsquarefree proper divisors also number 3. The
[3-4-3 / 4-2-4 fold](../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md)
changes the meaning of “fixed” when a tile coordinate is
forgotten. The Fano two-sheet count provides a second
concrete instance: intersection within a sheet obstructs
packing, while the other sheet supplies a unique mate.
The numerical equalities do not identify the group actions.

The [Berggren parent circuit](../../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md)
uses a typed three-coordinate state and integer maps
preserving the Pythagorean equation, not mere counting
of triples. The direct counterpart here is to retain
all three marked reads and preserve their valuation
balance under composition. The analogy ends at selection:
Berggren supplies a total inverse on primitive triples;
our signed lattice finds **one** physical repair, not a
target-dependent inverse for every prime. That is the
remaining Chan–Pak obligation.

Run `python replay.py` next to the JSON. Normal and optimized
Python give the same checks.
