# Dyadic two-cover amplifiers and their exact boundary

**PROVED identities and all-size composition; FINITE-EXACT input classification
from the accepted October 31h census. No full spectrum conclusion.**

## 1. Actual marked state and an empty-cut block

For a parent tournament A and an external marked vertex x with out-neighbors
S in A, retain

`F=H(A+x), h=H(A), g=# full paths of A ending in S,`

`q=# ordered spanning two-path covers (P,Q), Q nonempty, first(Q) in S.`

P may be empty. For a block B, E(B) counts all ordered spanning two-path covers
with either slot allowed empty. Append B after A in the parent ordinal sum;
thus every old parent vertex beats every B vertex. The **empty cut** means
every B vertex beats x. Directly from the inherited positive block action,

`F' = H F + (E-H)q + H g,  q'=E q,  g'=0,  h'=H h.`             (A1)

For the deficit c=q-F-g this becomes

`(c,q,h) -> (Hc, Eq, Hh).`                                     (A2)

This diagonal law is scalar response compression. It does not identify
different block orders as tournaments, or license later arbitrary re-marking.
Keep the actual ordered block word and all selected endpoints for those tasks.

## 2. Explicit family from the recovered nineteen atom

Adjacency rows are outgoing-neighbor bitsets, indexed from zero:

| Atom | Rows | vertices | H | E |
|---|---|---:|---:|---:|
| D | 8,1,3,6 | 4 | 5 | 32 |
| B19 | 38,60,56,49,1,16 | 6 | 19 | 256 |
| singleton | 0 | 1 | 1 | 2 |

Both nontrivial atoms are strong. Literal permutation and cut enumeration in
`verify.py` confirms H and E. The atoms themselves and D's response formulas
are inherited data; their joint amplification is the present application.

The inherited optimal deficit-nine source has parent rows

`12,121,2,100,45,69,21`

and x beats the vertices selected by bitmask 67. Its actual marked state is

`(F,q,g,h)=(449,512,54,117), c=9.`                              (A3)

It uses eight vertices. Append a copies of D, b copies of B19 and r singleton
blocks, always with empty cut, in any chosen order. Assume a,b,r>=0 and
a+b+r>=1. By (A2), after at least one block g=0, and

`q=2^(9+5a+8b+r), c=9*5^a*19^b, h=117*5^a*19^b.`

Therefore, with k=9+5a+8b+r,

`H(T)=2^k-9*5^a*19^b,    |V(T)|=8+4a+6b+r=k-a-2b-1.`          (A4)

This proves the family for every parameter choice, not merely for sampled
graphs. When a+b>0, every k>=9+5a+8b is allowed; when a=b=0, k>=10.
The omitted zero-operation case is essential: the seed has F=449, not
512-9=503, because its g=54 has not yet been reset. No minimum-order claim
for the constructed counts is made.

Examples: a=1,b=r=0 gives H=16339 on12 vertices; a=r=0,b=1 gives
H=130901 on14 vertices. Seventeen parameter controls replay the formulas;
the nine with at most16 vertices also use independent Hamiltonian subset DP
on the full adjacency matrices.

## 3. Only the empty cut preserves a dyadic slope

**Proposition.** Suppose an actual input has q>0 a power of two, and B has
E a power of two. Among cuts in B, q' is a power of two iff the cut is empty.

**Proof.** Cutting a Hamilton path immediately before each v in S injects
|S|h objects into q covers. Thus q>=|S|h, and q>0 implies h<=q. For a
nonempty block cut, its q digit Q is positive by the same injection. Covers
with second slot empty contribute H and are excluded from Q; hence
0<Q<=E-H<E. The general action q'=Eq+Qh gives

`Eq < q' < 2Eq.`

There is no power of two strictly between these successive powers. For the
empty cut Q=0 and q'=Eq is dyadic. QED. This is an all-size result for actual
inputs; formal states lacking h<=q need not satisfy it.

## 4. Complete bounded-atom classification, unbounded composition

The accepted 31h profile census covers **every strong tournament with H<=127**.
Busch's order-eleven minimum225 closes the size boundary. Two independent
augmentation routes give the same profile set. The retained filter Q<=256,
F_B<=512 cannot remove an empty cut: it has Q=0 and F_B=H<=127.

Scanning both files, whose hashes are in `verification.json`, yields precisely
these dyadic-E profile pairs:

`(H,E) = (1,2), (5,32), (19,256).`                             (A5)

This is a classification of numerical profile types, not a claim of one
isomorphism class for each pair. The nontrivial types occur at orders4 and6.
Their sizes are forced: the strong floor at order7 is25, the order-five
maximum is15 in the inherited small census, and the only nontrivial strong
count below9 other than3 is D's5.
The receipt reuses the accepted complete census; it does not claim a new
nauty enumeration.

For any ordinal sum A=A1 -> ... -> At, both H and E multiply. For E, restrict
each of the two ordered paths to each component. Conversely, concatenate the
first slots and second slots in component order; every connecting arc points
forward. These are inverse bijections, including empty slots.

Consequently, **for every parent at every size whose nontrivial SCCs each
have H<=127**, E is dyadic iff every SCC has one of the profile types (A5).
Writing s for its singleton count and a,b for the two nontrivial types gives

`H=5^a19^b,  E=2^(s+5a+8b),  |V(A)|=s+4a+6b.`                 (A6)

Thus an entire infinite attachment bank is classified from the bounded atom
census. Larger-H strong atoms are outside this conclusion.

## 5. Why this does not prove the remaining gaps theorem

For (A4),

`c/q=(9/512)(5/32)^a(19/256)^b 2^(-r) <= 9/512.`

Hence H(T)>=(503/512)2^k. If H(T)<=X, k=O(log X), and a,b,r each have
O(log X) choices. There are at most O((log X)^3) distinct outputs. The
same counting argument applies to any fixed finite bank of dyadic empty-cut
amplifiers and finitely many actual positive-dyadic-slope seeds, with a fixed
polynomial exponent: their
H/E ratios are at most1/2 because E>=2H.

The amplifier lane is therefore a source of explicit infinite families, not
a possible cofinite-coverage proof on its own. To cover the missing prime
lanes one needs unbounded source diversity, nonempty-cut excursions with
certified return, or genuinely different operations. The new inverse in
[DESCENT.md](DESCENT.md) retains those excursions within its declared grammar;
it still needs a uniform target-entry theorem.
