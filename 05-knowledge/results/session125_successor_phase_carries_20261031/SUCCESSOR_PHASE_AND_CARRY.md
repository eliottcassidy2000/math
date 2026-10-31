# Successor pairs, finite hidden phase, and unavoidable radix-49 carry costs

**Status:** PROVED structural consequences, conditional only on the stated
original-port Q49 constructor assumptions; FINITE-EXACT short-order input
from the independently enumerated 41/42 catalogue. Exact constants and sample
counts for Gerbicz's original tables are VERIFIED separately. This is not a
new proof of the classical all-order square-sum-cycle existence theorem.

## 1. Three meanings of successor must stay separate

Write `Pair(n)=(S_n,L_n)` for an original nice pair. Here `S_n` is a
Hamiltonian path on `[1,n]`, `L_n` one on `[1,n+1]`; both start at 1, finish
at 3 for odd order and 8 for even order, and every common label has the same
position parity in the two paths. Positions start at 1. Closing either path
is permitted because `1+3=4` and `1+8=9`.

1. **Independent neighboring realizers:** for each n, choose the two members
   of `Pair(n)`. This is enough to exhibit cheap pairs of consecutive cycles.
2. **A single-valued canonical sequence:** choose `C_n=close(S_n)` for every n.
   There is no reason that `L_n=S_(n+1)`; the difference is a boundary carry.
3. **A globally phase-coherent sequence:** choose paths `P_n` so every pair
   `P_n,P_(n+1)` shares position parity. This is much stronger than item 2,
   and is already impossible from order 41 through order 44.

The results below make these distinctions quantitative. They do **not**
exclude cheap consecutive edits for every conceivable cycle family.

## 2. Phase coherence forces finite balanced memory

**Theorem (finite phase defects).** Suppose an original-port Hamiltonian
path `P_n` is chosen for every `n>=n0`, and every adjacent pair has the same
position parity on common labels. Then there is a fixed phase function
`phi(v)` with:

- `phi(v)=v mod 2` for every `v>n0`;
- the defect set `D={v:phi(v) != v mod 2}` lies in `[1,n0]`;
- D contains equally many odd and even labels;
- `D` excludes 1, 3, and 8;
- `|D|>=6`.

**Proof of finite memory and balance.** Common-label agreement makes each
label's phase permanent. Among positions `1,...,n`, exactly `ceil(n/2)`
are odd. Subtract the odd-position count at n from the count at n+1: the
new label's phase must be its natural parity. At n0, the number of even
labels flipped into the odd-position class equals the number of odd labels
flipped out. The fixed endpoints force the three port phases to be natural.

**Proof of the six-defect lower bound.** Two adjacent nondefects have opposite
natural parity, so their sum is an odd square, hence 1 modulo 8. The outside
vertices therefore have four disconnected residue channels:

```
A={0,1}, B={2,7}, C={3,6}, D8={4,5}  (modulo 8).
```

A defect adjacent to a nondefect has the **same** natural parity as its
neighbor, so their sum is an even square, 0 or 4 modulo 8. A defect `d`
can attach only to residues `-d` and `4-d` modulo 8. Consequently it attaches
only to the superpair `A union D8` when `d mod 4` is 0 or 3, and only to
`B union C` when `d mod 4` is 1 or 2. No one defect can connect the two
superpairs.

Take a sufficiently large even n that each outside residue channel is
nonempty. The path starts at 1 and ends at 8, both nondefects in A. After
deleting all defect vertices, its surviving contiguous runs each lie in
one residue channel. There must be at least five runs: A at both ends,
plus at least one each in B, C, D8. Thus there are at least four intervening
nonempty defect runs. Moreover the path must leave superpair `A union D8`
to visit `B union C` and later return. The two crossing defect runs each
need at least two vertices. Therefore the total number of defects is at
least `4+2=6`. Endpoints are not discarded or counted as defects. For odd
orders, endpoints A and C give only the weaker lower bound 4. This is a
necessary bound, not an existence or sharpness claim for six defects. QED.

## 3. The first nice pair does not start a coherent phase family

The independent onset census in `../session125_origin32_numeration_20261031/origin_pairs/check_independent.py`
proves that a nice pair at 41/42 has exactly three possible common odd-position
masks. Bit `v-1` records label v. The three masks and their first impossible
fixed-phase continuation are:

| Mask | First failure | Short certificate |
|---:|---:|---|
| 1187208845135 | 44 | Degree saturation at 7 and 20 blocks internal vertex 29 |
| 1453217412461 | 44 | New internal vertex 44 has only allowed neighbor 37 |
| 1487448155501 | 43 | New internal vertex 43 has only allowed neighbor 38 |

The extension rule is forced: labels 42, 43, 44 must have natural parity.
For the first mask at n=44, the phase-compatible neighborhoods are

```
N(18)={7,31}; N(42)={7,22};
N(16)={9,20}; N(44)={20,37};
N(29)={7,20,35}.
```

All five displayed vertices are internal: the endpoints are 1 and 8.
The first two degree-two vertices force both incident path edges at 7;
the next two force both incident path edges at 20. Vertex 29 can then use
only 35 and cannot have internal path degree two. This requires no search
at order 44. The other two failures are immediate internal degree-one
obstructions. The checker includes valid fixed-phase paths at 43 for the
first two masks, so their stated first failures are exact.

**Corollary.** No original-port, common-position-parity path family exists
on the whole interval 41 through 44, hence none exists for every n>=41.
This does not rule out a different later starting order, relaxed ports, or
a single-valued family with phase changes.

## 4. Cheap independent pairs: exact within-pair recurrence

For any two paths define the old-edge deletion distance
`delta(P,Q)=|E(P) minus E(Q)|`, with undirected edges. If orders differ by
one, the number of added edges is one more than the deleted count.

The original Q49 lift uses 49 complete path blocks indexed by
`c=-24,...,24`, together with 24 singleton atoms. In a block, input label x
at position j is sent to

```
49*x + (-1)^(j+1)*c.
```

For output orders `n=49m+d` and `n+1`, where `24<=d<=72`, block c chooses
the longer input member precisely when

```
49 + (-1)^m*c <= d+ext,   ext=0 or 1.
```

Thus only the block `c=(-1)^m*(d+1-49)` changes input member between the two
outputs. Reversing any other block changes no undirected internal edge.
The common input phase makes the label map agree on shared vertices in
the changed block. Also the partition of all old labels into atoms is
unchanged, so an old internal edge cannot turn into a connector or vice
versa.

Let `d_pair(n)=delta(S_n,L_n)`. The 73 atoms have 72 connector edges. If
`K(d,m mod 2)` is the number of old connector edges absent from the new
word, then the following recurrence is **exact**:

```
d_pair(n)=d_pair(floor((n-24)/49))+K(d,m mod 2),  0<=K<=72.
```

There are finitely many base pairs at 41,...,2032. Hence
`d_pair(n)=O(log n)`. Closing the two paths changes the upper bound by at
most one old edge, so the resulting cycles also have `O(log n)` deleted
edges. This proves an existential neighboring-pair statement, not a cheap
single canonical family.

For the pinned original tables, the independent edit audit finds:

- maximum base path deletion count 1699, at n=2012;
- connector deletion constants from 9 through 44, over all 98 states;
- therefore `d_pair(n)<=1699+44*recursion_depth(n)`.

## 5. Full-offset lifting forgets phase, but amplifies disagreement

Put

```
H(t)=delta(L_(t-1),S_t).
```

Both paths have order t and the same endpoints. This is a disagreement
between two realizers of **the same** integer, not the within-pair distance.
At a noncarry boundary in the recursive range, the two words select the
same member at every offset. Only connectors can change, so `H(t)<=72`.

At a carry boundary `t=49r+24`, with `r>=42`, the longer output of
`Pair(t-1)` uses the longer member of `Pair(r-1)` at **all** offsets; the
shorter output of `Pair(t)` uses the shorter member of `Pair(r)` at all
offsets. These two parent paths may have different phases.

For any undirected parent edge `{x,y}`, the union of its 49 lifted edges
is nevertheless

```
{{49*x+c,49*y-c}: -24<=c<=24},
```

independent of which endpoint has odd position: replacing c by -c leaves
the set unchanged. Different parent edges lift to disjoint edge sets,
because nearest-center decoding recovers x and y uniquely. Thus the
full-offset lift multiplies internal-edge deletion distance by exactly 49.

There is also **no collision between either internal lift and either
connector set**. A connector incident to a singleton has an endpoint at
most 24, absent from internal lifts. A connector between two block endpoints
decodes to parent labels from `{1,b}`, where b is the common terminal port
3 or 8. Equal decoded labels cannot form a parent edge; distinct decoded
labels give `{1,b}`, which is absent from a Hamiltonian path of order
r>2 having endpoints 1 and b. Thus the two kinds of edges remain disjoint
even when the parent phases differ. If `K_carry` counts old connector
edges absent from the new connector word, the exact recurrence is

```
H(49r+24) = 49*H(r)+K_carry,  0<=K_carry<=72.         (1)
```

For fixed original transition words, `K_carry` depends only on `r mod 2`.
The independently parsed original tables give `K_even=40`, `K_odd=52`.
The weaker `49H-72 <= H_new <=49H+72` is valid but unnecessarily loses
information; the endpoint decoder proves the sharper nonnegative correction.

This argument explicitly allows the parent phases and atom partitions to
differ. Pretending that individual offset blocks remain aligned at a carry
would be incorrect; the **union over all signed offsets** is the needed
invariant.

## 6. Universal linear worst-case carries for this constructor

This conclusion is stronger than a defect in Gerbicz's particular stored
base choices. It holds for **any** original-port nice base pairs at
41,...,2032 and any valid Q49 whole-block words with the forced member
choices above.

If `H(42)=H(43)=0`, then `L_41=S_42` and `L_42=S_43` as oriented paths:
undirected equality plus the fixed starting vertex 1 determines the
orientation. These identifications glue Pair(41), Pair(42), Pair(43) into
a globally phase-coherent family at 41,42,43,44, contradicting Section 3.
So at least one `t0` in `{42,43}` has `H(t0)>0`. Two distinct Hamiltonian
paths with the same endpoints cannot differ by just one deleted edge:
the equal vertex-degree vectors would force its replacement to have the
same two endpoints. Hence `H(t0)>=2`.

Define the carry ray

```
t_k=49^k*(t0+1/2)-1/2.
```

Iterating (1) gives

```
H(t_k) >= 49^k*H(t0).
```

Now choose the single-valued canonical cycles `C_n=close(S_n)`, and write
`D_can(n)=|E(C_n) minus E(C_(n+1))|`. Let `D_pair(n)` be the analogous
within-pair cycle distance. The common closing edge at order n+1 implies
that the same-order cycle representation gap is exactly `H(n+1)`. The
set-difference triangle inequalities and the one-edge size difference give

```
H(n+1)-D_pair(n)-1 <= D_can(n) <= H(n+1)+D_pair(n).    (2)
```

Since `D_pair(n)=O(log n)`, equations (1)-(2) imply

```
limsup_(n->infinity) D_can(n)/n >= 4/87.
```

Indeed along `n=t_k-1`, the limiting lower ratio is at least
`2/(43+1/2)=4/87`. Thus this canonical short-member family has
**linear worst-case successor edit cost**, despite its cheap independent
nice pairs. The claim is about the specified original-port Q49 whole-block
construction, not every possible Hamiltonian-cycle sequence and not a
universal impossibility of logarithmic repairs.

For the original pinned tables, `H(42)=18` and the even carry correction
is always 40. Since the carry map preserves parity, this gives the exact
formula

```
H(t_k)=(113*49^k-5)/6,   t_k=(85*49^k-1)/2.
```

The two-sided bound (2) then proves the actual limit along this ray:

```
lim_(k->infinity) D_can(t_k-1)/(t_k-1)=113/255.
```

Thus about 44.3% of old cycle edges are deleted at these canonical carries.
This exact limit concerns the original pinned tables, unlike the universal
`4/87` lower bound. Exact samples:

| n | Within-pair path deletions | Canonical cycle deletions n to n+1 |
|---:|---:|---:|
| 2033 | 46 | 68 |
| 2080 | 28 | 36 |
| 2081 | 27 | 921 |
| 2082 | 35 | 44 |
| 98170 | 1697 | 90936 |
| 102041 | 48 | 45219 |

The corresponding same-order representation gaps are
`H(42)=18`, `H(2082)=922`, `H(102042)=45218`.

### 6a. Choosing Shorter or Longer at each order does not remove the obstruction

Let a **two-representative selector** choose, independently at every order n,
either `close(S_n)` or `close(L_(n-1))`. This obvious repair retains the
full universal lower bound, not merely half of it:

```
limsup_(n->infinity) deletion_cost(n)/n >= 4/87.
```

For the pinned original tables, the 42 carry ray gives at least `113/255`.
These are bank-restricted statements, not lower bounds for every conceivable
cycle family.

Let r be a sufficiently large carry, write `B_n=D_pair(n)`, and call the
left and right representatives at order r `U=close(L_(r-1))` and
`V=close(S_r)`. Their distance is H(r). The neighboring orders r-1 and r+1
are not carries, so their two bank representatives differ by at most 72.
For any neighboring bank choices X at order r-1 and Y at order r+1:

- If the middle choice is V, then
  `delta(X,V)>=H(r)-B_(r-1)-73`.
- If the middle choice is U, then
  `delta(U,Y)>=H(r)-B_r-73`.

For the first inequality, use
`H(r)<=delta(U,X)+delta(X,V)` and
`delta(U,X)<=B_(r-1)+1+72`. For the second, use
`H(r)<=delta(U,Y)+delta(Y,V)` and
`delta(Y,V)<=72+B_r+1`. Thus one adjacent step costs `H(r)-O(log r)`
regardless of the middle bank choice. This proves the full lower bound.

There is also a uniform upper bound for **every** bank-selector step:

```
delta(C_n,C_(n+1)) <= H(n)+B_n+H(n+1).
```

At most one of H(n),H(n+1) can be a large carry; the other is at most 72
outside the finite initial boundary. Consequently all bank selectors share
the same leading worst-carry rate when the original fixed tables are used;
Section 6b identifies it exactly.

**More freedom only at the middle order.** If the middle cycle W at a carry
is allowed to be completely arbitrary, while its two neighbors remain bank
choices, one still obtains a half-strength bound. Put
`A=close(S_(r-1))`, `Z=close(L_r)`. Directed triangles give

```
delta(A,Z) >= H(r)-B_(r-1)-B_r-2,
delta(X,Y) >= H(r)-B_(r-1)-B_r-146.
```

Since `delta(X,W)+delta(W,Y)>=delta(X,Y)`, one of these two steps costs
at least `[H(r)-B_(r-1)-B_r-146]/2`. This yields universal `2/87` and
original-42-ray `113/510` lower ratios under this local relaxation, but is
not a statement about unrestricted choices at all orders.

The edit checker includes a bounded finite control: for five orders around
each of carries 2082 and 102042, a two-state minimax dynamic program chooses
the cheapest worst-step bank sequence. It is checked against all 32 possible
selector words. The optimum maximum deletions are 916 and 45214,
respectively; the corresponding H values are 922 and 45218. This is an
exact finite bank-selection control, not a search over arbitrary cycles.

### 6b. The original constructor's exact global worst-carry rate

For the pinned original tables, **every** two-representative bank selector,
including the canonical short-member choice, has the same global limit
superior:

```
limsup_(n->infinity) deletion_cost(n)/n = 22273/24042.
```

This is approximately 0.92642043: along some carry subsequence, over 92.6%
of the old cycle edges must be deleted inside this constructor bank.
The value is not asserted for a different set of base certificates or an
arbitrary family of square-sum cycles.

Here is why the exact global answer reduces to a finite calculation.
Repeatedly undo `t ->49t+24` until reaching either a root s in `[42,2033]`
or a larger noncarry root. The mixed boundary root 2033 must be included:
its previous representation comes from a base pair while its new one is
recursive, so the generic noncarry bound does not cover it. For a root s,
the limiting carry ratio is exactly

```
R(s)=[H(s)+K_(s mod 2)/48]/[s+1/2].
```

Every noncarry root `s>=2034` has `H(s)<=72`, and therefore
`R(s)<=877/24414`. The finite scan of all 1992 roots 42 through 2033
finds a larger maximum at s=2003:

```
H(2003)=1855, K_odd=52, R(2003)=22273/24042.
```

The special mixed boundary has `H(2033)=1861`. Exact rational arithmetic
compares every root, not just a plotted sample. Along the maximizing ray,

```
t_k=(4007*49^k-1)/2,
H(t_k)=(22273*49^k-13)/12.
```

The exact recurrence bounds every other carry by the root maximum up to a
uniform `O(1/t)` ratio error. The `O(log n)` within-pair bound and Section
6a's lower/upper bank inequalities then prove the stated global equality
for every bank selector. The independent auditor and the producer both
compute the same maximizing root, rational value, and outside-root bound.

## 7. Numeration interpretation and limits

The centered carry map `t -> 49t+24` becomes simple multiplication after
adding one half: `2t+1 -> 49(2t+1)`. Carry depth is therefore controlled by
the power of 49 dividing `2t+1`, subject to stopping at the finite base
range. This is a concrete hidden digit coordinate, not merely an analogy.
The order labels 32 onward form the visible existence range, while the
port/phase data and the disagreement between duplicate realizers determine
how costly successor can be.

This resembles carry propagation in constrained numeration, including the
motivation from Zeckendorf decomposition. No map to Fibonacci digits,
Zeckendorf uniqueness theorem, or Fibonacci dynamics has been established
here; the proved coordinate in this packet is centered radix 49.

## 8. Replay, external data, and independent checks

```
python check_fixed_phase.py --origin-dir ../session125_origin32_numeration_20261031/origin_pairs
python -O check_fixed_phase.py --origin-dir ../session125_origin32_numeration_20261031/origin_pairs
python q49_edit_audit.py /path/to/gerbicz_squares.c
python -O q49_edit_audit.py /path/to/gerbicz_squares.c
```

The phase checker authenticates the sibling compressed reference through
its published `read_reference()` helper; no unpacked origin-certificate JSON
is required. It then reruns the independent finite degree-subset enumeration
at 41 and 42 (49 and 275 admissible paths), recovers the three masks,
validates the positive 43 paths and short negative neighborhoods, checks
all residue attachments, and checks hostile phase/endpoint mutations.
The weighted four-channel check returns lower bounds 6 for A-to-A and 4
for A-to-C, in agreement with the proof.

The Q49 checker reads only inert integer initializer data from the source
whose SHA256 is
`602c1919ee51101519a2f19f2c1817407006824f6f09c85a3b772a818a0b736e`.
It never imports, compiles or executes external C. It checks every base
edge distance, all 98 one-block changes and connector constants, selected
expanded paths/cycles and parity, exact recurrence samples, and canonical
triangle bounds, the 32-choice selector controls, and the finite exact global-limsup maximum. Both normal/-O receipts agree byte-for-byte.

The original source and its arrays are **not redistributed**. Source
acquisition and earlier full validity audits are documented in
[SOURCE.json](SOURCE.json) and the prior
[Gerbicz reconstruction](../session124_gerbicz_reconstruction_20261031/GERBICZ_RECONSTRUCTION.md). Thus reproducing the original-table numerical
samples requires that separate pinned source. The universal carry theorem
uses only the stated Q49 schema and the independent finite 41/42 catalogue,
not the numerical sample receipt or any particular base choices.
