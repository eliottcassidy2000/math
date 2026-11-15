# Internal modules are a mixed-radix address, not an arbitrary phase

An exact specialization of the inherited ordinal painted-block
transfer occurs when every ordinal component is uniformly marked
against a single apex. This gives a second place to change the
mathematical construction, beyond replacing the final outlet.

Let the nonempty blocks be B_1,...,B_s in ordinal order. An apex u
beats every vertex of B_i if w_i=1 and loses to every vertex of B_i
if w_i=0. Write h_i=H(B_i), P_i for its two-slot count with empties,
r_i=P_i/h_i, and D=∏h_i. The rooted source W is this collection plus u.

Define two rational addresses
```
A = Σ_(i: w_i=1) (r_i−1) ∏_(k>i) r_k,
L = Σ_(i<j: w_i=1,w_j=0) (r_i−1)(r_j−1) ∏_(i<k<j) r_k
    + #{ i: w_i=1,w_(i+1)=0 }.
```
The exact spectrum readouts are
```
H(W) = D(1+L),
2E_u(W) = D(A+w_s),
H(W)+2E_u(W) = D(1+A+L+w_s).
```
Here E comprises even-size paths starting at u, weighted by Hamilton
counts of their complements. The formula includes the terminal-mark
boundary; leaving off w_s is wrong even for a single marked block.

## Direct proof of the count polynomial

To count paths through u, split into left and right runs around u.
Each run traverses the selected pieces of B_i in ordinal order. For
nonempty left and right runs, the left finishes on an unmarked block
j and the right starts on a marked block i. Equality i=j is impossible.

If j<i, one needs adjacency j=i−1, since no run can use the interposed
blocks. Every block is then used entirely in its one run and contributes
h_i. This gives D times the number of 0→1 transitions. If i<j, the
right run's piece of B_i and the left run's piece of B_j are nonempty,
giving P_i−h_i and P_j−h_j, respectively. Intermediate blocks may
split arbitrarily into two optional pieces, with P_k choices; outside
blocks have only one possible slot, with h_k path choices. This gives
the corresponding term of L. If either run is empty, the contribution
is D times the appropriate indicator [w_1=1] or [w_s=0]. These two
indicators plus the number of 0→1 transitions equal 1 plus the number
of 1→0 transitions, proving the first formula.

For the even started readout, the part after u must have odd positive
size. Its initial block i is marked and the chosen first piece there
is nonempty. Later blocks may contribute empty or nonempty pieces.
For each nonempty later block the signed sum weighting the chosen
piece by (−1)^size is zero; this is the inherited two-slot parity
cancellation. One proof expands both ordinary Hamilton path counts
by the odd-cycle packing formula; each odd cycle or singleton
chooses either of two slots, canceling independently because its
order is odd. Thus for i<s, the even/odd split is exactly half of the
unsigned total giving the i-term of A. For i=s the same cancellation
includes the empty piece, whose weight is h_s and whose sign is +1.
Removing it gives signed mass −h_s, and hence the extra D in 2E.

## Integer transfer and deterministic aliasing

Update the word one block (h,P,w) at a time. Let C=D A, let N denote
H+2E, and let z be the previous mark. The empty start has
(D,C,N,z)=(1,0,1,0). The complete physical readout updates are
```
D' = h D
C' = P C + w(P−h)D
N' = h N + (2−w)(P−h)C + w(P−h z)D
z' = w.
```
These follow either by subtracting successive explicit sums or by
restricting the inherited four-coordinate transfer to uniform marks.
The selected last reader is part of the state. In particular a block
with w=0 is diagonal on the usual retained-root coordinates.

If r_i=2^ell is a dyadic power, replace B_i by ell transitive
vertices of exactly its mark. This divides both H and the rooted
readout by h_i, with no other change: a transitive ell-block has
(h,P)=(1,2^ell) and therefore the same normalized formulas. Several
such substitutions commute. In particular arbitrary C3 pieces,
whose profile is (3,12), are exact aliases for pairs of equally
marked transitive vertices, up to multiplication of the whole
result by three.

Even for nondyadic ratios, every integer coefficient of the output
as a polynomial in h_i,P_i is separately linear in the two
coordinates of each block. Thus gcd(h_i,P_i) divides the output.
With at least one marked block, N>h_i for every i (H≥D, E>0:
take the path u→v with v in any marked block, and a Hamilton path
exists in its complement). Prime targets therefore need every
block profile primitive, eliminating all inserted C3 modules and
other nonprimitive patches.

An independent exact control is a pair with the same three
coordinates (D,C,z) and therefore reusable under every common
uniform-mark continuation, although its difference is not a
unit. Use one strong-four core K=(h,P)=(5,32) and five singleton
blocks. Order K last and use marks 110110; order K first and use
marks 100000. The two states have (D,C,z)=(5,864,0), and N=1623
and 1711. Every further common uniform-mark suffix multiplies
the difference 88 by its product of ordinary counts, and completely
preserves the agreement of the other three coordinates. Merely
adding periodic changes that use only uniform-mark dyadic modules
cannot supply a new primitive ratio, even if their physical sizes
and commonly observed symmetries differ.

These are exact return and exclusion statements for a restricted
module interface. They do not show that any choice of mixed-radix
digits admits an absolute target; that is a physical inverse
problem, and the finite modulus tests need not settle it.

## A unit difference that really preserves continuation

Changing the non-dyadic core itself, not just moving a fixed core, gives
this exact paired construction. Each word has sixteen blocks, with
fifteen singletons and a single strong six-vertex block.

| final H+2E | special block index (0-based) | special outgoing rows | (h,P) | marks mask, least-significant bit first |
| ---: | ---: | --- | --- | ---: |
| 9,971,529 | 10 | (34,8,3,5,15,30) | (23,272) | 30,557 |
| 9,971,531 | 14 | (48,1,3,7,14,30) | (23,296) | 19,413 |

Both states have (D,C,z)=(23, 6,510,178, 0). The two ordinary
rooted remainders have H=3,461,351 and 3,461,353, respectively,
and both have E=3,255,089. Both are strong because their first
mark is 1 and final mark is 0. The small-core profiles and the
full counts have independent checks in the replay: it exhausts
all literal partitions of the 21 base vertices on the two sides
of the apex, selecting only the feasible endpoint orientations
and reading each induced order-six block by exact subset DP.

**Any shared continuation by uniform marked blocks preserves the
agreement of (D,C,z).** In particular any binary word of singleton
suffixes gives another physical pair of adjacent odd counts.
The underlying unrooted tournaments stay strong when the common
suffix ends in zero. For example, t added unmarked vertices give
the exact pair
```
13,020,356 * 2^t − 3,048,827,
13,020,356 * 2^t − 3,048,825,  t≥0.
```
This time it is the *full outgoing state for this context*, rather
than merely the affine slope in the reported count, that is locked.

The pair is initialized, not repeatable after arbitrary input.
Each word induces an affine map on the incoming C coordinate with
coefficient equal to the product of its slot counts. These two
coefficients are 272·2^15 and 296·2^15. The block words match C at
the empty prefix, but no common prefix with C>0 can preserve that
same lock. Thus it would be invalid to concatenate these as if
they independently supplied a unit interval; a family matching
at varying incoming C/D, or with equal products of slot counts,
is still needed.

## Length and reciprocal length are separate continuation moments

Let T be the product of the normalized block ratios, and put
t_0=1, t_i=r_1···r_i. Let M denote the union of intervals
[t_(i−1),t_i] having mark 1, disregarding boundaries of measure zero.
Write |M| for its ordinary length. The two potentially different
affine continuation coefficients have simple integral readings:
```
A = T ∫_M dx/x²,
coefficient of incoming C in outgoing N
  = D [2(T−1)−|M|],
L = ∫_(1≤x<y≤T) 1_M(x) (1−1_M(y)) dx/x² dy + #{1→0 transitions}.
```
The A formula follows by integrating on one marked block.
For the C coefficient, incoming C is multiplied by the product
of the P's before each block. The subsequent propagation into N
uses the product of the h's after that block. Divide their
weighted sum by D to obtain
Σ_i (2−w_i)(r_i−1)t_(i−1) = 2(T−1)−|M|.
For the L formula, a marked interval i contributes reciprocal
length (r_i−1)/t_i and an unmarked later j contributes length
(r_j−1)t_(j−1), giving exactly the pair interaction. Therefore,
at fixed D,T and common initial/final marks, a repeatable unit
patch asks for two physically admissible arrangements with equal
length and reciprocal length, but a difference 2/D in the ordered
interaction L. The initial mark is needed to match the dependence
on the incoming previous-mark bit.

This distinguishes ordered interaction from its two marginal
moments. It also rules out a tempting solution by plain reversal:
under reversing the block order, normalized times transform by
x ↦ T/x and exchange the two weighted measures in the pair term.
The difference of the pair terms is
∫_(x<y) [1_M(x)−1_M(y)] dx/x² dy = A−|M|.
The difference in transition counts is w_1−w_s. Thus if reversal
matches both moments and preserves the boundary marks, it cannot
supply a unit difference: A=|M| implies L is unchanged. Reversing
and complementing jointly preserves L unconditionally by pairing
each marked→unmarked term with its reversed complemented term.
If that operation matches the continuation moments and final mark
it also cannot change the count. A repeatable unit needs genuinely
different ordered correlations rather than just one of these
symmetries; the initialized unit above does not claim to supply it.
