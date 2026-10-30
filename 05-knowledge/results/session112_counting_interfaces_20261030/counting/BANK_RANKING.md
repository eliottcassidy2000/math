# Canonical ranking for the shared-residue bank action

**PROVED:** a constructive exact rank/unrank algorithm, polynomial in the number
D of bank coordinates and in `log Q`, for the specified named permutation
action below. **FINITE-EXACT:** independent necklace and literal-group controls
in `bank_rank.py`. This is a bank-orbit algorithm, not a full tournament
isomorphism algorithm. It uses neither factoring nor the new paper's counting
engine, and makes no historical-priority claim for necklace ranking.

## 1. Fixed interface, order and inherited action

Use the action from Session 111's
[`BANK_RECURRENCE.md`](../../session111_scale_connections_20261030/BANK_RECURRENCE.md),
not an extrapolation to arbitrary graph symmetries. Let p be an odd prime,
`k≥2`, `Q≥2`, `n=p^(k−1)` and `h=(p−1)/2`. The coordinate permutation group is

    G = C_n × C_p,       elements (a,b).

Its named coordinate cycles are:

- h top A-cycles of length n, shifted by a;
- h top B-cycles of length n, shifted by `a+(n/p)b`;
- p−1 lower A-cycles of length p^j for each `0≤j≤k−2`, shifted by a.

Their total length is `D=p^k−1`. Fix **this block order**: top A-cycles, top
B-cycles, then lower cycles in increasing j; fix the named order of equal
copies and a cyclic base point in every cycle. Inside a block use positions
0 through L−1. Color digits in `{0,...,Q−1}` are read most significant first.
The canonical representative is the lexicographically least whole vector
under G. Equal-sized named cycles are not permuted by this action.

For the ramified arithmetic interpretation, use the source's hypotheses
`p≥5` and simple metallic ramification. In binary orientation coordinates the
odd affine action becomes unsigned only after choosing an invariant origin.
This document orders coordinates **after that fixed-origin change**, not in
an unrecorded original orientation gauge. Transporting answers back requires
retaining that origin and the cycle bases. The algorithm and checker act on
the displayed explicit model; they do not silently reconstruct a raw ring
labeling. General Q is an ordinary coloring alphabet on this model.

The smaller p=3 experiment below is only a permutation-model control. It is
not an extension of the arithmetic ramification theorem to p=3.

## 2. Only k+2 block-boundary stabilizers are needed

After fixing any sequence of whole named blocks, its stabilizer belongs to:

    C_r = {(a,b): p^r divides a},       0≤r≤k−2;
    K_A = {(a,b): a=0};
    K_B = {(a,b): a+(n/p)b=0 mod n};
    I   = {(0,0)}.                                            (1)

These are k+2 subgroups. This closure is a statement at **completed block
boundaries**; it does not say that arbitrary incomplete prefixes have only
k+2 states.

**Proof.** A word on a cycle of length p^j has minimal rotation period p^v
for some `0≤v≤j`. Fixing it adds the congruence `p^v | a` for an A-block or
`p^v | a+(n/p)b` for a B-block. A B-block is top-sized. Below exponent k−1,
the latter condition is just `p^v | a`, since n/p is divisible by p^v.
All such conditions combine by their maximum exponent and give C_r.
A full-exponent A condition gives K_A, a full-exponent B condition gives K_B,
and imposing both gives I. Any lower condition is already automatic in K_A
and K_B. This proves both closure and the list.

For H from (1), the shift image on a block of length `L=p^j` is
`p^q Z/LZ`. Its q value is:

- H=C_r: `q=min(r,j)`;
- H=I: q=j;
- H=K_A on A, or H=K_B on B: q=j;
- the opposite top kernel: `q=min(k−2,j)`.

The last case on a lower A-block gives q=j and hence trivial action. These
formulas retain the shared residue between the top A and B shifts.

## 3. A necklace-period transfer refines the total recurrence

Suppose the current block has image shifts `p^q Z/LZ` under H. Group its
consecutive digits into packets of size p^q. The action becomes ordinary
cyclic rotation of

    m=p^(j−q) packets, over alphabet A=Q^(p^q).

Packet order agrees exactly with the original word order. Each orbit has a
unique lexicographically least packet word, a necklace representative.
Its least packet period is `ell=p^s` for `0≤s≤j−q`.

Its stabilizer in H is

    H_ell = H ∩ {g: shift(g)=0 mod p^(q+s)},                    (2)

which is another subgroup in (1). Crucially H_ell depends on the period, not
on the individual necklace. The number of such necklaces is

    P_A(1)=A,
    P_A(p^s)=(A^(p^s)-A^(p^(s−1)))/p^s  for s≥1.             (3)

Indeed words of length p^s whose least period is smaller are exactly the
words repeated from length p^(s−1). A primitive word has p^s distinct rotations.

Let `T_i(H)` count orbits of blocks i onward under H. Then

    T_end(H)=1,
    T_i(H)=sum_{ell|m} P_A(ell) T_(i+1)(H_ell).                (4)

**Proof.** First select an H-orbit of the current block and put it into its
least form. Whole tuples with that first block are then identified exactly
by its stabilizer in H. Any two transporters bringing the block to that same
least form differ by this stabilizer. Thus there is neither a missing copy
nor an extra transporter multiplicity. Group first-block necklaces by their
period and use (2)–(3).

This is a stabilizer-refined version of the earlier scalar Burnside recurrence.
It has O(pk) block positions, k+2 subgroup states and at most k period choices
per transition. It is already an efficient total-count algorithm, but ranking
needs period-resolved counts **below a specified word**. The next section
supplies them without enumerating necklaces or saturating a full threshold
set under G.

## 4. An exact threshold automaton for necklaces

Fix an alphabet size A, word length m and threshold word
`t=t_0...t_(m−1)`. A bi-infinite periodic word has every length-m window at
least t exactly when it avoids every finite forbidden pattern

    t_0...t_(j−1) a,        0≤j<m and a<t_j.                   (5)

This follows by the first differing digit. Patterns in (5) include short
first-difference witnesses, not merely length-m words.

Construct a deterministic prefix-suffix automaton with states `0,...,m−1`.
State s records the length of the longest suffix of the scanned history
which is a **proper** prefix of t. Its active shorter suffixes are exactly

    J_s = {0} ∪ {j≤s: t_(s−j)...t_(s−1)=t_0...t_(j−1)}.

The next symbol a is allowed iff

    a ≥ max_{j∈J_s} t_j.                                     (6)

Otherwise it finishes a forbidden pattern. After an allowed symbol, the
new state is the longest proper prefix of t that is a suffix of `t[:s]a`.
A complete equality with t is allowed; then use its longest proper border.
This transition retains every relevant suffix because all shorter active
prefixes are suffixes of the longest one.

Let M_t be the m×m matrix counting allowed alphabet symbols for each state
transition. The alphabet may be large: only the at most m distinct symbols
appearing in t require individual transitions. Every other allowed symbol
goes to state 0. Thus M_t is constructed using integer counts and comparisons,
not A separate branches.

### Trace lemma, including the circular-state obligation

For every positive ell,

    tr(M_t^ell)

is exactly the number of length-ell words whose infinite periodic repetition
avoids (5).

**Proof.** Matrix powers count labelled closed state walks. Every such walk
can be repeated arbitrarily many times. The state after m−1 symbols is fully
determined by those last m−1 symbols, because it records only proper prefixes
of length at most m−1. Thus an admissible periodic label sequence has **one**
periodic state sequence, not several; its state at each cut is determined by
the infinite past's last m−1 symbols. Conversely every label sequence avoiding
(5) supplies that unique state sequence and all its transitions are allowed.
A distinguished start position gives precisely a length-ell word. Hence the
trace counts words once. Determinism alone would not have established this;
the bounded-memory synchronization is necessary.

For `ell=p^s` dividing m, let `R_ell(t)=tr(M_t^ell)`. Primitive-period
necklace counts with their least length-m repetition at least t are

    N_1(t)=R_1(t),
    N_(p^s)(t)=[R_(p^s)(t)-R_(p^(s−1))(t)]/p^s,  s≥1.        (7)

The same forbidden-pattern property applies to the infinite repetition
independently of its chosen period. Thus subtraction removes precisely the
smaller-period words, and every remaining rotation orbit has size ell.
For threshold 0 all necklaces qualify; for threshold A^m none qualify.

## 5. Weighted canonical ranking and selection

At block i in state H, each exact-period necklace has continuation weight
`w_ell=T_(i+1)(H_ell)`. The number of complete canonical tuples whose current
block is strictly below a numeric threshold t is

    L_i(H,t)=sum_{ell|m} [P_A(ell)-N_ell(t)] w_ell.             (8)

This is an exact count of distinct tuples, because (4) gives every necklace
its true number of continuation orbits. Unlike arbitrary auxiliary-witness
multiplicity, these weights describe the target objects requested.

To unrank r, find the least block value x such that `L_i(H,x+1)>r` by binary
search over its `Q^L` possible values. Replace r by `r-L_i(H,x)`, output x,
compute its least packet period, and replace H by H_ell. Repeat for every
block. To rank an already canonical tuple, add `L_i(H,x)` at each block and
perform the same state updates. Equations (4) and (8) prove both algorithms
by induction. Values skipped by the search have no canonical necklace and
therefore contribute no jump. The first tuple and the last tuple are valid
boundary cases, not special assumptions.

The subgroup state alone does **not** process an unfinished prefix. The
threshold automaton does that work and can have m states. Closure to k+2
states resumes only when the whole block has been selected.

## 6. Complexity and scope

All matrices have dimension at most n≤D, and powers have exponent at most n.
Repeated squaring computes their traces with a polynomial number of integer
operations. Matrix entries and intermediate counts have polynomial bit length
in `D log Q`; the total orbit count is at most `Q^D`. There are polynomially
many block/state/period entries in (4), and binary search needs at most
`ceil(L log_2 Q)` threshold queries per block. The total is therefore polynomial
in D and `log Q`, using ordinary exact integer arithmetic. The implementation
uses uncomplicated dense matrix powers and bottom-up stabilizer tables, not an
optimized necklace-ranking library. Its rank operation requires a canonical
input tuple; it validates block lengths and color digits and checks each
block against its current stabilizer. It does not first canonicalize an
arbitrary input.

This is **not** polynomial in k: D is p^k−1, and writing one full bank state
already requires `D log Q` bits. It does not count prime bank orbits or full
tournament isomorphism classes. Enlarging G, permuting named coordinate cycles,
forgetting the fixed-origin gauge or imposing extra global predicates can
invalidate the small stabilizer family and its continuation weights.

The stopping point of this lane is an implemented, proof-backed rank/unrank
interface for the specified shared-residue action, with literal small-depth
controls. A full graph canonicalizer is not implied or attempted.

## 7. Independent finite controls

Run `python3 bank_rank.py` and `python3 -O bank_rank.py`; their outputs match.
All checks use explicit exceptions, not optimizable assertions.

- 73 exhaustive necklace thresholds cover alphabet/length pairs (3,1), (3,2),
  (2,4), (2,3) and (2,5). Literal rotation minima and exact periods independently
  audit the automaton trace and period subtraction.
- For every symbolic subgroup and block at the tested banks, literal pairs
  (a,b) audit both its shift image and each stabilizer transition.
- Period-transfer totals are checked against direct Burnside sums over every
  element of G.
- The p=3,k=2 permutation model has 8 coordinates, 256 labelled colorings and
  64 orbits. Every orbit, canonical rank and unrank is checked exhaustively.
- For p=5,k=2, the total is 692,224 over 24 coordinates and four subgroup states.
  Random ranks, both endpoints and random literal canonical forms are checked.
- For p=5,k=3, the total is
  170141183461707171770974935790600192 over 124 coordinates and five subgroup
  states. The same controls use every one of the 125 literal group actions.

These tests support the all-depth proof; they do not replace it. The all-depth
claim follows from the subgroup classification, stabilizer decomposition,
threshold trace lemma and ordinary exact-arithmetic cost bounds above.
