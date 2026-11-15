# Physically choosing a two-tape host at every macroscopic scale

Let U be the interval seven-circulant. Choose the strong tournament `C3[U,A,1]`,
where A is nonempty. Put Q=H(C3[U,A,1])/7. The block-run count gives
```
Q=259 P_A(2)-209 H(A),
```
where P_A(k) counts ordered covers by k paths with empties allowed. For example the
3- and 4- and 5-run color words contribute
`3h_U h_A+c2_U h_A+h_U c2_A+c2_U c2_A`;
substitute `h_U=175, c2_U=1638, c2_A=P_A(2)-2h_A`. This readout is inherited;
the two-tape asymptotics and inverse below are the new parts.

Take A on an m-vertex transitive spine plus apices a,b, using 0/1 rows a_i,b_i
for the arcs between the apices and spine (1=apex beats the site).
The arc between apices is arbitrary. Define
```
L_m(a,b) = sum_(i<j) (a_i b_j+b_i a_j) 2^(-(i-1)) (3/4)^(j-i-1).
R_m(a,b) = L_m(reverse(1-a), reverse(1-b)).
```
**Uniform two-end phase law**
```
Q / 4^m = (777/128) L_m R_m + O(m^4 (3/4)^m).
```
The leading product is independent of the arc between the apices.
Both endpoint functionals are bounded between 0 and 16.

## Proof by disjoint odd-cycle packings

For an odd-cycle packing, each uncovered vertex gives a factor k in P_A(k),
since it may be assigned to any slot; each cycle gives a factor 2k, since
it must all occupy one slot and has weight two in the odd-cycle expansion.
The spine itself contains no directed cycle, so there are at most two cycles.
The possibilities are (1) one cycle through a single apex, (2) two cycles through
the separate apices, (3) one cycle through both apices. In (3) the portions
between the two apices are two optional increasing spine runs. Formulas
below need neither an assumption about the mutual arc nor any probability model.

First count endpoint quadruples with two distinct start sites i<j followed by
two distinct end sites k<l, with j<k. The two starts attach to different apex
out-arcs, and the two ends attach to different apex in-arcs. For each such
choice there are two structural possibilities: pair by apex color to give
two separate cycles; pair crosswise to give one cycle through both apices.
All edge conditions are now paid. At a non-endpoint spine position the
remaining choices in either pairing are two uncovered slots, plus one
insertion into each of the assigned increasing runs whose endpoints enclose it.
The cover multiplicities (0,1,2 active intervals) are identical for the two
ways to pair the starts and ends. Write
```
n0=(i-1)+(m-l), n1=(j-i-1)+(l-k-1), n2=k-j-1.
```
The unsigned assignment mass is M=2^n0 3^n1 4^n2 for each pairing.
For the separate cycles, both internal-run cardinalities must be even.
The parity projection has unsigned coefficient 1/4 and three signed products.
In each signed product, every non-endpoint site has weight at most three:
on double coverage the factors are 2+1-1, 2-1+1, 2-1-1.
Thus this pairing contributes `16*(M/4+O(3^m))=4M+O(3^m)` to P2, because
two cycles have weight 4 each. For the single cycle through both apices,
there are four endpoint vertices plus the two apices, so the sum of the two
internal sizes must be odd. Its projection has unsigned coefficient 1/2,
and the other term is bounded by 3^m since both run signs are negative.
It contributes `4*(M/2+O(3^m))=2M+O(3^m)`.
Thus the total leading coefficient for each eligible quadruple is 6M.

The other possibilities have O(m^4) endpoint specifications and O(3^m)
contribution each. In the one-apex cycle or a mixed cycle with a short run
(length zero or one), at most one spine run has free internal sites. With
distinct endpoints in the wrong order, the two nonempty intervals can be
valid only if disjoint, again leaving at most one free run at a site.
This also bounds the empty packing and its contribution with only one apex.
All factors involving apices and reserved endpoints are absolute constants.

Divide 6M by 4^m. Since four sites are reserved, the answer factors as
`(6/4^4)·2^(-(i-1))·(3/4)^(j-i-1)·2^(-(m-l))·(3/4)^(l-k-1)`.
These are exactly the terms of `(6/4^4)L_m R_m` except that that product
also includes combinations with j>=k. Their outer and fringe exponents
sum to m+j-k-3>=m-3, so the total fake mass is O(m^4(3/4)^m).
Finally H(A)=O(3^m): any packing allocates each spine site to one of the
at most two increasing runs or to neither, and there are a bounded number
of apex-use patterns. This proves the formula for Q.

## Target-first placement at every large scale

On the infinite left row b=111..., the weight of a_j in L is exactly
`4·(3/4)^(j-1)`. Its use as the earlier of the two endpoints has weight
`4·2^(-(j-1))`. Its use as the later one gives the finite geometric sum
`4·(3/4)^(j-1)-4·2^(-(j-1))`. Hence
```
L(a,111...) = 4 sum_(j>=1) a_j (3/4)^(j-1).
```
This fills [0,16], since the two first-digit rest intervals overlap
(equivalently 3/4>=1/2), giving a greedy decoder in order of positions.
The all-zero right end makes R tend to 16. Truncating an endpoint functional
beyond position d changes it by O((3/4)^d): sum its geometric weights with
second endpoint j>d. Therefore physically splice the greedy a digits and
the b row of ones over the left floor(m/2) positions, with both rows zero
on the right half. For each real x in [0,1554], this constructs a host with
```
abs(Q - x·4^m) = O((2 sqrt(3))^m).
```
The constants are uniform, including x at the endpoints (the error bound
does not mean Q=0 occurs). Choosing a→b makes A itself strong for m>=2:
b beats the first spine point; both apices lose to the last, and a reaches b.
Choose m as a function of a prescribed large positive integer q so that
q/4^m lies in a fixed positive interval [c,16c] contained in [0,1554],
taking m even. This yields relative error O((3/4)^(m/2)) and absolute
error O(q^log_4(2sqrt(3))). This is not an exact integer selector or a
smaller error than the unrelated one-tape row reductions.

The even-order choice matters for prime cofactors. The slot polynomial
has parity P_A(-z)=(-1)^|A| P_A(z), since every odd cycle removes an even
number of slot choices versus all singles. Reducing z=2 modulo three
gives P_A(2)≡(-1)^|A|H(A); thus Q is divisible by three whenever |A|
is odd. For even m, Q≡-H(A) mod 3. A macroscopic window does not force
a given residue or solve this remaining joint constraint.

## Location of the critical entropy

For a spine with at most one apex, the first 1 and last 0 bound all
cyclicity. Deleting exterior singleton padding keeps H and multiplies P2
by powers of two. The extreme pair contributes 2(3^(L-2)+1) to P2 for
core length L, whereas there are at most 2^(L-2) core words. The number
of possible q through X is O(X^log_3(2)). This is the specialization of
the repository's earlier 3-state two-slot bound; it is not presented
as a fresh general capacity theorem.

For two apices, for each mark row take the closed trimmed interval
[third 1, third-last 0], treating missing or crossed positions as empty.
Outside each nonempty interval there are at most two exceptions to the
zero prefix and two to the one suffix. A row with empty interval has
a separator with at most two 1's before and two 0's after (if either
bit appears at most twice, take an end separator). A singleton interval
cannot occur: its left endpoint is a 1 and right is a 0.

Let t,a,b count spine positions covered by zero, one, two trimmed intervals.
For each nonempty interval, pick a start at or before its left endpoint
and an end at or after its right endpoint. There are three candidates
of each kind; after the first row is handled, the second can avoid both
reserved endpoints. Thus the chosen odd cycles may use disjoint runs with
the much larger actual enlarged intervals. Ignoring parity, for each
nonreserved position there are at least 2,3,4 choices according to its
trimmed cover number. Projection to even counts in each run costs at
most four: all other Fourier products are nonnegative (local weights
2±1±1). The loss from reserving at most four endpoints is a bounded
factor, as are the unused-apex factors. Thus P2≫2^t 3^a 4^b.

For fixed t,a,b the intervals have only poly(a+t+1) possible shapes:
their four endpoint events have finitely many types, intervals with
double overlap have that central overlap as one segment of fixed
length b, and all other run lengths sum to a+t. Exceptions and separators
are again encodable with poly(a+t+1) choices, independently of b:
every exception lies outside its own trimmed interval, hence outside
the double overlap. If a row is empty then b=0 and all spine positions
already count in a+t. Given the fixed outside data the row bits on
one- and two-covered positions have at most 2^a 4^b possibilities.

Because P2≥2H, Q≥(309/2)P2. Sum the presentation bound first over b
with 4^b≪X/(2^t3^a); the geometric sum is at most a constant times
that upper cutoff. The ending bound is
```
O(X) sum_(a,t>=0) poly(a+t+1) (2/3)^a 2^(-t).
```
It converges, and restricting a+t>=D gains exp(-cD) for an absolute
c>0. Thus in this family a positive-density residue or exact-selector
program can focus on two long overlapping tapes with bounded fringes;
merely forcing a big gap between the two free intervals cannot help.
The upper bound does not assert positive density for any subfamily.

## Relation to twisted collision and what remains

At macroscopic precision the leading two-end phase is an overlapping
digit system of ratio 3/4 with a product join. The small-error term
contains overlapping run boundaries and the Hamiltonian value itself.
Equality at an integer target would require a new lattice-scale
argument for this remainder, not merely tighter greedy digit choices
or a general continuum limit. In particular these statements neither
prove the earlier twisted-collision lemma nor imply that positive
density of q would prove Chan–Pak.

## A fixed-width repair matches the exact residue modulo three

The mod-three issue can be repaired jointly with the macroscopic placement,
not just solved in a separate family. Read a spine left-to-right. For each
increasing run keep one of five states: unstarted, open odd, open even,
closed odd, closed even. At each spine site either leave it uncovered
(weight z), or put it in exactly one run (weight one). To start a run
require the outgoing arc from its apex; one may close the run at its
assigned current vertex only if that vertex beats the destination apex.
An open run may skip intermediate sites. An unstarted run left unused is
accepted only in the mixed-apex cycle with the corresponding direct
apex arc. Each single-apex cycle needs a nonempty even run, the two
separate cycles each need a nonempty even run, and the mixed-apex cycle
needs odd total run size.

The four possible run sets (a-loop, b-loop, both loops, and both cross
runs) thus give vectors of sizes 5,5,25,25. At z=1, combine their final
observations with factors 2,2,4,2 and add 1 for the empty packing, giving
H(A) exactly. At z=2 the factors and uncovered-site weight become
4·2, 4·2, 4²,4 and 2, respectively. (The extra factor 2 in a single loop
is the uncovered other apex.)

Here is a finite all-length mod-three certificate (see verify_host_residue.py).
Use z=1 and coefficient vectors over F3. The only left prefix symbols
are 01,11 since b is identically 1 there. Starting with four unstarted
run sets, their reachable-state closure has size 73: successive new
layers have sizes 1,2,4,6,10,12,12,10,8,6,2,0. The pair of acceptance
functionals for the two choices of mutual arc, pulled backwards through
a positive all-zero right suffix, has period two: length 3 agrees with
length 1, and lengths 0,1,2 are distinct. Across every one of the 73
states and both suffix parities, there is a choice of three more left
digits a and one orientation of the mutual arc to meet each residue
of -H(A) mod 3. For an even right suffix and the fixed arc a→b, four
repair bits suffice for every state. All statements are checked on
literal finite closure and observation vectors, so they hold for
prefixes and suffixes of any length, not just hosts of bounded size.

For a strong physical template in this argument, take m
divisible by 4. Use the greedy phase encoder on only the first m/2−4
left spots, and fill the final four spots of the left half with the
repair word. The right half is all zero in both rows. The right suffix
length is even and |A| is even, so the cofactor Q≡-H(A) mod 3
holds. Choose among the repair words for the residue of the prescribed
target. This keeps the same asymptotic error order, as only a constant
number of greedy bits were set aside. At m divisible by 4 the scales
are in ratio 256, so a fixed scaled interval [4,1024] suffices to place
arbitrary large targets. The output cofactor now matches the target
exactly modulo three as well as falling in the macroscopic interval.
No claim is made for finer residue lattices or for equality.

Run `python3 replay_host.py` to compare the slot DP with independently
computed subset path counts at 802 hosts. Run `python3 verify_host_residue.py`
to replay the finite all-order certificate, or `python3 replay_host.py N`
to generate the two rows for a prescribed cofactor N. Every count
returned for these rows is an exact integer, computed with a number
of DP operations linear in their length. The error to N itself is
asymptotic and generally nonzero.

### Modulo-nine strengthening

The same fixed-arc physical template can hit every residue modulo nine
with eight repair bits instead of four. This time use both weights z=1
and z=2 in the run machine and reduce their direct sum modulo 9,
giving 120 coordinates. The final observation has the 259 and -209
weights built in; the empty packing omitted from this linear part is
the same for all repair words of the given length. The complete left
01/11 reachable closure has 209,955 states. Its coordinate for the
unstarted a-loop at z=2 is exactly 2^l modulo 9 at left length l,
detecting the parity of l. There are 104,977 even-prefix states.
The fixed-arc accepting row pulled back through zero symbols has
preperiod 2 and period 6, so all positive even suffixes use one
of the observations at lengths 2,4,6. For each observation initialize
each prefix state with the singleton set of its observed residues.
An eight-step backwards union along the two left-symbol transitions
gives the entire nine-residue set at each even-prefix state. These are
exact finite graph assertions in `verify_host_residue9.cpp`, compilable
using C++17 with assertions enabled.

Take m divisible by 4 with m>=16, so both the right suffix and the
left prefix before the eight repair bits are even. Reserve those eight
bits and choose the nearest of the completions in the target residue.
The splice keeps the same big-O absolute error, only changing its
constant. The command `python3 replay_host.py N` uses this stronger
mod-nine selector by default. This is still not a lattice-scale
approximation: the width grows unboundedly with the target.
