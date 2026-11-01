# Executable Hamiltonian counts, pointed inverse fibres and a square-shell law

**PROVED (identities and program semantics); FINITE-EXACT (explicit replay
ranges).** All graph counts are ordinary directed Hamiltonian-path counts,
not paths modulo reversal or `H/|Aut|`. The digital inverse is inherited from
Session124, and the weighted word formula from the October20 feedback-word
packet. Their actual finite FRACTRAN realizations and the fixed-H converse
refinement are the contribution here. No unrestricted spectrum-completeness
or shortest-program claim is made.

## 1. Exact carrier: pointed one-apex tournaments

A pointed one-apex tournament is `(T,x)` such that deleting the distinguished
vertex `x` leaves a transitive tournament. Its chain has its unique order
`v_0 -> v_1 -> ... -> v_(m-1)`; write `b_i=1` iff `x -> v_i`. Thus binary
words are in bijection with **pointed isomorphism classes**, not unpointed
classes and not Hamilton-path presentations. An isomorphism fixing `x` must
preserve the unique chain order, so different words are not identified.

Leading zeroes are source vertices and trailing ones are sinks. Their deletion
preserves H. Strip both. Either the empty word remains and H=1, or the reduced
word starts in 1 and ends in 0. Those words are exactly the nontrivial strong
pointed one-apex tournaments: the chain together with `x` gives a directed
spanning cycle, while a leading zero or trailing one prevents strongness.
Singletons are excluded from the word **nontrivial** below.

For `d>=1` put `w(1)=1`, `w(d)=2^(d-2)` for `d>=2`. Define

```
c(b) = sum_(i<j) b_i(1-b_j) w(j-i).
```

Every directed cycle contains `x`. Its minimum and maximum chain indices
`i<j` must have `b_i=1,b_j=0`, and its interior vertices form an increasing
subsequence. The cycle is odd exactly when it chooses an even number of the
`j-i-1` possible interior vertices, giving the weight `w(j-i)`. Therefore
`c(b)` counts ALL directed odd cycles, with different cyclic orders treated
as different cycles. The conflict graph is the clique on these `c(b)` cycles,
so its entire independence polynomial is `1+c(b)z`.

The Hamiltonian identity is

```
H(T(b)) = 1 + 2c(b).                                      (1)
```

Here is a proof independent of a black-box cycle identity. Removing `x` from
a path leaves two increasing blocks `A,x,B`. Empty `A` contributes `b_0`,
empty `B` contributes `1-b_(m-1)`. If every A-index precedes every B-index,
the cut is an adjacent `01`; otherwise put `i=min B<j=max A`. Its endpoints
must be an inversion `10`; the indices between them may be assigned freely
to either block, yielding `2^(j-i-1)` paths. Finally
`N_01-N_10=b_(m-1)-b_0`. Thus adjacent inversions receive their second copy
and every nonadjacent inversion contributes `2w(j-i)`, proving (1). The empty
word is the one-vertex tournament and has H=1.

### Faithful transfer, and its inherited digital normal form

Start with `(c,t,l)=(0,0,0)`. Appending a letter acts by

```
bit 0: (c,t,l) -> (c+t+l, 2t+l, 0),
bit 1: (c,t,l) -> (c,     2t+l, 1).                       (2)
```

Here `l` is the previous letter; `t` is the geometrically weighted sum of
all earlier ones. Equivalently, `Q=2t+l` is the ordinary binary value of the
word, with the first chronological letter most significant. Hence `Q'=2Q+b`.
The H-coordinate is precisely the inherited digital function

```
f(0)=f(1)=1,
f(2n+1)=f(n)                  (n>=0),
f(2n)=f(n)+n+(n mod 2)        (n>=1),
H(T(b))=f(Q).                                               (3)
```

This is the `h=1` specialization of Session124's actual-seed singleton
constructor. A new name for f would not constitute a new sequence theorem.
The sidecar Q matters: `f(4)=f(6)=5`, but `f(8)=9` and `f(12)=11`. Thus the
same scalar H can have different next-zero response.

## 2. A fixed 32-fraction machine for the whole word family

The actual ordered list is `one_apex.fractions.txt`, reproduced numerically
below. For a word b of length m, set

```
E(b)=2^m+sum_(i=0)^(m-1) b_i 2^i.
```

Unlike Q, the earliest letter is **least significant** in E. The leading 1
is a termination sentinel. Start the FRACTRAN program at

```
38 * 5^E = 2 * 19 * 5^E.
```

The final integer is exactly

```
107 * 2^H * 3^c * 11^floor(Q/2) * 13^(Q mod 2),             (4)
```

where H and c are from (1). Readout is by **prime exponent**, not the whole
integer. In particular the program retains the Q-sidecar, not only H.

```
161/475, 19/23, 29/95, 31/19,
185/203, 107/29, 215/217, 109/31,
205/259, 37/41, 53/37,
235/301, 43/47, 79/43,
17051/583, 53/59, 61/53,
1139/793, 61/67, 71/61,
803/1207, 71/73, 247/71,
287844/869, 79/83, 89/79,
19788/1157, 89/97, 101/89,
1133/1717, 101/103, 19/101.
```

The data primes are `H:2,c:3,N:5,Q_work:7,t:11,l:13,U:17`. Controls are
`19,23,...,103`, with successful halt107 and error109. The documented input
has only control19, H=1, N=E and all other registers zero.

**First-eligible proof.** Every numerator and denominator has exactly one
control prime, and these differ. Every reduced denominator has exactly its
source control; distinct-state fractions cannot fire. In each state the
nonzero/paired rule occurs before its zero fallback. Guard states make loop
returns a separate fraction, avoiding cancellation of a control factor.
Thus exactly one control remains throughout and literal first-eligible
semantics implements the following operations:

* At control19, consume N in pairs and put one unit per pair into register7.
  An odd residual selects control29; an even one selects31.
* If the quotient is zero, odd means the sentinel `N=1`, giving halt107.
  Even `N=0` gives error109 and is outside the permitted E>=1 inputs.
* Otherwise transfer the quotient back to N and enter the selected bit's
  transfer phases. Move each t-unit to two U-units; move l to U; transfer U
  back to t. On bit0 each old t- or l-unit also adds two H-units and one
  c-unit; bit1 sets l=1. This is exactly (2).

N becomes `floor(N/2)` after each processed bit; every transfer loop consumes
one nonnegative exponent. Therefore every permitted input terminates and
(4) follows by induction. All arithmetic conditions are proof obligations,
not inferred from finite replay. Swapping the paired rule with the odd
fallback makes the word10 halt with H=1 instead of3: fraction order matters.

**Exact runtime refinement.** Let s_2(Q) be its binary digit sum and take
`Q mod2=0` for the empty word. The number of actual fractions applied is

```
steps(b)=4E+3Q-6s_2(Q)-(Q mod2)+4m-2.                       (5)
```

A processed bit with current N contributes `4 floor(N/2)+1` parser steps
and `6t+4l+3` transfer steps; the sentinel costs two. Sum the identities
`sum_(j>=1) floor(E/2^j)=E-s_2(E)` and the corresponding prefix identity
for Q, using `s_2(E)=1+s_2(Q)`. Equation (5) also gives two steps on the empty
word. Unary exponent input can be exponentially expensive in word length;
the construction makes no fast binary-computation claim.

## 3. Exact bounded-H fibres, with no hidden order cutoff

For odd `h>=3`, let R_h be the number of **pointed isomorphism classes of
nontrivial strong one-apex tournaments** with exactly h Hamiltonian paths,
summed over every order. Put R_1=0 and R_h=0 for even h.

Every reduced word has a unique positive even address `Q=2n`, n>=1.
Since `f(n)>=1`, equation (3) gives `h>=n+1`; and `h<=2n+1` follows from
`f(n)<=n+1`, with `f(n)<=n` when n is odd. Consequently

```
R_h = #{n : (h-1)/2 <= n <= h-1,
             f(n)+n+(n mod2)=h},             h odd >=3.     (6)
```

The n>=1 condition is automatic here; omitting h>=3 would incorrectly
count the empty/transitive case at h=1. Equation (6) is the exact inherited
finite inverse, now interpreted as a class fibre rather than just a yes/no
constructor. Every solution gives the actual word `binary(2n)`.

To determine ALL R_h for h<=B, compute f(n) for 1<=n<=B-1 by its halving
recurrence and add each value `f(n)+n+(n mod2)<=B` to its target fibre. This
is an O(B) arithmetic-step algorithm, not polynomial in log B. It has no
chosen graph-order cutoff: `n<=B-1` is a proved necessary bound.

For h=3,5,...,31 the R_h sequence is

```
1,2,0,3,1,0,2,2,2,0,2,1,2,0,1.
```

The exact bounded replay at B=1,000,000 has 693,169 pointed strong classes,
249,076 distinct H-values, and maximum fibre size14 (first at H=176077).
These are restricted-family counts, not a new million-prefix unrestricted
coverage result. `FIBRES.json` keeps a digest of all coefficients and explicit
small witnesses; the algorithm regenerates them without external data.

## 4. Converse refines the inverse: exactly one odd fibre per square shell

Converse acts on words by

```
b -> complement(reverse(b)).                              (7)
```

This preserves reducedness and H. An isomorphism of pointed tournaments
must reverse the unique chain order when it identifies T with its converse.
Thus its fixed words are exactly

```
b = u · complement(reverse(u)),                           (8)
```

where u is a nonempty binary word beginning in1. Write q>=1 for the binary
value of u. There is exactly one u per q; no arbitrary leading zeroes are
allowed, since reduced b begins in1.

**Reflected-path formula.** The H-value of (8) is

```
G(q)=q^2+(q mod2)+2f(q)-1.                                 (9)
```

Indeed the two half-words contribute `2c(u)`. For cross inversions, a one
at position i of u meets a zero in the reflected half exactly when the
corresponding position j of u is one. Apart from the middle adjacent pair,
the geometric weights sum to q^2/2. The middle pair has weight1 instead of
1/2, giving the correction `(q mod2)/2`. Hence
`c(b)=2c(u)+(q^2+(q mod2))/2`, and (1) gives (9).

The first values, for q=1,2,..., are

```
3,9,11,25,31,45,51,81,91,117,127,165,179,213,227,
289,307,353,371,437,459,513,535,621,647.
```

**Square-shell theorem.** For every q>=1,

```
q^2 < G(q) <= (q+1)^2.                                   (10)
```

The lower bound uses f(q)>=1. If q is even, f(q)<=q+1 gives the upper
bound. If q is odd, `f(q)=f((q-1)/2)<=(q+1)/2`, yielding the strictly
stronger `G(q)<=q^2+q+1<(q+1)^2`. In particular
`G(q)<G(q+1)`, so no H-fibre has more than one pointed self-converse class.

Equality at the upper endpoint in (10) holds iff `q=2^a`, **a>=1**.
For even q=2n it is equivalent to f(q)=q+1, then to
`f(n)=n+1-(n mod2)`. An odd n forces n=1 by the sharper odd bound; an even
n recurses. Conversely every positive power2 with q>=2 attains equality.
The endpoint q=1 gives G(1)=3, not4.

Converse pairs all nonfixed elements of the finite fibre (6). Thus

```
R_h mod2 = 1  iff  h=G(q) for some q>=1.                   (11)
```

Together (10)--(11) say that **each consecutive square shell
(q^2,(q+1)^2] contains exactly one H-value whose pointed strong one-apex
fibre has odd size**. Equivalently, for each integer N>=1,

```
#{h <= N^2 : R_h is odd} = N-1.                           (12)
```

More generally, write J(X) for the number of fixed reduced words with
H<=X. For real X>=1 and k=floor(sqrt X), the exact formula is

```
J(X)=k-1+1_(G(k)<=X).                                     (12a)
```

All q<k lie below the cutoff, q>k lie above it, and only q=k is undecided.
The injectivity proved here makes J(X) also the number of odd R_h fibres
at or below X.

**Closest predecessor.** Section6 of the October20 feedback-word addendum
`BINARY_PRODUCT_DENSITY.md` already uses the same complement-reversal
involution and bounds this fixed-word count by J(X)=O(sqrt X). Its main
result is W(X)=X log2+O(sqrt X), where W counts all reduced words. Neither
the involution nor that density law is new here. Equations (9)--(12a)
refine its fixed-point sidecar to an exact all-cutoff formula and prove
fixed-H injectivity and parity; the million-bound total693169 is a finite
control consistent with the inherited ln2 asymptotic.

This is an all-order counting theorem derived from the actual converse
involution and the quadratic cross-inversion contribution. It is not a
numerological association between squares and primes.

### Scope hostile and inherited root-forgetting repair

Do NOT delete the apex mark in (11). At H=5 the two reduced words100 and110
are interchanged by converse and give R_5=2. Both are the same unpointed
strong four-tournament with row masks `(2,12,9,1)`. The anti-automorphism
`(0 1)(2 3)` exchanges its two valid apices, so the unpointed class is
self-converse although neither pointed class is fixed. Its unpointed H=5
fibre has size1, not2.

The inherited Session26 deletion-apex classification supplies an exact
repair, rather than guessing a quotient parity. Its needed special case
has the following short direct proof. In a reduced word of length m>=3,
deleting an internal chain vertex leaves the triangle consisting of x and
the first/last chain vertices. Deleting the first chain vertex leaves a
transitive tournament iff every remaining bit is0; deleting the last does
so iff every remaining bit is1. Thus a second apex exists exactly for the
two words `10^(m-1)` and `1^(m-1)0`. These represent the same unpointed
class `C3(1,1,n-2)`, n=m+1, with its long block transitive and
H=`2^(n-2)+1`. An anti-isomorphism exchanges its two singleton apices and
reverses the long chain. No automorphism can exchange the two apices, since
it would reverse their arc, and fixing one fixes the unique remaining
chain order. Every other class has its unique apex fixed by every
automorphism. Hence every strong one-apex class of order>=4 is rigid.
At m=2, C3 has three apices in one automorphism orbit and is counted once
both pointed and unpointed.

Let U_h and S_h count, respectively, unpointed nontrivial strong one-apex
classes and their self-converse classes over all orders. With
`e_h=1` iff `h=2^k+1` for some integer k>=2, the cited classification gives

```
U_h = R_h-e_h,
S_h = 1_(h in image G)+e_h.                               (13)
```

Thus the correct unpointed parity is the square-shell indicator **XOR** the
explicit double-apex exception. The first unpointed hostile is H=5; H=9
has two self-converse unpointed classes. The classification is inherited, and its needed special case is proved
above; no classification input is needed for the pointed theorem (10)--(12).

## 5. Actual square-shell sequence program

`square_shell_source.py` evaluates (9) with the exact inherited f recursion,
unrolled by successively halving q. `square_shell.fractions.00.txt` contains
its **1,371 fixed ordered fractions** (19 data registers,984 controls).
For every q>=1:

```
input: 71 * 5^q;   successful halt control:7937;
answer: exponent of3 = G(q); error exponent at2 =0.
```

Input q=0 is explicitly rejected with output0, not assigned a fictitious
strong class. The program uses no table of G-values. Its source terminates
because each q-loop strictly decreases q>1; multiplication and divisions
are expanded into finite counter operations. Correctness is (3),(9).

## 6. A universal graph-input H(T) program (separate, larger carrier)

For a labeled tournament on vertices0,...,v-1 encode its orientation by

```
b=sum_(0<=i<j<v) 1_(i->j) 2^(j(j-1)/2+i),
N=2^v(2b+1)-1.
```

Every permitted graph has exactly one code N; extract v as `v_2(N+1)` and
b from the remaining odd factor. Reject b>=`2^binom(v,2)` with output0.
The empty tournament has code0 and H=1. The one-vertex tournament has code1.

`hamilton_source.py` enumerates all v^v base-v words of length v, least
significant digit first, and keeps exactly those with distinct digits and
all forward path arcs. Its seen mask only sets a previously absent bit;
repetitions set an irreversible invalid flag. Every permutation appears
once, and every accepted permutation is exactly one Hamiltonian path.
For v=0 the sole empty word contributes one; no division by zero is entered.
For v>=1 every loop has a fixed finite bound. This is an explicit counting
algorithm with lossless graph input, not an invocation of universality.

The emitted list has **4,581 fixed ordered fractions**,31 registers and
3,279 controls. For every permitted code N,

```
input: 131 * 5^N; successful halt control:30689;
answer: exponent of3 = H(T); error exponent at2 =0.
```

Work registers may remain at halt. Invalid codes return0; legal tournaments
never do. Counter input is unary as a prime exponent, so even small
adjacency words can require large numbers of elementary steps.

### Common compiled first-eligible proof

`compiler.py` is a byte-preserving copy of the independently audited
Session132 compiler. It emits only INC, JUMP and DECJZ plus HALT. If p is
a data prime and q_i the current control, these become respectively
`p q_j/q_i`, `q_j/q_i`, and the ordered pair
`q_j/(p q_i), q_k/q_i`. No source and target control coincide; generated
loop guards prevent cancellation. Each fraction is reduced. On one-hot
inputs only the current state's fractions can apply; the nonzero DECJZ
fraction precedes the zero fallback. Induction proves literal fraction
steps coincide exactly with primitive counter steps.

The arithmetic compiler clears reused scratch registers, preserves named
operands when required, and implements truncating subtraction, repeated-add
multiplication and repeated-subtraction division. Every arithmetic loop
has a decreasing natural counter. In these two sources all intended
subtractions are nonnegative, except `v-1` at v=0 inside `v(v-1)`, where
both ordinary and truncating evaluation give zero. The source algorithms
above therefore prove all-input termination and output. Emitted-list
replays supplement, rather than replace, this argument.

## 7. What this contributes to the {7,21} question

For this word grammar, normalizing the word makes the first1 and last0
alone contribute `2^(m-3)` odd cycles for m>=3. If H<=21 then c<=10,
so m<=6. The 31 reduced words of lengths2..6 give the exact spectra

```
m=2: 3
m=3: 5
m=4: 9,11
m=5: 15,17,19,23
m=6: 25,27,31,33,37,39,45,47,51.
```

Thus7 and21 are absent at **every input length** in this grammar. So is13,
but an unrestricted five-tournament with row masks `(22,12,24,17,2)` has
H=13. Closing this one-apex machine's range therefore cannot prove full
spectrum completeness. The missing coordinate is not an extra scalar
counter: it is the availability of cycles avoiding the marked apex and the
actual cross-cycle incidence when feedback number grows.

Even unrestricted odd-cycle independence-polynomial arithmetic alone cannot
exclude7 or21: abstract cliques K3 and K10 have `I(G,2)=7,21`. A proof must
retain realizability as the **entire** directed-odd-cycle conflict graph of
a tournament. Our machine retains that condition for one-apex inputs only.

The all-order unrestricted exclusions are already inherited: THM-1370
reduces to strong components, uses nondecreasing strong floors3,5,9,15,25
and the finite order<=6 gap, and credits the pre-existing Moon/Busch route.
The new fraction lists do not replace those structural ingredients. Nor do
they prove the still-open claim that no further unrestricted H-values are
missing. The separate `tc=H/|Aut|` exclusion is a different theorem.

## 8. Other explicit sequence specializations

The same32-fraction list evaluates these actual tournament families:

* Alternating word `(10)^r`, r>=0:
  `H=(4^(r+1)+6r+5)/9 = 1,3,9,31,117,459,1825,...`.
  Its prefix state obeys `t_r=4t_(r-1)+1` and
  `c_r=c_(r-1)+2t_(r-1)+1`, proving the formula. It satisfies
  `H_(r+3)=6H_(r+2)-9H_(r+1)+4H_r`.
  Exact FRACTRAN time is `(22*4^r-16)/3+2r`, including r=0.
* Balanced word `1^r0^r`, r>=1:
  `H=(2^r-1)^2+2 = 3,11,51,227,963,3971,...`.
  This inherited run formula is not extended to r=0 (whose H is1, not2).
* Reflected words (8) give G(q) and the square-shell fibre-parity theorem.
* Equations (6),(13) give pointed and unpointed fixed-H multiplicity
  sequences, rather than ordinary all-tournament A/S/M counts.

## 9. Reproduction and inherited sources

From any directory, with Python3 standard library:

```
python3 /path/to/path_counter/check.py
python3 -O /path/to/path_counter/check.py
python3 /path/to/path_counter/fibres.py --bound 1000000
python3 /path/to/path_counter/build.py
```

`check.py` uses explicit exceptions, compares source-to-IR and every emitted
fraction byte, tests literal integer and exact exponent runs, direct graph
Hamiltonian DP and actual odd-cycle lists, malformed inputs, priority and
carry mutations, the H5 quotient hostile, the H13 grammar boundary, exact
runtime, inverse witnesses and square-shell endpoints. The receipt states
every finite range. It does not imply million-target FRACTRAN execution.

Inherited source pins (repository commit335577ae6d58):

* `05-knowledge/results/feedback_word_hamilton_spectrum_20261020/proofs.md`,
  Sections1--4: pointed words, weighted inversions, transfers and restricted
  H=13 obstruction. These are recovered, not re-claimed as new here.
* `05-knowledge/results/feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md`,
  especially Section6: the closest predecessor uses the same involution and
  proves J(X)=O(sqrt X); its W(X)=X log2+O(sqrt X) law is inherited. Our exact
  J(X), fixed-H injectivity and square-shell parity refine that sidecar.
* `05-knowledge/results/session124_digital_constructor_inverse_20261031/PROOF.md`,
  Sections1--4: f, actual-seed singleton inverse, q ancestry and fibre sidecar.
* `05-knowledge/results/session124_binary_bank_frontier_20261031/FINITE_BANK_PROOF.md`:
  the independent finite binary-extension obstruction; it does not assert
  this one-apex word language is unrestricted.
* `05-knowledge/results/session26_random_bridges_20261022/session26-blue/RESULT.md`,
  Section2, and its referenced apex-classification note: root-forgetting
  correction used only for (13).
* `01-canon/theorems/THM-002-ocf.md` and
  `THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md`:
  the all-tournament odd-cycle formula and explicit H9 polynomial collision.
* `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`:
  existing unrestricted omissions and priority correction.
* `05-knowledge/results/session132_fractran_tournament_recursion_20261031/fractran/`:
  inherited compiler and its audited first-eligible control semantics.

The independent Session133 audit is a separate sibling packet; producer
normal/-O equality is not described as independent verification.
