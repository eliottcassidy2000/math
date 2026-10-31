# Session 125 — the 1, 3, 8 bridge, with its proof boundary

## What is established

The triple `1,3,8` has two compatible classical meanings: it begins the
even-indexed Fibonacci sequence and is a `D(1)` Diophantine triple. The former
supports an exact positional arithmetic, the **Chung–Graham representation**.
This is not a newly discovered numeral system. Its existence, uniqueness and
connection to Zeckendorf arithmetic are established literature; see
[Burns, *Chung-Graham and Zeckendorf representations* (2025), Theorem 2 and
Section 2](https://arxiv.org/html/2502.18870v1).

Below are self-contained derivations of the particular identities and
algorithms needed for this project. No claim here says that a numeral, a
Diophantine triple, or an automaton constructs a Hamiltonian path/cycle in a
square-sum graph.

## 1. Uniform notation and the Diophantine triple

Use the conventional Fibonacci indexing

\[
F_0=0,\quad F_1=1,\quad F_{r+2}=F_{r+1}+F_r.
\]

Set

\[
A_{-1}=0,\qquad A_i=F_{2i+2}\quad(i\geq0).
\]

Thus `A_0,A_1,A_2,... = 1,3,8,21,55,144,...`, and

\[
A_{i+1}=3A_i-A_{i-1}\quad(i\geq0).
\tag{1}
\]

Indeed, `F_(r+4)=3F_(r+2)-F_r`, obtained by expanding the Fibonacci
recurrence twice, gives (1).

For successive terms `a=A_i`, `b=A_(i+1)`, `c=A_(i+2)`, the invariant

\[
a^2-3ab+b^2=1
\tag{2}
\]

holds: it is `1-9+9=1` initially and is unchanged by
`(a,b) -> (b,3b-a)`. Consequently

\[
ab+1=(b-a)^2,\qquad ac+1=b^2,\qquad bc+1=(c-b)^2.
\tag{3}
\]

In particular,

\[
1\cdot3+1=2^2,\quad1\cdot8+1=3^2,\quad3\cdot8+1=5^2.
\]

The positive regular extension of a Diophantine pair satisfies

\[
c=a+b+2\sqrt{ab+1}=a+b+2(b-a)=3b-a.
\tag{4}
\]

Equivalently, every consecutive triple
`(F_(2k),F_(2k+2),F_(2k+4))`, `k≥1`, is a `D(1)` triple.
This is classical, not a novelty claim. An author-hosted mathematical source
places these triples in the Hoggatt–Bergum quadruple family:
[Dujella, *Diophantine quadruples and Fibonacci numbers*, Section 2](https://web.math.pmf.unizg.hr/~duje/pdf/pondi.pdf).

### Exact overlap with square-sum edges, and its limit

At the unit vertex only, `1*b+1 = 1+b`. Therefore the two `D(1)` relations
involving `1` also certify the square-sum edges `1—3` and `1—8`.
But `3+8=11` is not a square; `3*8+1=25` uses a different operation.
So `1,3,8` is **not** a triangle of square-sum edges.

Nor are all terms of the chain pairwise `D(1)` compatible: `1*21+1=22`
already fails. The claim is about every three consecutive terms, not the
whole infinite set.

## 2. The canonical even-Fibonacci numeral system

A finite word `d_(k-1)...d_0`, written most-significant digit first, has value

\[
V(d)=\sum_{i=0}^{k-1}d_i A_i.
\tag{5}
\]

Its canonical constraints are:

1. Each digit belongs to `{0,1,2}`.
2. Between any two digits `2`, at least one digit is `0`.

Equivalently the forbidden factors are `2 1^r 2`, for all `r≥0`.
Initial zeros are harmless; the empty word has value zero. For display the
implementation uses the one-character word `0` for zero and removes other
initial zeros.

The two live legality states are:

| Current state | Read 0 | Read 1 | Read 2 |
|---|---|---|---|
| Free | Free | Free | Locked |
| Locked | Free | Locked | Reject |

Both live states accept at end of input. Start in Free. The same legality
condition works in either reading direction because every forbidden factor
is its own reversal.

### Interval theorem: existence, uniqueness, order and exact counts

Put

\[
B_k=A_k-A_{k-1}=F_{2k+1},\qquad k\geq0.
\]

For every `k≥0`:

* Legal length-`k` words starting Free represent each integer in
  `[0,A_k-1]` exactly once.
* Legal length-`k` words starting Locked represent each integer in
  `[0,B_k-1]` exactly once.
* In each list, lexicographic order is numerical order.

**Proof.** At length zero each list is the single empty word, of value zero,
and `A_0=B_0=1`. Suppose the statement holds at `k-1`, and write
`a=A_(k-1)`, `b=B_(k-1)`. A leading digit in a length-`k` word has weight `a`.
The Free possibilities are disjoint successive intervals:

\[
0+\text{Free}: [0,a-1],\qquad
a+\text{Free}: [a,2a-1],\qquad
2a+\text{Locked}: [2a,2a+b-1].
\]

The Locked possibilities are

\[
0+\text{Free}: [0,a-1],\qquad
a+\text{Locked}: [a,a+b-1].
\]

Since `A_k=2a+b` and `B_k=a+b`, these are exactly the required intervals.
Inductive uniqueness inside each block and disjointness between blocks
prove uniqueness; their order proves the lexicographic assertion. ∎

Consequences:

* Every nonnegative integer has one canonical word, modulo leading zeros.
* There are exactly `A_k=F_(2k+2)` legal length-`k` words.
* For `k≥1`, the largest length-`k` word is `2 1^(k-1)`, of value `A_k-1`.
* The largest word starting Locked is `1^k`, of value `B_k-1`.
* The numerical value is exactly the zero-based lexicographic rank among
  legal fixed-length words.

### Greedy construction is genuinely valid

Choose the largest weight not exceeding the integer, divide by that weight,
then repeat on the remainder and lower weights. No separate search is needed.
In the Free interval, the leading quotient is 0, 1 or 2, and the three
remainder ranges are exactly those in the proof. In the Locked interval it
is 0 or 1, and the two remainder ranges are exactly those above. Thus greedy
division follows the legality automaton even without explicitly storing its
state.

### Why Fibonacci appears in the automaton

The numbers of completions satisfy

\[
\begin{pmatrix}A_k\\ B_k\end{pmatrix}
=\begin{pmatrix}2&1\\1&1\end{pmatrix}
 \begin{pmatrix}A_{k-1}\\ B_{k-1}\end{pmatrix},
\qquad
\begin{pmatrix}2&1\\1&1\end{pmatrix}
=\begin{pmatrix}1&1\\1&0\end{pmatrix}^{\!2}.
\tag{6}
\]

This is an exact two-step Fibonacci counting relation. It does **not** mean
that taking consecutive pairs of ordinary Zeckendorf digits and reinterpreting
them as ternary digits preserves numerical value.

## 3. Positive carry algebra and termination

For nonnegative coefficient arrays, there is a direct local normalization,
not merely evaluation followed by greedy re-encoding. The value-preserving
rules are

\[
3A_0=A_1,
\tag{7}
\]

\[
3A_i=A_{i+1}+A_{i-1}\qquad(i\geq1),
\tag{8}
\]

and

\[
2A_i+\sum_{j=l+1}^{i-1}A_j+2A_l
=A_{i+1}+A_{l-1}\qquad(i>l\geq0).
\tag{9}
\]

The `A_-1=0` convention handles the least-significant boundary in (9).
Rules (7)–(8) are the recurrence. To prove (9), note that

\[
\sum_{j=l+1}^{i-1}A_j
= (A_i-A_{i-1})-(A_{l+1}-A_l).
\]

Substitution gives `3A_i-A_(i-1) + 3A_l-A_(l+1)`, which is the right-hand
side of (9). This also works when the interior sum is empty.

**Algorithm.** Whenever a coefficient is at least 3, subtract 3 there and
apply (7) or (8). Otherwise, if a forbidden block `2 1* 2` occurs, apply
(9), removing that block and incrementing its immediately higher place and,
unless `l=0`, its immediately lower place. Repeat.

**Termination proof.** Every coefficient stays nonnegative. The integer
potential `Σ d_i` strictly decreases at every rewrite: (7) decreases it by 2,
(8) by 1, and (9), with `r=i-l-1` interior ones, decreases it by `r+3` at
`l=0`, or by `r+2` otherwise. Consequently there are at most the initial
digit sum many rewrites. When none applies, all digits lie in `{0,1,2}` and
no forbidden block remains. The interval theorem then makes the result the
unique canonical word. This proves strategy-independent output for any
choice of available rewrites; the implementation chooses a deterministic
order.

**Minimum-summand corollary.** Among all finite nonnegative coefficient
representations of an integer using the weights `1,3,8,21,...`, the canonical
Chung–Graham representation has the smallest possible coefficient sum and
is the unique minimizer, modulo leading zeros. Indeed, any noncanonical
coefficient array admits at least one rewrite, and its normalization has
strictly smaller coefficient sum. This is a consequence of the proof here,
not an assertion of historical novelty. "Summand count" means multiplicity
is counted: a coefficient 2 contributes two summands.

Examples of noncanonical to canonical conversion:

* `22 -> 100`, because `2*3+2=8`.
* `212 -> 1000`, because `2*8+3+2=21`.
* `2120 -> 10001`, because `2*21+8+2*3=55+1`.

These moves are arithmetic identities. They do not prescribe a legal
square-sum graph extension or a local modification of a Hamiltonian object.

## 4. Successor and parity are exact

### Successor without decoding

Prepend a zero. For each position store whether its prefix leaves the
automaton Free or Locked. In a Free position the largest permitted digit is
2; in a Locked position it is 1. Find the rightmost position below its
permitted maximum, increment that digit, and replace every following digit
by zero. Drop initial zeros.

The resulting word is legal. Every lexicographically greater legal word
must first increase a digit somewhere, and the rightmost available increase
followed by the least possible suffix is the immediate lexicographic
successor. By the interval theorem this increases the integer by exactly 1.
The prepended zero guarantees that an incrementable position exists.

For example `2 -> 10`, `21 -> 100`, `121 -> 200`, and `211 -> 1000`.

### Parity is not the final digit

The weights have period-three parities

\[
A_i\bmod2=1,1,0,1,1,0,\ldots.
\tag{10}
\]

This follows from (1) modulo 2 and the initial pair `(1,1)`.
Thus for LSD-indexed digits,

\[
V(d)\equiv\sum_{i\not\equiv2\pmod3}d_i\pmod2.
\tag{11}
\]

For an MSB streaming formula, let `N` be the value of a prefix, and `P` the
value obtained by shifting all its weights upward one place. Appending digit
`d` gives the exact integer update

\[
(N,P)\longmapsto(P+d,\;3P-N+3d).
\tag{12}
\]

Start at `(0,0)`. This follows by shifting (5) and using
`A_(i+2)=3A_(i+1)-A_i`, with `A_1=3`. Modulo 2 it is

\[
(N,P)\longmapsto(P+d,\;N+P+d)\pmod2.
\tag{13}
\]

The product of these four register states and the two legality states has
eight live states; all are reachable. Witnesses `(locked,N,P): word` are

```
(False,0,0): empty    (False,1,1): 1
(True, 0,0): 2        (False,1,0): 10
(False,0,1): 11       (True, 1,0): 12
(True, 1,1): 21       (True, 0,1): 102
```

A total recognizer also needs its rejecting sink. No minimality claim is
needed. The parity output is the first register.

In particular `10` represents the odd integer 3, whereas `11` represents
the even integer 4. Odd/even status is not just the last digit. The two
successive words for any consecutive integer pair always have opposite
parity; that fact does not by itself specify which member of a geometric
pair is the odd one.

## 5. Actual relation to Zeckendorf, not just a resemblance

Zeckendorf uses all weights `F_2,F_3,F_4,...=1,2,3,5,...`, binary digits, and
no consecutive 1s. Chung–Graham uses only the even-indexed weights and the
different digit constraint proved above. They are distinct canonical
representations of the same integers.

A concrete conversion bridge is to split each Chung–Graham digit into
`d=min(d,1)+max(d-1,0)`. Insert zero coefficients at all odd Fibonacci
indices. Each resulting binary array is a valid Zeckendorf word, because
its nonzero coefficients occur only at even indices. Their numerical sum is
the original integer; ordinary Zeckendorf addition yields its canonical
Zeckendorf representation. This splitting construction is explicitly used
by [Burns, Section 2](https://arxiv.org/html/2502.18870v1).

For instance `17=201_CG=8+8+1`; the split is `(8+1)+8`, while its single
Zeckendorf expansion is `13+3+1`. The split's two addends are canonical,
but their unnormalized sum is not automatically canonical.

Conversely one may evaluate a Zeckendorf word and apply the proved
Chung–Graham greedy algorithm. A direct algebraic explanation for moving
between the two weight sets is

\[
F_{2j+1}=F_{2j+2}-F_{2j}=A_j-A_{j-1}\quad(j\geq1).
\tag{14}
\]

Thus if `z_r` are Zeckendorf coefficients (`z_1=0`), their signed even-weight
coefficients are

\[
c_j=z_{2j+2}+z_{2j+1}-z_{2j+3}\quad(j\geq0),
\tag{15}
\]

with omitted coefficients taken as zero. Because adjacent Zeckendorf 1s are
forbidden, each `c_j` lies in `{-1,0,1}`. Formula (15) is not yet a canonical
positive Chung–Graham expansion: negative coefficients still need conversion.
The positive-only rewrite algorithm above does not pretend to handle signed
input.

The common connection to Fibonacci/Ostrowski ideas is real, but the present
weights are not being asserted to be a standard continued-fraction
denominator sequence: their recurrence has the minus sign in (1). The
precise established name is Chung–Graham, and it is unnecessary to relabel it
as an Ostrowski system.

## 6. How to attach this safely to an origin at G_32

Let `G_m` retain its ordinary vertex labels `1,...,m` and ordinary square-sum
edge predicate. A coordinate choice can label graph **orders**, without
changing any vertex arithmetic:

* Zero-based: `G_m` has growth coordinate `t=m-32`; `G_32` is coordinate 0.
* One-based: `G_m` has natural-number coordinate `n=m-31`; `G_32` is number 1.

Either coordinate can be written in Chung–Graham numerals, and adding one
vertex corresponds numerically to the proved successor. In the first
convention order parity equals coordinate parity; in the second it is
opposite. The convention must be stated before discussing parity.

For a pair of orders `(m,m+1)`, the coordinate words are consecutive and
have opposite parity. This remains a necessary arithmetic observation, not
a sufficient certificate that both graphs admit any proposed Hamiltonian
objects. The hidden prefix `m<32` is not reconstructed or characterized by
the numeral system alone. That requires separate graph-theoretic evidence.

Examples in the unshifted integer coordinate:

| Integer | Chung–Graham | Zeckendorf |
|---:|---:|---:|
| 0 | 0 | 0 |
| 1 | 1 | 1 |
| 2 | 2 | 10 |
| 3 | 10 | 100 |
| 7 | 21 | 1010 |
| 8 | 100 | 10000 |
| 17 | 201 | 100101 |
| 20 | 211 | 101010 |
| 21 | 1000 | 1000000 |
| 31 | 1102 | 1010010 |
| 32 | 1110 | 1010100 |
| 33 | 1111 | 1010101 |
| 34 | 1112 | 10000000 |

## 7. Verification and provenance

`port_numeration.py` implements the independent arithmetic and normalization;
`test_port_numeration.py` uses explicit exceptions rather than `assert`, so
optimization cannot silently remove checks. Reproduce from this directory:

```sh
python test_port_numeration.py verification.normal.json
python -O test_port_numeration.py verification.optimized.json
```

Both runs passed 1,359,569 explicit checks. Coverage: all 88,573 ternary words
of lengths 0 through 10; 100,001 consecutive integer roundtrips and
successors; 1,000 `D(1)` triples; 10,841 local normalizations accounting for
118,646 rewrites; 1,000 random 512-bit cases; all eight live parity/legality
states. An independent unbounded-coin dynamic program checks the minimum
summand count and its unique minimizing coefficient array for all 100,001
integers from 0 through 100,000. These are finite regression checks supporting, not replacing, the
universal proofs above.

Primary research sources checked on 2026-10-31:

1. [Rob Burns, *Chung-Graham and Zeckendorf representations*, arXiv:2502.18870v1](https://arxiv.org/html/2502.18870v1).
   Directly inspected theorem and splitting/automata treatment.
2. [Sungkon Chang, *The Chung-Graham Expansion*, arXiv:2504.00314v1](https://arxiv.org/html/2504.00314v1).
   Independently confirms the historical expansion and generalizes equally
   spaced Fibonacci indices; no generalization claim is made here.
3. [Andrej Dujella, *Diophantine quadruples and Fibonacci numbers*](https://web.math.pmf.unizg.hr/~duje/pdf/pondi.pdf).
   Directly inspected identities and consecutive-triple discussion; no
   present-day status claim is inferred from this older survey.
4. Original Chung–Graham historical reference: *On irregularities of
   distribution*, in *Finite and Infinite Sets* (Eger 1981), Colloq. Math.
   Soc. János Bolyai 37, North-Holland, 1984, pp. 181–222.
   [Author-hosted scan](https://mathweb.ucsd.edu/~ronspubs/81_12_irregularities.pdf).
   The scan resolves as a 42-page PDF but the available browser returned no
   extracted text; theorem attribution above therefore rests on the
   directly inspected modern primary papers, not a claimed full reading of
   this scan. The conference year 1981 and publication year 1984 explain
   differing short dates in the literature.

No source worktree or remote repository was edited in this lane. No
publication result is implied by these scratch files.
