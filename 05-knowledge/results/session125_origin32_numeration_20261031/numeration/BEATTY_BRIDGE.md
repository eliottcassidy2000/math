# Exact Beatty division for the Chung–Graham numerals

This supplement leaves the previously verified files unchanged. It concerns
canonical arithmetic on `A_i=F_(2i+2)=1,3,8,21,...`, not a theorem that
square-sum graph growth is mechanical/Sturmian or Hamiltonian.

## 1. Main statement

Write

\[
\phi=\frac{1+\sqrt5}{2},\qquad
\beta=\phi^2=\frac{3+\sqrt5}{2},\qquad
\tau=\beta-2=\phi^{-1}=\frac1{\beta-1}.
\]

If the canonical Chung–Graham word for `n` is `w`, appending a zero has value

\[
P(n):=V(w0)=\lceil\beta n\rceil\qquad(n\geq0).
\tag{B1}
\]

The terminal legality state of `w` is Free exactly when
`P(n+1)-P(n)=3`; it is Locked exactly when that gap is 2. Hence appending a
digit `d` is the exact operation

\[
n\longmapsto P(n)+d,\qquad 0\leq d<P(n+1)-P(n).
\tag{B2}
\]

Every nonnegative integer `N` has the unique division

\[
p=\left\lfloor\frac N\beta\right\rfloor,
\qquad d=N-\lceil\beta p\rceil,
\qquad 0\leq d<P(p+1)-P(p).
\tag{B3}
\]

Thus the quotient drops the last canonical digit, and the remainder is that
digit. This is exact division based on an irrational slope, not an
approximation based on the limiting ratio of consecutive weights or an
assertion that these are ordinary power-of-beta place values.

## 2. Strict fractional bounds prove the shift formula

The even-Fibonacci identity is

\[
A_{i+1}-\beta A_i=\beta^{-i-1}\qquad(i\geq0).
\tag{B4}
\]

It follows immediately from
`A_i=(beta^(i+1)-beta^(-i-1))/sqrt(5)` and
`beta-beta^(-1)=sqrt(5)`. Thus if `n=Σ d_i A_i`, the shifted value satisfies

\[
P=\sum_i d_iA_{i+1}=\beta n+\rho,
\qquad \rho=\sum_i d_i\beta^{-i-1}.
\tag{B5}
\]

The reversed finite canonical word remains legal: the forbidden factors
`2 1^r 2` are unchanged by reversal. Regard the reversed digits as a
fractional expansion, read most significant first. For every finite legal
fractional word its value is:

* strictly less than 1 if its initial legality state is Free;
* strictly less than `tau` if its initial state is Locked.

**Proof by induction on length.** The empty fraction is zero, strictly below
both bounds. If the next digit is `d` and its remaining fraction is `T`,
the full fraction is `(d+T)/beta`.

* Free, digit 0 or 1: it is less than `1/beta` or `2/beta`, both below 1.
* Free, digit 2: the remaining state is Locked, so the value is less than
  `(2+tau)/beta=1`.
* Locked, digit 0: the remaining state is Free, giving a value less than
  `1/beta<tau`.
* Locked, digit 1: it is less than `(1+tau)/beta=tau`.

These are all permitted transitions. ∎

Therefore `0≤rho<1`, with `rho>0` whenever `n>0`. Since `P` is an integer,
(B5) proves (B1), including `n=0` separately by `P(0)=rho=0`.
The strictness is essential: the infinite legal fractional word
`2,1,1,1,...` has value 1, but no finite canonical word does.

## 3. The terminal state is exactly a fractional-phase threshold

The following proof does not infer the state merely from numerical tests
or from matching word counts. The infinite fractional word `1,1,1,...`
has value `tau`.

In an MSB-first canonical word:

* Locked means it ends in `2 1^r` for some `r≥0`.
* Free means it ends in `0 1^r`, or consists entirely of 1s (the empty word
  and zero are also Free).

When read in reverse, the first `r` fractional digits agree with the
infinite all-1 word. After them, a Locked word has residual fraction at
least `2/beta>tau`; a Free word with next digit 0 has residual less than
`1/beta<tau`, by the preceding strict bound. A finite all-1 word is also
strictly below the infinite all-1 word. Hence

\[
\text{Free}\iff \rho<\tau,
\qquad
\text{Locked}\iff \rho>\tau.
\tag{B6}
\]

Equality cannot occur in these finite words. Finally,

\[
P(n+1)-P(n)=\lceil\beta-\rho\rceil.
\]

Since `beta-tau=2` and `0≤rho<1`, this ceiling is 3 in the Free case and 2
in the Locked case. This proves the state/gap statement universally.

## 4. Division, reconstruction, and termination

Take `p=floor(N/beta)`. Then

\[
\beta p\leq N<\beta(p+1).
\]

Because `N` is an integer, this is equivalent to

\[
P(p)\leq N<P(p+1).
\]

Subtracting `P(p)` gives exactly the remainder range in (B3). By the proven
state/gap correspondence, that remainder is a legal digit after the
canonical word of `p`. The half-open intervals `[P(p),P(p+1))` are disjoint
and tile all nonnegative integers, so this parent and digit are unique.

Since `beta+1/beta=3`, the quotient is also the exact inverse companion
coordinate:

\[
\operatorname{parent}(N)=\lfloor(3-\beta)N\rfloor=3N-P(N).
\tag{B3a}
\]

For `N>0`, `p<N` because `beta>1`. Repeating the division therefore reaches
zero after finitely many steps. Record the remainders and reverse them:
the result is precisely the canonical Chung–Graham word. Conversely,
reading its digits left to right with `n <- ceil(beta*n)+d` recovers its
value. Zero has `(p,d)=(0,0)` and is handled as the terminating special case.

Examples:

| N | Parent `floor(N/beta)` | Remainder | Canonical word |
|---:|---:|---:|---:|
| 1 | 0 | 1 | 1 |
| 2 | 0 | 2 | 2 |
| 3 | 1 | 0 | 10 |
| 5 | 1 | 2 | 12 |
| 6 | 2 | 0 | 20 |
| 7 | 2 | 1 | 21 |
| 8 | 3 | 0 | 100 |
| 21 | 8 | 0 | 1000 |
| 32 | 12 | 0 | 1110 |

## 5. Exact phase dynamics, with its proper scope

Since `beta=2+tau`,

\[
P(n)=2n+\lceil\tau n\rceil.
\]

Consequently the indicator of the Free state is the upper mechanical word

\[
f(n)=\lceil\tau(n+1)\rceil-\lceil\tau n\rceil\in\{0,1\}.
\tag{B7}
\]

For every integer interval `[m,m+k)`, its exact count of Free prefixes is

\[
\sum_{n=m}^{m+k-1}f(n)
=\lceil\tau(m+k)\rceil-\lceil\tau m\rceil.
\tag{B8}
\]

This count differs from `k*tau` by strictly less than 1 and belongs to
`{floor(k*tau),ceil(k*tau)}`. Thus two intervals of equal length have counts
differing by at most one. The limiting Free frequency is the irrational
number `tau`; the sequence cannot be periodic or eventually periodic,
because such a binary sequence has rational limiting frequency. These
properties give the familiar irrational mechanical/Sturmian coding here.

This is a precise sense in which an initial fractional phase controls future
**numeral branching**. It does not say that square-sum graph existence,
Hamiltonian witnesses, or ninja-pair occurrence is a Sturmian sequence.

There is also an exact **lattice-strip/internal-coordinate model**. The
integer pair `(n,P(n))` lies in the strip `0≤rho=P-beta*n<1`. For a legal
appended digit `d`, the pair and internal coordinate evolve as

\[
(n,P)\longmapsto(P+d,\;3P-n+3d),
\qquad \rho\longmapsto\frac{\rho+d}{\beta}.
\tag{B8a}
\]

Digits 0 and 1 always remain in the strip; digit 2 does so exactly when
`rho<tau`. Digit 0 resets to Free, digit 1 preserves the current side of
the threshold `tau`, and digit 2 switches Free to Locked. Thus this window
recovers the full legality automaton. The internal coordinate is determined
by `n`, not an additional unconstrained hidden variable. This is a model of
canonical numeral construction, not ordinary pure-base-beta evaluation,
not a claim of a full cut-and-project scheme, and not a new graph inference.

The Wythoff relation is also exact. With the usual lower and upper Wythoff
sequences `a(n)=floor(phi*n)` and `b(n)=floor(phi²*n)`,

\[
P(n)=b(n)+1=n+a(n)+1\qquad(n>0),
\tag{B9}
\]

while `P(0)=0`. The exception at zero must not be dropped.

## 6. The 1,3,8 triple is literally a shift orbit

Appending zeros preserves canonical legality. The weight recurrence gives

\[
P(P(n))=3P(n)-n\qquad(n\geq0).
\tag{B10}
\]

Thus the orbit of 1 is

\[
1\longmapsto3\longmapsto8\longmapsto21\longmapsto55\longmapsto\cdots,
\]

which is exactly the even-Fibonacci weight chain. The earlier `D(1)`
invariant and these numeration/Beatty operations therefore act on the same
explicit sequence, not just similarly growing sequences.

### 6.1. Every positive shift orbit has a Diophantine parameter

Define the integer

\[
C(n)=n^2-3nP(n)+P(n)^2.
\tag{N1}
\]

For `n>0`, substitution of `P(n)=beta*n+rho` and
`beta²-3beta+1=0` gives

\[
C(n)=\rho(\sqrt5\,n+\rho)>0.
\]

The recurrence (B10) shows `C(P(n))=C(n)`: the quadratic form
`a²-3ab+b²` is unchanged by `(a,b) -> (b,3b-a)`.
Thus, for `a=n`, `b=P(n)`, `c=P(P(n))=3b-a`,

\[
ab+C(n)=(b-a)^2,\quad
ac+C(n)=b^2,\quad
bc+C(n)=(c-b)^2.
\tag{N2}
\]

Every three consecutive elements of a positive shift orbit therefore form
a `D(C(n))` triple. The orbit `1,3,8,21,...` has parameter 1, recovering the
user's triple; `2,6,16,42,...` has parameter 4.

### 6.2. Exact orbit classification, and a non-classifying invariant

For a positive integer, remove all trailing zero digits from its canonical
word. Let the remaining word represent `r`, and let `k` be the number of
removed zeros. Then

\[
n=P^k(r),
\]

where `r` is primitive: its last digit is 1 or 2, so it has no predecessor
under the positive shift map. This decomposition is unique, since
canonical words are unique and `P` appends exactly one zero. Two positive
integers lie in the same primitive shift chain exactly when their stripped
words agree. Their depths locate them in that chain.

The norm `C` is preserved on each chain but does not classify chains:

\[
C(5)=C(7)=11,
\qquad 5=12_{CG},\quad7=21_{CG}.
\]

These are distinct primitive cores. Also, norm preservation does not apply
to an arbitrary nonzero-digit quotient: `parent(5)=1`, but their norms are
11 and 1. It applies to moving backward only when removing a zero digit.

Zero is excluded from unique positive core/depth classification: `P(0)=0`
would make every depth possible. The code explicitly assigns the convention
`orbit_core(0)=(0,0)` and `C(0)=0`.

### 6.3. Norm-one characterization and the quadratic-unit identity

For a legal appended digit `d`, put `m=P(n)+d`. Because `parent(m)=n`,
(B3a) gives `P(m)=3m-n`. Expanding the norm yields

\[
C(m)=C(n)+d\bigl(2P(n)-3n\bigr)+d^2.
\tag{N3}
\]

Digit 0 preserves the norm. For `n>0`,
`2P(n)-3n>sqrt(5)*n>2`, so this integer is at least 3. Consequently each
later nonzero digit increases the norm by at least 4. The first nonzero
digit, appended to zero, gives norm 1 for digit 1 and norm 4 for digit 2.
It follows exactly that

\[
C(n)=1\iff (n)_{CG}=10^k\iff n=F_{2k+2},
\tag{N4}
\]

and

\[
C(n)=4\iff (n)_{CG}=20^k\iff n=2F_{2k+2},\qquad k\geq0.
\tag{N5}
\]

The lower bound of 4 matters for (N5): mere strict increase would not by
itself exclude reaching norm 4 from a word beginning with 1.

Finally let `D=2P(n)-3n`. Direct expansion gives

\[
D^2-5n^2=4C(n),\qquad D\equiv n\pmod2.
\tag{N6}
\]

The quadratic integer `alpha=(D+n*sqrt(5))/2=P(n)-n/beta` therefore has norm
`C(n)`. Under one canonical shift, it is multiplied by `beta`, whose
conjugate is `1/beta` and whose norm is 1. This is an elementary
Pell-type/quadratic-unit explanation of the invariant. It identifies the
norm-one chain precisely, without classifying arbitrary Pell solutions or
claiming any square-sum Hamiltonicity consequence.

Two important limitations:

1. **Repeated ceiling is not ceiling of a power.** `P(P(1))=8`, but
   `ceil(beta²)=7`, since `beta²=(7+3*sqrt(5))/2` lies strictly between 6
   and 7.
2. **Raw shifting and normalization do not commute.** The noncanonical word
   `22` has value 8 and normalizes to `100`. Yet the raw shifted word `220`
   has value 22, while the canonical shifted word `1000` has value 21.
   The formula `P(n)=ceil(beta*n)` requires the canonical representation.
   The least-place identity `3*A_0=A_1` does not shift to `3*A_1=A_2`;
   rather, `3*A_1=A_2+A_0`. Least-place carries explain the boundary effect.

## 7. Integer-only implementation

For `n>0`, put `q=isqrt(5*n*n)`. Since `n*sqrt(5)` is irrational,
`q<n*sqrt(5)<q+1`, giving the exact formulas

\[
\lceil\beta n\rceil=\left\lfloor\frac{3n+q}{2}\right\rfloor+1,
\qquad
\left\lfloor\frac n\beta\right\rfloor
=\left\lfloor\frac{3n-q-1}{2}\right\rfloor.
\tag{B11}
\]

Both functions are defined to return zero at zero. No floating-point
approximation or guessed epsilon is used, even for thousands of bits.

`beatty_bridge.py` provides the shift, parent, division, encoding, decoding,
state, and interval-count functions. `test_beatty_bridge.py` checks these
against the previously verified weight arithmetic and, independently,
against exact squared inequalities certifying each ceiling/floor. It also
verifies all five hashes in the old manifest before testing; the original
proof/code/reports remain immutable.

Run:

```sh
python test_beatty_bridge.py verification.beatty.normal.json
python -O test_beatty_bridge.py verification.beatty.optimized.json
```

Coverage includes integers 0 through 100,000, all legal words through length
8, 512 random 512-bit inputs, 32 random 4096-bit inputs, 3,000 large
Fibonacci-weight boundary cases, and 10,201 interval-count/balance cases.
The suite additionally checks `D(C)` identities, primitive-core/depth
classification, norm-one/norm-four classification, and Pell-type identities
for positive integers through 25,000, plus every legal child of those
integers. Both normal and optimized runs passed **1,983,501 explicit checks**
each; the JSON reports record these outcomes.
The JSON reports record actual outcomes and explicit check counts; finite
tests support rather than replace the universal proofs above.

## 8. Prior-art boundary

These are self-contained derivations for this project, not claims that
Chung–Graham arithmetic, Beatty/Wythoff sequences, or their connection is new.
[Rob Burns, *Shifting Zeckendorf and Chung-Graham representations*
(2025)](https://arxiv.org/html/2507.06162v1) directly studies these shift
operators, gives the even-Fibonacci/golden-ratio identity used in (B4), and
relates Chung–Graham and Zeckendorf shifts in Lemma 3. The introduction
discusses Wythoff and generalized Beatty descriptions. His convention
retains odd Fibonacci positions as zeros, so one compressed digit shift
here corresponds to two Fibonacci-position shifts there.

The underlying canonical representation and its automata were checked in
[Burns, *Chung-Graham and Zeckendorf representations*
(2025)](https://arxiv.org/html/2502.18870v1).
We do not claim that the exact presentation (B1)–(B11) is absent from prior
literature merely because it was derived afresh here.
