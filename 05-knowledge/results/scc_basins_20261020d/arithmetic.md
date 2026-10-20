# Golden SCC words, deficit powers, and the marked 105 boundary map

**PROVED:** the SCC/Fibonacci bijection, the triangle-potential identity,
fixed-insertion formula and asymptotics, restricted terminal count, and
marked golden-character incidence below. **FINITE-EXACT:** literal small
tournament inventories and the ranks of the specified boundary banks.
The existing arithmetic105 dictionary is inherited and replayed here.
No all-order black connectivity or Petersen-forcing theorem is inferred.

## 1. Recovering the meaning before using a symbol

The new request describes `D(z)` as a series of strong centers with
deficit two. This notation is not defined in either of the two cited
incoming packets. In particular, it must not be replaced silently by
the earlier, explicitly chosen metagraph degree polynomial.

Three actual parameters need separate names:

| Parameter | Definition | Four-vertex strong tournament R4 |
|---|---|---:|
| Black presentation deficit | `beta=h-b`, with h Hamiltonian-path orbits and b blue orbits | `5-3=2` |
| Symmetry slack | `delta=3f-s`, feedback f and common moved support s | `3*1-0=3` |
| Heavy-gap deficiency | `D=2k+1-h_heavy` in a capped weighted-core presentation | Depends on the declared budget/core |

[Width and symmetry, Theorems4–6](../continuation_width_symmetry_20261020/README.md)
prove that symmetry slacks one and two are impossible; its zero-slack
alphabet is `{1,C3}`. At slack three its strong exceptional component is
`C3[1,U,1]`, for any singleton/triangle word U on at least two vertices.
Thus R4 is only the smallest of that *different* family. The
[deficiency equation](../defect_transitions_20261020/strata.md) concerns
heavy gaps and backward endpoints, and is a third grading.

The compatible interpretation of the new `D(z)` is therefore

```text
D_2(z) = sum_[strong center T with beta(T)=2] z^|T|.
```

This is a reconstruction from the user's description, not recovered old
notation. Its all-order classification belongs to the companion center
analysis. The results below require only the actual atom R4; they do not
assume the assertion `D_2(z)=z^4` in order to prove it.

Inheritance: the closest mechanism is unique ordered SCC decomposition;
the corrected near miss is confusing the three deficiencies above; the
hostile is an order join of two R4 atoms; the least-used sidecar is the
terminal symbol in a finite prefix code. Our concept board is SCC words,
triangle potential, black deficit, golden recurrence, and marked closure
incidence. These objects are related by specified maps below.

## 2. A genuine golden-ratio SCC language

Let `a`, `b`, and `c` denote, respectively, the singleton, directed
triangle C3, and unique strong four-vertex tournament R4. Consider only
ordered SCC words in these three atoms. Their vertex weights are1,3,4.
Unique SCC decomposition makes these words distinct ordinary tournament
types. If `a_n` counts the words on n vertices, including the empty word,

```text
F(z)=sum a_n z^n = 1/(1-z-z^3-z^4)
                    =1/[(1+z^2)(1-z-z^2)].             (1)
```

This is a restricted family. Its first terms are

```text
1,1,1,2,4,6,9,15,25,40,64,104,169,273,441,714,1156.
```

At order five it has six ordinary types, while there are twelve tournament
types in total. It is not a golden formula for unrestricted tournament
enumeration. The neighboring value104 must not be changed to105 to fit
the arithmetic analogy.

### A bijection, not just a shared denominator

A Fibonacci composition is a word on symbols1 and2, of total weight n.
Parse it from the left using the prefix tokens

```text
1  <-> a,       21 <-> b,       22 <-> c.                (2)
```

If the last unparsed symbol is2, retain it as a terminal flag. The token
set is prefix-free and the parser consumes every word, so this gives a
weight-preserving bijection

```text
Fibonacci words <-> SCC words x {empty, terminal2}.
```

The flag contributes two vertices to the word weight; it is additional
data, not a tournament atom. Its necessity is visible at total weight2:
the Fibonacci word `2` has no token, while `11` gives two singleton SCCs.
Forgetting the flag is not a weight-preserving identification.

Writing `Fib_0=0,Fib_1=1`, equation(2) proves directly

```text
Fib_(n+1)=a_n+a_(n-2),        a_j=0 for j<0.            (3)
```

If `L_n` are Lucas numbers and `c_n=[z^n](1+z^2)^(-1)`, then

```text
a_n = [L_(n+2)+2c_n+c_(n-1)]/5,                       (4)
c_(2j)=(-1)^j, c_(2j+1)=0, c_(-1)=0.
```

Indeed the partial fractions of(1) are

```text
(z+3)/(5(1-z-z^2)) + (z+2)/(5(1+z^2)).
```

Consequently `a_n~phi^(n+2)/5`. The number5 here comes from the exact
partial fractions; no action of a five-element group on these words is
being asserted.

### The map preserves the descent potential

A directed triangle contributes one directed3-cycle; R4 contributes two.
No directed cycle crosses ordered SCCs. Hence a word with t triangles
and r R4 atoms has

```text
c3=t+2r.
```

Under(2), this is precisely the number of domino symbols2 in the decoded
Fibonacci word. With the terminal flag epsilon, the count is `c3+epsilon`.
Thus the bivariate identity is

```text
F(z,y)=1/(1-z-yz^3-y^2z^4)
      =1/[(1+yz^2)(1-z-yz^2)],                         (5)
```

where y marks directed triangles. This is the same potential used in
[all-order black descent](../black_descent_20261020c/cycles.md), not a new
analogy chosen only by dimension. The full word and terminal flag are
retained by the inverse map.

**Operation boundary.** A Hamiltonian-path reversal can leave this SCC
sublanguage. Equation(5) therefore does not conjugate the black metagraph
to a graph of local domino moves, and does not prove connectivity of
the terminal basins.

One operation nevertheless already has an authenticated image. The
[four-black SCC interchange](../black_descent_20261020c/scc.md), section5,
changes `A ->1 ->C3` to `A ->C3 ->1` when A begins with two dominating
singletons. Under(2) this is the rewrite `code(A)121 <-> code(A)211`,
with the guard that `code(A)` begins11. The four actual Hamiltonian paths
are its certificate, including any intermediate tournaments outside the
language. This is a valid guarded rewrite, not an assertion that every
unrestricted adjacent swap has a black implementation.

## 3. What squares, cubics and quartics count

Fix exactly r occurrences of R4. A word decomposes uniquely as

```text
W_0 c W_1 c ... c W_r,
```

where each W_i is a word on singleton and triangle atoms. Put
`A(z)=1/(1-z-z^3)`. The exact series is

```text
F_r(z)=z^(4r) A(z)^(r+1).                              (6)
```

Thus the powers are convolution powers of independent gap words.
`D(z)^r` alone, even if `D(z)=z^4`, only records the r distinguished atoms;
it omits their r+1 contexts. Explicitly,

```text
[z^n]F_r(z)
 = sum_(t>=0, 3t+4r<=n)
     (n-3r-2t)! / [r! t! (n-4r-3t)!].                 (7)
```

Let rho be the positive solution of `rho+rho^3=1`, and
`lambda=1/rho`, so `lambda^3=lambda^2+1`. At each fixed r,

```text
[z^n]F_r(z)
 ~ rho^(3r-1) / [r! (1+3rho^2)^(r+1)] * n^r lambda^n. (8)
```

Proof: the positive zero of `1-z-z^3` is simple and is the unique zero
of smallest modulus. Near it the denominator is
`rho(1+3rho^2)(1-z/rho)+O((z-rho)^2)`. Substitution into(6) gives a pole
of order r+1; the binomial coefficient of that pole has leading term
`n^r/r!`. Other poles are exponentially smaller. The argument and
constant apply for every fixed nonnegative integer r.

This gives a precise sense in which two, three and four inserted atoms
produce quadratic, cubic and quartic prefactors. It does not say the
coefficients are literal squares, cubes or fourth powers of a sequence.
At fixed r the exponential rate remains Narayana's lambda, about1.4656;
summing over unbounded r changes it to phi, about1.6180. This distinction
also agrees with the incoming fixed-symmetry-slack growth theorem:
in this restricted family, `f=t+r`, `s=3t`, and `delta=3r`.

### Black deficit is not additive under these insertions

For any such SCC word, Hamiltonian paths concatenate and automorphism
groups multiply. Since C3 has three paths and three automorphisms, while
R4 has five paths and trivial automorphism group,

```text
h=5^r.
```

The word is self-converse exactly when it is a palindrome. If it is not,
`b=0`. For a palindrome, an outer SCC pair contributes the path-orbit
count of its left member, and the central atom contributes its own blue
count. Therefore

```text
b=5^floor(r/2) * (3 if the central atom is R4, else1).
```

In particular the joins of one through four copies of R4 have

| r | h | b | beta=h-b |
|---:|---:|---:|---:|
|1|5|3|2|
|2|25|5|20|
|3|125|15|110|
|4|625|25|600|

Already the second row rejects additive beta and a naive interpretation
of `D^2` as a deficit-four series. Formula(6) counts inserted atoms,
which is a different and well-defined grading.

## 4. A counted family of actual terminal boundary states

All three atoms are self-converse, so converse reverses their SCC word.
The palindrome series and the converse-merged class series are

```text
P(z)=(1+z+z^3+z^4)F(z^2),
M(z)=(F(z)+P(z))/2.                                   (9)
```

Within this restricted language, beta is zero precisely for palindromes
with no R4 atom. Their series is

```text
P_0(z)=(1+z+z^3)/(1-z^2-z^6).
```

The only strong black-active member of this language is R4 itself.
Consequently the restricted nonstrong black-active merged types have
the exact series

```text
T(z)=M(z)-P_0(z)-z^4
    =z^4+2z^5+4z^6+6z^7+12z^8+18z^9+32z^10
      +50z^11+85z^12+134z^13+222z^14+... .             (10)
```

These are actual strict triangle minima by
[the SCC basin theorem](../black_descent_20261020c/scc.md). They are a
subset of all terminals, not a count of distinct connected components.
In particular `T_n~phi^(n+2)/10`. Equations(5) and(10) provide an explicit
language of terminal boundary states with a preserved potential, while
the black routes between those states still require their path witnesses.

## 5. The old 105/golden map reaches a precise closure bank

The [105 research atlas](../105_research_atlas_20261018f.md), section3,
already gives the arithmetic character construction

```text
105 -> Q(sqrt3,sqrt5,sqrt7) -> eight quadratic characters
    -> H8 plus a repeated pair ->105 matching atlases.
```

Its field/class-group statement is inherited, not reproved here. The
following part uses only its explicit character cube. Address characters
by `F2^3`, with radical labels in binary order

```text
1,3,5,15,7,21,35,105.
```

Mark the trivial character0. The golden character is the vector2,
corresponding to sqrt5, since it detects whether a field automorphism
sends `phi` to `1-phi`. Translation by that character gives the parallel
matching

```text
{1,5}, {3,15}, {7,35}, {21,105}.                       (11)
```

Now apply the actual incidence R from
[the two105 carriers](../black_descent_20261020c/105_17.md), section2:
delete the pair containing the marked0; choose one of the other three
pairs; combine that pair with the free point to make a cubic junction;
retain the other two pairs as equality wires.

For(11), its three resulting junction closures are exactly

| Junction triple | Equality wires |
|---|---|
|3,5,15|7–35 and21–105|
|5,7,35|3–15 and21–105|
|5,21,105|3–15 and7–35|

The three triples are the three Fano lines through the golden character.
For any nonzero direction d, the same construction gives the three Fano
lines through d, with complementary wires parallel to d. All seven
directions therefore specify21 distinct junctions inside the full105 bank.

This preserves labelled containment, the chosen character, and the exact
closure kernel. The action on characters is multiplication by a chosen
character; it is not a field automorphism permuting quadratic subfields.
The trivial character and the three radical-coordinate choices remain
part of the decoder.

### The bank has a precise, substantial loss

Let J be the21-by91 matrix of these junction responses on the odd
seven-port colour orbits. Let Q sum each direction's three rows, giving
seven parallel responses. The exact calculation gives

```text
rank J=21,
Q Q^T=20 I+7 all-ones,
rank Q=7, eigenvalues69 once and20 six times.           (12)
```

An independent fingerprint is

```text
char(J J^T)=(t-23)(t-5)^8(t^2-21t+94)^6.
```

Every row is an actual junction probe, so every one annihilates all
seven inherited Kempe blind vectors `w_i`. The golden mark selects a
well-defined piece of the closure bank; it does not recover the missing
Kempe-path coordinates or force a Petersen minor. In particular it must
not be identified with the visible seven or invisible seven merely
because each has seven coordinates.

## 6. Which other golden105 inheritance survives

The earlier atlas already distinguishes three arithmetic clocks:

| Arithmetic carrier | Order of the golden unit |
|---|---:|
|`Z[phi]/(105)`|80|
|`Z[phi]/(Phi_105(phi))`|105|
|`Z[phi]/(16271615641)`|210|

Its exact evaluation is

```text
Phi_105(phi)=5224431949+8453308464 phi,
Norm=16271615641=21211*767131.
```

The checker independently replays the quotient root, its exact order105,
the conjugate root's order210, and the scalar-modulus period80. These
different clocks cannot be substituted for each other or for a path
count of105. Nor can one repeated eight-label permutation have order105
or80; its possible orders are too small. The marked incidence(11) is the
lawful boundary map, rather than an alleged shared105-step dynamics.

The [strong-tournament/golden orbit paper trail](../strong_tournament_balance_golden_rigidity_20261014.md)
uses still another observer: Hamiltonian-path *pair balance*, and golden
unit orbits at prime powers. Its THM-4655 rotation argument is related to
choosing cycle cuts, but neither that balance statistic nor its Pisano
clock counts the black deficit beta. Its explicit correction of the
fixed-cycle balance claim is consistent with the later need to change
cycles in black descent. We do not import its open balance conjecture.

Finally, the inherited Fibonacci Euclid map sends parameters13,8 to
`(105,208,233)`. It preserves the Pythagorean identity and its actual
parameters. Consecutive Fibonacci ratios tend to phi, for which

```text
(phi^2-1)/(2phi)=1/2.
```

Thus those two legs approach ratio1:2, not the isosceles ratio1:1.
Equal Euclid legs instead require parameter ratio `1+sqrt2`. The golden
SCC growth in(1) and the geometric angle are distinct quantities; the
boundary-support pairing remains the faithful obstruction coordinate.

## 7. Verification and transfer contracts

Run from the repository root:

```text
python3 05-knowledge/results/scc_basins_20261020d/arithmetic_probe.py
python3 -O 05-knowledge/results/scc_basins_20261020d/arithmetic_probe.py
```

The receipt records4,180 literal Fibonacci-word codec tests through
weight16;440 literal SCC tournament path/blue inventories through order12;
the separate four-R4 hostile; exact recurrence coefficients through80 for
r=0,...,4; all91 seven-port colour orbits and21 marked junctions; the
golden clock checks; and the normalized source hash. Conditions use
explicit exceptions rather than optimization-disabled assertions.
Normal and optimized runs passed with identical stdout and byte-identical
JSON receipts.

| Connection | Preserved | Required retained data / failure boundary |
|---|---|---|
|SCC words -> Fibonacci words|Vertex count and directed-triangle potential with flag correction|Full token word and terminal flag; black moves need not stay in the language|
|r R4 insertions -> convolution power|Exact class count and insertion number|r+1 gap words; insertion number is not additive black deficit|
|Converse -> word reversal|Actual class involution and palindrome fixed types|SCC word; counts do not determine a basin's graph neighbors|
|Golden character -> three junctions|Marked Fano incidence and closure tests|Character/address mark; all seven Kempe coordinates remain invisible|
|Fibonacci parameters -> Euclid triple|Primitive-triple parameters and leg ratio|The ratio is not the golden SCC growth rate|

The useful next operation is to transport an authenticated black route
between terminal SCC words through(2), retaining every intermediate
tournament that leaves the language. A purely scalar recurrence, even
one with an exact golden factor, does not supply that missing route.
