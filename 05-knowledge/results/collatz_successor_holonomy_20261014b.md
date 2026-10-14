# Native successor blocks: a finite-group return has no intrinsic descent sign

**Status.** PROVED: source-checked concatenation, both all-height modular-return families, the exact paid native cells and their density, their Mersenne complement, and an explicit infinite ROOT-certified subfamily. FINITE-EXACT: the declared controls below. CITED / INHERITED: AGP's finite-group construction principle and the repository's native-word and terminal-phase mechanisms. OPEN: universal source coverage and Collatz convergence. The infinite ROOT family is an explicit instance of the existing terminal compiler, not a new general basin theorem.

## 1. What transfers from a zero-product construction

[Alford–Granville–Pomerance, *There are infinitely many Carmichael numbers*](https://math.dartmouth.edu/~carlp/PDF/paper95.pdf)
uses specially constructed primes and a finite-group product condition to
meet Korselt's divisibility requirements. Its selected product is a new
integer under construction. The mechanism is useful here as a question:
which finite summary of several steps is sufficient for the target? The
analytic supply of primes is not transferred to Collatz.

The closest exact Collatz inheritance is
[signed affine cycles, §7](arithmetic_braids2_20260917_signed_cycles.md): a
finite-ring annihilator condition is necessary, but the integer carry still
matters. The [clock-holonomy compiler](collatz_clock_holonomy_20261007.md)
already distinguishes a checked common-future loop with nonzero integer
charge from a zero-charge loop. The
[native terminal compiler](collatz_grounded_child_ports_20261007f.md)
constructs ROOT families by a full ternary phase, retaining the source.
The recent [log-gap packet](collatz_log_gap_20261014a.md), and incoming
[THM-4633](../../01-canon/theorems/THM-4633-collatz-cycles-as-rational-approximations-of-log2-3-with-the-admissible-minimum-budget.md),
separate a cycle's clock from its carry. The real/parity-number dictionaries
in THM-4632/4634 are not premises for a native successor rule here.

Board: **ordered blocks / integer seam / native residue / modular affine
return / strict source payment / supplied ROOT terminal**. A finite quotient
retains a congruence, not the size of the lifted integer. Constructing a
source that matches a word is also different from routing an unchanged
prescribed source.

We use the positive odd map `U(n)=(3n+1)/2^a`, with the actual exponent
`a=v_2(3n+1)`. A valuation word `w` has carrier

\[
F_w(x)=\frac{P_wx+B_w}{Q_w},\qquad P_w=3^{|w|},\quad Q_w=2^{\sum w}.
\tag{1}
\]

Its exact odd native cell is

\[
x\equiv (Q_w-B_w)P_w^{-1}\pmod{2Q_w}.
\tag{2}
\]

The extra bit in `2Q_w` is essential: it makes the final value odd and fixes
every prescribed valuation exactly. This is an inherited native-word fact.

## 2. The actual linked-list operation

A receipt stores `(source,target,word)` and replays the word at that source.
Two receipts may be concatenated only when the first target is the second
source, as integers. Their composite keeps the original source and the
concatenated word. The verifier rejects a nonempty suffix after ROOT.
A supplied strict ROOT word at the final child may then be appended.

For example,

\[
43\xrightarrow{(1)}65\xrightarrow{(2,2)}37
\]

is an authenticated composition. An unrelated receipt `49 -> 37` cannot
replace the second block: it has the same eventual target but the wrong
input seam. The multiplier also does not permit reordering:

\[
F_{(1,2)}(x)=(9x+5)/8,\qquad
F_{(2,1)}(x)=(9x+7)/8.
\]

Their native cells are respectively `11 mod16` and `9 mod16`. The word
`(2,1)` is illegal at source 11. This is the elementary order/carry barrier
to treating arbitrary Collatz words as a commutative bag of factors.

For `gcd(M,6)=1`, (1) does induce an invertible affine transformation of
`Z/MZ`. Keeping its translation is stronger than retaining only its
multiplier. Even this full quotient has the obstruction below.

## 3. Full affine identity can grow at every native source

Let `M>=2` be coprime to 6 and choose positive `k` with
`(3/2)^k=1 modM`. Such a `k` exists in the finite unit group. Then

\[
F_{1^k}(x)=-1+(3/2)^k(x+1)\equiv x\pmod M.
\tag{3}
\]

On the entire positive native cell

\[
n=2^{k+1}t-1,\quad t\ge1,
\qquad U^k(n)=2\cdot3^kt-1>n.
\tag{4}
\]

Every letter is the actual valuation 1, and every intermediate state grows.
Thus even an identity in the full affine group, including modular carry,
does not imply a cycle, a descent, or ROOT. For `M=5,k=2`, the least example
is `7 -> 11 -> 17`, returning modulo 5 while increasing by 10.

This works simultaneously for any finite list of moduli coprime to 6:
apply it to their least common multiple. Increasing a fixed finite list of
modular observations does not repair the missing integer-height comparison.

## 4. A controlled concatenation that always pays

For every `k>=1`, set

\[
W_k=1^k\,2^{2k},\qquad
P=3^{3k},\quad Q=2^{5k},\quad
B=Q+P-2^{k+1}3^{2k}.
\tag{5}
\]

Here superscripts on letters mean repetition. This is a word of `3k` odd
edges and cost `5k`, not a formal iteration at an unspecified source.
Its exact native cell is

\[
n+1=2^{k+1}u,\qquad
3^ku\equiv1\pmod{2^{4k}},\qquad u>0.
\tag{6}
\]

Indeed the first block ends at `2*3^k*u-1`. For the following `2k`
valuations to be exactly 2, its difference from 1 must be divisible by
`2^(4k+1)`, which is precisely (6). Equivalently (2) gives the same single
residue modulo `2^(5k+1)`. This also proves `v_2(n+1)=k+1` exactly, because
`u` is odd.

**All-height payment theorem.** Every positive native source of `W_k`
has an actual endpoint `h` with `1<h<n`.

For the proof, let `a=(3/2)^k`, `b=(9/16)^k` and
`lambda=ab=(27/32)^k<1`. The first block expands about -1 and the second
contracts about 1, giving

\[
F_{W_k}(x)=\lambda x+1+\lambda-2b.
\]

Its fixed point `z_k` satisfies

\[
z_k=\frac{1+\lambda-2b}{1-\lambda}\le\frac{23}{5}.
\tag{7}
\]

To check (7), put `f_k=14*(27/32)^k-5*(9/16)^k`. Then `f_1=9` and
`f_(k+1)-f_k=(35/16)(b-lambda)<0`; this is exactly (7), with equality
only at `k=1`. For `k>=2`, the native source is at least `2^(k+1)-1>=7`.
For `k=1`, (6) gives `n=43 mod64`, whose least positive member is 43.
Thus every native source exceeds `z_k`, proving strict payment.
Finally `h=1+2*3^(2k)*(3^k*u-1)/2^(4k)>1`: the displayed quotient is a
positive integer. No earlier ROOT is possible, since `U(1)=1`.

The first two rays are

\[
43+64t\longmapsto37+54t,\qquad
455+2048t\longmapsto325+1458t,\qquad t\ge0.
\tag{8}
\]

The source grows for its first `k` odd steps before the paid endpoint.
This is a proved unbounded-depth family, but not a claim that every initial
growth pattern can be followed by the prescribed twos.

### Native coverage and its exact complement

The source determines its only possible index: `k=v_2(n+1)-1`. Therefore
the cells (6) are disjoint. Their relative natural densities among positive
odd integers are `2^(-5k)=32^(-k)`, giving

\[
\operatorname{dens}_{\mathrm{odd}}\left(\bigcup_{k\ge1}\text{cell}(W_k)\right)
=\frac1{31}.
\tag{9}
\]

This countable union has a genuine natural density: for `k>=K`, all its
members lie in `n=-1 mod2^(K+1)`, whose relative odd density is `2^(-K)`.
This uniform tail bound justifies passage from finite sums. The cells are
disjoint from the immediate-descent class `n=1 mod4,n>1`. With ROOT 1
handled separately, that declared baseline plus this bank has relative
**paid-dependency** density `1/2+1/31=33/62`. This is not density relative
to the repository's much larger existing banks, and not ROOT density.

Every positive Mersenne `n=2^E-1,E>=2` misses the entire bank. It would
require `k=E-1,u=1` and `v_2(3^k-1)>=4k`, whereas

\[
v_2(3^k-1)=\begin{cases}1&k\text{ odd},\\2+v_2(k)&k\text{ even},\end{cases}
\]

which is always less than `4k`. This preserves the previously difficult
Mersenne source boundary instead of hiding it inside a density statement.

## 5. The full finite-group return has both signs

Define

\[
d(M)=\operatorname{lcm}\left(
\operatorname{ord}_M(3/2),\operatorname{ord}_M(9/16)\right).
\]

The full affine map of `W_k` is the identity modulo `M` **iff** `d(M)|k`.
Indeed its multiplier is `ab` and its translation is `1+ab-2b`.
Identity forces `ab=1` and `2(1-b)=0`; since 2 is invertible, `b=1,a=1`.
The reverse implication is immediate. Thus the complete neutral subbank
has relative odd density `1/(32^d-1)` by the same tail proof.

For every such `k`, the growing word in (3) and the paid word in (5) both
represent the identity in the same full finite affine group. The integer
carry, native source and threshold decide the opposite signs. This is an
obstruction to the proposed transfer, not to AGP: Korselt's requested target
is a divisibility statement, whereas Collatz payment is an inequality at a
prescribed integer.

The initial quotient state does not fix this loss. For any chosen
`rho modM`, both native cells are single dyadic residue classes; their
periods are coprime to `M`. CRT therefore supplies infinitely many positive
sources in each cell with the same prescribed residue `rho`. Both full
affine quotient maps are the identity, yet their integer payment signs
are opposite. For `M=5,k=2,rho=2`, the exact examples are
`7 -> 17` via `(1,1)` and `8647 -> 6157` via `(1,1,2,2,2,2)`.

## 6. Grounded subfamilies, with an exact source test

The paid endpoint in (5) has no automatic ROOT proof. To ground an infinite
subfamily, require the one-step terminal endpoint

\[
h_j=(4^j-1)/3,\quad j\ge2,
\qquad
n_j=\frac{Q(4^j-1)-3B}{3P}.
\tag{10}
\]

It is integral exactly when

\[
4^j\equiv1+3B Q^{-1}\pmod{3P=3^{3k+1}}.
\tag{11}
\]

The right side is a principal unit. The order of 4 modulo `3^(3k+1)` is
`3^(3k)`, by `v_3(4^t-1)=1+v_3(t)`. Hence (11) specifies one exponent
phase, obtained by testing three digits at each ternary lift. Advancing
within that phase gives infinitely many `j>=2`.

Since `0<B<2Q`, every such `j` gives a positive numerator in (10).
That numerator is odd, as is its divisor `3P`, so `n_j` is positive odd.
The exact affine equality with odd endpoint `h_j` supplies (2), hence every
valuation in the head is native. The endpoint is above ROOT, and the final
valuation is `2j>=4`, yielding a strict first-hit certificate

\[
W_k\,(2j):\quad n_j\longrightarrow1.
\tag{12}
\]

For `k=1`, the phase is `j=10 mod27`, and (10) becomes the explicit family

\[
\boxed{n_s=\frac{32\cdot4^{10+27s}-101}{81},\qquad
\text{ROOT word }(1,2,2,20+54s),\quad s\ge0.}
\tag{13}
\]

The least displayed source is 414251, its head endpoint is 349525, and the
full word has four odd edges and shortcut cost 25. For example source 43
is paid by the same head but does **not** satisfy this one-step terminal
condition. It has a separate supplied strict child certificate from 37,
which the code verifies rather than discovering inside the production API.

`root_word_at_source(n,k)` preserves the given integer. It tests the exact
head endpoint and whether `3h+1` is a power of four with exponent at least 2,
then strictly replays the candidate. It never replaces `n` by a nearby
constructed family member. Symbolic plans store `k,j` and query `n_j mod m`
by computing the numerator modulo `3P*m` **before** exact division by `3P`.
This retains the full denominator for arbitrary even or ternary moduli.
The explicit materializer has a size cap; symbolic phase proofs are not
claims of constant-cost expansion of gigantic sources.

## 7. Exact verification and boundary of the result

The companion script checks complete native residues for `k=1,2`; eight
members for each `1<=k<=32`; all 33 moduli `5<=M<=101` coprime to 6, with
both growth signs and the complete tested neutral-index interval; exhaustive
principal-unit phases through depth 7; six literal ROOT sources for `k=1,2`;
and symbolic phase readers through `k=16`, including index `10^20+7`.
The native union contains exactly 132 odd sources at most 8191, independently
counted by its disjoint arithmetic progressions. Mersenne controls cover
`2<=E<=128`, in addition to the all-exponent proof above.

Hostiles reject a mismatched successor seam, reordered words at an unchanged
source, ROOT padding, adjacent-source borrowing, boolean/float aliases,
corrupted receipts and plans, wrong native cells, and insufficient explicit
materialization budgets.

```text
python -B 04-computation/experiments/collatz_successor_holonomy_20261014b.py
python -B -O 04-computation/experiments/collatz_successor_holonomy_20261014b.py
```

Both runs give **6,226 checks** and the saved
[JSON](collatz_successor_holonomy_20261014b.json), LF SHA256
`eafd419754f827a6883564520575aff59092f5eaeab3390624fc88363c6166ff`.
The constructive quantifier is: for each permitted modulus there are native
paid cells and explicit ROOT-certified subfamilies. The unproved target
would be a source-preserving selector that covers every input and grounds
every emitted child. These statements are not interchanged.
