# Atomic prefix mass, multifractal addresses, and discounted Collatz receipts

2026-10-05. **PROVED, scoped:** the source measure and cylinder formula;
its lower local-dimension spectrum; the atomic coverage criterion; the
discounted-flow equivalence; the height-profile and sibling-mass formulas.
**CITED:** the normalized algorithmic-randomness criterion and the previously
recovered semigroup theorem. **FINITE-EXACT:** the accompanying program.
**OPEN:** decay of the mass without root certificates, construction of the
discounted flow without first knowing convergence, and universal Collatz.
These are elementary author-audited results, not a canon promotion or a
claim of literature priority.

[Program](../../04-computation/experiments/collatz_effective_prefix_mass_20261005.py)
and [exact output](collatz_effective_prefix_mass_20261005.json).

## 1. Inheritance and the changed target

The anchor is universal entry into sound root certificates. The niche is
computable atomic measures on a binary address tree. The wildcard is the
finite-observer interpretation of the repository's old Vitali/tournament
analogy. The concept board is **prefix machines / source height / labelled
defects / sign-sensitive survivors / discounted positive flow / quotient loss**.

Closest proved mechanisms:

- [Universal weak receipts, W1 and G1](collatz_universal_weak_receipts_20261005.md)
  give total weak coverage but retain a labelled endpoint defect.
- [Four helpers, T2](collatz_fusion_helpers_20261005.md) preserve that entire
  defect under a stock-supported common-future replacement.
- [THM-4512, coefficient-descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
  and [exact ports](creative_descent_20260925.md) distinguish a cylinder's
  contracting coefficient from its exceptional least representative.
- [THM-2966, critical-run fair extractors](../../01-canon/theorems/THM-2966-spine-normal-form-for-critical-run-fair-extractors.md)
  preserve undecided spine terms rather than discarding them prematurely.

Canonical hostile: identical parity words are realized on the two Collatz
signs, while the minus sign has the known cycles with minima 1, 5 and 17.
Corrected near miss: endpoint coefficient balance is not all-prefix safety;
see the [Pascal prefix correction](collatz_pascal_prefix_gate_20261005.json).
Least-used sidecar: the ordinary size of the **canonical source representative**,
not just the count of legal words or their binary/ternary grades.

The new target is a mass statement in which **one exceptional positive integer
cannot have zero weight**. This removes one quantifier mismatch. It does not
establish the mass estimate required by that target.

## 2. What the paper's randomness criterion does and does not supply

[Reimann, Information vs Dimension, sections 4.3–4.4 and 5](https://arxiv.org/html/2408.05121v1)
connects prefix complexity, randomness for a computable measure, and pointwise
dimension. Use the normalized criterion

\[
K(x\upharpoonright m)\ge -\log_2\mu[x\upharpoonright m]-O(1)
\]

for a measure-random point. The minus sign and additive constant matter; the
quoted/rendered passage is not usable literally without them. Effective
dimension is the liminf of complexity divided by prefix length. A universal
lower semicomputable semimeasure is a different object from a computable
probability measure. Those are the literature inputs used here; the following
construction and proofs are given explicitly.

Write

\[
T_\sigma(n)=\begin{cases}n/2&n\text{ even},\\(3n+\sigma)/2&n\text{ odd},\end{cases}
\qquad \sigma\in\{+1,-1\}.
\]

**P1 — computable-itinerary obstruction.** For every fixed integer n, the
infinite parity itinerary Q_sigma(n) is computable, regardless of whether
the orbit reaches a cycle. A program given n and m simply executes m steps.
Consequently

\[
K(Q_\sigma(n)\upharpoonright m)\le K(n)+K(m)+O(1)=O(\log m),
\qquad \dim_{\rm eff}Q_\sigma(n)=0.
\]

A positive entropy randomness lower bound therefore cannot apply to all
these itineraries. Conversely, zero effective dimension does not imply
eventual periodicity: computable nonperiodic sequences exist. Description
length also does not bound the time needed to discover a terminating route.
This recovers, in a sharper form, the warning in the
[mod-six compression audit](collatz_mod6_20260922_compression_and_lean_audit.md).

The reroute is to change the measure, retaining the integers as atoms.

## 3. A complete prefix machine for choosing integer sources

Choose h>=1 with probability 2^-h, then uniformly choose an odd integer
1<=n<2^h. This is implemented by the prefix code

\[
1^{h-1}0\;\operatorname{bin}_{h-1}(j),\qquad n=2j+1,
\quad0\le j<2^{h-1}.
\]

Each codeword has length 2h-1. Its Kraft sum is
sum_h 2^(h-1) 2^(-(2h-1))=1. The same integer has different descriptions
at every sufficiently large h; injectivity is required for descriptions of
(h,j), not for their decoded integer values. The all-one input stream never
finishes the unary header and has fair-coin measure zero.

Let ell(n)=bit_length(n). The induced probability measure on odd 2-adics is

\[
\boxed{w_n=\mu\{n\}=\frac8{3\,4^{\ell(n)}}},\qquad
\boxed{\mu\{n\ge2^H\}=\frac2{3\,2^H}}\quad(H\ge1).\tag{1}
\]

Indeed n appears at every h>=ell(n), with weight 2^(1-2h).
The geometric sum gives (1). In particular w_1=2/3, w_3=1/6,
w_5=w_7=1/24, and w_17=w_27=1/384.

**P2 — exact cylinder formula.** For K>=1 and canonical odd
1<=r<2^K,

\[
\boxed{\mu(r+2^K\mathbb Z_2)=w_r+\frac4{3\,4^K}.}\tag{2}
\]

For mixture levels h<K, a cylinder contains its canonical integer r exactly
when h>=ell(r). For h>=K, it occupies the fraction 2^(1-K) of the finite
uniform head. Summing those two geometric series proves (2). Thus mu is
computable, assigns positive weight to every positive odd integer, and has
full topological support despite being concentrated on a countable set.

Each parity word of length K determines exactly one source residue modulo
2^K. To append a parity bit, one of the two source lifts has that parity,
since the current affine numerator has odd coefficient. This proves by
induction that Q_sigma is an isometry of the binary address trees. Pushing
mu through Q_sigma gives a computable measure on itineraries with the same
cylinder formula, evaluated at the word's canonical source r_sigma.

For a port n=r+2^K k, the extra coordinate k=(n-r)/2^K is the actual lift;
the bit length of r measures the mass of the whole cylinder, not the height
of every lift. Neither coordinate may be silently substituted for the other.

## 4. The fractal is in the distribution, not only in its support

Put t_K=4/(3*4^K). Refinement of a source cylinder gives

\[
\mu[r]_{K+1}=w_r+t_K/4,\qquad
\mu[r+2^K]_{K+1}=4^{-K}=3t_K/4.\tag{3}
\]

Since w_r>=2t_K, a newly revealed high binary one has conditional probability
at most 1/4. A high zero can have conditional probability arbitrarily close
to one. Hence there is no fixed positive information cost per refinement.

At a fixed positive integer n, the cylinder mass decreases to w_n>0.
At a negative odd integer -m, for sufficiently large K the representative
is 2^K-m with bit length K, so its cylinder mass is exactly 4^(1-K).
The latter has local dimension 2, though its parity sequence is computable
and has effective dimension zero. It is not random for this mu.

Every atom is mu-Martin-Lof random: a test component of measure at most
2^-j cannot contain an atom whose mass exceeds 2^-j. This is consistent
with P1, because cylinder information remains bounded along the atoms.
Randomness is relative to a measure; it is not an intrinsic assertion of
unpredictability here.

**P3 — explicit lower local-dimension spectrum.** For an odd 2-adic x, let
r_K be its canonical representative and ell_K=ell(r_K). In the usual
binary ultrametric,

\[
\underline d_\mu(x)=2\liminf_{K\to\infty}\frac{\ell_K}{K},\qquad
f_\mu(\alpha):=\dim_H\{x:\underline d_\mu(x)=\alpha\}=\frac\alpha2
\quad(0\le\alpha\le2).\tag{4}
\]

The first equality follows from (2), uniformly bounded between
(8/3)4^-ell_K and 4*4^-ell_K. For the spectrum's upper bound, set
beta=alpha/2. Infinitely many K have ell_K<=(beta+epsilon)K; at most
2^((beta+epsilon)K) cylinders do so. Summing their s-dimensional costs over
all K gives Hausdorff dimension at most beta+epsilon, then beta.

For 0<beta<1, construct a Cantor set in stages. Take a free block ending at
t_j, with t_j much larger than all preceding lengths; fix a one at its end
and at least every j-th position inside it. Then force zeros until length
floor(t_j/beta). Choose t_j growing fast enough that these intervening
fixed ones have asymptotic density zero. The longest new zero block makes
ell_K/K approach beta; elsewhere the bounded gaps between the forced ones
prevent a smaller limiting ratio. The lower density of free positions is
beta. Giving their bits independent fair weights assigns each K-cylinder
mass 2^(-F(K)), where F(K) is its number of free positions and
liminf F(K)/K=beta. The elementary mass distribution estimate for every
s<beta gives dimension at least beta. For beta=1, forcing ones at square
positions has zero density and ell_K/K->1, giving dimension one. For beta=0,
positive integers supply a nonempty level set and the upper bound is zero.

There are 2^(ell-2) odd integers of bit length ell>=2, and one of length 1.
The same comparison proves, for every real q,

\[
\tau(q):=\lim_{K\to\infty}-\frac1K\log_2
  \sum_{r\ {m odd}<2^K}\mu([r]_K)^q=\min(0,2q-1).\tag{5}
\]

The shell sum is comparable to sum_(ell<=K) 2^((1-2q)ell);
at q=1/2 it grows linearly with K. Its Legendre infimum is alpha/2 on
[0,2], agreeing with (4). All of this coexists with **dim_H(mu)=0**:
a countable set carries its entire mass. The support's dimension is one.
These three dimensions, and effective dimension, must remain distinct.
The breakpoint 1/2 comes from our chosen squared-height weighting; it is
not evidence of a distinguished arithmetic prime or a Collatz exponent.

## 5. Universal coverage becomes an exact mass-one question

**P4 — atomic coverage criterion.** For any set C of positive odd sources,

\[
\mu(C)=1\quad\Longleftrightarrow\quad C\text{ contains every positive odd integer}.
\tag{6}
\]

Omitting n loses at least w_n. No randomness assumption is needed.

Let a prefix machine parse the code in section 3, then search for a checked
root certificate for n. Its halting probability is

\[
\Omega_{\rm ROOT}=\sum_{n\text{ certified}}w_n.
\]

With a complete enumeration of actual finite routes, Omega_ROOT=1 is
equivalent to positive Collatz. This particular machine is **not a universal
Turing machine**. A universal prefix machine does not make a particular
partial search total. Its undecidability phenomena cannot be silently
transferred to this fixed arithmetic machine either.

In contrast, the inherited semigroup theorem gives a total search for
weak receipts, so **Omega_WEAK=1 is already proved**. Universal weak coverage
and root coverage now occupy the same measure space, with the endpoint
defect as the still-unpaid difference. Changing the measure has not repaired
that defect.

For actual shortcut root time tau(n), with infinity allowed, put

\[
S_t=\sum_n w_n\,\mathbf1_{\{\tau(n)>t\}}.
\]

Every S_t is computable: simulate a finite head for exactly t steps and
bound its remaining mass by (1). Monotone convergence and (6) give
S_t->0 iff every source reaches 1. Any proved upper bound S_t<=b(t)->0
would be decisive. Once b(t)<w_n, source n cannot remain unresolved.
For negative integers use their positive magnitudes under T_minus and a
specified set of cycle targets. Taking {1,5,17} tests entry into the three
known cycles; the assertion that this catalogue is exhaustive remains open.

A universal lower semicomputable source semimeasure can give larger weight
to simply described huge integers, but is not a drop-in computable prior.
No computable probability on a countably infinite domain dominates all
lower semicomputable semimeasures within a fixed constant: computably pick
distinct n_j with mu(n_j)<2^(-2j), then place mass 2^-j at n_j. A universal
semimeasure dominates this new one, contradicting such domination by mu.
Program-based priors therefore need an explicit time budget or an
effectivity sidecar. Simplicity and fast certification are separate resources.

## 6. What the two signs reveal about prefix mass

Let W_D consist of length-D parity words starting odd whose coefficient
3^(number of odd steps)/2^(prefix length) is >=1 at **every** nonempty prefix.
Let N_(D,ell,sigma) count their canonical source representatives of bit
length ell on sign sigma. The exact prefix mass is

\[
\boxed{\mu(S_{D,\sigma})=
\frac{8\sum_{\ell=1}^D N_{D,\ell,\sigma}4^{D-\ell}+4|W_D|}
     {3\,4^D}.}\tag{7}
\]

Reflection n->-n conjugates the signs, so the two representatives of any
word satisfy r_minus=2^D-r_plus. The word count and Haar mass agree exactly;
the height histogram need not. In the extreme all-odd word,

\[
\mu(S_{\text{all-odd},+})=4^{1-D},\qquad
\mu(S_{\text{all-odd},-})=\frac23+\frac4{3\,4^D}.\tag{8}
\]

This is an exact sign-sensitive distinction made entirely from finite
prefixes and ordinary representatives. It answers a different question
from whether a bounded residue observer can distinguish every future orbit.

The known minus cycles with minima 1,5,17 stay above their minima.
If a minus prefix had coefficient <1, its negative carry would send a
positive source strictly below itself. Therefore these three minima survive
every coefficient prefix, giving the rigorous lower bound

\[
\mu(S_{D,-})\ge w_1+w_5+w_{17}=91/128.
\tag{9}
\]

| Shortcut depth D | Survivor cylinders | Haar mass, conditional odd | mu plus | mu minus |
|---:|---:|---:|---:|---:|
| 8 | 19 | 19/128 | 0.008117676 | 0.716247559 |
| 16 | 2114 | 1057/16384 | 0.007556504 | 0.711350212 |
| 24 | 286581 | 286581/8388608 | 0.007379534 | 0.711110131 |
| 26 | 1037374 | 518687/16777216 | 0.007374011 | 0.711107644 |

Displayed decimals are rounded presentations of exact rationals in the JSON.
The plateau is informative: 27 alone contributes 1/384 through depth 58;
its first coefficient exit is at shortcut depth 59. A small collection of
low ordinary representatives can dominate a very large word tree.

For the plus sign a coefficient exit is not automatically actual descent.
At an exact port n=r+Qk with endpoint s+Pk and P<Q, actual descent is
equivalent to k>(s-r)/(Q-P). The program retains these exceptional lifts.
Through D=26 the entire accumulated exceptional set is {1}. This is a
finite-depth, all-lifts statement, not an all-depth claim. Coefficient
survival and actual root survival remain different events.

Separate direct root computations certify all 65,536 odd sources below
2^17. The maximum shortcut root time in that head is 223, at 106239.
For the minus sign, entry into {1,5,17} takes at most 230 in the same head,
at 126465. Thus in both cases S_256<=1/196608. This is a rigorous tail bound
from a finite census, not a bound tending to zero with t for a fixed proof.

## 7. Incoming work: the same sibling rule has a different mass

Incoming commit `9eb19b215d` develops the
[two-sheet sibling automaton](collatz_two_sheet_receipts_20261005.md).
Its exact quarter-child rule U(4s+1)=U(s) is useful independently of its
probabilistic interpretation. The following derivation checks the deep cell.

Let x=t*2^k-1 with odd t, k>=4, and U=oddpart(3n+1). Put

\[
y=9t2^{k-3}-1,\qquad X=U^2(x)=2y+1.
\]

For k-4 simultaneous U steps the pair remains X_i=2y_i+1.
At the next step the valuation is

\[
b=1+v_2(3^{k-1}t-1)\ge2.
\]

The first structural merge is certified when b=2, equivalently

\[
\boxed{t\equiv3^k\pmod4.}\tag{10}
\]

At simultaneous time k-3 the pair is a quarter-child pair; equality follows
at time k-2. If b=3, the stated automaton next enters a negative-sign state
with exponent at least two and stops certifying; if b>=4 it stops directly.
These are failures of that structural rule, not proofs that the two orbits
never meet. For x=31, t=1, k=5, the rule fails, but x's orbit reaches 35
and then 53, while its shadow y=35 reaches 53 in one step.

**P5 — change-of-measure test.** In the exact cell E_k={v_2(x+1)=k}, the
Haar conditional probability of (10) is 1/2. Formula (2), however, gives

\[
\mu(E_k)=3/4^k,\quad
\mu(t\equiv1\bmod4\mid E_k)=11/12,\quad
\mu(t\equiv3\bmod4\mid E_k)=1/12.
\]

Thus the structural merge probability is **11/12 for even k and 1/12 for
odd k** under mu. The sidecar is the ordinary size of the first lift, not a
new dynamical assumption. Neither probability is a universal pathwise
guarantee, and the incoming automaton's success is sufficient for a common
future, not necessary. The targeted audit does not audit every census or
claim in that incoming note.

## 8. A stronger structure: summable discounted edge flow

The atomic measure translates coverage exactly, but proving its decay still
needs dynamics. Here is a second representation which makes positivity of
labelled receipts the central condition.

Let V be the positive odd integers other than 1. Kill the U edge when it
reaches 1. For positive weights w on V define the incoming operator

\[
(\mathcal K w)(m)=\sum_{n\in V:U(n)=m}w(n),\qquad m\in V.
\]

**P6 — discounted-flow equivalence.** Fix any rational 0<rho<1. Universal
positive Collatz is equivalent to the existence of strictly positive,
summable weights w on V satisfying

\[
\boxed{\mathcal K w\le\rho w.}\tag{11}
\]

If Collatz holds, such weights can be computable with a computable tail
bound. This converse is a representation theorem, not a construction that
avoids the unknown convergence assumption.

**Proof, sufficient direction.** Every surviving edge n->U(n) obeys
w(n)<=rho*w(U(n)), since its contribution is one term of the incoming sum.
An orbit remaining in V would have w(U^t(n))>=rho^(-t)w(n). This eventually
exceeds the finite total W=sum w, a contradiction. In particular the first
root time is at most 1+floor(log(W/w(n))/log(1/rho)). This also excludes
nonroot cycles, not only unbounded orbits.

**Proof, converse.** Assume all odd sources reach 1. Enumerate n_j=2j+1,
j>=1, and let tau_j be the actual odd root time. Give source n_j injection
a_j=2^-j rho^(tau_j). Send the discounted positive path receipt

\[
w=\sum_{j\ge1}a_j\sum_{t=0}^{\tau_j-1}\rho^{-t}[U^t(n_j)].\tag{12}
\]

Its total mass is at most rho/(1-rho), with tail over j>J at most
rho*2^-J/(1-rho). It is positive at every source, and exact cancellation of
the labelled path terms gives K w=rho(w-a)<=rho w. Finite heads can be
computed by the assumed terminating simulations and the displayed tail
bound makes w computable. This proves the equivalence. Initial powers of
two recover the full positive-integer statement. The identical argument
works for T_minus or U_minus killed on a specified known cycle set.

The representation is richer than an integer: **source description,
labelled flow, discount clock, incoming mass, and terminal certificate**.
Ordinary multiplication does not preserve (11). Weak receipt products may
erase exact labels, and common-future replacements preserve the old endpoint
defect but generally change discount lengths. Retain both sidecars.

**P7 — why the simplest weights cannot prove (11).** Taking w from (1)
already fails at n=3->5: w_3/w_5=4. More strongly, for every t>=1,

\[
n_t=2^{t+1}-1,\quad U^t(n_t)=m_t=2\,3^t-1,
\quad \frac{w_{n_t}}{w_{m_t}}
=4^{\ell(m_t)-\ell(n_t)}\longrightarrow\infty.\tag{13}
\]

All t intermediate valuations are one and avoid the root. Thus no fixed
block length makes this source measure a strictly contracting incoming
flow. Nor can a reweighting bounded above and below by fixed multiples of
these weights have strict one-step contraction: iterating that inequality
would contradict (13). A repair must permit an unbounded address-dependent
correction or keep a separate climb-time coordinate. This is a proved
obstruction to a particular norm, not to all summable flows.

**P8 — a credit correction that pays the climb, and its exact refuel bill.**
The next candidate was

\[
h(n)=v_2(n+1),\qquad \widetilde w(n)=w_n\,8^{-h(n)}.
\]

For every n=3 mod4, h(U(n))=h(n)-1 and the bit length increases by at most
one. Therefore w_n/w_(U(n))<=4 and
tilde_w(n)/tilde_w(U(n))<=1/2. This proves strict edge payment on **every
valuation-one climb**, with a positive summable explicit weight. It is a
successful local repair of (13).

But for every odd h>=3 the actual valuation-two edge

\[
n_h=(2^{h+2}-5)/3\longrightarrow m_h=2^h-1
\]

has h(n_h)=1, h(m_h)=h, and

\[
\widetilde w(n_h)/\widetilde w(m_h)=8^{h-1}/4.
\tag{14}
\]

Thus a single halving event can introduce arbitrarily much new climb fuel.
The first edges are 9->7, 41->31 and 169->127; the family recurs by
n->4n+5 and m->4m+3, preserving 3n+1=4m. This identifies the next missing
coordinate: **the cost of replenishing the climb register**, not merely
the register's current value. Suppression must propagate to the actual
incoming sources as well. Infinite lookahead defined using unknown root
times would make that repair circular; a finite symbolic family of incoming
inequalities with a summability bound is the concrete remaining target.

## 9. Recovering the useful part of the Vitali analogy

The historical [Vitali-wall reflection](../../07-reflections/the-vitali-wall-measure-vs-set-and-why-the-core-needs-construction-not-measure-s551o.md)
contains two overstatements now explicitly corrected: positive measure
does certify nonemptiness, and the AP speeds 1,...,n-1 have closed lonely
set {a/n:gcd(a,n)=1}, not all n polygon vertices. To see the latter, the n
points 0,t,...,(n-1)t must have all pairwise circular separations at least
1/n. Circle packing forces equal spacing, and distinctness forces t to be
a primitive n-th fraction. At n=6 only 1/6 and 5/6 qualify.

The survivor sets here are Borel; the integer target is countable. There is
no actual nonmeasurable set or choice obstruction. The useful transfer is
**an observable can forget exactly the information the target needs**.
A fixed atomless measure forgets isolated integer exceptions; our atomic
measure retains them. It does not thereby prove that they are absent.

The old [THM-168 lambda note](../../01-canon/theorems/THM-168-lambda-completeness.md)
has a genuinely finite version. Let lambda_uv count cyclic triples
containing u and v. Exhaustive deterministic search finds labelled
seven-vertex tournaments with arc masks 82468 and 84707, obtained by
reversing all arcs on {0,1,2,3}, with identical **labelled** lambda vectors
but Hamiltonian-path counts 141 and 143. Both counts are checked by subset
DP and all 7! permutations. The JSON specifies the arc encoding completely.

The correct phrase is “H is not measurable with respect to the finite
sigma-algebra generated by lambda”: H is not constant on a lambda fibre.
All subsets of this finite discrete tournament space are Borel. There is
no Lebesgue-versus-Borel hierarchy or nonconstructive tournament choice.
Also, the directed triangle has three Hamiltonian paths, not two; the
old [THM-169 sketch](../../01-canon/theorems/THM-169-vitali-atom-characterization.md)
must not be treated as a proved general classification on the basis of
its sampled census and that erroneous step.

## 10. Typed transfers and the next helper questions

| Source -> target | Map and preserved predicate | Lost information / necessary sidecar | Decisive test |
|---|---|---|---|
| AMM12592 spines -> source codes | Stop on a unary header; retain residual spine mass | Fair-coin null ray versus a specified arithmetic atom | All-odd positive/negative word, equation (8) |
| Word counts -> Collatz coverage | Attach canonical r and its atomic cylinder weight | Coefficient exit versus actual endpoint order | Exact exceptional lift k>(s-r)/(Q-P) |
| Two-sheet automaton -> probabilistic rule | Keep the same quarter-child guard | Distribution of the first lift t | 1/2 versus 11/12 or 1/12 |
| Weak receipts -> discounted flow | Keep labelled positive edges and a time factor | Prime cancellation loses labels; path replacement changes clocks | Check every incoming inequality (11) |
| Tournament lambda -> orbit statistics | Quotient objects by a finite observer | Fibre variation in the target | Masks 82468/84707 and sign-reflected prefixes |

Four sharpened helper questions:

1. **Address-adapted flow.** Can a computable correction to (1), organized
   by climb length and actual carry, satisfy (11) with a summability proof?
   Any correction uniformly comparable to (1) is ruled out by (13). Begin
   with an unbounded but summable climb register; test complete cylinders
   and the negative-cycle controls before asserting a global inequality.
2. **Defect-to-mass payment.** Can universally available weak receipts be
   repaired into positive discounted flow, charging every fusion defect to
   a summable bank of already rooted sources? Endpoint cancellation without
   the discount clock is insufficient. The inherited supported-replacement
   rule is the safe local move; its two path lengths must now be retained.
3. **Low-height survivor exclusion.** Can the height histogram in (7) be
   bounded uniformly, rather than just its total count? A bound on the
   contribution of each ordinary-height shell is the missing arithmetic
   input. The 27 plateau and minus-cycle mass (9) are mandatory controls.
4. **A weaker total machine with a strengthening budget.** The weak-receipt
   machine already halts with mass one. Can its unresolved *repair* mass
   have an effective tail bound independent of brute-force root search?
   A finite table of small tails is insufficient; the desired bound must
   tend to zero and preserve actual source identity.

The reusable move is to **change the observable while keeping the target's
quantifiers**. It succeeds here in distinguishing the signs and in exposing
the right positivity condition. It has not yet supplied a universal decay
law or a flow satisfying that condition.

## 11. Reproduction and limits

```bash
python3 04-computation/experiments/collatz_effective_prefix_mass_20261005.py > /tmp/prefix-mass.json
python3 -O 04-computation/experiments/collatz_effective_prefix_mass_20261005.py > /tmp/prefix-mass-O.json
cmp /tmp/prefix-mass.json /tmp/prefix-mass-O.json
```

All decisions use integers or Fraction; decimals only display exact masses.
Controls include independent mixture versus cylinder sums, binary refinement,
direct parity iteration for both signs through depth 14, exact word-tree
enumeration through 26, direct root-time replays, all odd sources below
2^17 for both target sets, complete t mod 256 sibling tests for k=4,...,24,
exact q=0,1,2 partition sums, and independent tournament path counts.
No finite computation proves (4), (6), or (11); their proofs are above.
No finite mass estimate is promoted to universal Collatz.
