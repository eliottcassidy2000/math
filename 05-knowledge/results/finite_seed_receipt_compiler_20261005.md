# Finite seed obligations and reusable infinite receipt families

2026-10-05. **PROVED:** sound finite conditional receipt DAGs; exact seed
discharge; and an infinite constructed source family for every admissible
controller program and fixed odd seed coprime to three. **FINITE-EXACT:**
the declared checker and symbolic/literal controls. **OPEN:** coverage of
arbitrary positive odd sources. The construction repackages inherited
inverse closure into a reusable obligation compiler; it does not claim a
new Collatz basin or a proof of universal convergence.

Artifacts: [program](../../04-computation/experiments/finite_seed_receipt_compiler_20261005.py)
and [saved output](finite_seed_receipt_compiler_20261005.out).

## 1. Inheritance and the logical object

The closest mechanisms are the checked common-future transport and grounded
closure in [fair frontier extension](fair_frontier_extension_20261004.md),
the exact first-hit inverse AST in
[inverse ray addresses](inverse_ray_ternary_addresses_20261004.md), and the
native five-letter controller in
[translation receipts](translation_phase_decoder_20261005.md). Its payment
constants come from [paid guard budget](paid_guard_budget_20261005.md).
The canonical hostile is an ungrounded cycle of mutually equivalent
convergence claims. The corrected near miss is mistaking a checked affine
relation for a legal actual word, or a conditional receipt for a completed
root certificate. The least-used sidecar is the **set of named, still-open
constant seed obligations**.

The concept board is **finite constants / quantified schemas / exact
receipts / source identity / well-founded dependency / grounded ROOT**.
Write \(U(n)=(3n+1)/2^{v_2(3n+1)}\) on positive odd integers, with the
certificate convention that a word stops at its first visit to one.
\(\operatorname{Root}(n)\) means such a finite first-hit word exists.

A checked receipt

\[
U_u(n)=U_v(d)
\]

gives \(\operatorname{Root}(d)\Rightarrow\operatorname{Root}(n)\).
The converse is also true, but the DAG records the chosen proof direction.
Both words are independently replayed, with exact valuations and no
root padding. An actual transition is the special case \(v=\varnothing\).
A `smaller` receipt additionally requires \(d<n\), making it admissible
as a natural-number induction step. A general finite join does not require
that inequality; it must instead participate in an acyclic derivation
with explicit leaves.

The checker has four types of inference nodes plus the checked ROOT leaf:

- `seed(s)` has conclusion \(\operatorname{Root}(s)\) and obligation
  exactly \(\{s\}\), where \(s\) belongs to a supplied finite set.
- `root(1)` has no obligation.
- `actual`, `join`, and `smaller` inherit their dependency's obligations
  after checking the corresponding exact receipt.

A bundle of conclusions takes the union of its obligation sets, so shared
seeds and subproofs are stored once. A depth-first active stack rejects
circular dependencies before they can be mistaken for an induction proof.
This is structural induction on a finite DAG. If a universally quantified
rule were supplied instead, its domain and decreasing rank would require
a separate proof; checking a finite list of instances would not suffice.

The example requests \(7,9,11\). The receipts are

\[
7\xrightarrow{1123}5\xleftarrow{1}3,\qquad
9\xrightarrow{2}7,\qquad11\xrightarrow{123}5.
\]

Their exact obligation set is \(\{3,5\}\). The separately supplied words
\(3\xrightarrow{14}1\) and \(5\xrightarrow{4}1\) discharge all three.
The reverse pair \(\operatorname{Root}(7)\Leftarrow\operatorname{Root}(3)
\Leftarrow\operatorname{Root}(7)\), with no seed or ROOT leaf, is rejected.
Merely declaring \(\operatorname{Root}(7)\) as an assumption leaves that
obligation visible; it does not discharge itself.

After discharge the checker cuts the supplied child first-hit word at its
checked prefix and prepends the source receipt. Determinism makes the
prefix agreement mandatory. It replays the completed word and independently
constructs the inherited inverse AST, retaining the original source.

## 2. A finite assumption supporting an infinite family

Use the inherited alphabet and its native source guards:

| Letter | Affine dependency | Source guard | Conservative credit |
|---|---|---|---:|
| H | \((729n+669)/1024\) | \(155\bmod2048\) | 2 |
| G | \((9n+5)/8\), actual word 12 | \(11\bmod16\) | -1 |
| A | \((81n+85)/128\), actual word 1213 | \(187\bmod256\) | 1 |
| B | \((81n+73)/128\), actual word 1123 | \(7\bmod256\) | 1 |
| L | \((9n-3)/16\) | \(219\bmod256\) | 4 |

This uses one credit for both A and B, including B's exceptional source
seven. It does not silently apply the stronger source-sensitive two-credit
rule there. For a marked word \(w\), require:

1. it is nonempty, avoids the pair LB, and does not end in L;
2. every prefix has nonnegative sum of the displayed credits.

These requirements define an unbounded language, not a finite list of
checked programs. Let

\[
F_w(n)=\frac{Pn+B}{Q},\qquad P=3^R,\quad Q=2^A.
\]

Here the coefficient \(B\) is unrelated to the letter named B. By the
inherited native-language theorem, a positive odd source is legal for the
whole word iff

\[
Pn+B\equiv Q\pmod{2Q}.
\]

Ending in L would instead require the stricter target condition
\(F_w(n)\equiv11\pmod{16}\). We retain the ending restriction explicitly.

Fix a positive odd **constant** \(s\) with \(3\nmid s\), without assuming
its orbit has already been checked. We construct sources conditional on
the single proposition \(\operatorname{Root}(s)\).
Choose \(a\ge4\) satisfying

\[
2^a s\equiv1+3B Q^{-1}\pmod{3^{R+1}}. \tag{1}
\]

The right side is one modulo three. The parity of \(a\) is therefore even
for \(s\equiv1\pmod3\) and odd for \(s\equiv2\pmod3\). After factoring
out that parity, powers of four give exactly one exponent class modulo
\(3^R\). Consequently (1) has exactly one class modulo \(2\cdot3^R\).

For completeness, \(v_3(4^{3^j}-1)=j+1\): the base case is immediate,
and cubing \(1+3^{j+1}u\), with \(3\nmid u\), raises the valuation by
exactly one. Each existing exponent residue thus has exactly one of its
three lifts meeting the next ternary digit. This proves existence,
uniqueness, and the finite digit-lifting algorithm used by the compiler.

Take a sufficiently large representative \(a_0\) and put

\[
a_t=a_0+2\cdot3^R t,\qquad
y_t=\frac{s2^{a_t}-1}{3},\qquad
n_t=\frac{Qs2^{a_t}-Q-3B}{3P},\qquad t\ge0. \tag{2}
\]

The implementation uses the explicit sufficient cutoff
\(a_0\ge\max(4,\operatorname{bitlength}(\max(1,Q+3B)))\), within the
required residue class. It is a sufficient cutoff, not a claimed least
positive exponent.

**Theorem.** For every admissible \(w\), every such seed \(s\), and every
\(t\ge0\), (2) is a distinct positive odd integral source with an exact
first-hit-safe actual receipt to \(s\). Its source-to-seed word has
\(R+1\) odd edges and total halving cost \(A+a_t\).
Its controller endpoint satisfies \(s<y_t<n_t\).
Thus every member is conditional on the same single constant seed
obligation, and one supplied ROOT certificate for \(s\) completes the
entire infinite family.

**Proof.** Equation (1) gives exact divisibility by \(3P\) in (2).
The source numerator is odd because \(Q\) is even and the reduced carry
\(B\) is odd. The cutoff gives positivity. The terminal value \(y_t\)
is positive odd, exceeds \(s\), and has exact edge
\(y_t\xrightarrow{a_t}s\). Direct substitution gives
\(F_w(n_t)=y_t\), so the whole native guard holds.

Conservative funding makes the inherited potential
\(E(x,k)=(x+5)(9/8)^k\) decrease strictly at funding actions and remain
unchanged at G. Every prefix has \(k\ge0\), so the complete nonempty
word gives \(y_t<n_t\), measured against the original source. The
convergence receipt itself does not need this inequality; it additionally
supplies a legitimate smaller dependency.

To obtain the actual word, start with the known terminal word \((a_t)\)
and read \(w\) backwards. G/A/B prepend their actual words. H changes
\((b,\text{tail})\) to
\((1,2,1,1,1,2,b+2,\text{tail})\). L changes
\((1,2,b,\text{tail})\) to
\((1,2,1,1,b+2,\text{tail})\).
These are the inherited common-future substitutions. Ending non-L gives
enough terminal data without querying the still-assumed seed's orbit;
the no-LB native rule retains each L prefix. More explicitly, prepending
A or G starts 12; prepending B starts 11 and cannot follow L; H and L
substitutions start 12. The single initial terminal word is never
consumed by a terminal L.

The resulting word has length \(R+1\), cost \(A+a_t\), and affine carry
\(Q+3B\). Its sole potentially enormous exponent stays last and is at
least \(a_t\). Equation (2) is precisely its cleared affine equation
with odd endpoint \(s\). Final oddness forces each prescribed valuation
exactly: after removing a first exponent, the remaining carry is odd,
so the first divided value is odd, and induction applies. Positivity is
preserved at each actual edge.

For \(s>1\), an earlier visit to one is impossible because the ordinary
odd map then stays at one. For \(s=1\), an earlier visit would force
all remaining valuations to be two, contradicting the last valuation
being at least four. Thus the receipt contains no root padding. Appending
a supplied first-hit seed word proves the result. Increasing \(t\)
strictly increases \(n_t\), giving infinitely many sources. \(\square\)

This proof is uniform in program length, proved by structural induction
through exact substitutions. A finite checked seed set can be used with
all programs and all parameters; no finite test is being promoted into
that induction theorem.

## 3. A concrete discharged family

Choose \(w=LG\). It is funded and legal, with
\((P,Q,B)=(81,128,53)\). Choose the single constant seed seven.
Equation (1) gives

\[
a_t=46+162t,\qquad
n_t=\frac{896\cdot2^{46+162t}-287}{243}.
\]

For every \(t\ge0\), the actual source word

\[
(1,2,1,1,a_t+2)
\]

reaches seven. The supplied, independently replayed seed word
\((1,1,2,3,4)\) then reaches one. The family therefore has the explicit
first-hit certificate

\[
(1,2,1,1,48+162t,1,1,2,3,4).
\]

There are ten odd edges and \(74+162t\) ordinary steps. The source remains
in the exact native cylinder \(219\pmod{256}\). This is inherited inverse
construction and certificate reuse, not evidence that every number in
that cylinder belongs to this thin exponential family.

The same compiler gives \(a=116+162t\) for seed one and
\(a=93+162t\) for seed five. Seed one has no outstanding obligation;
the checked ROOT object already supplies it.

## 4. What is stored, proved, and still required

The reusable record is \((w,s,t)\), together with one checked seed
certificate shared by every family member. The program's exact carry
and exponent phase are reconstructed. The source may be vastly larger
than its record. Reading it modulo \(M\) uses

\[
n_t\bmod M=
\frac{(Qs2^{a_t}-Q-3B)\bmod(3PM)}{3P},
\]

which retains the precision needed before division. This is exact modular
exponentiation, not an approximation or a source-residue substitution.
The inverse AST stores the large valuation as an integer label; it need
not expand the source or perform that many individual halvings.

The production APIs never search for an unknown seed's route. They either
return an explicit obligation, or consume and check a supplied first-hit
word for that exact seed. A wrong source, padded root, absent proof, or
conditional proof offered as a completed seed is not accepted.

A finite list \(\operatorname{Root}(s_1),\ldots,\operatorname{Root}(s_k)\)
is different from a finite list of sentences such as
\(\forall t\,\operatorname{Root}(f_i(t))\). The latter have infinitely
many instances and need their own uniform arguments. The theorem above
supplies one such argument for the stated constructed families. It does
not turn an arbitrary family label or an unproved domain-cover statement
into a discharged finite seed.

In particular, no assertion here says every positive odd integer enters
one of these families. A universal proof by finite seed reduction would
still require a proved all-source reduction with a well-founded rank,
or another complete rule that terminates at the finite seed set.

## 5. Reproducible scope and adversarial controls

Run from the repository root:

```text
python -B 04-computation/experiments/finite_seed_receipt_compiler_20261005.py
python -B -O 04-computation/experiments/finite_seed_receipt_compiler_20261005.py
```

The exact universe is the 88 admissible conservative-funded programs of
length one through three, five seeds \(1,5,7,11,13\), and parameters zero
and one: 880 completed compact ASTs. The independent modular source
readers are compared at binary and ternary depths \(1,2,5,11,31\), making
8,800 comparisons. All 203 members under the declared 20,000-bit expansion
budget also receive literal receipt replay, original-source payment
checks, and expanded inverse-AST identity checks. Larger members remain
symbolic; the largest terminal exponent tested is 1,439,843,110, with at
most nineteen source-to-seed odd edges.

Six short seed words are explicitly supplied and checked; no convergence
search is used to manufacture these inputs. Nineteen hostile controls
cover a genuine common-future cycle, missing and wrong seed proofs,
undeclared assumptions, corrupted receipts, numeric type aliases,
root-padding, invalid LB/terminal-L/unfunded programs, seeds divisible
by three, parameter types, and expansion limits.

The useful next operation is to match an actual supplied source to a
proved family guard or receipt while preserving its identity, then expose
the exact residual seed or induction obligation. Merely increasing the
number of checked constants does not resolve an unproved universal
domain-cover statement.
