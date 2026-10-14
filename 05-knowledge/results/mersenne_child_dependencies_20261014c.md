# Mersenne child dependencies: retained roots and the boundary of exponent recursion

**Status:** PROVED exact fixed-word test, guarded reset-pair transport and
certificate criteria; FINITE-EXACT declared word universes and retained
ROOT transports. Universal Collatz convergence, coverage of every exponent,
and ROOT for arbitrary members of the inherited phases remain OPEN.

Companions: [checker](../../04-computation/experiments/mersenne_child_dependencies_20261014c.py)
and [saved JSON](mersenne_child_dependencies_20261014c.json). No orbit search
is performed by the production APIs or the default experiment.

## 1. What was already completed, and what was only a dependency

Write $M_E=2^E-1$ and $U(n)=(3n+1)/2^{v_2(3n+1)}$ on positive odd
integers. A ROOT word stops at its first visit to 1. A paid common-future
receipt connects an original source to a smaller child using two actual
valuation prefixes with the same endpoint. It is not a ROOT word by itself.

The closest mechanisms are the
[two-twos grounding packet](collatz_two_twos_grounding_20261007c.md),
[reset-pair theorem, section 2](mersenne_switch_parity_f21_compression_20261006.md),
[context cancellation](collatz_context_cancellation_20261007d.md), and
[marked completion consumer](collatz_marked_completion_20261007c.md).
The canonical hostile is the legal inverse dependency $31\rightsquigarrow27$,
whose forward evidence is $27\xrightarrow{(1,2)}31$: it does not independently
ground either port. The corrected near miss is to read a type, a native
phase, or a positive-size drop as a completed ROOT proof. The useful sidecar
is the actual child type together with its source-owned certificate.

The live board is exponent phase, ordinary odd child, reset pairing,
ordered carry, certificate profile, and unresolved terminal. The anchor
is the still-unproved child obligations; the niche is a complete test of
whether a fixed inverse word even stays within Mersenne numbers; the
wildcard is transport of already retained ROOT evidence.

The following four children were **already grounded** before this turn:

| E | First-hit odd rank | Total halving cost | Already grounded parent |
|---:|---:|---:|---:|
| 6125 | 28762 | 51712 | 6129 |
| 18269 | 90248 | 161309 | 18273 |
| 46941 | 226498 | 405932 | 46945 |
| 98477 | 475425 | 852008 | 98481 |

Their compressed positive valuation words are retained in the inherited
[source-bound data file](collatz_two_twos_grounding_20261007c.cert.json),
whose checker performs forward and reverse authentication. This packet
replays the first two seeds; the other two rows are explicitly inherited
metadata, not a claim of fresh full replays in this run.

In contrast, the [eight-bit family](collatz_eight_bit_completion_20261007c.md)
and its later continuations retain smaller-child obligations. In particular
the actual residual parameter 23 has
$E=99708993705$ and the retained alternatives
$D\in\{3,4,7,8,\ldots,28\}$ in the
[source-specific exit](collatz_residual_t23_20261007f.md).
Neither an orbit of that giant Mersenne nor its ROOT certificate is supplied
here. The five additional transports below do not change that status.

## 2. A fixed contracting inverse word has at most one Mersenne pair

For a nonempty positive valuation word $w$ of length $L$ and cost $A$, let

\[
 F_w(x)=\frac{Px+B}{Q},\qquad P=3^L,\quad Q=2^A.
\]

We consider the contracting inverse slope $Q<P$, with potential child
$h=(Qn-B)/P$. Its native integrality guard still belongs to the original
source $n$. Asking whether it maps $M_K$ to a smaller nonroot Mersenne
$M_F$, $2\le F<K$, is a stronger restriction than asking for an arbitrary
positive odd child.

**One-candidate theorem.** Put $R=Q+B-P$. If such a pair exists, then

\[
 \boxed{F=v_2(|R|),\qquad 2^{K+A}=P2^F+R.}                 \tag{1}
\]

If $R=0$ there is no pair. Otherwise (1) gives at most one candidate:
test that $P2^F+R$ is a positive power of two, recover $K$, require
$2\le F<K$, and authenticate the supplied forward word from $M_F$ to
$M_K$. These tests are sufficient as well as necessary.

**Proof.** The endpoint equation is exactly
$2^{K+A}-P2^F=R$. Since $K+A>F$ and $P$ is odd, the difference has
2-adic valuation $F$. If it vanishes, unique prime factorization would
force $3^L$ to be a power of two, impossible for $L>0$. Conversely the
tests retain the exact affine equality. A final positive odd endpoint
forces every intermediate valuation to be native: a nonintegral dyadic
denominator cannot disappear under a later odd numerator, and an even
intermediate makes the next step nonintegral. Intermediate values are
positive. An earlier visit to ROOT cannot be followed by an actual endpoint
$M_K>1$. The checker also literally replays the candidate. ∎

There is no hidden negative-$R$ family for genuine positive valuation
words. Appending a letter $a$ gives

\[
 R_{w,a}=3R_w+(2^a-2)Q_w.
\]

Starting at $R=0$, this proves $R\ge0$, with equality exactly for an
all-ones word. The absolute value in (1) is harmless algebraic notation;
it is not a claim that negative carries are native positive words.

**Finite-bank consequence.** Every finite bank of such fixed words has
at most finitely many Mersenne-to-smaller-Mersenne pairs. It cannot by
itself be a cofinal recursion on the exponent. Variable-length grammars,
two-sided common-future receipts, and ordinary odd children are outside
this obstruction.

There are stronger complete results for two inherited finite universes:

| Fixed word universe | Retained words | Pairs at any exponent sizes |
|---|---:|---:|
| All positive valuation words $L\le8$, $2^A<3^L$ | 953 | 0 |
| All nonempty words of generator length at most 8 in $G_1,G_5,G_{17}$ | 9840 | 0 |

The first is the complete universe in the
[bounded inverse cover](collatz_bounded_inverse_cover_20261007e.md).
Its finiteness follows from the displayed cost inequality; the script
enumerates all compositions of each permitted cost. The second uses the
authenticated inverse words and chronology of
[child normal forms](collatz_child_normal_forms_20261007.md).
The new conclusion concerns closure in the Mersenne subcategory, not a
new inverse family or a failure to pay arbitrary odd children. The
one-candidate proof makes these **all-exponent** verdicts for the declared
finite word banks. No extrapolation in the word length is made.

For example $G_5(31)=27$ is positive, native and smaller, but is not a
Mersenne number. More generally the inherited phase $K=5\pmod {18}$
has first child $(2^{K+3}-13)/9$ divisible by 3 and terminal for that
three-generator router. Leaving the exponent coordinate is a real type
change, not permission to replace the child by a convenient Mersenne.

## 3. Reset pairs give a lossless quotient with a ROOT boundary

The following is an explicit specialization of the inherited reset switch,
not a new general switching theorem. For even $E\ge4$, write $s=v_2(E)$.
Then the strict actual words

\[
 \begin{array}{c|c}
 M_E&1^{E-1},\ 3+s\\
 M_{E-1}&1^{E-2},\ 2,\ 1+s
 \end{array}
\]

meet at

\[
                 \frac{3^E-1}{2^{2+s}}.                         \tag{2}
\]

After the initial one-run the source is $2\cdot3^{E-1}-1$.
For even $E$, the elementary lifting identity
$v_2(3^E-1)=2+v_2(E)$ proves its final valuation. For the child,
$E-1$ is odd, so its first reset is 2, after which the same identity
gives the displayed next valuation. Initial runs grow, and their
preterminal states exceed 1. Both words have length $E$ and their costs
differ by one. Equation (2) proves a common-future equivalence, including
transport of a supplied ROOT word in either direction.

The $E=2$ exception must not be padded: $M_2=3$ has the direct word
$(1,4)$ to ROOT, while $M_1=1$ has the empty word. The formal child word
$(2,2)$ would already start at ROOT and is rejected.

Thus the pairs $\{M_{2j-1},M_{2j}\}$, $j\ge2$, may be treated as one
dependency class **when this exact lift is retained**. Payment still
points from the larger number to the smaller one. Bidirectional transport
of an independently supplied proof is a different logical operation.

## 4. A two-bit phase continuation and ten source-owned ROOT recipes

The inherited paid exit at

\[
 E=18269\pmod {2^{48}}
\]

has source head $u$ and partner head $v$ recorded in the grounding packet,
with $|u|=24$, $|v|=25$, costs 50 and 48, and $v$ beginning with 5.
Its child has even exponent $E-1$, with $v_2(E-1)=2$ throughout the phase.
Replace that first partner letter 5 by $(2,3)$:

\[
              v'=(2,3),v_2,\ldots,v_{25}.                       \tag{3}
\]

The complete carrier equality is

\[
 F_{v'}(y)=4F_u(x)+1\quad\text{when }x+1=9(y+1).
\]

The checker retains both its slope and carry. This is exactly composition
with (2), so it retains the **same original-source guard**, now paying
$M_E\rightsquigarrow M_{E-2}$. If the actual terminal valuation after $u$
is $c$, the words are $1^{E-1},u,c$ and $1^{E-3},v',c+2$.
All positive phase members satisfy the inherited proper-prefix growth
cutoff, and both paths are strict first-hit paths. This adds no source
guard mass. With the inherited parent entrance it also gives
$M_K\rightsquigarrow M_{K-6}$ on $K=18273\pmod {2^{48}}$;
arbitrary members still need the final child's ROOT proof.

At the retained point $E=18269$, the existing ROOT word transports first
to 18268 and then to **18267**. The direct (3) splice is checked to be the
identical valuation word. An independent reverse replay reconstructs
$M_{18267}$ from ROOT. No orbit discovery is used.

The public `named_root` consumer reconstructs ten exact source-owned words:

| Exponents | Odd rank | Cost minus exponent |
|---|---:|---:|
| 6125, 6126, 6129, 6130 | 28762 | 45587 |
| 18267, 18268, 18269, 18270, 18273, 18274 | 90248 | 143040 |

The seeds are 6125 and 18269. Their old parent receipts give 6129 and
18273; the old exit gives 18268. Reset-pair transport supplies the other
five source receipts. This describes the new packet's exports, not a
priority claim that those integers' convergence was previously unknown.
The JSON retains source, rank, cost, peak height and a hash of each decoded
word; the recipe retains its original compressed seed dependency.

## 5. A terminal profile is an obstruction to a specified join, not to ROOT

For a supplied strict ROOT word at $M_E$ let $r_E$ be its odd rank and
$A_E$ its halving cost. Consider the fixed affine bit-deletion ansatz

\[
 F_w(X)=F_v\!\left((X+1)/2^{E-F}-1\right),\qquad E>F.    \tag{4}
\]

There are strict actual prefixes satisfying (4) iff

\[
              (r_E,A_E-E)=(r_F,A_F-F).                         \tag{5}
\]

Comparing slopes in (4) forces equal word lengths and a cost difference
$E-F$. Appending the common strict ROOT suffix preserves these equalities,
proving necessity. Conversely use the complete supplied ROOT words:
their slopes agree by (5), and equality at the original source fixes the
constant term. This is the arbitrary-drop version of the inherited
[half-child criterion](half_child_extension_obstruction_20261004.md).

All ten receipts above retain their displayed profile. The four inherited
seed ranks are distinct, so no pair of those four seeds can be merged by
this strict affine ansatz at any word length. They still all have ROOT
proofs, and their orbits can join at ROOT with different clocks. The
profile does not obstruct those broader joins. The smallest hostile
$M_3=7$, $M_2=3$ has strict words $(1,1,2,3,4)$ and $(1,4)$ with unequal
profiles; longer search cannot repair (4) for that pair.

This is the relevant limit of the 561 simultaneous-prime analogy.
Several compatible obligations can be compiled together, but a proof
dependency cannot be grounded by cycling among its assumptions. Here the
retained ROOT seeds are independently authenticated; the dependency DAG
only transports their exact suffixes. A modular native phase or an
off-category terminal such as 27 supplies no replacement for those seeds.

## 6. Verification and the next obligation

Run `python -B -X utf8 04-computation/experiments/mersenne_child_dependencies_20261014c.py`
and repeat with `-O`. Public APIs validate exact integer, word and receipt
types before source comparison. Materialized Mersenne sources are capped
at 100000 bits. Symbolic phase controls include parameter $10^{30}$ and
do not expand those sources. ROOT padding, wrong phase, substituted child,
altered valuation and a valid proof at another source are rejected.

The declared controls comprise even exponents 4 through 128; the complete
953 and 9840 word universes; independent direct substitution at exponents
2 through 32; the ten streaming ROOT words; the direct/two-stage 18267
word equality and reverse replay; and malformed/boundary receipts.
The large check count mostly consists of individual valuation validations,
not millions of independently explored sources.

The research moves used here are the repository's **attack a proposed
bound before extending it** and **verify the consequence, not the model's
own assumptions** cards. The finite inverse banks were tested for actual
Mersenne return, rather than merely for integral smaller children; the
positive consequence was checked at the original source by ROOT replay.
No new general meta-pattern is proposed.

The next requirement is a genuinely authenticated complementary rule or
ROOT proof for the supplied unresolved child, possibly outside Mersenne
numbers. A finite bank that necessarily leaves this category cannot close
that obligation merely by relabelling its endpoint as a smaller exponent.

Normal and optimized execution agree on **8,522,665 exact checks** and
the saved JSON. LF-normalized SHA256:
`a8b2f5f81dab81ba42a5c6b336b924c36e9dea8d0ca51bedb93f99fbf9b87ade`.
