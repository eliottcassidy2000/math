# Beyond 27: a certified growth family and the exact poset/arithmetic interface

**Status: PROVED, elementary + INDEPENDENTLY AUDITED; FINITE-EXACT controls.
Collatz remains OPEN.** These constructions do not depend on accepting either
of the two newly submitted poset papers. They specify precisely where those
papers' probability laws could, and could not, enter a Collatz argument.

Use the shortcut map `T(n)=n/2` for even n, `(3n+1)/2` for odd n. Stopping
times below use this convention, not the convention counting `3n+1` and the
following halving separately. We stop the count at the first visit to 1;
the algebraic map itself is not made absorbing there.

Canon: [THM-4502, certified 41-tail growth family](../../01-canon/theorems/THM-4502-collatz-41-tail-growth-family.md)
and [THM-4503, parity posets and height selection](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md).
Reproduce with `python -X utf8 04-computation/experiments/crossroads_poset_20260926_bridge.py`;
[script](../../04-computation/experiments/crossroads_poset_20260926_bridge.py),
[raw output](crossroads_poset_20260926_bridge.out).

## 1. Inheritance and what each small object actually carries

The closest proved mechanisms are the actual CRT orbit families of
[THM-4501, recursive motif families](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md)
and the exact finite-word carry of
[THM-4476, thin divergent orbits](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md).
The canonical hostile is the all-odd address: every finite prefix has positive
integer realizations, but their compatible least residues tend to the negative
2-adic fixed point -1. The corrected near miss is replacing the source by a
new representative whenever a word is rearranged. The useful underused
coordinate is the arithmetic carry attached to each linear extension.

The live board contains: a known convergent tail; a ternary phase clock;
two chains of parity events; the carry/residue of an extension; a height-selected
law; and finite versus escaping boundary witnesses. No cosmetic tournament is
introduced: the intrinsic relation is the partial order of events.

| Source and target | Map / preserved fact | Destroyed information / required sidecar |
|---|---|---|
| known tail 41 to new starts | inverse affine construction / actual convergence | occurrence frequency and initial growth must be proved separately |
| congruence modulo powers of 3 to family | unique phase lift / integrality | small phases are not random ternary digits |
| positive-slope words to two-chain extensions | successive odd/even event labels / every prefix slope | source height unless the exact carry is kept |
| adjacent extensions to each other | swap incomparable events / endpoint multiplier | the source changes at an exactly known binary digit |
| all residues to starts below X | arithmetic selection / actual positive integers | uniformity, convex support, and XYZ inequalities |
| finite family to dyadic closure | convergence of residues / local congruences | ordinary counting density and positive-integer realization |

## 2. A long-growth family containing 27, with termination built in

For `k>=1` put `p_k=2*3^(k-1)`. There is a unique phase `r_k` in `[0,p_k)`
satisfying

    41*2^(r_k) = -1 mod 3^k.                         (2.1)

For every `t>=0`, set `r=r_k+p_k t` and

    u=(41*2^r+1)/3^k,
    n_(k,t)=2^k u-1.                                (2.2)

These are positive integers. For `0<=j<=k`, their actual trajectory obeys

    T^j(n_(k,t))=3^j 2^(k-j)u-1.                    (2.3)

Indeed the right side is odd for j<k and the next shortcut step gives the
next displayed value. Consequently all first k iterates lie strictly above
the start, and

    T^k(n_(k,t))=41*2^r -> ... -> 41 -> ... -> 1.     (2.4)

The first arrow uses exactly r halvings. The tail from 41 has 69 shortcut
steps, checked directly in the retained program. No earlier value in (2.3)
or the halving tail is 1. Thus the total stopping time is exactly `69+k+r`.

The phase exists because 2 has multiplicative order `2*3^(k-1)` modulo `3^k`.
For completeness, `v_3(4^m-1)=1+v_3(m)`: the factorization at a factor 3
increases valuation by exactly one, while an exponent prime to 3 does not
increase it. This proves the asserted order; the unit group has that same
size, so the powers exhaust its units. The phases satisfy the explicit
three-choice recursion

    r_(k+1)=r_k+p_k d_k,  d_k in {0,1,2},            (2.5)

where exactly one digit makes (2.1) true at the next level. They are
nondecreasing and tend to infinity: a bounded sequence would eventually
be constant, making the fixed positive integer `41*2^r+1` divisible by all
powers of 3. No distributional assertion about the digits d_k is needed.

| k | least phase r_k | least-phase start n_(k,0) | shortcut stopping time |
|---:|---:|---:|---:|
| 1 | 0 | 27 | 70 |
| 2 | 4 | 291 | 75 |
| 3 | 10 | 12439 | 82 |
| 4 | 28 | 2173995791 | 101 |
| 5 | 28 | 1449330527 | 102 |
| 6 | 28 | 966220351 | 103 |
| 7 | 1000 | a 1002-bit integer | 1076 |
| 8 | 2458 | a 2459-bit integer | 2535 |
| 9 | 6832 | a 6833-bit integer | 6910 |

The starts need not increase with k: a phase plateau makes them decrease.
They are not claimed to be stopping-time record holders or atypically slow
relative to their bit length. What is proved is arbitrarily long initial
growth, a precise recursive family, and a certified common tail.

**The 27 exception.** For r>0, u is odd and `v_2(n+1)=k`, so the initial
odd run has exactly length k. For r=0, divisibility forces k=1, and gives
27, which actually has two initial odd steps. Formula (2.3) promises k
steps, not that no further odd step occurs. This exception must not be
discarded when classifying the family.

For each fixed k put `B_k=2^(p_k)`. The within-family recursion is

    n_(k,t+1)=B_k n_(k,t)+(B_k-1)(1-(2/3)^k).        (2.6)

The intercept is an integer since `3^k` divides `B_k-1`. At k=1 this is
`n -> 4n+1`, the old 27 comb. At k=2 it becomes `n -> 64n+35`. Thus this
construction extends that comb to arbitrarily long initial growth.

The integrality, prescribed growth, and landing construction works with
any fixed positive odd tail a coprime to 3, replacing 41 by a. Convergence
then requires a separately certified tail from a to 1. The exceptional
r=0 cases depend on v_3(a+1); the special 41 classification is not imported.
Nor is the first-hit time formula imported for a=1, where the construction
can start at 1 itself. This is a family-building construction, not a proof
that all integers belong to one of these families.

## 3. The exact occurrence order and a zero-dimensional closure

Let H be the union of (2.2) over k>=1,t>=0. The parametrization is injective.
For r>0, the valuation of n+1 determines k, then (2.2) determines r. The
exceptional n=27 cannot duplicate another pair: its valuation would force
k=2, which would require `41*2^r=62`.

**Counting theorem.** As X tends to infinity,

    #(H intersect [1,X]) = Theta(log X).             (3.1)

Proof: set `L=log_2(X+1)`. Since `n+1=2^k u` with u>=1, `k<=L`. Also

    41*2^r+1=3^k u <= (3/2)^k(X+1),
    r <= L log_2 3-log_2 41 = O(L).

At a fixed k, possible r form a progression of spacing p_k, so there are
at most `R/p_k+1` below any common nonnegative cap R=O(L). Summing uses
`sum_(k>=1) 1/p_k=3/4` and at most L values of k, giving O(L). Conversely
k=1 supplies `(82*4^t-1)/3`, giving Omega(log X). In particular H has
natural density zero. No exact leading constant in (3.1) is asserted.

The program counts the full union exactly by integer logarithms and checks
it against direct trajectories for every odd source at most 100000. Those
sources are `27,109,291,437,1749,6997,12439,18659,27989`.

**Closure theorem.** In the 2-adic integers,

    closure(H)=H union {q_k=(2/3)^k-1:k>=1} union {-1}. (3.2)

For fixed k and t tending to infinity, the formula

    n_(k,t)=q_k+41*2^(k+r)/3^k                       (3.3)

gives the limit q_k. When k tends to infinity all starts are congruent to
-1 modulo `2^k`, so they tend to -1, as do the q_k. Conversely any sequence
has either a subsequence with k tending to infinity or one with k fixed;
in the latter case bounded r stabilizes at a member of H, while unbounded
r gives q_k. This classifies every possible limit.

The closure has **both upper box and Hausdorff dimension zero**, and Haar
measure zero. Here is the stronger quantitative cover, not just a reliance
on countability. Modulo `2^J`, k>=J gives only -1. For each k<J, r>=J-k
gives q_k; the remaining progression contributes at most `J/p_k+1` values.
Thus fewer than `11J/4` dyadic cylinders cover the closure, with harmless
rounding at small J. A linear number of cylinders of radius `2^-J` gives
upper box dimension zero, hence Hausdorff dimension zero.

This separates three previously conflated notions. The family has exact
recursive structure; its occurrence rate is logarithmic; its dyadic closure
is zero-dimensional. The larger positive-slope language from THM-4495 has
a positive entropy dimension, but that is a different set. Neither type
of dimension determines convergence of an arbitrary prescribed integer.

## 4. Expanding parity words are linear extensions of a two-chain poset

Fix a odd letters and b even letters, a+b=ell. A word is called expanding
here if **every nonempty prefix** with i odd and j even letters satisfies

    3^i > 2^(i+j).                                  (4.1)

This condition implies actual first-ell-step no descent for every positive
integer realizing the word. It is not necessary for one small source:
the positive affine carry can compensate a contracting multiplier (source
1 with word 10 returns to 1 with multiplier 3/4; no example with source
greater than 1 is asserted here). For
the exact actual-no-descent carrier, see the integer note, section 4.

Define `m_j=min{i>=0:3^i>2^(i+j)}` for j>=1. Equivalently
`m_j=floor(j log_(3/2)2)+1`; the inequalities can be decided with integers.
When b=0 the language has just the all-odd word. When b>0 and a<m_b it is
empty. In the remaining cases define P_(a,b) by

    A1<...<Aa, B1<...<Bb, and A_(m_j)<B_j for 1<=j<=b. (4.2)

**Bijection.** Expanding words with these counts are exactly the linear
extensions of P_(a,b), with an A read as 1 and a B as 0. At a zero event
B_j, at least m_j odd events must already have occurred; that is exactly
(4.1) there. An odd event increases the logarithmic multiplier, so no
additional check is needed between zero events. This proves both directions.

P has width at most two. Exponentially many extensions do not make its width
large. The newly attached large-width theorem therefore does not directly
apply to this encoding. The older bounded-width/large-range result is a
closer comparison, provided the **uniform full extension law** is retained.
For example, with b=1 and a large, B1 has a-2 incomparable elements.

The exact program compares independent word enumeration with a minimal-element
poset recursion in all 119 fixed-count universes through ell=14. It also
checks the width-two 1/3 balance bound in every nontrivial finite universe;
that finite check is not offered as a new proof of the known general theorem.

## 5. The arithmetic information carried by an extension

For a word w of length ell and a odd letters,

    T_w(n)=(3^a n+C_w)/2^ell,
    C_w=sum_(j:w_j=1) 2^j 3^(number of ones after j),
    r_w=-C_w*3^(-a) mod 2^ell.                      (5.1)

Indices j start at zero. The last quantity is the unique starting residue
that realizes w. Conversely, following the parity of an integer in that
residue class gives exactly w; this follows inductively by extending the
modulus from `2^j` to `2^(j+1)`, where precisely one lift has the next
required parity. This is stronger than endpoint integrality alone.

An allowed adjacent swap `10 -> 01` at position j, with s odd letters
after the pair, changes the carry and starting residue by

    C_new-C_old=2^j 3^s,
    r_new-r_old=-2^j 3^(s-a) mod 2^ell.              (5.2)

Thus the difference of the two residues has exact 2-adic valuation j.
Endpoint multiplier and event counts are preserved, but the source changes
at binary digit j. This is an exact obstruction to treating extension
rearrangements as alternative executions of the same fixed integer.

For a word with one zero after s initial ones, (5.1) simplifies to
`r_s=-1-(2/3)^s mod 2^ell`. Its ordering in ordinary height need not agree
with the insertion position of the zero. The height-cut counterexample is
already visible with six events:

| word | least residue modulo 64 |
|---|---:|
| 110111 | 27 |
| 111011 | 39 |
| 111101 | 47 |
| 111110 | 31 |

These are all four extensions of `A1<...<A5, A2<B1`. Sampling starts at
most 31 selects the first and last insertion positions. No poset on these
same labels has exactly those two extensions: all relations shared by
them already permit the other two. Adding auxiliary variables is not ruled
out, but then its marginal law and arithmetic decoding require proofs.

The [width audit, section 4](crossroads_poset_20260926_width.md) computes a
stronger analytic witness. Take the order statistics of six independent
Uniform(0,1) draws, independently of the selected rank order; attach them
to that order, let F be the resulting coordinates, and set
Z=7F. Then

    Cov(Z_A2-Z_A3, Z_B1-Z_A2)=1/4,

whereas the full four-extension law gives `-5/32`. The selected chambers
also have nonconvex support. Thus neither nonpositive XYZ covariance nor
the uniform convex order-polytope representation survives this arithmetic
cut. This refutes a transfer, not the original poset inequality.

## 6. A quantitative repair where there are many full residue periods

There is a useful surviving theorem. Choose any m nonempty distinct
residue classes modulo N. For a uniform integer up to X conditioned on
these classes, write `X=qN+s`, `0<=s<N`, and let t of the classes occur
among `1,...,s`. The class weights are q+1 on t classes and q on the rest.
Their total variation from the uniform class law is exactly

    TV=t(m-t)/(m(qm+t)) <= 1/(4q),  q>=1.            (6.1)

This follows by summing the positive deviations from 1/m; it is checked
against direct integer sampling in the independent width probe. For a
fixed-count parity stratum N=`2^ell`, uniform extension balance at least
b therefore gives height-selected balance at least `b-TV` for the same
pair. At whole periods the laws agree exactly.

So the poset machinery can transfer quantitatively when `2^ell=o(X)`.
At source-sensitive horizons ell comparable to or larger than `log_2 X`,
the error need not vanish. For the six-event hostile q=0,t=2,m=4, TV=1/2.
For one fixed source and ell tending to infinity this theorem supplies
no equidistribution. This boundary is arithmetic, not a numerical weakness
in the poset papers' balance constants.

## 7. Consequences and precise next questions

The family supplies positive controls with arbitrarily long growth and
exact density/closure. Any proposed rule that excludes them solely by a
long expanding prefix is false. The bridge supplies a genuine small poset
encoding and the precise arithmetic sidecar that rearrangements lose.
The height-cut example supplies a minimal hostile in the probed universe,
and (6.1) supplies the strongest elementary uniform-law repair.

Three next constructions are well-posed:

1. Derive an inequality for the actual weighted extension law with a carry-
   dependent defect that remains useful when q=0; a TV estimate alone is
   then too large. Test it first on the four-word 27/31 example.
2. Randomize the order in which certificates about one fixed source are
   revealed, retaining its actual parity/carry, instead of randomizing its
   execution word. Prove that a balanced comparison has a descent-relevant
   consequence before invoking a balancing theorem.
3. Use a summable full-support measure on positive starting integers and
   exact height-capped cylinders. Unlike Haar density, vanishing mass for
   this law would exclude each positive integer separately. The integer
   note states the precise generating-function reduction and its missing
   decay estimate; no such estimate is proved here.

Audits: root derived the family and two-chain bridge; automata independently
derived the bridge and checked the family; geometry and flow independently
checked the counting and closure arguments, including the exceptional 27.
Normal and optimized Python replay, the actual trajectory checks, and the
independent uniform-poset/height-cut census are recorded in the session audit.
