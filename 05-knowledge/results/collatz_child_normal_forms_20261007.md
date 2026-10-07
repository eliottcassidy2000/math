# Guarded child normal forms and their terminal frontier

**Status: PROVED elementary all-length algebra and guarded reduction;
FINITE-EXACT implementation controls. Supplied-paper results are accepted
premises for mechanism comparison only. Universal ROOT coverage is not claimed.**

The three inverse-cycle child maps admit a closed, lossless affine normal
form. Their guards are disjoint, so the current source itself selects the
next rule. Every positive source reaches a terminal frontier in finitely
many such reductions. That frontier is not automatically ROOT-grounded:
the least Mersenne example is `31 -> 27`.

## 1. Inheritance and the precise object

The closest mechanism is the corrected negative-cycle branch in
[THM-4594, maximal class-decided sieve](../../01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md):
its sign barrier applies to **unrefined binary classes**, and the added
ternary guards license the inverse maps below. The all-depth inverse-12
operation is already implemented in
[the collision grammar](collatz_collision_dp_20261007.md), in its section
“Arbitrary ternary depth after the binary collision.”
[Affine guarded lifts](collatz_affine_guarded_lifts_20261004.md) supplies
the broader distinction between formal affine group words and legal
integer paths. This package specializes to a contracting three-letter
alphabet; it does not replace that general inverse grammar.

The canonical hostile is a legal smaller child whose known forward route
returns to its parent: `27 --(1,2)--> 31`. This is a dependency, not a
proof of either orbit's ROOT fate. The corrected near miss is erasing the
ternary guard when transporting a negative-cycle inverse. The least-used
sidecar is the **ordered affine carry**, which here supports an exact word
decoder. The live board is: native ternary guards; ordered carries;
Mersenne exponent phases; positive height; terminal obligations; compatible
paper-inspired reductions.

Write `U(n)=(3n+1)/2^v2(3n+1)` on positive odd integers. The generators,
read chronologically as child operations, are

| Label | Child map | Exact positive-odd input guard | Actual word from child back to input |
|---|---|---|---|
| 1 | `(2n-1)/3` | `n=2 mod3` | `(1)` |
| 5 | `(8n-5)/9` | `n=4 mod9` | `(1,2)` |
| 17 | `(2048n-2363)/2187` | `n=2170 mod2187` | `(1,1,1,2,1,1,4)` |

In uniform notation let `(a,b,c)` be respectively `(1,1,1)`, `(3,2,5)`,
or `(11,7,2363)`. Then

\[
G(n)=\frac{2^a n-c}{3^b},\qquad
c=\ell(3^b-2^a),\quad \ell\in\{1,5,17\}.
\tag{1}
\]

Thus `-ell` is its fixed point and `G(n)+ell=(2^a/3^b)(n+ell)`.
No arbitrary common-future or forward-prefix map is silently included in
this alphabet. Such operations can be composed with its authenticated
receipts, but require their own source guards and payment comparison.

## 2. An exact terminating router, with an explicit missing conclusion

The guards are pairwise disjoint: G1 needs `2 mod3`, while G5 and G17 need
`1 mod3`; G5 needs `4 mod9`, while G17 needs `1 mod9`.
The least positive odd native inputs are respectively `5,13,4357`, with
children `3,11,4079`. All later native odd inputs increase the child by
a positive amount. Thus every native child is positive and odd. Equation
(1), with `2^a<3^b`, makes it strictly smaller than its positive parent.

There is a common quantitative height bound. Put
`lambda=2048/2187`. For every legal step,

\[
G(n)+1\le\frac{2^a}{3^b}(n+1)\le\lambda(n+1).
\tag{2}
\]

The first inequality uses `ell>=1`; equality there holds for G1.
The integer inequality `2*2048^11<2187^11` implies that eleven legal
steps more than halve `n+1`. Therefore the router takes at most
`11*bitlength(n)` steps. The production API enforces this proved cap;
it does not use an orbit-discovery cap as a theorem premise.

The output is a terminal **obligation**. Indeed none of these maps has
the value 1 at a positive odd integer input: solving `G(n)=1` gives
respectively `n=2,7/4,2275/1024`. Consequently a nonroot source can never
reach ROOT by this router alone. A supplied proof at a terminal can still
ground the input, as described in section 5.

## 3. A closed normal form that retains chronology

For a chronological generator word `w`, define its carrier

\[
G_w(n)=\frac{2^A n-C}{3^B}.
\tag{3}
\]

The empty word has `(A,B,C)=(0,0,0)`. Appending `(a,b,c)` gives

\[
(A,B,C)\longmapsto
(A+a,B+b,2^a C+c3^B).
\tag{4}
\]

Every nonempty generated carry is positive, odd, and a 3-unit. The exact
native cylinder is

\[
n\equiv C2^{-A}\pmod{3^B},\qquad n>0\text{ odd}.
\tag{5}
\]

**Final integrality is sufficient for all prefix guards.** For a prefix
`u` and remaining tail `v`, composition gives
`C_w=2^(A_v) C_u+C_v 3^(B_u)`. Reducing the final numerator modulo
`3^(B_u)` proves the prefix numerator divisible by that power. Each
prefix is therefore integral, and induction with section 2 makes every
intermediate positive and odd. Necessity is immediate. This sufficiency
is special to the present positive-native alphabet; integrality alone
does not establish legal positivity for arbitrary signed inverse words.

**The full affine representation is free.** Equal rational affine maps
have equal slope exponents `A,B`, hence equal `C` and the same cylinder
(5). That cylinder has arbitrarily large positive odd members. Its first
legal generator is unique by the disjoint guards. Cancel that invertible
affine map and repeat. A nonempty residual word cannot be the identity
because its real slope is strictly below 1. Hence equal carriers imply
identical marked words, including their order.

The decoder is constructive. Determine the first generator from (5),
using only candidates with `b<=B`. If its label is `(a,b,c)`, the tail has

\[
A_t=A-a,\quad B_t=B-b,\quad
C_t=\frac{C-2^{A-a}c}{3^b}.
\tag{6}
\]

Peel until the identity, and re-encode to authenticate the result. This
uses exact rational-affine data, not a rounding-based identification.

The slope alone loses both count and order information:

| Chronological word | A | B | C |
|---|---:|---:|---:|
| `(1,17)` | 12 | 8 | 9137 |
| `(17,1)` | 12 | 8 | 6913 |
| `(5,5,5,5)` | 12 | 8 | 12325 |

Even the elementary swap matters:
`G5(G1(n))=(16n-23)/27`, whereas
`G1(G5(n))=(16n-19)/27`. Their difference is `4/27`, and their native
source cylinders are disjoint. The formal noncommutativity is not a
license to apply both orders to one supplied integer.

Thus the three-register form is nonvacuous: it retains the entire word,
computes the exact source cylinder, supports exact composition and future
guard queries, and has a common decreasing height. It is not merely an
encoding of all integers. It is also not a fixed-size memory claim:
`A,B,C` and the guard modulus grow with the word.

## 4. Every mixed word has an exact Mersenne phase or a proved empty domain

For `M_K=2^K-1`, equation (3) becomes

\[
G_w(M_K)=\frac{2^{K+A}-D}{3^B},\qquad D=2^A+C.
\tag{7}
\]

The ordered exponential normal form is closed under appending letters:
`D -> 2^a D+c3^B`. Its source guard is exactly

\[
2^{K+A}\equiv D\pmod{3^B},\qquad K\ge1.
\tag{8}
\]

For a nonempty word, `D` is divisible by 3 precisely when its first
letter is G1. At that first letter the respective values are `3,13,4411`;
later letters multiply the residue by a unit. Thus words starting with
G1 have no Mersenne sources.

Every other finite word has **exactly one** phase

\[
K\equiv K_w\pmod{2\,3^{B-1}}.
\tag{9}
\]

This uses the elementary fact that 2 has exact order `2*3^(B-1)` modulo
`3^B`. For completeness, the binomial expansion gives
`v3(4^(3^j)-1)=j+1` by induction. Multiplication by an exponent coprime
to 3 preserves that valuation; odd powers of 2 are `-1 mod3`. The order
therefore equals the number of units, proving that every unit has exactly
one base-two logarithm. The implementation lifts this logarithm one
ternary digit at a time, testing three candidates; it does not enumerate
the Mersenne source. Every resulting `K_w` is odd.

The single-letter phases are `G5: K=5 mod6` and
`G17: K=733 mod1458`. Arbitrarily mixed finite depth is permitted:
every word starting with either of these two letters has infinitely many
positive Mersenne sources. For a simple explicit depth witness, repeated
G5 is legal exactly while its shifted numerator permits division by 9:

\[
\#\text{ initial G5 steps at }M_K
=\left\lfloor\frac{1+v_3(K-2)}2\right\rfloor
\quad(K\ge1\text{ odd}).
\tag{10}
\]

Here `v3(M_K+5)=v3(2^K+4)=1+v3(K-2)`; the `K=1` case is checked
directly. Taking `K=2+3^(2r-1)` gives exactly `r` initial G5 steps.
No finite maximum word depth describes all these sources.

On the other hand, if `K=5 mod18`, its first G5 child

\[
h=\frac{2^{K+3}-13}{9}
\tag{11}
\]

is divisible by 3. No generator then applies. Its actual word `(1,2)`
returns to `M_K`; it supplies no independent ROOT proof. This gives an
infinite exact frontier family, beginning with `31 -> 27`.

**Exact finite-depth frequencies.** A specified word has relative
natural density `3^-B` among positive odd inputs. Distinct words of the
same length have disjoint native cylinders. Hence the sources with at
least `k` router steps have density

\[
\left(\frac13+\frac19+\frac1{2187}\right)^k
=\left(\frac{973}{2187}\right)^k.
\tag{12}
\]

These are finite unions of arithmetic progressions, so no countable
additivity or unproved density exchange is used. Among **odd exponent
indices K**, not among integer source values `2^K-1`, a word starting
with G5 or G17 has relative phase density `3^(1-B)`. Thus for `k>=1`
the corresponding exponent-index survival density is

\[
\frac{244}{729}\left(\frac{973}{2187}\right)^{k-1}.
\tag{13}
\]

Neither frequency is a probability of convergence to ROOT. It measures
how long this particular reduction alphabet remains applicable.

## 5. What a supplied terminal certificate actually transfers

If the router reads `g1,...,gk`, its terminal child has an actual forward
word back to the immutable source: concatenate the listed forward blocks
in **reverse generator order**. Final oddness implies their exact
valuations. No earlier ROOT occurs, because the endpoint is the positive
source greater than 1 and actual `U(1)=1`.

The API `discharge(reduction, terminal_root_word)` re-authenticates the
reduction, strictly replays the supplied terminal ROOT word, checks that
the above return word is its prefix, and extracts the remaining source
ROOT suffix. At source 1 both words are empty. This is a proof transfer
from a genuinely supplied obligation; the production APIs never discover
that terminal certificate. Merely knowing the return prefix cannot
discharge it.

Parent work on grounded seed families supplies a complementary closure
mechanism: an already authenticated seed lift can retain its ROOT proof
while successive inverse-word guards refine the parameter phase. The
present theorem concerns the nongrounded normal form and its frontier;
it does not assume that stronger premise or duplicate its construction.

## 6. Three supplied papers: exact mechanisms and transfer boundaries

The following papers are accepted supplied premises. We do not audit their
headline classifications. The elementary results above have independent
proofs and do not depend on those headlines.

| Source and inspected pages | Mechanism retained | Target here; missing data and cheapest hostile |
|---|---|---|
| [The Campana–Peternell conjecture in dimension six](C:/Users/Eliott/Downloads/main-3.pdf), OpenAI, pp.16–18 | Cubic slice relations are glued with explicit overlap kernels; modular ranks provide upper bounds, while rational identities and independent witnesses supply the missing direction. | The slope/count projection is only a relaxation of the marked affine carrier. Keep the carry and exact native cylinder. The three distinct `(A,B)=(12,8)` examples refute identification by slope. No Fano hypothesis is transported. |
| [The Deligne–Drinfeld conjecture](C:/Users/Eliott/Downloads/paper-21.pdf), OpenAI, pp.15–16,20 | A filtered deletion representation retains ordered data; injectivity of the leading projection needs a separate kernel argument. | Our first-guard decoder proves injectivity by cancellation instead of assuming a leading-coordinate projection is faithful. The `(1,5)/(5,1)` swap is the minimal chronology hostile. The paper's Lie-algebra representation is not identified with these affine maps. |
| [The Isoperimetric Conjecture for the Cubic Flat Three-Torus](C:/Users/Eliott/Downloads/article-1.pdf), OpenAI, pp.4–5,28 | Two attained endpoint bounds propagate through an interval only with a retained differential inequality; interval certificates preserve outward-rounded dependent quantities. | A common height inequality (2) propagates finite termination across all legal branches. Finite examples alone would not do so. The terminal `27` shows why terminating a reduction is weaker than grounding its terminal. No isoperimetric differential equation is assumed for Collatz. |

The concrete transfer is the same disciplined interface in all three:
keep the compatibility witness that makes a projected calculation legal.
Here that witness is an exact source cylinder plus the ordered carry and
immutable source. Positivity of a representation or closure of a family
alone does not replace it.

## 7. Reproduction and finite scope

Run

```text
python -B -X utf8 04-computation/experiments/collatz_child_normal_forms_20261007.py
python -B -O -X utf8 04-computation/experiments/collatz_child_normal_forms_20261007.py
```

The exact universe is all 1,093 generator words through length 6; three
positive source lifts per word; every word split for carrier composition;
all phase classes, with 44 literal Mersenne replays whose phase is at most
4096; all 8,192 odd sources through 16,383; and supplied-terminal proof
transfer for the 512 odd sources through 1023. The latter proofs are
discovered only inside the explicitly experimental test helper, separately
from the production router and transfer APIs. The largest router depth
in the finite source census is 9, first at 2429. There are 4,549 distinct
terminal states in that census.

Positive controls include arbitrary generator composition, exact final
integrality, literal forward receipts, phase lifting and depth examples.
Hostiles include empty Mersenne G1 domains, the all-height 3-divisible
terminal family, slope/order aliases, incompatible source guards,
forged reductions, earlier ROOT padding, and bool/float input aliases.
The code uses explicit checked exceptions, so optimization does not remove
its mathematical tests. The saved output records **87,614 exact checks**.

The proved result is a lossless, closed **guarded reduction grammar** with
an explicit terminal frontier. Closing that frontier by additional rules
or authentic terminal proofs remains a separate task.
