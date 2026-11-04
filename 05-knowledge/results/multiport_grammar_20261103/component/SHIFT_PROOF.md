# A bounded-intercept inverse and the complete order-five return boundary

2026-11-03. **PROVED all-order inverse; FINITE-EXACT complete return certificate.**
Independent analytic/source/reconstruction audit passed. The all-order result
decides every fixed-clock-shift return in any
finite actual strong-component bank. Its order-five specialization exhausts
the common-entry arithmetic return tree. This is not a native source-preserving
program, an unrestricted tournament obstruction, or a global spectrum proof.

## 1. Incoming work and the precise increment

Incoming [S184's order-four theorem](../../session184_proof_frontiers_20261104/spectrum_four/PROOF.md)
at main `0653d73390` proves the seven-state common-entry return tree for all
suffixes with SCCs of order at most four, even allowing arbitrary fixed
integer clock shifts and arbitrary replacement by a new marked suffix.
It uses an initial-zero-prefix quotient and a bounded 2-adic valuation of
the first marked component. That conclusion is credited, not claimed anew.

The full order-five bank breaks that valuation argument. Below we give an
explicit normalized family with unbounded dyadic marginal scale, even at the
same odd marginal parts needed by the return problem. The extra coordinate
that restores a finite inverse is the **joint defect**: fixing the intercept
`c=h-j` bounds both `h` and the marginal of the normalized suffix.

The count/defect formula and seven-coordinate product remain the inherited
[S182 law](../../session182_physical_returns_20261104/induction/PROOF.md).
The complete fixed-marginal inverse and audited 170-signature actual bank are
in [PROOF.md](PROOF.md), with [independent audit](AUDIT.md). The live board is
joint defect; normalized intercept; dyadic scale; actual prefix representatives;
scalar return states; named physical return. The hostile is an unbounded
valuation sequence whose intercept also changes, so valuation alone cannot
decide a fixed requested digit.

## 2. MPG-BOUNDED-INTERCEPT: a finite bound for every finite bank

Let the bank be finite, nonempty, and contain actual nonempty strong
tournaments with specified port subsets. A word is admissible when its last
component is marked; call it normalized when its first is also marked.
The bank must have at least one marked atom, or there are no admissible words.
Let `M=max_A h_A`.

For every actual nonempty component, `p_A>=2h_A`. For a marked component,
`f_A>=h_A` and `g_A>=h_A`: split a Hamilton path immediately after or before
a chosen marked vertex. These are injections into the respective two-slot
marginals. All other signature coordinates are nonnegative.

**Lemma.** For a normalized word with prescribed nonzero intercept `c=h-j`,

```text
h <= B(c) := max(M^2, |c|).                              (1)
```

Proof. With one or two atoms, the Hamilton count is a product of at most two
bank values, hence at most `M^2`. With at least three atoms, peel the first
marked atom `A`, leaving a tail `D` with at least two atoms and marked final
atom `Z`. The tail marginal satisfies

```text
f_D >= p_(tail before Z) f_Z >= 2h_(tail before Z) h_Z
    = 2h_D.
```

The joint product contains `g_A f_D`, so
`j >= g_A f_D >= 2h_A h_D = 2h`. Therefore `c=h-j<=-h`, giving
`h<=|c|`. This argument uses the actual joint cross term; omitting it would
destroy the bound. QED.

For a normalized word with at least two atoms, write its first signature
as `(H,P,F,G,J,E,S)` and its remaining tail in lowercase. Since the full
joint contains `G f_tail` and `G>=H`,

```text
f_tail <= (h_word-c)/H,
f_word = F h_tail + P f_tail
       <= ((F+P) B(c) - P c)/H.                         (2)
```

Thus all normalized words at this intercept obey the explicit finite bound

```text
f_word <= L(c) := max( max_A f_A,
  max_(marked A) floor(((f_A+p_A)B(c)-p_A c)/h_A) ).      (3)
```

The first term covers a one-atom word. There is no integral-radix assumption,
no restriction on `v2(f_A)` relative to `v2(p_A)`, and no chosen order cutoff.
The fixed-marginal inverse can now test any requested normalized `(c,f)`
inside a proved finite universe. This lemma is independent of the special
alphabet and valuation normalization in S184.

## 3. Complete prefix normalization for all clock shifts

For an actual source, write `H_t=c_s+2^t f_s`. A target suffix with a fixed
integer clock shift `q` realizes digit `d` at all legal clocks precisely when

```text
c_target = C := 3c_s+d,
2^q f_target = 3f_s.                                   (4)
```

For `d=2,4`, the integer `C` is nonzero. Put `R=odd(3f_s)`, where `odd`
removes all powers of two. Then (4) asks exactly for target intercept `C`
and odd marginal part `R`; the dyadic factor determines `q`. A finite shift
is legal at all sufficiently large input clocks. No assertion is made for
input clocks where the shifted nonempty prefix would have negative or zero
length.

Factor a target word uniquely as its maximal unmarked initial prefix `P`
followed by a normalized word `D`. The exact scalar laws are

```text
c_target=h_P c_D,          f_target=p_P f_D.             (5)
```

The empty prefix is the scalar identity `(h_P,p_P)=(1,1)`. Hence

```text
h_P divides |C|,          odd(p_P) divides R.            (6)
```

Generate the multiplicative monoid of pairs `(h_A,odd(p_A))` for unmarked
atoms, retaining only pairs satisfying (6). This is a finite subset of
`divisors(|C|) x divisors(R)`. Save one actual prefix word for every pair.
Multiplying by `(1,1)` does not change these two scalar coordinates and
only changes a dyadic clock factor. Other repeated products either reach
a new bounded pair or can be discarded. This proves termination without
an arbitrary prefix-length cutoff.

For each attained pair `(u,v)`, test all

```text
c_D=C/u,       f_D=2^a (R/v) <= L(C/u),       a>=0,      (7)
```

using the complete fixed-marginal inverse and keeping words with first
and last atoms marked and the prescribed joint intercept. Prepend the
recorded actual representative of `(u,v)` to each accepted word.

**Theorem (MPG-SHIFT-INVERSE).** Equations (1)–(7) give a terminating,
necessary-and-sufficient test for scalar acceptance with arbitrary fixed
integer clock shift, for every finite actual strong-component bank and
nonzero requested intercept.

Necessity follows from the literal prefix factorization and bounds. For
sufficiency, the recorded prefix and normalized word construct a literal
admissible suffix, with target `c` and odd marginal part by (5). Choose
`q=v2(3f_s)-v2(p_P f_D)`. This gives (4) and is legal for sufficiently large
clocks. The full attachment is strong by its terminal port. QED.

This prefix-pair quotient is only a quotient for existence of these scalar
count families. An unmarked singleton prefix may change `s` to zero even
though its pair is `(1,1)`. Full observer values, old names and native guards
are therefore not identified by this quotient. Every accepted representative
still retains its full seven coordinates and actual graph.

## 4. Why the old valuation bound fails, and the defect repairs it

Take the actual strong order-five atom with rows `(24,1,3,6,14)` and port
mask16 (only vertex4). Its audited signature is

```text
(h,p,f,g,j,e,s)=(11,100,12,27,0,1,4).                    (8)
```

Append singleton atoms whose port bits, from the first appended atom to
the last, are the binary digits of a positive integer `m`, least significant
first. The last bit is one, so the word remains admissible and normalized.
If `j_T(m)` is the actual joint of that transitive suffix and `b0=m mod2`,
the exact product gives

```text
f=12+100m,
j=27m+11j_T(m)+b0,
c=11-27m-11j_T(m)-b0 <=11-27m.                         (9)
```

For any fixed `R=7*3^l`, powers of two run through all units modulo25:
their order is20 (`2^10=-1 mod25`, and no proper divisor of20 gives1).
Thus there are arbitrarily large integers `a` with `2^a R=3 mod25`.
Set `m=(2^a R-3)/25>0`. Then

```text
f=2^(a+2) R.                                           (10)
```

Consequently **the dyadic valuation of a normalized suffix is unbounded
even for the return problem's fixed odd marginal parts**. This is a genuine
new five-vertex behavior; deleting initial zero-port singletons does not
remove it. However (9) sends `c` to negative infinity with `m`. It supplies
no fixed-intercept return. The bounded-intercept lemma captures exactly the
sidecar that this valuation experiment would otherwise lose.

Late incoming [S185 carriers](../../session185_proof_frontiers_20261104/carriers/PROOF.md)
contains the same marked-five signature and the explicit odd-part63 instance,
with an exact inverse for one fixed front and an arbitrary transitive tail.
This concrete family is therefore credited as concurrent overlap. The general
bounded-intercept inverse above covers arbitrary finite-bank component words
and prefixes; it is the additional mechanism used for our complete trees.

## 5. Complete SCC-order-five return tree

Apply MPG-SHIFT-INVERSE to all 170 audited signature representatives, covering
every actual marked SCC of order at most five. Start at the all-marked TT3
suffix: `(c,R)=(-5,7)`. For each accepted arithmetic state, test both digits
`2,4`. The resulting complete state tree is

```text
{empty, 2, 4, 22, 42, 422, 4222}.                       (11)
```

| Address | c | R | Normalized/prefix representatives accepted | Next digits |
|---|---:|---:|---:|---|
| empty | -5 | 7 | specified initial TT3 | 2,4 |
| 2 | -13 | 21 | 1 | 2 |
| 4 | -11 | 21 | 2 | 2 |
| 22 | -37 | 63 | 2 | none |
| 42 | -31 | 63 | 1 | 2 |
| 422 | -91 | 189 | 1 | 2 |
| 4222 | -271 | 567 | 2 | none |

There are14 successor decisions,87 exhaustive fixed-marginal calls, and
1,463 normalized candidate visits. The prefix-pair possibilities in this
tree are `(1,1)` and `(5,1)`; the latter is tested when divisibility permits,
including the easy-to-miss target `24` with `C=-35`.

The order-five bank supplies new literal representatives at `4`, `22` and
`4222`, but no new arithmetic state or edge. For example `V5_c_19` at `4`
has `(h,p,f,g,j,e,s)=(11,100,42,60,22,4,8)`, giving the same count family
as the old `f=21,c=-11` state after a clock shift. At `4222`, the word
`T1 T0 C1 V5_8_10` has `(27,4032,567,2313,298,3,27)`, with the required
intercept `-271` but a different physical defect. These are explicit
controls against treating scalar arithmetic states as full named states.

All future scalar target tests depend only on `(c,odd(f))` when arbitrary
replacement and arbitrary fixed clock shifts are permitted. Thus exhausting
both successors of all seven states proves the entire address language,
not just its first five sampled levels. It allows a strictly larger set of
operations than any native old-name-preserving return program in this
attachment architecture. Its negative conclusion therefore applies to
those more restrictive programs as well.

**Consequence.** A clock-uniform `+4` return from either original sibling
in this same root/transitive-prefix architecture requires a suffix SCC of
order at least six, or departure from the architecture or clock-uniform
notion. Arbitrary numbers of SCCs through order five, arbitrary changes of
old ports, and arbitrary fixed clock shifts do not suffice.

This is stronger than the fixed-clock tests in PROOF.md and extends the
incoming order-four boundary. It does not exclude a single-clock return,
a suffix chosen differently for every clock, a new attachment frame, or
larger SCCs. In particular it does not contradict the separate finite
scalar spectrum through2000; that census chooses a fresh witness for each
integer and has no uniform-family return requirement.

For a concrete boundary, the separate scalar receipt realizes91 using the
all-marked four-vertex atom at clock2, with `(f,c)=(27,-17)`, and97 using
`V5_29_13` at clock1, with `(f,c)=(55,-13)`. These match the first numerical
targets of `24` and `44`, but their count families have different slopes and
intercepts. The theorem does not assert that91 or97 is an absent count.

## 6. Reproduction and validity gate

The producer is
[`shifted_returns.py`](../../../../04-computation/multiport_grammar_20261103/shifted_returns.py).
It retains every normalization query and every accepted literal word, row,
port mask and signature in `shifted_returns.json`. Both left and right
inverses are compared at every query. All nine accepted non-initial
representative witnesses are independently recounted from their literal
graphs; no witness exceeds the declared recount threshold14.

```text
python3 04-computation/multiport_grammar_20261103/shifted_returns.py
python3 -O 04-computation/multiport_grammar_20261103/shifted_returns.py --output /tmp/shifted_returns_optimized.json
```

The producer receipt has SHA256
`b8851f11aa198697374a822b5a3400ff25be163542b14a7ddc8d2d43ee38b4a9`.
The [independent proof/source/reconstruction audit](SHIFT_AUDIT.md) rebuilds
the bank from permutation covers and uses the larger envelope `M^2+|c|`,
checking93 normalized queries and all14 successor decisions by a separate
full-signature-set inverse. Its53,411 normal/optimized gates agree; all nine
witnesses and18 full attachments are independently checked. Root promoted
the candidate after this audit. The global coverage and native two-branch
return obligations remain open.
