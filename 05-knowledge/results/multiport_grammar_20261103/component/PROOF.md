# A complete finite inverse beyond the integral component chart

2026-11-03. **PROVED:** complete inverse algorithms for any finite bank of
actual nonempty strong marked tournaments, first at fixed marginal and then
at fixed Hamilton count. **FINITE-EXACT:** the full strong-component banks
through orders four and five, including every port subset, in the finite
universes below. **OPEN:** unbounded integer coverage and physical return
closure. An accepted scalar target is not an operation on an old named source.

## Inherited mechanism and the new question

The exact seven-coordinate ordered-sum law is inherited from
[S182, terminal-port strongness and ordinal components](../../session182_physical_returns_20261104/induction/PROOF.md).
The earlier [singleton/triangle inverse](../../multiport_defect_grammar_20261103/COMPONENT_LANGUAGE.md)
already decides every word length in that narrower language. It identifies
the strong four-vertex atom as the next test: its `p/h=32/5` leaves the
integral mixed-radix representation used in that inverse. We extend the
inverse to arbitrary finite actual strong-atom banks. We do not reprove a
new general substitution theorem or infer physical equivalence from a
matching signature.

During this session, incoming main `0653d73390` supplied
[S184's complete SCC-order-four return theorem](../../session184_proof_frontiers_20261104/spectrum_four/PROOF.md).
It independently establishes the stronger four-bank conclusion that the
entire reachable address tree has seven states even with arbitrary fixed
integer clock shifts. Our four-bank computation is therefore an independent
overlap/control, not a new claim of that obstruction. The general finite-bank
inverse, full order-five bank, and fixed-target scalar decisions below are
the useful extensions here. Its valuation normalization also identifies
the next question for the order-five bank.

The inherited positive control is the common-clock entry with both digits
from `H_t=7*2^t-5`. The hostile is its missing second `+4` return, with
intercepts `-35` for `24` and `-29` for `44`. The corrected near miss is
equating target acceptance with a reusable source-preserving rewrite. The
least-used sidecar is the unrestricted two-slot count `p`, which supplies a
second, independent finite inverse even when `p/h` is not integral.

Our concept board is: actual atoms; paired marginal/joint data; two inverse
directions; terminal SCC guards; scalar coverage; named return compatibility.
The source is a literal ordered word of marked components. Its signature
map preserves the ordered-assembly evaluator and terminal-port guard. It
forgets vertex names and internal arc guards; those must remain in actual
representatives for any subsequent physical compiler.

## 1. Definitions and bounds

Write `sigma(B,W)=(h,p,f,g,j,e,s)` for the inherited signature: Hamilton
paths; unrestricted ordered two-slot covers with empty slots allowed;
first-slot port endpoint; second-slot port start; both conditions together;
Hamilton paths ending and starting at ports. For `C -> D`, with lowercase
and uppercase signatures respectively, the exact law is

```text
(hH, pP, fH+pF, gP+hG, jH+hJ+eS+gF, hE, sH).             (1)
```

Let a finite atom bank contain actual nonempty strong tournaments and
specified port subsets, possibly empty. A word is admissible exactly when
its final atom has a nonempty port set. Its full root attachment with any
nonempty transitive prefix is then strong, because its port set meets the
terminal SCC. Every atom and every word has `h>=1` and `p>=2`.

For any chosen port `w`, cut each Hamilton path immediately after `w`.
This injects all `h` paths into the first-end marginal `f`. Different
ports give disjoint subsets, so

```text
f >= |W| h; in particular f >= h for every admissible word. (2)
```

The joint covers are a subset of the first-end covers with nonempty
second slot. The excluded second-empty covers are exactly the `e`
Hamilton paths ending at a port. Thus

```text
0 <= j <= f-e.                                           (3)
```

Let `C=max p_A/f_A` over marked atoms in the bank. If there are no marked
atoms, there are no admissible words and both inverse problems are empty.
Otherwise `C` is a finite positive rational. If the last atom is `A`,
then `f_word >= p_prefix f_A` and `p_word=p_prefix p_A`, with the same
inequality also valid for a one-atom word. Consequently

```text
h_word <= f_word,          p_word <= C f_word.            (4)
```

These inequalities are restrictions on realized words. They do not assert
that all nonnegative signature vectors satisfying them are realizable.

## 2. MPG-FINITE-INVERSE: fixed marginal, every word length

**Theorem.** For a finite actual strong-atom bank and a positive integer
`f`, the following terminating algorithm enumerates every admissible word
with marginal `f`, and no others.

First enumerate the finite multiplicative monoid of atom Hamilton counts
up to `f`. Each possible total `h` must be in this set by (1) and (2).
For a state `(h,f)`, try every first atom `A`, with signature beginning
`(H,P,F)`. There are two alternatives:

1. Stop with the one-atom word exactly when `A` is marked and `(h,f)=(H,F)`.
2. For a nonempty tail, require integer values
   `h_tail=h/H` and `f_tail=(f-F*h_tail)/P`, both positive, then recurse.

Any nonempty tail is itself admissible. The displayed identities are
precisely the Hamilton and marginal coordinates of (1), so necessity and
sufficiency follow by induction. Moreover `f_tail<=f/2`, since `P>=2`
and `F>=0`. This proves termination and covers arbitrary numbers of
unmarked singleton atoms; no tournament-order or word-length cutoff was
chosen. Each physical atom word has a unique first-atom decomposition.
The admissible word length is at most `floor(log2 f)+1` as a consequence,
not an input limit. QED.

There is a genuinely separate right inverse. Enumerate multiplicative
pairs `(h,p)` of atom counts satisfying `h<=f, p<=Cf`. The pair set is
finite. Peel the last atom `A`, writing

```text
h_prefix=h/H, p_prefix=p/P,
f_prefix=(f-p_prefix*F)/H.                               (5)
```

Require integrality, `h_prefix,p_prefix>=1`, and `f_prefix>=0`. The prefix
need not be admissible: it can have an empty terminal port set or no ports
at all. Stop at the empty prefix `(h,p,f)=(1,1,0)`, and impose the marked
last-atom condition only on the complete word. Since `p` decreases by at
least two at each peel, this inverse also terminates and is complete.
Erroneously requiring `f_prefix>0` or a marked prefix would lose actual
words; the independent implementation retains these boundary cases.

The number of words can grow quickly. This is a finite decision theorem,
not a polynomial-time assertion in the bitlength of `f`.

## 3. MPG-FIXED-TARGET: complete scalar acceptance

For an admissible word and a common prefix clock `t>=1`, the actual
attachment has

```text
K = h-j+2^t f,              physical root defect = j.     (6)
```

**Theorem.** Given a positive target `K`, membership in this finite-bank
attachment language is decidable by finitely many of the preceding inverse
calls and literal arithmetic tests, with no order cutoff.

Indeed, (3) implies

```text
K >= h+e+(2^t-1)f >= 1+f,      K >= 2^t.                (7)
```

Thus it suffices to enumerate `1<=f<=K-1`, and
`1<=t<=floor(log2 K)`, and test (6) on each complete word. Conversely a
passing word is an actual strong attachment, so the test is sufficient.
`K=1` is supplied separately by a singleton tournament, outside this
nonempty-prefix/nonempty-suffix attachment grammar. QED.

This provides an exact selector that either returns a fresh actual witness
or certifies absence in the chosen language. It does not prove the answer
is always yes for the desired integers, nor does it retain an arbitrary
old source's names, size, arcs, or root.

## 4. The complete small-component banks

For SCC size at most four, use the singleton with either port bit, the
directed triangle with canonical masks `0,1,3,7`, and the strong four-vertex
tournament with rows `(8,1,3,6)` and all sixteen port masks. Independent
labelled enumeration gives strong counts `1,0,2,24` in orders one to four,
and marked isomorphism-class counts `2,0,4,16`. Hence these 22 atoms cover
every possible marked SCC of order at most four. Here `C=32/5`.

For order five the source enumerates all `2^10=1024` labelled tournaments,
retains all 544 strong ones, and computes all 32 port subsets on each.
The 17,408 cases have 148 different full seven-coordinate signatures,
disjoint from the earlier 22, giving a 170-element bank with `C=28/3`.
Each retained signature has literal adjacency and a literal port set.
Independent permutation-and-cut enumeration checks all 17,408 signatures.

Full-signature equality is a congruence for (1), and a strong component is
marked exactly when `f>0`. Thus this quotient preserves scalar word
acceptance and the terminal-port condition. It does not identify marked
isomorphism classes or physical names. Candidate counts in this larger
bank count representative-signature words, not all literal labelled words.
All finite decisions below cover **every suffix whose SCCs have order at
most five**, even though suffix order itself is unbounded.

## 5. Fixed-clock return tests

Start from `H_t=-5+7*2^t`. A digit word in `{2,4}` of length `l` requires
`f=7*3^l` and intercept `a_0=-5, a_i=3a_(i-1)+d_i`. The complete inverse
gives the following exact all-length acceptance tests:

| Level | f | Four-bank candidates | Four-bank accepted words | Five-bank candidates | Five-bank accepted words |
|---|---:|---:|---|---:|---|
| 0 | 7 | 2 | initial target | 2 | initial target |
| 1 | 21 | 5 | 2, 4 | 8 | 2, 4 |
| 2 | 63 | 9 | 22, 42 | 14 | 22, 42 |
| 3 | 189 | 18 | 422 | 31 | 422, 442, 444 |
| 4 | 567 | 29 | 4222 | 53 | 2242, 4222 |

Thus **neither bank repairs `24` or `44` at any word length** in the fixed
prefix clock. The order-five bank accepting `442` and `444` does not repair
the missing intermediate `44`: endpoint acceptance need not be closed
under taking the preceding arithmetic digit. This is an explicit hostile
to reading a successful deeper target as a returned physical induction.

The root defect is recomputed as the actual joint `j`, not held constant.
The receipts retain all actual candidate words, rows, masks and seven
coordinates. In the default four-bank run all 63 candidates have order at
most twelve and are independently recounted; in the five-bank run all 108
are recounted. Selected accepted attachments at two distinct clocks give
14 and 20 further count, strongness and physical-root-defect controls.

## 6. Bounded spectrum signal and exact remaining question

The complete fixed-target algorithm, plus the separate singleton seed,
was run for every positive odd integer at most 2000. There is no suffix
order cutoff. The finite missing sets have sizes:

| Largest allowed SCC | Missing odd values through 2000 |
|---|---|
| 1 | 456 values |
| 3 | 58 values |
| 4 | exactly 7, 21, 173 |
| 5 | exactly 7, 21 |

For the five-bank universe, both inverse directions agree for every
`f=1,...,1999` (72,190 candidate words). Every one of the 997 nontrivial
accepted targets has an actual strong representative independently checked
by literal subset-path DP, including its physical root defect. The
singleton accounts for the remaining accepted value 1.

This suggests the sharply stated **OPEN** sufficiency question: does the
same attachment language with SCCs of order at most five, together with
the singleton seed, realize every positive odd integer except 7 and 21?
The finite receipt proves only the displayed range. A proof would supply
the realization half of the global classification using this restricted
language; the exclusions and parity still require their inherited proofs.
A counterexample would instead identify a precise target requiring a
larger SCC, a different attachment model, or another construction.

Even affirmative scalar coverage would not itself yield a reusable named
two-branch program. That additionally requires a construction respecting
old names and exits, valid native guards, and a return condition. The
[ordered-end obstruction](../closure/PROOF.md) excludes a broad simple
choice of such programs. The [coloured forest grammar](../transport/PROOF.md)
can evaluate a richer supplied word and its changing defects, but does
not supply the missing accepting word. These three interfaces are kept
distinct.

## 7. Reproduction and independent audit

Run from the repository root:

```text
python3 04-computation/multiport_grammar_20261103/component_inverse.py
python3 -O 04-computation/multiport_grammar_20261103/component_inverse.py
python3 04-computation/multiport_grammar_20261103/component_inverse.py --max-atom-order 5 --output 04-computation/multiport_grammar_20261103/component_inverse_five.json
python3 -O 04-computation/multiport_grammar_20261103/component_inverse.py --max-atom-order 5 --output 04-computation/multiport_grammar_20261103/component_inverse_five.json
python3 04-computation/multiport_grammar_20261103/spectrum_acceptance.py
python3 -O 04-computation/multiport_grammar_20261103/spectrum_acceptance.py
```

All validity guards raise exceptions explicitly and survive optimization.
The inherited triangle-only inverse is a regression control. The two
independent inverse directions agree on the dense range `f=1,...,100` and
on every fixed-clock target tested. Exact hashes, independent permutation
controls and source pins are recorded in [AUDIT.md](AUDIT.md) and
[REPLAY.md](REPLAY.md). No external novelty claim beyond the targeted
repository comparison is made.
