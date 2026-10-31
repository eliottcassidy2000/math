# Session 125 independent audit: Chung–Graham / G_32 coordinate bridge

## Verdict

**PASS.** The producer's interval theorem, greedy encoder, positive carry
normalizer, unique minimum-summand corollary, successor, parity tracker,
`D(1)` triples and Zeckendorf identities are correct under their stated
hypotheses. No mathematical or implementation correction was required.
The automaton statement can be strengthened: the complete DFA for legal
**even-valued** words has exactly **nine** states, of which eight are live.
A full proof and machine-checkable certificate appear below.

This is an arithmetic encoding of a chosen graph-order coordinate. It is not
an extension theorem for square-sum Hamiltonian paths or cycles. The established
Chung–Graham representation is not a newly discovered number system. The
minimality result below is a direct verification for this specified language,
not a historical novelty claim.

## Independent proof checks and important boundaries

### 1. Intervals, rank and greedy versus normalization

Let `A_i=F_(2i+2)` and `A_-1=0`. For suffix length `k`, the free-state values
are exactly `[0,A_k-1]`; the locked-state values are exactly
`[0,A_k-A_(k-1)-1]`. At `k=0` both sets consist of the empty word and value
zero. Setting `a=A_(k-1)` and `b=A_(k-1)-A_(k-2)`, the free leading-digit
blocks have lengths `a,a,b`, and the locked blocks have lengths `a,b`.
They are adjacent, disjoint numerical intervals because
`A_k=2a+b` and `A_k-A_(k-1)=a+b`.

This simultaneously proves the producer's count, lexicographic rank and
uniqueness claims. In particular, quotient/remainder greedy stays inside the
appropriate next-state interval. It is not merely empirically observed to
avoid forbidden blocks. On the other hand, local normalization is a separate
algorithm that works directly on nonnegative coefficients; it is not defined
by evaluating and calling the greedy routine.

Edge cases checked explicitly include empty input, all leading zeros, the
least-significant boundary `A_-1=0`, maximal words, and overflow to one more
digit under successor. Uniqueness always means modulo leading zeros.

### 2. Carry identities and the actual unique minimum

The carry rules have nonnegative outputs. Their digit-sum decreases are:

- `3 A_0 -> A_1`: decrease 2;
- `3 A_i -> A_(i+1)+A_(i-1)`, `i>=1`: decrease 1;
- `2 A_i + A_(i-1)+...+A_(l+1)+2 A_l`
  `-> A_(i+1)+A_(l-1)`: decrease `r+2`, where `r=i-l-1`,
  except the `l=0` boundary has decrease `r+3`.

If an input is not canonical, some coefficient is at least 3, or all
coefficients are at most 2 and a forbidden `2 1* 2` factor exists. Thus a
rewrite is always available on a noncanonical input. The integer potential
is nonnegative and decreases strictly, so every rewrite strategy terminates.
Every irreducible word is legal and therefore is the unique interval normal
form. This proves strategy independence without requiring local confluence.

Crucially, any *noncanonical* coefficient array begins with a strict decrease
and ends at that canonical array. Hence the canonical array has strictly
smaller summand count than **every other nonnegative coefficient array** of
the same value. This proves uniqueness of the optimum, not just optimality
of greedy, and not just uniqueness among legal words.

The finite independent check counts unordered coefficient vectors. Its
coin-type recurrence partitions arrays into `d_j=0` and `d_j>0`; removing one
copy of coin `j` gives a bijection in the latter case. These disjoint cases
explain why its number of optimal arrays is not a count of ordered
compositions. A separate exhaustive enumeration of coefficient arrays for
values at most 500 verifies both DP costs and counts.

The statement does **not** permit negative coefficients, arbitrary real
coefficients, or additional odd-indexed Fibonacci denominations. For example,
allowing the ordinary Fibonacci denomination 2 already introduces other
minimum decompositions in related problems.

### 3. Exact Zeckendorf relation

Splitting each ternary digit into two binary digits and inserting zeros at
odd Fibonacci indices produces two individually legal Zeckendorf arrays.
Their numerical sum is correct; their unnormalized coefficientwise sum need
not itself be Zeckendorf-canonical. The producer makes this distinction.

For a Zeckendorf array `z_r`, with `z_1=0`, substituting
`F_(2j+1)=A_j-A_(j-1)` gives
`c_j=z_(2j+2)+z_(2j+1)-z_(2j+3)`. Adjacency exclusion puts these coefficients
in `{-1,0,1}`. This is an exact signed identity, not a license to feed negative
coefficients to the positive normalizer. The documented evaluation followed
by greedy re-encoding is a valid conversion in that direction.

### 4. The triple and the graph boundary

The invariant `a^2-3ab+b^2=1`, propagated by `(a,b)->(b,3b-a)`, gives
`ab+1=(b-a)^2`, `a(3b-a)+1=b^2`, and the corresponding third pair. Thus all
consecutive triples of even Fibonacci terms are `D(1)` triples.

At the unit vertex the equations `1*b+1` and `1+b` coincide, explaining the
square-sum edges `1--3` and `1--8`. They do not make `3--8` an edge:
`3+8=11`, whereas the different operation `3*8+1=25` is square. Nor do
nonconsecutive even Fibonacci terms form an infinite pairwise `D(1)` set:
`1*21+1=22` is already a counterexample.

For fixed ordinary graph labels, `t=m-32` and `n=m-31` are two different
coordinate conventions. They give, respectively,
`m mod 2 = t mod 2` and `m mod 2 != n mod 2`. Consecutive shorter/longer
orders have opposite parity in either convention, but changing convention
changes which *coordinate* is called odd. None of this proves graph
Hamiltonicity or reconstructs the hidden prefix below order 32.

## Minimal automaton theorem

Define

`E = { w in {0,1,2}* : w contains no factor 2 1* 2 and V(w) is even }`.

Leading zeros are allowed; the empty word has value zero and is accepted.
Read words most-significant digit first. A live state is `(g,N,P)`, with
`g` either Free or Locked, and `N,P` the current and one-place-shifted values
modulo 2. Start at `(Free,0,0)`. On digit `d`, update

`(N,P) -> (P+d, N+P+d) mod 2`.

The grammar transitions are `Free: 0->Free, 1->Free, 2->Locked` and
`Locked: 0->Free, 1->Locked, 2->sink`; the sink is absorbing. A live state is
accepting precisely when `N=0`. This recognizes `E` by the exact integer
shift identity and the forbidden-factor grammar.

**Theorem.** This nine-state complete DFA is minimal. In particular all
eight live states are necessary, not just reachable.

**Reachability.** The access words are:

| State | Access word |
|---|---|
| Free,0,0 | empty |
| Free,1,1 | `1` |
| Free,1,0 | `10` |
| Free,0,1 | `11` |
| Locked,0,0 | `2` |
| Locked,1,0 | `12` |
| Locked,1,1 | `21` |
| Locked,0,1 | `102` |
| sink | `22` |

**Pairwise distinguishability.** Consider two live states.

1. If their `N` registers differ, the empty suffix distinguishes them.
2. If their `N` registers agree but their `P` registers differ, suffix `0`
   distinguishes them: its new current register is the old `P`, and zero
   is legal from either grammar state.
3. Otherwise their registers agree and their grammar states differ. Suffix
   `2` sends the Locked state to the sink and the Free state to a live state.
   From that latter state append zero, one or two zeros until `N=0`.
   This always succeeds: zero acts on registers by
   `T(N,P)=(P,N+P)`, which fixes `00` and cycles the other three states
   `01 -> 11 -> 10 -> 01`. Thus each nonzero cycle contains a state with
   first coordinate zero. All these zeros are legal.

Finally every live state is distinguishable from the sink using zero, one
or two zeros by the same register-cycle argument. Therefore all nine
reachable states have different right languages. The Myhill–Nerode theorem
gives the lower bound nine, attained by the displayed DFA. QED.

The certificate `minimal_dfa_certificate.json` supplies the transition
table, access words, independent partition refinement, and all 36 pairwise
distinguishing suffixes. This theorem is for the stated leading-zero language;
changing display conventions, the reading direction, or the accepted parity
requires a separately specified automaton.

## Reproducible independent verification

Run from the containing `numeration` directory:

```sh
python audit_port_numeration/independent_audit.py \
  audit_port_numeration/independent.normal.json
python -O audit_port_numeration/independent_audit.py \
  audit_port_numeration/independent.optimized.json
```

All checks raise explicit exceptions rather than relying on `assert`.
The original producer suite was also replayed normally and under `python -O`;
those two reports each record 1,359,569 successful checks.
Both independent runs passed **3,639,641 explicit checks**. Independent
partition refinement yields class counts `2 -> 7 -> 9` and confirms that
all 36 pairs already have distinguishing suffixes of length at most two.
The independent reports cover 265,720 words, 250,001 coin minima with unique
minimizers, 1,365 exhaustive rewrite-strategy inputs, 1,000 randomly chosen
rewrite strategies, 100,001 Zeckendorf bridge values, 1,000 `D(1)` triples,
and all nine DFA states / 36 state pairs. These finite tests complement the
universal proofs rather than replace them.

`producer_input.sha256` records the producer snapshot audited. No source
worktree or remote repository edits or publication actions were performed.

## Primary-source check

Independently inspected on 2026-10-31:

- [Rob Burns, *Chung-Graham and Zeckendorf representations*, Theorem 2 and
  Section 2](https://arxiv.org/html/2502.18870v1): confirms the even-indexed
  weights, digit restrictions, uniqueness, and two-binary-array splitting
  bridge. Its digit convention is LSD first with odd-index zero slots, so
  the producer's compressed MSD words require the documented change of
  convention.
- [Sungkon Chang, *The Chung-Graham Expansion*, Section 1](https://arxiv.org/html/2504.00314v1):
  independently confirms the classical attribution and exact expansion
  rule. This modern primary paper also explicitly distinguishes general
  nonnegative coefficient sequences from the canonical rule.
- [The original Chung–Graham author-hosted scan](https://mathweb.ucsd.edu/~ronspubs/81_12_irregularities.pdf)
  resolves as a 42-page PDF, but this browser supplied no extracted text.
  This audit does not claim to have read its theorem directly. Historical
  identification above rests on the inspected primary research papers;
  all specific arithmetic and minimality claims were proved independently.

## Portable package checks

The audit script locates its producer module relative to its own file, not
the current working directory. `producer_input.sha256` uses paths relative
to this audit directory; verify it here with `sha256sum -c producer_input.sha256`.
`portability_verification.json` records normal and optimized replays after
relocating this package and invoking it from an unrelated working directory.
The original scratch audit is preserved; these copies differ only in package
paths, reproduction instructions and regenerated packaging metadata.
