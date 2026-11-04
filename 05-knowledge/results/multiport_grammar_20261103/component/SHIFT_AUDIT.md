# Independent audit: bounded-intercept normalization and the order-five tree

2026-11-03. **PASS — complete analytic/source/replay audit.** This audit accepts
the theorem and scope in [SHIFT_PROOF.md](SHIFT_PROOF.md), not a native return
program or an unrestricted spectrum claim. The producer owns status promotion.
No producer proof, source, or receipt was edited by this audit.

## 1. Analytic bound, including the empty and short cases

The theorem uses a finite bank of actual nonempty strong components with actual
port subsets, and admissible words end in a marked component. Empty-port atoms
are allowed internally. A normalized word starts and ends marked. If the bank
has no marked atom, no admissible word exists; otherwise all maxima are finite.

For every marked atom, splitting any Hamilton path before or after a chosen
port proves both `f_A >= h_A` and `g_A >= h_A`. The unrestricted two-slot count
has its two empty-slot families, so `p_A >= 2h_A`. These facts are actual
combinatorial inequalities, not assumptions about a coordinate cone.

For a normalized word of at least three atoms, the tail after its first atom
contains a nonempty prefix before its marked last atom. The exact product law
therefore gives `f_tail >= 2h_tail`. The first marked atom contributes the
nonnegative cross term `g_A f_tail` to the joint, giving `j_word >= 2h_word`.
Thus `c <= -h_word`. With one or two atoms, `h_word <= M^2` instead. These cases
prove exactly `h_word <= max(M^2,|c|)`.

For a normalized word of at least two atoms, `j_word >= h_A f_tail` and
`h_tail=h_word/h_A`; substitution into the marginal product proves the stated
bound

```text
f_word <= ((f_A+p_A)B(c)-p_A c)/h_A.
```

Taking the floor is valid because the left side is an integer. The separate
maximum over one-atom marginals pays that case. No assertion that `p_A/h_A`
is integral, or that `v2(f_A)<v2(p_A)`, occurs in this proof.

The nonzero-intercept qualification is sufficient and explicit. In fact every
actual suffix has odd intercept: attach a one-vertex transitive prefix and use
Rédei on `H_1=c+2f`. In the advertised return problem, `3c+2` and `3c+4` are
also automatically nonzero for integral `c`.

## 2. Prefix normalization and the necessary converse

The maximal initial unmarked prefix leaves a nonempty normalized word because
the original last atom is marked. For a nonempty unmarked prefix the scalar
laws are exactly `c_total=h_P c_D` and `f_total=p_P f_D`. The empty prefix is
handled as a scalar identity, not as a fictitious seven-coordinate actual atom.

Consequently the pair `(h_P,odd(p_P))` belongs to the finite product of divisor
sets of `|C|` and `R`. All generator coordinates are positive integers. Any
generator other than `(1,1)` increases at least one coordinate; a discovered
pair need not be expanded again. The stored path to a pair provides an actual
unmarked prefix representative. Products excluded by either divisibility test
cannot become divisible later, so pruning is complete.

For each attainable pair, the normalized marginal is `2^a R/odd(p_P)` and has
the proved finite upper bound. This gives finitely many fixed-marginal inverse
queries. The converse reconstructs a literal suffix and chooses

```text
q = v2(3f_source) - v2(p_P f_D).
```

The shift may be negative. The theorem correctly requires the identity only at
legal clocks, equivalently at all sufficiently large source clocks. The final
marked component supplies the exact terminal-SCC strongness guard.

The proof explicitly preserves the quotient's loss boundary. For example,

```text
marked singleton:                (1,2,1,1,0,1,1)
unmarked singleton -> marked:    (1,4,2,1,0,1,0)
```

have the same intercept and odd marginal, but different `s`. Thus the prefix
pair is complete for existence of scalar count families and is not a full
named-observer equivalence. Actual representatives, masks, and full signatures
are retained for witnesses. Additional native predicates require their own
sidecars. This qualification is essential and is present in the primary note.

## 3. Source audit and all-order termination

Read in full:

- [`shifted_returns.py`](../../../../04-computation/multiport_grammar_20261103/shifted_returns.py);
- its imported [`component_inverse.py`](../../../../04-computation/multiport_grammar_20261103/component_inverse.py);
- the actual `shifted_returns.json` receipt and the bank in
  `component_inverse_five.json`.

`f_bound` implements the proved floor bound, with a separate one-atom maximum.
`prefixes` stores actual representative words and a finite divisor-state queue.
`successors` tests every allowed dyadic normalized marginal, compares the left
and right fixed-marginal inverses, filters first marking and the full joint
intercept, and reconstructs literal target graphs. No order threshold limits
acceptance: the threshold14 only controls optional direct recounting.

The imported left inverse decreases positive `f` by at least a factor two in
every nonempty recursive tail. Its one-atom base case is marked. The right
inverse strictly decreases `p`; its finite `(h,p)` monoid uses the valid bounds
`h<=f` and `p<=C_bank f`. The latter follows by isolating the last marked atom
and comparing its `p_A/f_A` ratio. Empty prefixes are allowed only at the exact
right-inverse base case `(h,p,f)=(1,1,0)`.

The order-five bank enumerates all labelled five-vertex tournaments and every
actual port subset on each strong one. Quotienting by all seven coordinates is
sound for this ordered-product scalar acceptance language, and retains actual
representatives. Marked versus unmarked is preserved since `f=0` exactly for an
empty port set. It is not a quotient for arbitrary future native operations.

The outer breadth-first loop records whether pending states remain. In this
receipt it actually exhausts all successors at depth5. Thus its configured
depth limit did not terminate the computation prematurely. All fourteen
successor decisions are present, with both digits checked from every accepted
arithmetic state.

## 4. Independent reconstruction, not a replay-only audit

The new self-contained checker
[`shifted_audit.py`](../../../../04-computation/multiport_grammar_20261103/shifted_audit.py)
imports no producer module. It independently enumerates every labelled strong
tournament of orders1–5, computes path covers by literal permutations, and
rebuilds the entire 170-element full-signature bank. The labelled strong counts
are `1,0,2,24,544`; the independent bank equals the producer bank as a set.
Each producer atom's actual graph and mask is separately recounted.

For the complete successor decisions, this checker stores sets of full
signatures instead of the producer's literal atom words. It uses the strictly
no-smaller envelope

```text
B_audit(c)=M^2+|c|,
```

which produces **93** normalized inverse queries rather than the producer's87.
The all-order completeness argument is still the analytic bound and decreasing
inverse, not a sampled vertex-order cap. This separate reconstruction recovers
exactly the producer's full-signature answer sets on all fourteen targets.
It also independently checks each displayed prefix pair, normalized intercept,
and tight producer bound.

All nine accepted non-initial literal witnesses are recounted by an additional
end-point-removal recurrence, using the converse tournament for start counts.
For each witness, two actual full attachments are constructed and checked for
strongness, Hamilton count, and the physical root defect `F_root-H=j`. This is
18 independently checked full graphs. It includes every new five-vertex-atom
witness, not only the smallest old representative.

The complete arithmetic addresses are

```text
empty, 2, 4, 22, 42, 422, 4222.
```

In receipt order, the fourteen distinct full-signature witness counts are

```text
1,2,2,0,1,0,0,0,1,0,2,0,0,0.
```

The extra literal representatives occur at `4`, `22`, **and `4222`**. Arithmetic
state uniqueness must not be promoted to physical word uniqueness. The primary
now states this correctly.

## 5. The valuation hostile and its exact missing coordinate

The new atom `(24,1,3,6,14)` with port mask16 independently has signature
`(11,100,12,27,0,1,4)` and is strong. For its appended binary transitive word of
mask `m>0`, the marginal and joint are exactly

```text
f=12+100m,
j=11 j_T(m)+27m+(m mod2).
```

The final bit of the binary word is marked; the initial atom is marked. Thus
this is a normalized actual word. The `m mod2` term is the first-vertex start
count of the transitive word, multiplied by the atom's Hamilton endpoint count1.
Seven literal masks are checked independently in the audit source.

The powers of2 modulo25 have order20 and exhaust its unit group. Every
`R=7*3^l`, `l>=0`, is a unit, so the congruence `2^aR=3 mod25` has an unbounded
arithmetic progression of solutions in `a`. For large solutions,
`m=(2^aR-3)/25` is a positive integer and `f=2^(a+2)R`. The normalized marginal's
dyadic valuation is therefore genuinely unbounded at the relevant fixed odd
parts. The explicit joint forces `c<=11-27m`, so the intercept cannot remain
fixed along this construction.

The checker verifies the complete twenty-element unit orbit and63 exact
arithmetic instances, in addition to the literal graph controls. These controls
do not replace the elementary unbounded congruence argument. The conclusion is
exactly that the old valuation bound does not extend, while the new defect
bound still applies. It is not an example of a missed fixed-intercept return.

## 6. Accepted consequence and excluded overstatements

Every suffix whose SCCs have order at most five belongs to the declared actual
component language after replacement by signature representatives. The complete
arithmetic search permits arbitrary suffix replacement and arbitrary fixed
clock shifts. It is therefore more permissive than a native old-name-preserving
return in the same architecture.

The rejected addresses `24` and `44` prove the stated consequence: a
clock-uniform `+4` return from either original sibling requires a suffix SCC of
order at least six, or a departure from the stated root/transitive-prefix
architecture or uniform-clock requirement. This does not exclude a
single-clock witness, clock-dependent suffix choice, larger SCC, or a different
attachment frame. It does not prove or refute the global integer spectrum.

## 7. Replays and frozen pins

Fresh producer normal and optimized runs, with their receipts written only to
`/tmp`, are byte-identical to the frozen 22,977-byte producer receipt. The
independent checker also passes normal and optimized modes byte-identically:
**53,411 always-active audit gates**, with its 877-byte output frozen in
[`SHIFT_AUDIT.out`](SHIFT_AUDIT.out).

```text
python3 -B 04-computation/multiport_grammar_20261103/shifted_returns.py --output /tmp/shift-audit-normal.json
python3 -B -O 04-computation/multiport_grammar_20261103/shifted_returns.py --output /tmp/shift-audit-optimized.json
python3 -B 04-computation/multiport_grammar_20261103/shifted_audit.py
python3 -B -O 04-computation/multiport_grammar_20261103/shifted_audit.py
```

| Artifact | Bytes | SHA256 |
|---|---:|---|
| `shifted_returns.py` | 5,203 | `970794e4810090c17c6aa4ad324b7cb1a6ba0bfc656511ed2ff0800ac2867b34` |
| `component_inverse.py` | 12,220 | `a101041d142c2c5821eca7cd91b1dec145fbc71c671f4315e505e538f29a2dc5` |
| `shifted_returns.json` | 22,977 | `b8851f11aa198697374a822b5a3400ff25be163542b14a7ddc8d2d43ee38b4a9` |
| `component_inverse_five.json` | 101,800 | `225e228394ac04faea87ec663fb1391ecb1256a27eae61b44cbb4d5db99978f6` |
| `shifted_audit.py` | 13,327 | `d136f56df8a75b48f814c7016fa707b0ea479d5bb239cc7398076cf69b106df8` |
| `SHIFT_AUDIT.out` | 877 | `7a73ddedad68eb237d0fadb1a2451cbb5a8ac8f106e28cd0d1cf448619b4d86d` |

The audit requests no mathematical or source correction. The producer's final
prose now correctly says independent reconstruction rather than a full forward
search; that clarification was reread and accepted. The complete independent
algorithm is described precisely in §4 above. All finite controls retain their
declared finite scopes.
