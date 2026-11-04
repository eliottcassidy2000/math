# Independent audit: finite-bank inverse and scalar acceptance

2026-11-03. **INDEPENDENT ANALYTIC / SOURCE / FINITE-EXACT AUDIT PASS.**
The unbounded conclusion is a terminating inverse for a declared finite atom
bank. The bounded spectrum receipt does not prove unbounded coverage or a
source-preserving two-branch return grammar. No external novelty claim is made.

## 1. Audited interfaces and proof

The seven coordinates are those of the proved
[S182 ordered-component law](../../session182_physical_returns_20261104/induction/PROOF.md):
`(h,p,f,g,j,e,s)`. In particular the joint update retains both `eS` and `gF`.
An atom is a nonempty strong tournament with an actual specified port set.
Internal atoms may have no ports; the final atom must have at least one.
The attachment uses a nonempty transitive prefix and the same marked root
convention as S182. Its count is `h-j+2^t f`, and its actual root defect is `j`.

Fix any marked vertex in the final atom. Cutting each Hamilton path immediately
after this vertex injects the Hamilton paths into the covers counted by `f`.
Thus `1<=h<=f`. More generally `f>=|W|h`. Every atom has `p>=2h>=2`, by the
two possible empty slots. Consequently the possible `h` at fixed positive `f`
are exactly tested within the finite multiplicative monoid of the atom Hamilton
counts bounded by `f`; factors equal to one introduce no new monoid values.

For a nontrivial word with first atom `(H,P,F,...)`, the remaining word retains
the terminal-port guard and must satisfy

```
h_tail = h/H,
f_tail = (f-F*h_tail)/P.
```

Both values must be positive integers. Conversely each recursively realized
tail satisfying these equalities gives the claimed marginal. Since `F>=0` and
`P>=2`, `f_tail<=f/2`. The single-atom stopping rule is exactly
`(h,f)=(H,F)` with nonempty ports. Stopping and continuing are both tested;
one must not discard longer words merely because a singleton solution exists.
This proves completeness and termination at every positive `f`, without a
chosen length cutoff. As a consequence, a word has at most
`1+floor(log2(f))` atoms.

The independent right inverse retains `(h,p,f)`. For a last atom `(H,P,F,...)`,

```
h_prefix = h/H,
p_prefix = p/P,
f_prefix = (f-p_prefix*F)/H.
```

The first two must be positive integers and the last a nonnegative integer.
The prefix need not retain the final-port guard: its `f` may be zero. The formal
empty prefix is exactly `(1,1,0)`, an identity only for these three coordinates;
it is not silently used as an empty seven-coordinate signature. Every
nonempty step decreases `p` by a factor at least two. The last atom of the
original word has `F>0`, so

```
p_total <= C*f,    C=max_terminal_atoms(P/F).
```

This bounds the independently enumerated `(h,p)` monoid. The implementation
admits arbitrary prefix ports and filters the actual final atom after the
right recursion. It does not impose the old singleton/triangle condition
`h divides f`: the literal C4 singleton-port control has `h=5,f=7`.

For a fixed target `H0` and clock `t>=1`, the covers counted by `j` are a subset
of those counted by `f` after removing the `e` covers with the second slot
empty. Hence

```
j <= f-e,
H0 = h-j+2^t*f >= h+e+(2^t-1)*f >= 1+f,
H0 >= 2^t.
```

Therefore `1<=f<=H0-1` and `1<=t<=floor(log2(H0))` are exhaustive limits.
They prove decidability of exact scalar acceptance in this finite-bank
attachment grammar. The singleton count one is explicitly a separate seed.
These bounds do not concern a clock-zero construction or arbitrary physical
rewrites of a given source.

## 2. Independent finite universes and quotient coverage

The separate verifier is
[`audit_component.py`](../../../../04-computation/multiport_grammar_20261103/audit_component.py).
It enumerates all directed tournament orientations, tests strongness by its
own vertex-reachability traversal, and counts ordered two-path covers by
literal full permutations and cuts. It does not derive those counts from
the supplier's subset DP or composition law. A full permutation and its cut
uniquely specify the two ordered path slots, including either empty slot.
Endpoint vectors and the joint endpoint matrix give all port observables.

The complete labelled strong counts at orders 1–5 are `1,0,2,24,544`.
Independent relabelling of every tournament and every port mask through order
four gives strong shape counts `1,0,1,1` and marked isomorphism class counts
`2,0,4,16`. Thus the 22-atom bank is complete for all ordered SCC words with
each component of order at most four. In particular `(8,1,3,6)` represents
the unique strong four-vertex shape, its 16 literal port subsets are distinct
marked classes, and the four C3 masks are the four marked cyclic orbits.

For order five the verifier retains all 544 labelled strong tournaments and
all 32 port masks: **17,408 raw marked cases**. Every seven-coordinate value
agrees with the supplier's independently implemented subset DP. There are
148 distinct values, all new relative to the old 22, giving 170 representatives.
The full sorted signature set has SHA256
`68076e4ba43ee493e70cb5c9fd8b6daff62b83d707cec3c7764f65dd15d4bfe1`
when encoded with the verifier's compact JSON convention.

Full-signature replacement is valid for scalar ordered assembly because the
seven-coordinate law is closed and every representative is realized by an
actual strong marked tournament. Nonempty ports are preserved (`f>0`). Fresh
disjoint copies realize every representative word. This quotient does not
identify the physical names, port positions, local guards or native programs
of its preimages. Its word counts are representative-word counts, not counts
of all physically labelled words or all marked isomorphism classes at order
five.

All **484 ordered pairs of the 22 base atoms**, with up to eight actual
vertices, were independently assembled and counted by permutations/cuts.
Every coordinate of the product, including both joint cross terms, agrees.
All guarded pairs satisfy the stated fixed-f and fixed-target bounds.
The terminal constants are exactly `C=32/5` and `C=28/3` for the four- and
five-vertex banks, respectively.

## 3. Replays and accepted scope

Normal and optimized runs of the independent verifier reproduce **18,770
always-active gates** and byte-identical JSON. They include all inverse
comparisons for `f=1,...,100`, and the four-vertex queries `f=7*3^level` at
levels 0–4. The representative-word counts there are `2,5,9,18,29`; accepted
digit words are respectively the empty word, `2/4`, `22/42`, `422`, and `4222`.
Thus adjoining all four-vertex SCCs does not repair the first missing `24/44`
clock-uniform targets. This is an all-length statement within that declared
bank, not an obstruction to unrestricted actual rewrites.

The supplier's normal and optimized component runs were independently repeated
for both the order-four and order-five banks and compared with their frozen
receipts. The bounded scalar program was independently rerun in optimized
mode through 2000. It checks both inverses at every `f=1,...,1999`, evaluates
72,190 representative candidates, and directly recounts all 997 selected
non-singleton actual witnesses, their strongness, and root defects. The only
missing positive odd targets through 2000 in the order-five grammar are 7 and
21. A larger exploratory run is not part of this audit or its proof claims.

The one cosmetic minimum-order comparison was corrected before this final
audit: the candidate and stored orders now both include the added root.
The earlier expression merely overwrote ties; it never changed acceptance.

Reproduction commands, from the worktree root:

```
python3 04-computation/multiport_grammar_20261103/audit_component.py
python3 -O 04-computation/multiport_grammar_20261103/audit_component.py --output /tmp/audit_component_optimized.json
python3 04-computation/multiport_grammar_20261103/component_inverse.py --max-atom-order 4 --output /tmp/component_four.json
python3 04-computation/multiport_grammar_20261103/component_inverse.py --max-atom-order 5 --output /tmp/component_five.json
python3 -O 04-computation/multiport_grammar_20261103/spectrum_acceptance.py --bound 2000 --output /tmp/spectrum_audit.json
```

Frozen source/receipt pins at this audit:

| Artifact | SHA256 |
|---|---|
| `audit_component.py` | `d96549dccd67aba2c98f788f8ba3a9d0dfc81fc28b278172cc923666f04ed353` |
| `audit_component.json` | `b1719d7a7e787ef2947dfe8cddc4efb34fb2eb1eb5add80b1a78e2f9d3581d93` |
| `component_inverse.py` | `a101041d142c2c5821eca7cd91b1dec145fbc71c671f4315e505e538f29a2dc5` |
| `component_inverse.json` | `a4415b55139346e6facc7e611f5a0070966357510687fc34c3c433acceab4145` |
| `component_inverse_five.json` | `225e228394ac04faea87ec663fb1391ecb1256a27eae61b44cbb4d5db99978f6` |
| `spectrum_acceptance.py` | `87922b09f6883c259ec9c58c022a71302644afcc454a9c105044503f309cccd3` |
| `spectrum_acceptance.json` | `cb2c54346f46f719ee2d3fa946873a1c59c18c0d9f3c3fcdcc97f5c84c77027e` |

No mathematical correction remains in the audited inverse or implementation.
The missing future obligation is construction and closure of a physical
source-preserving return class, followed separately by any unbounded spectrum
surjectivity theorem.
