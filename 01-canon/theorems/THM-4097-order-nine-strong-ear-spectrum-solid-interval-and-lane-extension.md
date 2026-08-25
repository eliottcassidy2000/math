---
id: THM-4097
title: "Order-nine strong-ear spectrum, solid interval, cut-field carrier, and prime-lane extension"
status: >
  PROVED + CITED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. Moon's
  strong-subtournament theorem, already imported by THM-1370, makes every
  strong tournament of order at least four a nonconstant one-vertex ear of a
  strong parent. An exact subset-DP census of all 1,526,032 such ears from the
  6,008 strong order-eight classes independently reproduces the full
  1,482-value order-nine strong H-spectrum. It contains every odd value from
  85 through 2881. Balanced cut weight four misses exactly
  89,93,105,125; weight three supplies all four, so {3,4} is an exact
  two-weight basis and no singleton weight is exact. The exact response splits
  into a symmetric cut energy plus a zero-sum orientation field. The ordinary
  strong-prime lane is now furnished through 2879 and the seven-prime lane
  through 7*409. THM-4102/4104 later move the current targets to 80407 and
  7*11527=80689. Global completeness and all-order tiling remain OPEN.
source: >
  codex-frontier-synthesis-creative-20260825g;
  codex-padic-zeta-tournament-20260825
depends_on:
  - THM-001-redei
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-1862-order-join-reduction-principle
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
related:
  - THM-012b-insertion-decomposition
  - THM-1862-order-join-reduction-principle
  - THM-4051-tournament-order-seven-strong-base-exact-frontier
  - THM-4093-rational-edge-diagonal-gauge-and-padic-tournament-zeta-tangent
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - THM-4102-selected-order-ten-strong-ear-solid-interval
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - HYP-2879-strong-ear-atom-calculus
  - HYP-9029-strong-interval-tiling-law
  - MISTAKE-402
  - MISTAKE-505
references:
  - >
    J. W. Moon, "On Subtournaments of a Tournament," Canadian Mathematical
    Bulletin 9(3) (1966), 297--301, Theorem 2,
    https://doi.org/10.4153/CMB-1966-038-7.
script: 04-computation/tournament_order9_strong_ear_spectrum_thm4097.py
output: 05-knowledge/results/tournament_order9_strong_ear_spectrum_thm4097.out
data: 05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out
independent_audit_script: 04-computation/tournament_order9_strong_ear_spectrum_thm4097_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_order9_strong_ear_spectrum_thm4097_independent_audit.out
cut_field_audit_script: 04-computation/order_nine_strong_ear_cut_field_thm4097_independent_audit.py
cut_field_audit_output: 05-knowledge/results/order_nine_strong_ear_cut_field_thm4097_independent_audit.out
script_sha256: 610ca5850b272e0e75c574f2c1a710a0b96c75cc7191b1e1f1a03dfbdd1378d6
output_sha256: 0c3c9690ad5877d86480693af7ce97d8936c90e21b65677eefec72234c933dc0
independent_audit_script_sha256: b58de51efa200374e6014c8aeace4086fc7df1b2b73db36a01c0fcc4d2dd7943
independent_audit_output_sha256: 71058bc6b31ba26a59a58bb7f1e5366e767ed43db1f4bb2ddb80642b02acbfa6
certificate_atlas_script: 04-computation/order9_strong_ear_interval_certificate_codex_20260825.py
certificate_atlas_audit_script: 04-computation/order9_strong_ear_interval_independent_audit_codex_20260825.py
certificate_atlas: 05-knowledge/results/order9_strong_ear_interval_certificates_codex_20260825.tsv
certificate_atlas_output: 05-knowledge/results/order9_strong_ear_interval_certificate_codex_20260825.out
certificate_atlas_audit_output: 05-knowledge/results/order9_strong_ear_interval_independent_audit_codex_20260825.out
certificate_atlas_script_sha256: 8687aa30d3e282afe69cb80a5df67f146c5ffd98e0e87845a06969c589727877
certificate_atlas_audit_script_sha256: f8347ff0fb3e3412ab0b87802074da68e85c7feb589a676374f2c0592f81fd8d
certificate_atlas_sha256: 4aa1304da5c70e401c4953707820d377aeaf2418638dce747246c383309b3136
certificate_atlas_output_sha256: 1fb90b6e3c622c46cd99dfff7f8ced8312b740a6ac74d7b23e93d99082a3db12
certificate_atlas_audit_output_sha256: 7fe896cebf2e5641e07be6de56a0d25a0b952c9dd752fd5c7a93279aa625f975
engine: 04-computation/strong_H_spectrum_m8_isoclass_monad_s5.py
engine_sha256: 6ab922de4a8b6f6c15ee0ca7e0b036c3821b3e800dbdf961de72194e73346419
cut_field_audit_script_sha256: 1f0bd30d9a2e48f05243a6cbbcb2a3a203367a6e37cd456874c82cb631bb3775
cut_field_audit_output_sha256: 56efda1aa77391ba42530799dea95605b91fe22797cf855b00ac56c22b0443f2
historical_histogram: 05-knowledge/results/h_spectrum_n9_histogram_monad_s6.tsv
historical_histogram_sha256: e7d5594879d4c3af739cb94ca8cfd944879c4d586747d993dd6687e60126552f
data_sha256: 27fbef5b06fcf0369eeb602e513c3802ea171492e1292a3f6afa3efeadef9f55
semantic_sha256: b72cfe0e2c2ca2302f54933060299fc302cc9f241ef79154aead9524b426f0d5
independent_semantic_sha256: 7c2f9f5c9a586d1eed5b0329a294b7541c1c492682cb9ecb76c2eefc2d14b167
cut_field_semantic_sha256: de2b67a2c9bc0a349d33f7c9e996508a53116fbfe1e4764edc2e918005acd736
hash_basis: >
  raw LF bytes for files; canonical compact JSON for the compiler ledgers;
  the cut-field semantic hash is the ordered ALL,W1,...,W7 text ledger
audit: >
  PASS. The hash-pinned primary compiler evaluates all 1,526,032 ears and
  directly rechecks a retained child for every one of the 1,482 values; its
  independent companion cross-checks both older order-nine universes and
  literal boundary witnesses. A certificate atlas independently freezes one
  parent/cut/child witness for each of the 1,399 values in the solid interval.
  A third all-parent audit verifies all cut-field identities and integrality
  gates, every E_r spectrum and the exact {3,4} basis, the 623/735 hostile,
  target multiplicities, and both triangle/H firewalls. Normal/-O streams
  match the frozen LF artifacts after CRLF-to-LF normalization (MISTAKE-402).
---

# THM-4097 -- order-nine strong ears, the solid interval, and the cut field

**PROVED + CITED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.**

This theorem promotes the finite order-nine content of HYP-9029, closes the
strong-ear reducibility subproblem of HYP-2879, and advances both strong-prime
lanes isolated by THM-4094. It does **not** prove the all-order interval tiling
law or the global H-spectrum conjecture.

The main structural addition is an exact loss ledger. A one-vertex insertion
is not controlled by cut size, or even by an unoriented weighted cut alone.
It is controlled by a symmetric cut energy together with a zero-sum
orientation field.

## 1. Strong parents and nonconstant ears

For a tournament `T`, let `H(T)` be its number of directed Hamiltonian paths.
Let

```text
S_n^str={H(T):T is a strong tournament on n vertices}.       (1)
```

If `x` is a new vertex and `S` is a subset of `V(T)`, write `T+x_S` for the
child with

```text
x -> v  iff  v in S.                                      (2)
```

When `T` is strong, `T+x_S` is strong exactly when `S` is nonempty and proper.
Indeed those conditions give `x` both an out-neighbor and an in-neighbor in
`T`, and the strong parent supplies every remaining reachability. The two
constant cuts make `x` a source or sink and are not strong.

Moon's exact Theorem 2 says that the minimum number `s(n,k)` of induced strong
`k`-subtournaments in a strong `n`-tournament is `n-k+1`. This is the
classical strong-subtournament theorem already imported in the proof of
[THM-1370 -- H-spectrum omissions 7 and 21](THM-1370-h-spectrum-omits-7-21-all-n.md),
now with the source pinned to Moon 1966, Theorem 2. At `k=n-1` it supplies at
least two induced strong parents. Therefore:

> **Theorem 1.1 (all-order strong-ear reducibility; cited input).** Every
> strong tournament of order at least four is `T+x_S` for a strong parent `T`
> and a nonconstant cut `S`. Iterating reaches `C_3` through nested strong
> induced subtournaments.

Thus enumerating all strong parents and all nonconstant cuts is an object-
complete construction, not merely a sampled value generator. Only one of the
two guaranteed deletions is needed below.

## 2. The exposed-slot matrix by subset path DP

For `a,b in V(T)`, define:

```text
s_b = # Hamiltonian paths of T starting at b,
e_a = # Hamiltonian paths of T ending at a,
Q_ab = # old-vertex orders with a immediately before b
       and every adjacency except the a,b slot valid in T.   (3)
```

Then the insertion formula from HYP-2879 is

```text
H(T+x_S)
 = sum_(b in S) s_b
   + sum_(a notin S) e_a
   + sum_(a notin S,b in S) Q_ab.                          (4)
```

The third term includes both old valid slots and slots repaired by `x`.

There is a faster exact construction of `Q`. For `A subseteq V(T)`, let
`E_A(a)` count directed spanning paths of `T[A]` ending at `a`, and let
`B_A(a)` count those starting at `a`. Splitting an exposed order immediately
after `a` gives the bijection

```text
Q_ab = sum_(A: a in A, b notin A) E_A(a) B_(V\A)(b).       (5)
```

Both families are ordinary subset path DPs. Equation `(5)` computes the full
boundary state in `O(n^2 2^n)` arithmetic operations per parent instead of
enumerating `n!` orders.

## 3. Symmetric cut energy plus orientation field

Put

```text
w_ij = (Q_ij+Q_ji)/2,
h_i  = s_i-e_i+(col_i(Q)-row_i(Q))/2.                      (6)
```

The definitions a priori put the weights and field in `(1/2)Z`, but Rédei
parity makes every `w_ij` and `h_i` an integer. Put `F(S)=H(T+x_S)`. Every
value of `F` is odd, and its Boolean second difference is

```text
F({i,j})-F({i})-F({j})+F(empty)=-2w_ij.                  (6a)
```

The left side is even, so `w_ij` is integral. Then

```text
F({i})-F(empty)=h_i+sum_(j!=i)w_ij                       (6b)
```

shows that `h_i` is integral as well. Finally,
`sum_i s_i=sum_i e_i=H(T)` and total row sum equals total column sum, so

```text
sum_i h_i=0.                                               (7)
```

For a cut `S`, define

```text
cut_w(S)=sum_(i in S,j notin S) w_ij.                      (8)
```

> **Theorem 3.1 (cut-field identity).** For every tournament `T` and every
> insertion cut `S`,
>
> ```text
> H(T+x_S)=H(T)+cut_w(S)+sum_(i in S)h_i,                  (9)
> H(T+x_(V\S))=H(T)+cut_w(S)-sum_(i in S)h_i.             (10)
> ```

### Proof

Let

```text
C_in =sum_(a notin S,b in S)Q_ab,
C_out=sum_(a in S,b notin S)Q_ab.                          (11)
```

Then `cut_w(S)=(C_in+C_out)/2`, while internal row/column terms cancel and

```text
sum_(i in S)(col_i(Q)-row_i(Q))/2=(C_in-C_out)/2.          (12)
```

Subtract `H(T)=sum_i e_i` from `(4)` and combine `(11)`--`(12)` to get `(9)`.
Equation `(7)` and complement symmetry give `(10)`. QED.

The quotient

```text
(starts,ends,Q,oriented S) -> (H(T),w,{S,V\S})             (13)
```

preserves the parent, labelled symmetric weights, and the unoriented cut, but
destroys which complementary orientation has which response. To recover an
oriented response one must retain the orientation and the field sum
`sum_(i in S)h_i`; equivalently, retain the antisymmetric part of `Q` together
with `s-e`.

## 4. Exact order-nine universe and independent reconstruction

The parent tower has the exact nonisomorphic tournament counts

```text
1,1,2,4,12,56,456,6880          (orders 1 through 8),     (14)
```

matching `A000568`. Exactly `6008` order-eight classes are strong. For each,
the referee evaluates all `2^8-2=254` nonconstant cuts, hence

```text
6008*254=1,526,032 strong ears.                            (15)
```

The resulting value image agrees exactly with the frozen order-nine spectrum
from the independent `191536`-class canonical generator:

```text
|S_9^str|=1482,       min S_9^str=75,       max S_9^str=3357.   (16)
```

There are zero missing and zero extra values. The exact spectrum is the
ordered data file in the frontmatter; compactly,

```text
S_9^str
 ={75,81}
  union {85,87,...,2881}
  union U_9,                                                  (17)
|U_9|=81,       U_9 subset {2885,2887,...,3357}.              (18)
```

Here `(18)` specifies an ambient range, not full coverage; the exact maximal
odd intervals are printed by the referee. In particular the lower holes are
`77,79,83`, and the central band consists of exactly `1399` consecutive odd
values.

The three exact paths divide the verification burden as follows.

1. The primary compiler evaluates all `1,526,032` ears, retains one labelled
   child for every one of the `1,482` values, and rechecks each retained child
   by direct Held--Karp DP and strongness.
2. Its first independent companion compares both older complete order-nine
   universes and checks eight boundary witnesses by Held--Karp DP and literal
   permutation enumeration.
3. The separate cut-field audit recomputes all `1,526,032` ear rows, checks
   one direct child per parent, and verifies every cut-field identity. It also
   checks all `168224` symmetric-pair and `48064` field-entry integrality gates,
   every cut-energy integrality gate, all seven `E_r` spectra, and the semantic
   digest
   `de2b67a2c9bc0a349d33f7c9e996508a53116fbfe1e4764edc2e918005acd736`.

## 5. A sharp two-weight basis

For `1<=r<=7`, let

```text
E_r={H(T+x_S):T strong of order 8, |S|=r}.                (19)
```

Taking converses reverses every Hamiltonian path and replaces `S` by its
complement. Hence

```text
E_r=E_(8-r).                                               (20)
```

The exact cut-weight census is

| `r` | `|E_r|` | values of `S_9^str` missed |
|---:|---:|---:|
| 1,7 | 962 | 520 |
| 2,6 | 1294 | 188 |
| 3,5 | 1467 | 15 |
| 4 | 1478 | 4 |

The balanced self-converse weight is sharp:

```text
S_9^str \ E_4={89,93,105,125}.                            (21)
```

Every value in `(21)` belongs to `E_3`, and therefore

```text
S_9^str=E_3 union E_4.                                    (22)
```

No single cut weight has the full spectrum. Thus `{3,4}` is an exact
two-orbit basis, while weight four is the unique best singleton by image
cardinality. Direct positive controls for the four repairs are

| child `H` | parent `H` | weight-three cut mask |
|---:|---:|---:|
| 89 | 49 | `0x70` |
| 93 | 45 | `0x68` |
| 105 | 65 | `0x68` |
| 125 | 69 | `0x38` |

The labelled parent adjacency masks are frozen in the referee, and direct DP
checks every displayed child.

This clarifies the parity pattern behind the older `7->8` result. On an odd-
order parent the two middle weights are the same converse orbit and cannot
repair one another; at `7->8`, weights three and four both miss `49,75`. On
the even order-eight parent, the self-dual weight four and the neighboring
orbit `{3,5}` are genuinely distinct, and together they are exact.

## 6. The 623/735 orientation-field hostile

Take the strong order-eight parent with adjacency bitmasks

```text
[126,60,248,112,32,64,146,59],       H(T)=91.             (23)
```

For `S` with mask `105`, namely `S={0,3,5,6}`, the exact field is

```text
h=[288,84,121,-50,-282,-263,-31,133],
cut_w(S)=588,                 sum_(i in S)h_i=-56.         (24)
```

Both cuts are nonconstant and hence give strong children, but

```text
H(T+x_S)       =91+588-56=623,
H(T+x_(V\S))  =91+588+56=735.                             (25)
```

Thus the two orientations of the same parent and same unoriented bipartition
have different responses. Cut size and symmetric cut energy do not recover
the oriented response; the field sign is load-bearing.

## 7. Solid interval and the two strong-prime lanes

Equation `(17)` immediately promotes the HYP-9029 order-nine observation:

> **Corollary 7.1 (solid strong interval).** Every odd integer from `85`
> through `2881` is the Hamiltonian-path count of a strong order-nine
> tournament.

Combining exact strong spectra through order nine gives

```text
p in S_str  for every odd prime p!=7 with p<=2879,         (26)
7p in S_str for every odd prime p!=3 with p<=409.          (27)
```

These are exactly the two lanes of THM-4094 in a much longer finite prefix.
In particular, its former frontier targets are now explicit:

- `613` already has strong order-eight witnesses; two labelled adjacency rows
  are stored in the referee, and a balanced order-nine ear gives an additional
  control;
- `623=7*89` is the strong child in `(23)`--`(25)`.

The first ordinary-lane value not furnished by the exact spectra through
order nine is the prime `2887`. The first seven-prime-lane value not furnished
is

```text
7*419=2933.                                                (28)
```

These are **frontier targets, not asserted gaps**: either may occur at a
higher order. The upper order-nine tail also furnishes the sporadic ordinary
primes

```text
2903,2917,2963,2971,2999,3001,3011,3019,3023,
3037,3041,3049,3061,3067,3079,3083,3109,3299              (29)
```

and the sporadic seven-prime targets for `p=431,433,439,443`.

**Successor note.** THM-4102 realized `2887,2933`; THM-4104 then realized
`14657,14777` inside a selected order-eleven image. The current finite-data
lane targets are `80407` and `7*11527=80689`. Equations `(26)--(29)` remain
the exact order-at-most-nine boundary proved here.

Using THM-1862 order-join multiplicativity, the exact finite support now gives
every odd `m<=2885` except `7,21`. The next value not forced by these witnesses
and multiplication is the prime `2887`. This extends the proved coverage
prefix; it does not prove that `2887` is absent.

## 8. Triangle and p-adic graph-zeta firewall

The prime-lane valuation is not the triangle valuation controlling the first
term of the graph-zeta tangent in THM-4093. Two exact strong hostiles give both
failed directions:

```text
(H,c3)=(623,16):   7 divides H but not c3,
(H,c3)=(33,7):     7 divides c3 but not H.                (30)
```

The second witness has order-six adjacency masks
`[14,60,24,48,33,5]`. Therefore neither `7|H => 7|c3` nor its converse holds.
The lawful carrier for this open-path insertion problem is `(w,h)`; closed-
walk traces and the Bowen--Lanford zeta retain different information. No
Kubota--Leopoldt, irrationality, or LRC consequence is claimed.

## 9. What closes and what remains open

This theorem changes three frontier labels.

1. **PROVED from a cited classical input:** strong-ear reducibility holds at
   every order, so that subproblem of HYP-2879 is closed.
2. **PROVED + FINITE-EXACT:** the insertion formula, subset-DP `Q` compiler,
   cut-field decomposition, exact order-nine spectrum, central two-weight
   basis, solid interval, and lane prefixes above.
3. **OPEN:** prove overlapping solid intervals at every higher order, or
   otherwise construct every ordinary-prime and seven-prime strong atom.
   THM-4102/4104 compute selected order-ten/eleven images and move the current
   lane tests to `80407,80689`; the complete images remain uncomputed.

The observed central basis is a finite-order theorem, not an all-order
conjecture smuggled into the proof. A next computation should retain the full
`(w,h)` carrier and test central cut-weight **orbits**, not only raw weights.

## 10. Reproduction

Run from the repository root:

```text
python 04-computation/tournament_order9_strong_ear_spectrum_thm4097.py
python -O 04-computation/tournament_order9_strong_ear_spectrum_thm4097.py
python 04-computation/tournament_order9_strong_ear_spectrum_thm4097_independent_audit.py
python -O 04-computation/tournament_order9_strong_ear_spectrum_thm4097_independent_audit.py
python 04-computation/order_nine_strong_ear_cut_field_thm4097_independent_audit.py
python -O 04-computation/order_nine_strong_ear_cut_field_thm4097_independent_audit.py
```

After CRLF-to-LF transcript normalization (MISTAKE-402), each normal/optimized
pair must match its corresponding frozen output. The primary compiler matches
`05-knowledge/results/tournament_order9_strong_ear_spectrum_thm4097.out`
byte for byte. The older independent universe controls are
`04-computation/strong_H_spectrum_m9_isoclass_monad_s6.py` and
`05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out`.

The additional cut-field referee
`04-computation/order_nine_strong_ear_cut_field_thm4097_independent_audit.py`
does not import the primary THM-4097 program. It independently recomputes all
`1,526,032` ear rows, carrier integrality, the sorted support, target
multiplicities, every cut-cardinality spectrum `E_1,...,E_7`, the exact
`{3,4}` basis, complement hostile, and triangle firewall; LF-normalized normal
and optimized modes match its frozen output, including the semantic digest
`de2b67a2c9bc0a349d33f7c9e996508a53116fbfe1e4764edc2e918005acd736`.
