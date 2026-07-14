---
id: THM-773
title: Prime-seven sheet monodromy and the exact tournament-atlas fibre
status: PROVED (elementary token/moment/holonomy theorem; tournament-atlas bridge finite-exact over all 5,040 owner assignments)
source: codex-2026-07-14-S6
depends_on:
  - THM-771   # exact seven-owner sheet defect and strict endpoint convention
  - HYP-6825  # canonical merged-metagraph addresses and inverse tiling fibres
related: [THM-754, HYP-3802, HYP-6835, HYP-6840]
verification: 04-computation/lrc14_prime7_sheet_monodromy_metagraph_codex_S6.py
  (+ 05-knowledge/results/lrc14_prime7_sheet_monodromy_metagraph_codex_S6.out and .json)
---

# THM-773 — Prime-seven sheet monodromy and tournament fibre

## 1. Exact prime-sheet tokenization

Let `W={w_0,...,w_6}` be seven positive integers, none divisible by `7`.
At `c=7` and threshold `1/14`, define the strict bad-sheet sets

```text
B_a(x)={k in F_7 : ||w_a(x+k)/7||<1/14}.
```

Put

```text
E_a={x in R/Z : w_a x in Z+1/2}.
```

If `x` is outside every `E_a`, let `n_a(x)=round(w_a x)`, with the unique
nearest integer.  Then `B_a(x)` is the singleton

```text
k_a(x)=-w_a^(-1)n_a(x)  (mod 7).                         (1)
```

If `x in E_a`, then `B_a(x)` is empty.  Thus away from the endpoint set the
seven-owner incidence deck is a labelled seven-token configuration

```text
k(x)=(k_0(x),...,k_6(x)) in F_7^7.                         (2)
```

The strict sheet cover is exact if and only if `a -> k_a(x)` is a bijection
onto `F_7`.

## 2. A complete finite-field moment certificate

For a token configuration put

```text
P_k(X)=product_a (X-k_a) in F_7[X],
S_m(k)=sum_a k_a^m       in F_7,       1<=m<=6.
```

The following are equivalent:

1. the owner tokens form an exact sheet partition;
2. `P_k(X)=X^7-X`;
3. the six moments are

   ```text
   (S_1,S_2,S_3,S_4,S_5,S_6)=(0,0,0,0,0,-1).             (3)
   ```

Consequently the six-moment vector decides the instantaneous exact-cover
predicate on the prime sheet.  The companion exhaustive audit checks all
`7^7=823,543` token configurations: exactly `7!=5,040` satisfy (3), with no
false positives or negatives.

## 3. Arbitrary owner counts and the eight-owner wall factorization

The prime-sheet polynomial also gives an exact carrier beyond seven owners.
For any number `r` of unramified owners at a non-event point, form

```text
D_k(X)=product_(a=1)^r (X-k_a).
```

Because `X^7-X` is the square-free product of all seven linear factors over
`F_7`, the sheet deck is fully covered if and only if

```text
                         X^7-X divides D_k(X).              (4)
```

At an endpoint, omit the event owners, whose strict bad sets are empty, and
apply (4) to the remaining-owner polynomial.  Therefore, with exactly eight
owners, a covered wall must be a **simple** event: two event owners would
leave at most six tokens.  At a simple event, coverage survives if and only if
the other seven owners form an exact permutation state.  Every covered
eight-owner wall consequently factors as

```text
(absent event owner) + (seven-owner heptagon stalk).         (5)
```

The seven-owner stalk maps to the same `n7-a267` fibre described below.  The
exact occupancy-profile audit checks all `C(14,6)=3,003` multisets of eight
tokens: precisely seven cover (one choice of duplicated residue), and these
are precisely the profiles whose polynomial is divisible by `X^7-X`.

This factorization resolves the structure of HYP-6835's published event
survivor

```text
W=(108,169,143,213,206,197,30,162),       x=19/216.
```

Owner `108` is the unique endpoint owner.  The remaining tokens are
`(6,5,3,1,4,2,0)`, an exact `F_7` state whose least-path atlas mask is `32153`
at node `n7-a267`.  Thus the `r=8` example is not outside the heptagon atlas;
it contains that atlas as its codimension-one event stalk.

## 4. Endpoint monodromy and holonomy

On crossing an event of owner `a` in the positive `x` direction,
`n_a` increases by one and hence

```text
k_a -> k_a-w_a^(-1)  (mod 7).                              (6)
```

For a simultaneous event, apply (6) to every event owner.  More generally,
after crossing `m_a` events of owner `a`,

```text
k_a' = k_a-m_a w_a^(-1)  (mod 7).                          (7)
```

If the initial and final chambers are both exact, their first moments vanish,
so necessarily

```text
sum_a m_a w_a^(-1)=0  (mod 7).                             (8)
```

This is the **prime-sheet first-moment holonomy**.  It is necessary but not
sufficient; after applying (7), the full six-moment test (3) is complete.
For every `r=1,...,6`, the exact higher update is

```text
Delta S_r=sum_a ((k_a-m_a w_a^(-1))^r-k_a^r).              (9)
```

At a wall, each event owner has an empty strict bad set.  If the chamber on
one side was exact, removing those owner tokens exposes their distinct old
sheets.  This is the algebraic form of THM-771's chamber-locking tear.

There is also a global carry which a circle-valued `x` coordinate hides:

```text
k_a(x+1)=k_a(x)-1  for every a.                           (10)
```

Thus one circuit of the base translates the whole sheet fibre.  Exactness is
invariant, but a marked sheet is not.  Any transport quotient which resets
`x=1` to `x=0` without (10) loses genuine monodromy.

## 5. Exact map to the merged tournament metagraph

Suppose the token assignment is exact.  There are two natural owner
tournaments.

**Marked-cut gauge.** Orient `a -> b` when `k_a<k_b` in the marked order
`0<...<6`.  This is always transitive, so all `5,040` assignments map to the
HYP-6825 root node

```text
n7-a000.
```

**Circular gauge.** Orient

```text
a -> b  iff  k_b-k_a in {1,2,3}.                          (11)
```

Relabelling owner `a` by its sheet `k_a` identifies (11) with the rotational
heptagon tournament `R_7=C_7(1,2,3)`.  Thus every assignment maps to the one
merged node

```text
n7-a267.
```

After putting owners in sheet order, the two gauges disagree precisely on
the long chords of lengths `4,5,6`: `3+2+1=6` edge flips.  The exact atlas
accordingly places `n7-a267` at local depth six, with local address word
`SSRRSS`, blue/black root distance six, and blue/black word `BBBBBB`.

There are `7!=5,040` labelled owner assignments, `6!=720` rotation orbits,
and `360` dihedral orbits.  The circular tournament has `175` Hamiltonian
paths and automorphism group of order `7`, so its fixed-path staircase fibre
has

```text
175/7=25                                                      (12)
```

masks.  The exact audit chooses the lexicographically least Hamiltonian path
of each labelled circular tournament, relabels it to the explorer's fixed
path, and encodes its staircase mask.  The resulting map from all `5,040`
assignments is surjective onto **exactly the 25 masks stored at `n7-a267`**.
The full assignment-to-mask/fibre-index map is in the compact JSON; its SHA256
is `ed10c0914976a7f47b6ff4ef103052e2874964fd51967b17ecf8a852e3d1344e`.

This answers the finite mapping question exactly, but it also proves a
limitation: an ordinary tournament isomorphism-class node is not an index of
owner-labelled sheet tilings.  It is a base point with a nontrivial stalk.

## 6. The continuation-equivalence separation

Even `(node, fixed-path mask, instantaneous owner assignment)` is not a
continuation-complete state.  The exact audit finds the shared assignment

```text
(0,3,2,5,1,4,6),       least-path mask 27833,
```

in both rows

```text
W_1=(1,2,3,4,5,8,10),   chamber (5/16,7/20),
W_2=(1,2,3,4,5,8,11),   chamber (7/22,3/8).
```

The next future of `W_1` is `(event owner 10, free sheet 6)`, whereas the next
future of `W_2` is `(event owner 4, free sheet 5)`.  They have the same two
isomorphism-class nodes, the same six moments, the same owner assignment, and
the same canonical fibre mask, but different legal next outputs.  Hence they
are separated by the continuation/Nerode test proposed in HYP-6825.

The smallest faithful state for this prime-sheet movie must therefore retain,
at minimum,

```text
(owner-to-sheet assignment, inverse steps w_a^(-1), endpoint order/phases,
 global sheet carry),                                           (13)
```

over the core-safe metric base.  Moments and tournament nodes are exact
checksums of parts of (11), not replacements for it.

## Proof

For (1), write `n=round(wx)`.  The choice
`k=-w^(-1)n mod 7` makes `wk+n` divisible by seven, and

```text
|w(x+k)-7q|=|wx-n|<1/2
```

for the corresponding integer `q`.  Every other residue differs from a
multiple of seven by a nonzero integer, so its distance is greater than
`1/2`.  At `wx=n+1/2`, the two adjacent candidates lie at equality `1/2` and
the strict inequality excludes both.  This proves the token formula and the
endpoint statement.  Seven singleton bad sets cover seven sheets exactly iff
their labels are a permutation.

The polynomial equivalence is the factorization
`product_(u in F_7)(X-u)=X^7-X`.  A permutation gives (3) by the standard
finite-field power sums.  Conversely, Newton's identities are invertible in
degrees `1,...,6` over `F_7`.  From (3) they give

```text
e_1=...=e_5=0,       e_6=-1.
```

Hence `P_k(X)=X^7-X-e_7`.  Every token `k_a` lies in `F_7`, so
`k_a^7-k_a=0`; evaluating `P_k(k_a)=0` forces `e_7=0`.  Therefore
`P_k=X^7-X`, including multiplicities, and the tokens are exactly `F_7`.

For an arbitrary number of owners, `X^7-X` is square-free with roots exactly
the elements of `F_7`.  It divides `D_k` exactly when every one of those roots
occurs among the owner tokens, which is exactly full sheet coverage.  At a
wall the event-owner factors are absent.  The eight-owner conclusions follow
by counting the remaining factors; the displayed survivor and mask are
verified exactly by the companion computation.

Equations (6)--(10) follow directly from
`n_a(x)=floor(w_a x+1/2)`.  Summing (7) proves (8), and binomial expansion
gives (9).  The tournament identifications follow by relabelling owners with
their sheet positions.  The six-flip count is the number of unordered pairs
at linear separations `4,5,6`.  Formula (12) is Hamiltonian paths modulo the
free automorphism action; the remaining atlas and continuation statements are
the complete finite enumerations recorded by the verification artifact. ∎

## Tournament Analysis and preservation audit

- **Vertices:** owners for the two displayed gauges.  Sheets, endpoint events,
  fibre masks, and proof obligations were also considered; none alone is a
  faithful continuation carrier.
- **Pairwise observable:** owner precedence in the marked or circular sheet
  order.
- **Switch/gauge:** cut precedence versus three-step circular precedence.
- **Tie Hamiltonian path:** owner order `0<...<6`; unused because both relations
  are tournaments.
- **Cut fingerprint:** score histogram `{0:1,...,6:1}`, zero directed
  triangles, seven singleton SCCs, one Hamiltonian path.
- **Circular fingerprint:** score histogram `{3:7}`, fourteen directed
  triangles, one seven-vertex SCC, 175 Hamiltonian paths.
- **Edge flips:** six.
- **Preserved by ordinary isomorphism class:** only the unlabelled transitive
  or heptagon shape.
- **Destroyed:** owner assignment, inverse winding, endpoint phase/order,
  next free sheet, and global carry—the fields that decide future LRC
  witnesses.

This challenges the convention that the merged-metagraph node is the whole
object.  In the LRC application it is the base of an owner-labelled transport
stalk.
