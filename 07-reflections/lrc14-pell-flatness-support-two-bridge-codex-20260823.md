# Support-two flatness deck versus Pell/conductor selectors

**Status:** hostile-audited research synthesis.  **PROVED** for the exact
integer compiler and its loss ledger below; **FINITE-EXACT** for the declared
ratio-deck, semigroup, mod-13, and Pell-prefix universes; **OPEN / RESERVED**
for every unpromoted THM-3742--3745 statement and for LRC(14).

## Verdict

There is one nontrivial exact bridge, but it is an address compiler rather
than an LRC implication.

If a square-triangular/conductor Pell depth has

```text
delta_d=binom(d,2)=q^2,
x=2d-1,                         s=2q,
```

then its half-Hadamard pair

```text
(a,b)=(x-s,x+s)=(P_(2k-1),P_(2k+1))                 (1)
```

is a primitive odd support-two ratio whose relation mass is

```text
a+b=2x=4d-2.                                          (2)
```

Consequently the `l1<=356` flatness ratio deck contains exactly the first
three positive compiler rows:

```text
k   d   q    (a,b)       a+b
1   2   1    (1,5)         6
2   9   6    (5,29)       34
3  50  35    (29,169)    198.                         (3)
```

The next row `(169,985)` has mass `1154` and is outside the deck.

This bridge does **not** force the flatness relation to be Pell-selected, add
rank to THM-2052, or imply bad loneliness.  The decisive destroyed
coordinates are:

- the support-owner labels and common speed scale when one keeps only the
  primitive ratio deck;
- midpoint parity, exact height, and Pell central sign/depth after reduction
  to the mod-13 projective observer; and
- the exact ratio itself if one keeps only an `l1` or conductor-degree scalar.

THM-778 restores the first loss only after common scale and labels are
reattached.  No present sidecar restores semantic owner/root/arrival.

## 1. Truth-status audit of the requested files

At `origin/main@84bca3ba8` all four incoming namespaces are still explicitly
outside the proof graph:

- `THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle` is
  `RESERVED / UNPROVED EMPTY STUB`;
- `THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction` is
  now `CITED + PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`;
- `THM-3744-pell-prefix-loneliness-constant-carry-exact-formula` is
  `RESERVED / UNPROVED EMPTY STUB`; and
- `THM-3745-monomial-plane-branch-conductor-triangular-pell-selector` is
  `RESERVED / UNPROVED EMPTY STUB`.

Therefore none is used as a proved dependency here.  The adjacent tracked
flatness script/reflection prove an audited candidate reduction but explicitly
leave THM-3743 reserved.  Its normal and optimized replays agree byte-for-byte
with the frozen output.  The inherited proved canon used here is:

- THM-3335, `square-triangular-pell-markov-pythagorean-selector`, for the
  integral Pell/half-Hadamard identities;
- THM-778, `centered-christoffel-endpoint-skew-product`, for exact pair-word
  scale/parity semantics; and
- THM-2052, `finite-height-forces-high-rank-bounded-relation-code`, for the
  rank-eleven relation space and rank-twelve terminal.

THM-3743 gives the following **cited + proved algebraic** necessary
condition: a hypothetical LRC(14) counterexample has some nonzero relation
`r in Z^13`, `r dot n=0`, with `||r||_1<=356`.  Choosing an `l1`-minimal row
makes it Graver.  If its support is two, positivity forces opposite signs and
its primitive absolute coefficients form

```text
D_356={(p,q): 1<=p<q, gcd(p,q)=1, p+q<=356}.           (4)
```

The exact cardinality is `19,314`.  If the minimal relation has support at
least three, there is no pair ratio and the bridge in this report stops at
once.

Typing is literal here: if the two physical speeds are `(gp,gq)`, their
primitive relation is `(q,-p)`.  After forgetting the two owner positions and
ordering, the absolute coefficient pair is the same unordered pair `{p,q}`
as the reduced speed ratio used in `(4)`.

## 2. Exact conductor-to-ratio compiler

Use the proved square-triangular Pell state from THM-3335:

```text
[x_(k+1)]   [3 4][x_k]
[s_(k+1)] = [2 3][s_k],             (x_0,s_0)=(1,0),   (5)
```

where `s_k=2q_k`,

```text
x_k^2-2s_k^2=1,
d_k=(x_k+1)/2,
binom(d_k,2)=q_k^2.                                    (6)
```

Apply the half-Hadamard transform

```text
a_k=x_k-s_k,                       b_k=x_k+s_k.         (7)
```

Then `(5)` becomes

```text
(a_k,b_k) -> (b_k,6b_k-a_k),                           (8)
```

so induction from `(1,1)` gives exactly

```text
(a_k,b_k)=(P_(2k-1),P_(2k+1)).                         (9)
```

Pell Cassini gives `gcd(a_k,b_k)=1`; both entries are odd.  Equations
`(6)--(7)` give `(2)` immediately.  This proves a lossless injection from a
positive conductor-square depth to an exact primitive ratio.  On its image,
the inverse scalar is still visible:

```text
d=(a+b+2)/4.                                           (10)
```

Thus the exact integer bridge does not lose conductor depth.  It loses LRC
placement only when the ratio is detached from its coordinate owners and
common speed scale.

Since `d_k` is strictly increasing and begins

```text
2,9,50,289,...,
```

the cap `(2)` gives precisely `(3)`.  This is the strongest positive result
of the session.

### Monomial conductor boundary

The elementary monomial slice of the reserved THM-3745 proposal can be proved
without the general `F` theorem.  For

```text
A_d=k[b^d,b^(d+1)] subset k[b],                         (11)
```

the value semigroup is `<d,d+1>`.  Its conductor is `d(d-1)`, its Frobenius
number is `d(d-1)-1`, and the number of gaps is

```text
delta(A_d)=d(d-1)/2.                                   (12)
```

Hence `delta(A_d)` is square exactly on the square-triangular degrees.  The
probe checks the exact gap set and sharp conductor for every `2<=d<=30`, and
checks

```text
{d<=356: binom(d,2) is square}={2,9,50,289}.            (13)
```

This does not prove the reserved general claim for `k[F(b),bF(b)]`, its
module decomposition, conductor equality, or squarefree multi-branch
interpretation.

## 3. The mod-13 central/projective quotient and its failure boundary

Reduce the matrix in `(5)` modulo `13`:

```text
H=[3 4;2 3] in SL_2(F_13).                              (14)
```

Direct multiplication gives

```text
H^7=-I,                         H^14=I,                 (15)
```

with no earlier occurrence.  Therefore the vector state has period `14` and
the half-angle projective point `[a_k:b_k]` has period `7`.  In the affine
chart `a/b`, its cycle from depth zero is

```text
1, 8, 6, infinity, 0, 11, 5.                           (16)
```

The central sign is not visible projectively.  More sharply, since
`d=(x+1)/2`, the lift at depth `k+7` satisfies

```text
d_(k+7)=1-d_k mod 13.                                  (17)
```

Thus the same projective half-angle point has two different conductor-degree
residues unless at a fixed residue of this affine involution.  The exact
fourteen residues are

```text
1,2,9,11,3,5,12,0,12,5,3,11,9,2.                      (18)
```

This identifies the central-sign coordinate destroyed by the THM-3742-style
projectivization.

The quotient is also far too coarse as a flatness selector.  Reducing all
`19,314` ratios in `(4)` to `P^1(F_13)`, the seven classes in `(16)` contain
`9,651` ratios and the complementary seven classes contain `9,663`.  Projective
Pell-orbit membership retains essentially half the deck, not the three exact
compiler rows `(3)`.

### Minimal parity hostile

The two exact deck ratios

```text
(1,5),                         (1,18)                   (19)
```

both map to projective class `8` modulo `13`.  THM-778 distinguishes them:
the first reduced pair is odd/odd and has a simultaneous midpoint block;
the second is odd/even and has none.  Therefore the mod-13 projective Pell
observer does not preserve even the pair tie predicate, before any semantic
arrival question is asked.

## 4. What the primitive flatness ratio itself forgets

The exact deck ratio is much better than its mod-13 shadow: it retains
coprime parity and the full centered-Christoffel word shape.  But THM-778 says
that a physical speed pair is

```text
(u,v)=g(p,q),                                           (20)
```

and its word is `g` repetitions of the reduced word.  For an odd/odd reduced
pair it has exactly `g` simultaneous blocks, at

```text
(2ell+1)/(2g),                    0<=ell<g.             (21)
```

The primitive ratio deck discards `g`.  The hostile

```text
(1,5):  one tie at 1/2,
(2,10): two ties at 1/4 and 3/4                         (22)
```

has identical reduced ratio and different physical wall schedules.  Owner
labels are an additional independent loss: an unordered ratio does not say
which two of the thirteen runner coordinates carry it.

The correctly typed connection contract is therefore

```text
source:       square-triangular/conductor Pell depth k
target:       primitive support-two flatness ratio deck D_356
map:          k -> {P_(2k-1),P_(2k+1)}
preserved:    primitive exact ratio, l1 mass 4d_k-2,
              odd/odd reduced tie existence, conductor depth on the image
destroyed:    coordinate owners and physical common scale g
sidecar:      ordered owner pair and g; then THM-778 reconstructs the word
cheap test:   cap intersection (3), followed by scale hostile (22).
```

After the further map to `P^1(F_13)`, exact ratio, parity, height, central
sign, and depth are also destroyed; `(19)` is the cheapest decisive test.

## 5. THM-2052: why this cannot be the missing rank increment

Let `W=W_(91^6,3)(n)` be THM-2052's span of all support-at-most-three,
height-at-most-`91^6` relations.  Every support-two flatness row has

```text
support=2,                       height<=||r||_1<=356<91^6.
```

Therefore it already belongs to `W` by definition.  In the rank-eleven
branch it cannot be the `a notin W` row that raises the relation code to rank
twelve.  A Pell-selected pair merely identifies a particularly structured
edge inside the existing code/cluster graph.  To use it on a two-anchor star,
one must retain the owner pair; the unlabeled ratio deck cannot say whether
the relation locks two leaves, a leaf to an anchor, or the two anchors.

This is a decisive no-bridge to the THM-2052 terminal:

```text
support-two flatness/Pell selector => inside-W refinement,
not a new independent relation.                                  (23)
```

A support-three row of height at most `91^6` is likewise one of the generators
and hence lies inside `W`.  A support-at-least-four row can lie either inside
or outside the **span** `W`; support alone does not decide span membership, so
that branch needs an exact code-membership test.  Neither higher-support branch
has a Pell pair ratio.

## 6. Pell-prefix hostile and THM-3744 boundary

The origin THM-3744 file is an empty reservation, so its proposed all-length
formula is not inherited.  The independent probe computes the exact maximum
for every prefix

```text
S_m={P_1,...,P_m},                         1<=m<=13.    (24)
```

The complete candidate set is rigorous: on each sawtooth cell the lower
envelope is piecewise linear, so a maximum occurs at a corner or pairwise
intersection, hence at `t=j/D` with

```text
D in {2P_i, P_i+P_j, |P_i-P_j|}.                       (25)
```

The exact finite values are

```text
m : M(S_m)
1 : 1/2       2 : 1/3       3 : 1/3
4 : 4/13      5 : 3/10      6 : 21/71
7 : 5/17      8 : 120/409   9 : 17/58
10: 697/2379  11: 29/99     12: 4060/13861
13: 99/338.                                                (26)
```

Each maximizer in `[0,1/2]` is uniquely `t=M(S_m)`.  These rows agree through
length thirteen with the candidate three-case formula

```text
M(S_m)=
  [P_m-P_(m-1)+1]/[2(P_m+1)]              if m is even;
  x_r/[2P_(2r+1)]                          if m=4r+1;
  P_(2r+1)/x_(r+1)                         if m=4r+3.   (27)
```

Equation `(27)` is **FINITE-EXACT through `m=13` only** in this report.  The
constant-carry all-length proof and endpoint uniqueness for every length
remain reserved THM-3744 obligations.

The thirteen-speed prefix is the decisive hostile:

```text
M({P_1,...,P_13})=99/338>1/14.                         (28)
```

It contains every pair `(a_k,b_k)=(P_(2k-1),P_(2k+1))` for `1<=k<=6`,
including all three compiler rows in `(3)`.  Hence a genuine thirteen-speed
row can carry the conductor/Pell ratios, their mod-13 orbit addresses, and
their THM-778 tie blocks while being very safely lonely.  Selector alignment
is not a sufficient LRC predicate.

This is not a claim that a hypothetical counterexample's minimum-width row is
one of those prefix pairs, nor that Pell prefixes classify any LRC stratum.

## 7. Rejected scalar bridge

A tempting but untyped identification is

```text
flatness relation mass p+q = conductor degree d.       (29)
```

It is not the lawful compiler `(2)`.  Even as a scalar filter it loses the
ratio: on the square-triangular shells at most `356`, the fibres in `(4)` have

```text
d=2: 0,             d=9: 3,             d=50: 10,
d=289: 136.                                             (30)
```

The zero at `d=2` comes from the excluded equal pair `(1,1)`.  Thus `(29)`
does not select a unique endpoint word or speed ratio.  The repaired bridge
is the half-Hadamard compiler with mass `4d-2`, not equality of two unrelated
scalar labels.

## 8. Reproduction and scope

Run

```powershell
python3 -B 04-computation/lrc14_pell_flatness_support2_bridge_20260823.py
python3 -B -O 04-computation/lrc14_pell_flatness_support2_bridge_20260823.py
```

The normal and optimized outputs agree with the frozen transcript.  The probe
checks:

- fifteen exact Pell/conductor compiler depths;
- the full `19,314`-ratio cap-`356` deck and all fourteen projective classes;
- the mod-13 central/projective periods and all fourteen conductor residues;
- sharp monomial conductors for `2<=d<=30` and every square delta degree
  through `356`;
- the parity, scale, and scalar-shell hostiles; and
- exact Pell-prefix loneliness for every length `1..13`.

Raw SHA-256 hashes are

```text
7691b9fa3857c0cd39adaf6a45e6bfeb58afe40919d8181eb13947789e4a64c5
  lrc14_pell_flatness_support2_bridge_20260823.py
29c2bba3455339912c6c591faf49028b2a34e5479675249022a05c9ff9d64095
  lrc14_pell_flatness_support2_bridge_20260823.out
```

**Strongest survivor:** the injective compiler `(1)--(3)` plus the typed
THM-778 scale/owner sidecar.

**Stopping reason:** flatness does not force its minimum relation into the
three-row Pell image; projective reduction loses even pair-tie parity; and any
support-two survivor already lies inside THM-2052's rank-eleven code.  A next
test would need the actual labelled two-anchor star spaces and ask whether a
specific owner pair has one of the three exact ratios.  Even a positive hit
would be an address refinement, not LRC(14).

No theorem namespace is promoted here, and no LRC(14) conclusion is claimed.
