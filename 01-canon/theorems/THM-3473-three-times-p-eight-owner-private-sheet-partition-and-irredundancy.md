---
id: THM-3473
title: "Three-times-p eight-owner private-sheet partition and irredundancy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For p=14k-1, the
  THM-3469 eight-owner cover has exactly eleven Boolean sheet-support types:
  eight positive singleton atoms, one
  two-owner repair packet, and two four-owner backbone packets.  Exact
  private counts prove representation-relative irredundancy for every k>=1;
  the residue-deficit automaton has period three and its Berggren U-spine
  pullback has minimal period 63.  Global rank minimality and LRC(14) do not
  follow.
source: codex-2026-08-15-private-sheet-partition
audit: >
  self-contained nearest-p chart, strict support classification, residue
  counts, Boolean deletion derivatives, multiplicity spectrum, sparse marker
  polynomial, bidirected two-section, harmonic lanes, and period-63 U-spine
  pullback; exact 5259000-sheet / 42072000 owner-incidence and rational-route
  gates; independent clean-room derivation of all eleven supports, strict
  repair endpoints, private and multiplicity formulae, graph quotient, and
  minimal deficit word; independent 21018000-sheet private-word and
  1379964-cell endpoint audits; dependency, script, output, semantic, ID,
  security, and normal/optimized/stored replay gates
depends_on:
  - THM-3469-three-times-p-half-twist-eight-owner-cover-boundary
related:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum
  - THM-3472-odd-modulus-zero-half-conjugacy-and-global-zmc-rank-equality
script: 04-computation/lrc_three_p_private_sheet_boolean_atlas_thm3473.py
output: 05-knowledge/results/lrc_three_p_private_sheet_boolean_atlas_thm3473.out
script_sha256: 8831f932deb61075e754887f570390388f09ed6bd84092bd4c47dc50d7a588be
output_sha256: 7951215e40dfd8af18984beba39fb5cce6c1c25053894ecbbacbde7158af7fab
semantic_sha256: 1a269572f813dcfafa4d99b9dda39e6fd36514005ff671ef7eec3f6416467fcc
hash_basis: LF-normalized bytes
---

# THM-3473 -- three-times-p private sheets and the Boolean support atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The universal support classification, immutable package, and independent
hostile rederivations pass their stated gates.

## 1. The affine chart and three kinds of cells

Let `k>=1`, and put

```text
p=14k-1,        q=3p=42k-3,
h=6k-1,         A=7k-1.                                (1)
```

Use the THM-3469 owner order

```text
U=(1,p-1,p+1,2p-1,2p,2p+1,3p-6,3p-1)=(u_1,...,u_8).  (2)
```

For a sheet `ell in Z/qZ`, write `x=2ell+1` modulo `2q=6p`.  Since `p` is
odd, there is a unique nearest-`p` chart

```text
x=mp+y  (mod 6p),       m in Z/6Z,       |y|<=A.        (3)
```

Indeed the six intervals of `p=2A+1` consecutive integers centred at `mp`
partition `Z/6pZ`.  There is no tie.

For `sigma in {+1,-1}`, interpreted modulo six, define

```text
C_(m,sigma)={ell=(mp+y-1)/2 mod q:
              |y|<=h, mp+y==sigma (mod 6)},

D_m={ell=(mp+y-1)/2 mod q:
       |y|<=h, mp+y==3 (mod 6)},

G_(m,sigma)={ell=(mp+y-1)/2 mod q:
              h<|y|<=A, mp+y==sigma (mod 6)}.          (4)
```

The congruences in `(4)` make every displayed numerator odd.  The `C` cells
are the affine core, the `D` cells are the order-three backbone core, and the
`G` cells are the septimal repair gap.

For owner `u_i`, let `B_i` be its literal half-twist mask and let

```text
P_i=B_i minus union_(j!=i) B_j                          (5)
```

be its private-sheet set.

## 2. Exact private-sheet theorem

The private sets are exactly

```text
P_1 =C_(0,+1) union C_(0,-1),
P_2 =C_(1,+1) union C_(5,-1),
P_3 =C_(1,-1) union C_(5,+1),
P_4 =C_(2,+1) union C_(4,-1),
P_5 =D_1 union D_2 union D_4 union D_5,
P_6 =C_(2,-1) union C_(4,+1),
P_7 =union_(m,sigma) G_(m,sigma),
P_8 =C_(3,+1) union C_(3,-1).                          (6)
```

Write `e_r(k)=1` when `k==r (mod 3)` and zero otherwise.  In the owner order
`(2)`, their exact cardinalities are

```text
(|P_1|,...,|P_8|)=
 (4k,
  4k-2e_1,
  4k-2e_0,
  4k,
  8k-2e_2,
  4k,
  4k-2e_2,
  4k).                                                  (7)
```

Every entry in `(7)` is positive for every `k>=1`.  Therefore deleting any
one owner from `(2)` uncovers exactly `P_i`, and the eight-owner presentation
is irredundant for every member of the family.

This is presentation-relative, not a statement that the global cover rank is
eight.  In fact THM-3469 proves global rank four whenever `k==2 (mod 3)`,
while `(7)` still makes all eight owners essential inside this particular
presentation.  A different four-owner cover exists; it is not obtained by
deleting four owners from `(2)`.  Thus numerical grade forgets cover ancestry
and the cover system has no basis-exchange property of the kind a matroid
would supply.

## 3. Proof of the support classification

For `u=ap+b` with `b in {+1,-1}`, `(3)` gives

```text
u x == (a sigma+b m)p+b y  (mod 6p),                  (8)
```

where `sigma=x mod 6`.  When `sigma=+1` or `-1`, the twelve THM-3469
channels make the coefficient of `p` vanish for exactly one owner:

| `m` | `sigma=+1` | `sigma=-1` |
|---:|---|---|
| 0 | `u_1` | `u_1` |
| 1 | `u_2` | `u_3` |
| 2 | `u_4` | `u_6` |
| 3 | `u_8` | `u_8` |
| 4 | `u_6` | `u_4` |
| 5 | `u_3` | `u_2` |

That owner's cyclic numerator distance is `|y|`.  Every nonselected
`ap+/-1` owner has distance at least

```text
p-|y| >= p-A=7k>h.                                    (9)
```

The backbone owner `2p` covers exactly `sigma=3`.  For every odd `x`, the
repair owner satisfies

```text
dist_(6p)((3p-6)x,0)=3p-6|y|.                         (10)
```

At the core boundary this is `6k+3=h+4`, so the repair is absent.  At the
first gap point it is `6k-3<=h`, and it remains at most `h` through `A`.
Consequently every nonbackbone core/gap cell has singleton support exactly as
listed in `(6)`.

It remains to classify `sigma=3`.  In the core, `(8)` gives

```text
m=0: support {u_1,u_4,u_5,u_6},
m=3: support {u_2,u_3,u_5,u_8},
m in {1,2,4,5}: support {u_5}.                         (11)
```

In the gap, `(9)` excludes every `ap+/-1` owner while `(10)` and the
backbone both fire, giving support `{u_5,u_7}`.  Equations `(6)` and `(11)`
therefore list the complete support atlas; in particular there are no hidden
triple, quintuple, or higher intersections.

## 4. Counts and the eleven Boolean atoms

The elementary residue lemma behind all boundary anomalies is

```text
#{y in [-h,h]:y==0 (mod 6)}=2k-1,
#{y in [-h,h]:y==r (mod 6)}=2k  for r!=0.             (12)
```

Also

```text
p mod 6 = -1,+1,3  as k mod 3 = 0,1,2.               (13)
```

Applying `(12)` to the twelve channel cells proves the first, second, third,
fourth, sixth, and eighth entries of `(7)`.  Applying it to
`D_1,D_2,D_4,D_5` proves the fifth.  There are `2p=28k-2` odd nonmultiples of
three modulo `6p`; subtracting the affine-core census gives the seventh.

Let `I(ell) in {0,1}^8` be the owner-incidence vector.  The only eleven
supports and their counts are

| support | count |
|---|---:|
| `{u_i}` | `|P_i|`, for `1<=i<=8` |
| `{u_1,u_4,u_5,u_6}` | `2k` |
| `{u_2,u_3,u_5,u_8}` | `2k-1` |
| `{u_5,u_7}` | `2k+2e_2` |

Equivalently the complete multiplicity word has alphabet `{1,2,4}` and

```text
#{ell:|I(ell)|=1}=36k-2-2e_2,
#{ell:|I(ell)|=2}= 2k+2e_2,
#{ell:|I(ell)|=4}= 4k-1.                              (14)
```

These sum to `q=42k-3`.  The Boolean deletion derivative

```text
partial_i(ell)=I_i(ell) product_(j!=i)(1-I_j(ell))     (15)
```

is exactly the indicator of `P_i`.  Thus `(7)` is a full Boolean
realization of irredundancy, not merely eight selected witnesses.

Equivalently, the squarefree support enumerator

```text
H_k(z_1,...,z_8)=sum_ell product_(i:I_i(ell)=1) z_i   (15a)
```

has the exact sparse form

```text
H_k=sum_i |P_i|z_i
    +2k z_1z_4z_5z_6
    +(2k-1)z_2z_3z_5z_8
    +(2k+2e_2)z_5z_7.                                (15b)
```

Thus the private counts are the first homogeneous part, while the entire
multiple-cover locus is carried by only three higher monomials.  Under the
common marker `z_i=t`, `(15b)` becomes

```text
(36k-2-2e_2)t+(2k+2e_2)t^2+(4k-1)t^4.               (15c)
```

This marker polynomial is a compact loss-controlled carrier: its first jet
at zero sees every private sheet, unlike the pairwise quotient below.

As `k` tends to infinity, singleton sheets have limiting proportion `6/7`.
Each of the three multiple-support packets has limiting proportion `1/21`,
so the entire overlap locus has limiting proportion `1/7`.  The same
septimal threshold that built the cover therefore reappears exactly as its
asymptotic overlap mass.

## 5. The generalized tournament is two tetrahedra and one edge

There is an intrinsic pairwise observable, but no intrinsic orientation.
Take the eight owners as vertices and put both arcs `i->j` and `j->i` exactly
when the owners co-occur on some sheet.  Missing pairs receive no arc.  The
resulting bidirected two-section is

```text
K_4({1,4,5,6}) union K_4({2,3,5,8}) union K_2({5,7}). (16)
```

The two tetrahedra and the repair edge meet only at the backbone vertex
`5`, corresponding to owner `2p`.  There are `13` undirected edges, `26`
directed arcs, and degree sequence

```text
(3,3,3,3,7,3,1,3).                                   (17)
```

Thus a tournament with missing and both-way edges is an exact quotient of
this family.  What it preserves is pairwise coactivity existence.  It loses
sheet locations, all counts, the eight singleton atoms, and the `k mod 3`
state.  In particular `(16)` alone cannot prove irredundancy because private
supports are invisible to every edge.  The needed sidecar is the eleven-row
support/count table above.  No arbitrary orientation is introduced.

## 6. Ternary and harmonic pullbacks

The deficit pattern in `(7)` is a three-state automaton:

```text
k==0: the p+1 channel loses two private sheets;
k==1: the p-1 channel loses two private sheets;
k==2: the backbone and repair each lose two private sheets.       (18)
```

For each state `r`, the subset `{k in N:k==r (mod 3)}` has both natural and
harmonic coefficient `1/3`:

```text
sum_(k<=N,k==r mod 3) 1/k = (1/3) log N+O(1).          (19)
```

This is the precise subset-of-the-harmonic-series realization attached to
the Boolean atlas.

On the Berggren U-spine `q_t=(2t+1)^2+2`, THM-3469 puts the affine family on

```text
t mod 21 in {5,8,12,15}.                              (20)
```

Colour a lane point by `k mod 3` and every off-lane point by a fourth symbol.
The resulting word has minimal period `63`.  In one period the off-lane
symbol occurs `51` times and each ternary colour occurs `4` times.  A shift
by `21` changes `(t,state)=(5,0)` to `(26,1)`, while a shift by `9` changes
`(3,off)` to `(12,0)`; these exclude both maximal proper divisors of `63`.
Hence each colour has ambient natural/harmonic coefficient `4/63` and
conditioned lane coefficient `1/3`.

The period-63 colouring records a residue sidecar on a Fibonacci/Berggren
branch.  It does not turn the branch recursion into a three-child tree, and
it does not identify cover ancestry with the global rank word.

## 7. Scope boundary and exact companion

This theorem classifies one literal half-twist presentation.  It produces no
endpoint current, H^1 class, bispectrum, physical row, decrement, Jacobian
counterexample, or LRC(14) conclusion.  LRC(14) remains open.

Run from the repository root:

```bash
python 04-computation/lrc_three_p_private_sheet_boolean_atlas_thm3473.py
python -O 04-computation/lrc_three_p_private_sheet_boolean_atlas_thm3473.py
```

The standard-library companion compares the nearest-`p` support classifier
with direct strict cyclic-distance masks on all `5,259,000` sheets for
`1<=k<=500`, making `42,072,000` owner-incidence tests.  It checks the residue
lemma through `k=10000`, an independent rational-mask route through `k=60`,
all eleven support/count identities, every deletion deficit, the bidirected
two-section, the period-63 U-spine word, security and dependency gates, and a
frozen semantic digest.  It uses explicit exceptions under `-O` and performs
no file write, dynamic evaluation, subprocess, or network action.
