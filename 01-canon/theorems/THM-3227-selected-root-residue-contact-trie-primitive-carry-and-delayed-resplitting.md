---
id: THM-3227
title: "Selected-root residue contact trie, primitive carry, and delayed resplitting"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Derivative-normalized jets of a finite family at one selected common simple
  root form a coordinate-invariant residue contact ultrametric.  Its prefix
  partitions are a rooted q-ary trie over a residue field of size q; every
  split carries an affine weighted-cotangent label, and the exact label census
  is leaves minus one.  Each nonzero edge has additive order p, the residue
  characteristic.  Over an unramified p-adic lift, p-divisible lower jets do
  not contaminate the first residue-visible divided carry.  Reduction can
  postpone and recreate splits, and ramified pi-small lower jets can destroy
  p-divisibility; both boundaries are explicit.  No global root/owner selector
  or physical carrier is supplied.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent exact companion pins promoted THM-3221; checks
  44,229 contact-ultrametric triples, 40 sharp full-trie nodes over alphabets
  of sizes 2,3,4,5,9, the leaves-minus-one affine census, 18 extension-field
  additive orders, 60 nonlinear coordinate changes, 273 primitive divided
  carries with p-divisible lower jets, 15 exact ramified-tail formulas, and
  the delayed depth-two to depth-three/four resplitting hostile.  Normal,
  optimized, and stored transcripts agree with the LF-normalized hashes below.
  An independent immutable audit rederived the contact ultrametric and affine
  leading-difference law, checked the q-branch versus p-order distinction,
  reconstructed the unramified divided-carry factorization, and verified the
  delayed-resplitting and ramified counterexamples.  It replayed both normal
  and optimized companions against the stored transcript and accepted the
  pinned LF-normalized hashes.
depends_on:
  - THM-3221-selected-root-osculating-separation-and-minimal-jet-prime-carry
related:
  - THM-2020-gmc2-finite-place-channel-separation
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
  - THM-3086-arbitrary-cluster-composition-chambers-and-alternant-clutch-holotopy
script: 04-computation/selected_root_residue_contact_trie_carry_thm3227.py
output: 05-knowledge/results/selected_root_residue_contact_trie_carry_thm3227.out
script_sha256: b3768ccc2e43b241feea7405869b1e6cf111ffedde6fbe59f0c2fe6be4c79e63
output_sha256: f80673e593adbd6bdb5efed3e96bf634652f9234bcdd0fcc8f5a3316b98e1b46
hash_basis: LF-normalized bytes
---

# THM-3227 -- selected-root residue contact trie, primitive carry, and delayed resplitting

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3221 canonically isolates the first live osculating coefficient of one
selected-root transition.  A finite family has more structure than a bag of
pairwise coefficients: agreement through successive jet depths gives a
rooted contact tree.  At a finite place this tree is a literal coefficient
prefix trie.  Its branching is controlled by the residue-field size `q`,
while every charged edge resets with the residue characteristic `p`.

The distinction `q` versus `p`, and the failure of naive tree contraction
under reduction, are both load bearing.

## 1. Residue coefficient words

Let `R` be a DVR with maximal ideal `m_R`, finite residue field

```text
k=R/m_R,       |k|=q,       char(k)=p.                    (1)
```

Fix a source point `a`, and let `f_i` be a finite family of polynomial or
formal germs with

```text
f_i(a)=0,       f_i'(a) in R^*.                           (2)
```

Normalize as in THM-3221:

```text
phi_i(u)=f_i(a+u)/f_i'(a)
        =u+sum_(n>=2) a_(i,n)u^n.                         (3)
```

For a fixed depth `D`, retain the reduced coefficient word

```text
w_i=(abar_(i,2),...,abar_(i,D)) in k^(D-1).               (4)
```

Coincident words are one leaf.  For distinct leaves define their contact
depth

```text
nu(i,j)=min{n in {2,...,D}: abar_(i,n)!=abar_(j,n)},       (5)
```

and put `nu(i,j)=infinity` for equal words.

## 2. Contact is an ultrametric

For every triple,

```text
nu(i,k)>=min(nu(i,j),nu(j,k)).                            (6)
```

If the two depths on the right are unequal, equality holds.  This is the
coefficientwise nonarchimedean law: below the minimum both differences
vanish, and at a uniquely smaller depth its nonzero leading coefficient
cannot cancel against a difference which is still zero.

Equivalently, for every `n` the relation

```text
i ~_n j  iff  abar_(i,r)=abar_(j,r) for 2<=r<=n           (7)
```

is an equivalence relation, and the partitions refine as `n` increases.
Their Hasse diagram after deleting unary repetitions is a rooted contact
trie.

## 3. Nonlinear coordinate invariance

Let `v=psi(u)` be an integral source-coordinate change with unit linear term
`rho=psi'(0)`.  THM-3221 gives

```text
phi_i^v=L_rho composed phi_i^u composed chi,
chi=psi^(-1),       L_rho(u)=rho u.                       (8)
```

Suppose a block agrees through depth `n-1`.  At depth `n`, every transformed
coefficient has the form

```text
abar_(i,n)^v = common_n + rhobar^(1-n) abar_(i,n)^u.      (9)
```

The term `common_n` depends on the common lower prefix and on the chart, but
not on `i`.  Therefore contact depths are unchanged, and child-label
differences transform by the exact weighted-cotangent factor

```text
delta_n^v=rhobar^(1-n)delta_n^u.                          (10)
```

Thus the rooted trie and its affine labels are intrinsic to the supplied
simple-root branch.

## 4. Sharp q-ary branching and the affine-label census

A block which first splits at depth `n` has one child for each occurring
value of `abar_(i,n)`.  Hence

```text
number of children <=q.                                  (11)
```

The bound is sharp: all `q` coefficient values occur among the normalized
polynomials `u+alpha u^n`, `alpha in k`.  More globally,

```text
number of distinct normalized D-jets <=q^(D-1),          (12)
```

and the full coefficient-word bank attains equality.

At a split node `v` of depth `n`, choose one base child and record the
differences of every other child coefficient from it.  If `d(v)` is the
number of children, this uses `d(v)-1` elements of the weighted line
`(T_a^*)^(n-1)`.  For `N` distinct leaves,

```text
sum_(split nodes v)(d(v)-1)=N-1.                          (13)
```

This is the standard rooted-tree Euler count.  It is also an exact decoder:
for two leaves, take their lowest common ancestor and subtract the two child
labels there.  The result is their first transition coefficient.  Hence the
trie plus `N-1` affine labels recovers all pairwise first-separation defects.
Cycle sums inside any one split vanish; any nontrivial global holonomy must
come from branch gluing outside this rooted affine trie.

## 5. Edge order is p, not q

Let two children first separate at depth `n`, with nonzero label difference
`delta in k^*`.  By THM-3221 their reduced transition has the form

```text
Hbar(u)=u+delta u^n mod u^(n+1).                          (14)
```

The minimal quotient is the additive group of `k`.  Therefore

```text
Hbar^r !=identity for 1<=r<p,
Hbar^p =identity                    mod u^(n+1).           (15)
```

Every nonzero additive element of a characteristic-`p` field has order `p`,
even when `q=p^e` is larger.  Thus `q` controls the maximum number of sibling
branches, while `p` controls the reset length.

## 6. Primitive divided carry with p-divisible lower jets

Assume now that `R` is unramified over `p`, so `m_R=pR`.  More generally, the
argument below needs only that the displayed lower coefficients lie in
`pR` and that division by `p` is lawful.

Let an integral transition through depth `n` be

```text
H(u)=u+h_2u^2+...+h_n u^n mod u^(n+1),
h_j in pR for 2<=j<n,       h_n in R^*.                  (16)
```

Then the residue transition first separates at depth `n`, and

```text
(1/p)[u^n](H^p-u) == h_n mod pR.                         (17)
```

In particular the divided coefficient is a unit.

To prove `(17)`, work in the `n`-jet group.  Factor

```text
H=L composed G,
L=u+h_2u^2+...+h_(n-1)u^(n-1),
G=u+h_n u^n.                                              (18)
```

The top layer `F^n` is central in this jet group, so

```text
H^p=L^p composed G^p.                                    (19)
```

The coefficient `[u^n](L^p-u)` is an integral polynomial in the lower
`h_j`.  Its linear part is `p` times the `u^n` coefficient of `L`, hence
zero; every remaining monomial has degree at least two and lies in `p^2R`.
Meanwhile

```text
G^p=u+p h_n u^n mod u^(n+1).                             (20)
```

Dividing the `u^n` coefficient of `(19)` by `p` and reducing modulo `pR`
gives `(17)`.  Thus earlier p-divisible jets do not contaminate the first
residue-visible carry.

## 7. Two sharp reduction boundaries

### 7.1 Reduction is not naive graph contraction

Over `Z_3`, take the three normalized coefficient words in depths two through
four

```text
(0,0,0),       (3,0,1),       (6,1,0).                   (21)
```

Over the fraction field all three pairs first split at depth two.  Modulo
three they become

```text
(0,0,0),       (0,0,1),       (0,1,0).                   (22)
```

The last word separates at depth three, while the first two separate at
depth four.  Thus reduction can postpone an early split and create nested
internal nodes which did not occur in the generic contact tree.  The residue
trie must be rebuilt from the full reduced coefficient words; it is not
recoverable by merely contracting nonprimitive generic edges.

### 7.2 Ramified pi-smallness does not imply p-carry integrality

Let `p=3`, and let `R` be a totally ramified DVR with uniformizer `pi` and
`pi^5=3`.  In the fifth jet take

```text
H(u)=u+pi u^2+u^5.                                       (23)
```

Its reduction first separates at depth five, but direct composition gives

```text
[u^5](H^3-u)=3+10pi^4.                                   (24)
```

The second term is not divisible by `3=pi^5`.  Hence the divided expression
in `(17)` is not even integral.  Agreement modulo the maximal ideal is not
enough in a ramified DVR; the stronger p-divisibility of the lower jets, or a
separate ramified Nottingham/Witt correction, is load bearing.

## 8. Frontier connections and refusals

### Finite-place Gaussian separation

THM-2020 exposes a first nonzero finite-place **factorial-channel valuation**.
THM-3227 exposes a first nonzero finite-place **coefficient-jet digit**.  Both
produce ultrametric separation and a carry, but no proved map identifies
their addresses.  THM-2022 already closes NC2/GMC(2); this trie is structural
sidecar information, not a replacement proof.

### LRC root and carrier gates

For the prime residue fields `F_7` and `F_13`, the trie is respectively
7-ary or 13-ary and every charged edge has the matching prime reset.  This
does not put the trie on a lawful LRC physical packet.  THM-2820's Boolean
carrier rigidity and THM-2591's root-selector Cech boundary remain in force:
a supplied polynomial root branch is not a canonical owner/root chart.

### Scale-tree holotopy

THM-3086's rooted scale tree records ordered asymptotic cluster composition;
the present trie records coefficient-prefix contact.  Their tree identities
look alike, but there is no canonical functor from one to the other.  A
future bridge must specify which Newton cluster selects which simple root and
show that its first residue digit survives the whole physical layer.

## 9. Connection contract

```text
source:      finite family of normalized jets at one selected simple root;
map:         reduced coefficient prefix and first unequal residue digit;
target:      rooted q-ary affine trie with p-order charged edges;
preserved:   contact depth, weighted leading difference, residue branch;
destroyed:   generic-tree chronology, absolute scalar, global ownership;
sidecar:     root selector plus physical carrier/whole-layer noncancellation.
```

## 10. Exact companion

The assertion-independent companion

```text
04-computation/selected_root_residue_contact_trie_carry_thm3227.py
```

pins promoted THM-3221 and checks exactly:

```text
44,229 ultrametric triples;
40 sharp full-trie nodes over q=2,3,4,5,9;
the leaves-minus-one affine census;
18 extension-field nonzero additive orders;
60 nonlinear source-coordinate changes;
273 unramified primitive carries;
15 copies of the ramified formula 3b+10a^4;
the delayed depth-two to depth-three/four hostile.          (25)
```

Ordinary and optimized runs byte-match

```text
05-knowledge/results/selected_root_residue_contact_trie_carry_thm3227.out
```

and the LF-normalized hashes are pinned in the frontmatter.

QED.
