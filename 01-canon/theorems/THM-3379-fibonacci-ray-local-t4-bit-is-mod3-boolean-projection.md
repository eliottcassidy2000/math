---
id: THM-3379
title: "Fibonacci-ray local T4 bit as a mod-three Boolean projection"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On THM-3339's
  three distinguished normalized Fibonacci/Berggren rays, THM-3364's labelled
  local-T4 outer-order bit is
  exactly epsilon(k)=1 iff k=1 mod 3.  More strongly, the three ancestry rays
  are the three median-channel fibres of the six Farey order states:
  R0 corresponds to median a, R1 to median b, and R2 to median c; epsilon is
  the indicator of the middle fibre R1.  Cassini parity distinguishes the two
  reversed orders within each median fibre.  The flipped-T4 support is
  {4,7,10,...}, with eventual cyclotomic period three, density 1/3, scar 5/18,
  and an explicit harmonic constant.  This collapse fails on the full ternary
  tree and does not recover words, triples, owners, global arcs, LRC current,
  or JC flux.
source: codex-2026-08-14-fibonacci-ray-t4-projection
audit: independent comparator, word-order, Fourier, harmonic, and scope audit
depends_on:
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
  - THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase
related:
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
script: 04-computation/fibonacci_ray_t4_mod3_projection_thm3379.py
output: 05-knowledge/results/fibonacci_ray_t4_mod3_projection_thm3379.out
script_sha256: 803065a6ce419554b04403b13537889547da53bd1e30aa35edc5e5caca549ae8
output_sha256: d3cf7be59fe988171b73184b9cf26974b26bdfbcd1e50bbb3d40075a7d666778
hash_basis: LF-normalized bytes
---

# THM-3379 -- the Fibonacci-ray `T4` bit is a median-channel projection

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and connection contract

[THM-3339](THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md)
puts the normalized consecutive-Fibonacci triples on three exact Berggren
parameter rays and assigns every Fibonacci index one of six Farey channel
orders.  [THM-3364](THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase.md)
makes every Berggren parent and its three children a labelled transitive `T4`
with one reset/XOR outer-order bit.  Neither theorem alone identifies those two
finite-state objects.

Their exact intersection is:

| field | connection |
|---|---|
| source | normalized Fibonacci indices `k>=2`, their three ray words, and their six Farey channel-order states |
| target | the labelled local-`T4` bit `epsilon` and its support as a subset of `N` |
| map | ancestry residue `k mod 3` -> median channel -> `1_[median=b]` |
| preserved | ancestry-ray class and which of `P<A<C<B` or `P<C<A<B` is the local order |
| destroyed | Cassini parity, the order-versus-reversal choice, the exact word, parameter, triple, owner and current |
| sidecar | `k mod 2` splits each median fibre into its two reversed order states |
| hostile | equal depth, equal branch counts, or equal terminal branch do not determine `epsilon` on the full ternary tree |

Here `A,B,C` are THM-3339's parameter branches and are respectively
THM-3364's `L,M,R`.  No literature-priority claim is made for this elementary
combination of the two proved automata.

## 2. The comparator algebra

For a primitive positive Euclid parameter `u=(m,n)`, `0<m<n`, put

```text
(a,b,c)=(n^2-m^2,2mn,n^2+m^2),
Delta(u)=b-a=m^2+2mn-n^2.                                (1)
```

There is no tie: `Delta(u)=0` would make `n/m=1+sqrt(2)`.  Encode the local
outer order by

```text
epsilon(u)=0 if Delta(u)>0,
epsilon(u)=1 if Delta(u)<0.                              (2)
```

For

```text
A=[ 0 1;-1 2],       B=[0 1;1 2],       C=[1 0;2 1],    (3)
```

direct expansion gives

```text
Delta(Au)= n^2+2mn-m^2 > 0,
Delta(Bu)=-Delta(u),
Delta(Cu)=-(n^2+2mn-m^2) < 0.                            (4)
```

Thus, without a finite search,

```text
epsilon(Au)=0,
epsilon(Bu)=epsilon(u) xor 1,
epsilon(Cu)=1.                                          (5)
```

The two labelled local tournaments are

```text
epsilon=0:  P < A < C < B,
epsilon=1:  P < C < A < B.                              (6)
```

After forgetting the labels they are the same transitive `T4`; the bit is a
labelled-order coordinate, not an unlabelled tournament invariant.

## 3. Restriction to the three Fibonacci rays

THM-3339 proves that, from the root parameter `(1,2)`, the normalized triple at
index `k>=2` has the root-to-child word

```text
k=3r+2:  (BA)^r              = R2,
k=3r+3:  A(BA)^r             = R0,
k=3r+4:  C(BC)^r             = R1.                       (7)
```

For `r=0`, the first word is empty and the root has `epsilon=0`.  For `r>0`,
`(BA)^r` ends in `A`; every `A(BA)^r` also ends in `A`; every `C(BC)^r` ends in
`C`.  The reset laws in (5) therefore prove for every `k>=2` that

```text
epsilon(k)=1  iff  k=1 mod 3  iff the normalized word lies on R1. (8)
```

So the apparent six-state clock times an extra Boolean bit is not a
twelve-state object on this carrier.  The bit is already a function of its
mod-three ancestry coordinate.

## 4. The six orders factor through their median channel

THM-3339's six Farey channel orders, indexed by `t=k mod 6`, are

| `t` | order | median | `epsilon` | Cassini sign |
|---:|:---:|:---:|---:|---:|
| 0 | `cab` | `a` | 0 | + |
| 1 | `cba` | `b` | 1 | - |
| 2 | `bca` | `c` | 0 | + |
| 3 | `bac` | `a` | 0 | - |
| 4 | `abc` | `b` | 1 | + |
| 5 | `acb` | `c` | 0 | - |

States `t` and `t+3` are reversals.  Quotienting the six orders by reversal
retains exactly the median channel, and (7)--(8) sharpen to

```text
R0 <-> median a,       R1 <-> median b,       R2 <-> median c,
epsilon=1_[median=b].                                      (9)
```

Cassini parity is the orientation inside each reversal pair.  Consequently
the `R1` fibre is the disjoint union of states `S_1` and `S_4`.

This is a Boolean projection, not a character `S3 -> C2`: its fibres have
sizes two and four, whereas a nontrivial group character has fibres three and
three.  The word `xor` in (5) describes branch dynamics; it does not turn (9)
into a group homomorphism.  Also, the six objects in this table are the six
order states of three channels, not the six vertices of THM-3339's separate
edge-product `T6`.

## 5. Cyclotomic and harmonic profile

The actual flipped-`T4` index support is

```text
H_F={k>=2:epsilon(k)=1}={4,7,10,...}
   ={n>=1:n=1 mod 3} minus {1}.                          (10)
```

Modulo finite startup, its period-three word is `(0,1,0)`.  With
`omega=exp(2 pi i/3)`, its nonzero Fourier coordinates are

```text
Phi_H(0)=1/3,
Phi_H(1/3)=omega^(-1)/3,
Phi_H(2/3)=omega^(-2)/3.                                (11)
```

In the period-six chart, (10) is `S_1 union S_4`; its odd modes cancel and its
even modes are exactly (11).  The transform is a convolution idempotent as in
THM-3364.  Its harmonic residue and THM-3359 deletion scar are

```text
delta(H_F)=1/3,                 scar(H_F)=5/18.          (12)
```

Finite startup matters to the constant term.  The digamma identity

```text
psi(1/3)=-gamma-pi/(2 sqrt(3))-(3/2)log(3)              (13)
```

and deletion of the initial element `1` give

```text
sum_(k<=x,k in H_F) 1/k
 = (1/3)log x + C_F + O(1/x),
C_F=gamma/3 + pi/(6 sqrt(3)) + log(3)/6 - 1.            (14)
```

Thus harmonic density sees only that one third of indices fire; the two
nonzero modes recover which residue class fires, and the transient `-1`
remembers that the Fibonacci carrier starts at `k=2`.

## 6. Why the collapse is not a full-tree law

The reset/XOR automaton retains word order away from the selected rays.  From
the root `(1,2)`:

```text
word   parameter   epsilon
A      (2,3)          0
B      (2,5)          1
AB     (3,8)          1
BA     (5,8)          0
CB     (4,9)          0.                                (15)
```

The first two words refute depth-only determination.  `AB` and `BA` have the
same branch-count vector but opposite bits, so abelianizing the word fails.
`AB` and `CB` have the same terminal branch `B` but opposite bits, so a suffix
reset heuristic also fails.  The exact whole-level counts remain

```text
#epsilon_0(d)=(3^d+(-1)^d)/2,
#epsilon_1(d)=(3^d-(-1)^d)/2,                            (16)
```

which do not select the three rays.  Nor does (9) create a global tournament:
THM-3364's distinct parameters `(1,8)` and `(4,7)` still tie at hypotenuse 65.

## 7. Verification and stopping boundary

The companion pins both theorem dependencies and uses no floating point or
`assert`.  It checks (4)--(6) on `6,120` primitive branch instances; verifies
all `603` ray indices `2..604`, including normalized Fibonacci parameters,
local `T4` orders, six-state medians and (8); computes (11) exactly in
`Q(omega)` and checks period-six convolution idempotence; and checks `88,573`
full-tree nodes through depth ten together with the hostile words (15).
Ordinary and optimized outputs are byte-identical.

The theorem identifies one Boolean quotient on three distinguished rays.  It
does not classify the full Berggren tree, reconstruct a word or triple, choose
an affine owner, orient nonlocal ties, or supply LRC phase/current or JC flux.
