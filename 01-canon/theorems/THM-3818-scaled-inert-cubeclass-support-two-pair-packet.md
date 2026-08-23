---
id: THM-3818
title: "Scaled inert cube classes recover the support-two pair packet"
status: >
  PROVED ALGEBRA + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED.  The 5,855 THM-3793 primitive cube sums occupy distinct classes
  modulo rational cubes.  Hence their physical cube sums are injective across
  every positive common scale and recover the selected support-two flatness
  covector, exposed facet pair, and THM-778 two-clock packet after labelled
  placement.  If the entire physical pair sum is all-inert, the address is
  also unique among unrestricted positive two-cube representations.  Eleven
  ambient residues give an exact pair-sum-lattice covering reduction.  A
  taxicab hostile blocks unrestricted decoding, and two rows with the same
  unique minimum relation block every pair-local semantic-arrival claim.  No
  LRC(14) row is excluded.
source: root + lrc_reversible_address / incoming-signal extension session, 2026-08-23
audit: >
  CANONICAL SELF-AUDIT + INDEPENDENT HOSTILE AUDIT PASS.  An
  assertion-free companion enumerates the exact 94-sum, 5,855-ratio universe,
  computes every rational cube class independently by trial-factor exponent
  reduction and maximal-cube-divisor scanning, checks 1,376,002 synchronized
  phases, 11,710 facet identities, 46,840 scale controls, and both hostiles in
  1,775,955 active gates.  Normal and optimized runs byte-match the frozen
  output.  An independent multiplicative enumeration tested all 17,137,585
  unordered pairs by gcd reduction, found zero rational-cube collisions, and
  decoded 585,500 scaled controls.  The audit repaired the distinction between
  a pair merely occurring in a row and the selected minimal support-two
  relation.  The finite separation must not be upgraded to an unrestricted
  singleton theorem without the physical all-inert hypothesis.  A second
  independently written companion checks 456,690 labelled placements, an
  unrestricted sum-divisor decoder on 32 all-inert scale controls, the exact
  full residue schedule, and split/exponent/off-grid hostiles in 2,764,096
  active gates.
depends_on:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
  - THM-778-centered-christoffel-endpoint-skew-product
related:
  - THM-3718-lrc-complete-atom-orbit-defect-saturation-and-semantic-boundary
  - THM-3731-lrc14-valuation-owner-equivariance-and-semantic-packet-boundary
script: 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.py
output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.out
script_sha256: e37ae5e8adf3414081e76d56eb01607172a6364703dca76b41ac096c9f6d77c3
output_sha256: 6199147668a19c0a2ed403b14d889caf65078c675068f000856202cb68a861a8
independent_script: 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_independent_audit_thm3818.py
independent_output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_pair_packet_independent_audit_thm3818.out
independent_script_sha256: 0bb6b60238335695addce0010182f5222da4c4a86c6f59912bd0f36123d063ca
independent_output_sha256: 37aa867d653d7123435f9841cdc772c3cd71c289162006f80aeb6079a17ad7bc
independent_semantic_sha256: 86c481b03158ba3cb7024ef8739fde640be331d9856ce5d4fc0c1b3b4fcc06cb
hash_basis: raw LF bytes
---

# THM-3818 -- scaled inert cube classes recover the support-two pair packet

**PROVED ALGEBRA + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  This theorem strengthens THM-3793's finite LRC address from
raw injectivity to separation modulo rational cubes and then records the exact
scale, facet, word, and residue consequences.  It is a restricted address
theorem, not a sufficient badness certificate or an LRC(14) proof.

## 1. The finite primitive atlas

Let

```text
R={ (p,q) in Z_(>0)^2 :
      p<q, gcd(p,q)=1, p+q<=356,
      every prime ell|(p+q) is 2 mod 3,
      and v_ell(p+q)<=2 }.                             (1)
```

THM-3743 supplies the ambient support-two cap, while THM-3793 gives

```text
94 admissible pair sums,             |R|=5,855.       (2)
```

For a positive integer `m`, write

```text
kappa(m)=[m] in Q_(>0)^x/(Q_(>0)^x)^3.                (3)
```

### Finite-exact cube-class separation

The `5,855` classes

```text
kappa(p^3+q^3),                    (p,q) in R,         (4)
```

are pairwise distinct.  Equivalently, for two members of `R`,

```text
(p^3+q^3)/(r^3+s^3) is a positive rational cube
              => (p,q)=(r,s).                         (5)
```

This is **FINITE-EXACT** on the explicitly frozen universe (1).  It is not
claimed beyond `p+q<=356`.

## 2. Two exact cube-class algorithms

For the first computation, trial-factor

```text
m=p^3+q^3
```

and retain the finite exponent vector

```text
(v_ell(m) mod 3)_ell.                                 (6)
```

For the independent integer computation, put

```text
c(m)=max{c in Z_(>0):c^3|m},       K(m)=m/c(m)^3.     (7)
```

Since

```text
p^3+q^3<(p+q)^3<=356^3,
```

the second algorithm scans every possible `c` through `356`, without
factoring `m`.  Prime valuations prove

```text
kappa(m)=kappa(n) iff K(m)=K(n).                       (8)
```

The algorithms agree row by row, and the `5,855` integers `K(m)` are
distinct.  The canonical rows

```text
[p,q,m,c(m),K(m)]
```

in increasing `p+q`, then increasing `p`, have SHA-256

```text
feee0062ab6afb2d19f00e5397f0650f3738b860d1cdc3c1ab49193b94636a26. (9)
```

There are exactly `43` primitive values with `c(m)>1`.  Their cube-root
distribution is

```text
c=7:36,             c=13:4,       c=19:2,       c=31:1. (10)
```

The first is

```text
1^3+19^3=6860=7^3*20.                                  (11)
```

Thus the largest cube divisor of an address is not automatically its physical
gcd scale.  The primitive row must be decoded before the scale is read.

## 3. All-scale injectivity and the decoder

For `g>=1` and `(p,q) in R`, define the physical address

```text
A(g,p,q)=(gp)^3+(gq)^3=g^3(p^3+q^3).                  (12)
```

Then

```text
A:Z_(>0) x R -> Z_(>0)
```

is injective.  Indeed, equality of two values in (12), followed by passage to
(3), identifies the primitive pair by (5).  Cancelling its positive cube sum
then gives `g^3=h^3`, hence `g=h`.

The exact decoder on the image is:

1. compute `kappa(N)`;
2. look up the unique `(p,q) in R` with class (4);
3. verify that `N/(p^3+q^3)` is an integer cube;
4. take its positive cube root `g`.

The native dilation operation is lawful:

```text
(g,p,q)->(kg,p,q)                corresponds to A->k^3 A. (13)
```

This proves injectivity **inside the scaled atlas**.  It does not say that the
integer `A(g,p,q)` has no representation outside that atlas.

There is one stronger specialization.  If every prime divisor of the physical
pair sum

```text
D=g(p+q)                                                (13a)
```

is `2 mod 3`, then THM-3793 says that `A(g,p,q)` has exactly one positive
distinct two-cube representation, with no atlas restriction on a competitor.
The exponents in `g` may be arbitrary; the exponent-at-most-two condition is
only on the primitive sum `p+q` from (1).  This unrestricted representation
can be decoded without a cube-class table: for each candidate sum `d=x+y`
dividing `N=x^3+y^3`, compute

```text
xy=(d^3-N)/(3d)                                       (13b)
```

and retain the positive integral roots of `T^2-dT+xy`.  The independent
companion verifies this decoder on 32 all-inert scale controls.  The identities

```text
1729=1^3+12^3=9^3+10^3,
515375=15^3+80^3=54^3+71^3
```

are sharp split-prime and primitive-exponent-three hostiles.  For arbitrary
scale, only the within-atlas injectivity above is claimed; (21) is the
corresponding out-of-atlas warning.

## 4. Selected labelled covector and exposed facets

Suppose the chosen `l1`-minimal THM-3743 relation has support two and its
labelled speeds satisfy

```text
n_i=gp<n_j=gq.                                           (14)
```

The canonical primitive support-two covector is

```text
a=q e_i-p e_j,        a dot n=0,       ||a||_1=p+q.     (15)
```

Thus the address plus labelled placement recovers the selected relation up to
the unavoidable global sign.  If runner coordinates are not ordered by
increasing speed, one additional orientation bit records which owner carries
`gp`.

For THM-3743's projected cube, pairing with `a` commutes with projection.  The
two cube faces exposed by (15) are

```text
F_-: x_i=1/14,  x_j=13/14,
F_+: x_i=13/14, x_j=1/14,                              (16)
```

and the exact functional interval is

```text
[(q-13p)/14, (13q-p)/14].                              (17)
```

Its width is

```text
(13q-p-q+13p)/14=(6/7)(p+q),                           (18)
```

as required by THM-3743.  Therefore the selected covector, sign partition,
exposed facet pair, and width are recovered.  Minimality comes from the
selection hypothesis, not from the address: for example, the row `{1,...,13}`
contains the admissible pair `(1,4)` of covector norm five, while `(1,2)` gives
a norm-three relation.  If an atlas pair merely occurs in an ambient row, the
address recovers that primitive pair covector and its facets, but cannot call
it the first flatness pull.  The address does **not** recover which lattice
points occupy the affine sections of (17), nor the quotient gauge and lifted
coefficients needed for a second flatness pull.

## 5. The THM-778 two-clock packet

THM-778 now gives, without another sidecar,

```text
W(gp,gq)=W(p,q)^g,                                     (19)
```

the continued fraction of `q/p`, its midpoint phase cocycle beginning at
`s=1`, and the tie count

```text
g, if p and q are both odd;
0, otherwise.                                          (20)
```

Thus scale, the complete two-owner endpoint word, and its internal phase/tie
packet are recoverable.  This pair phase is not a global LRC observation phase,
and a runner owner label is not THM-3718/3731's semantic root owner without a
typed common-base map.

## 6. The unrestricted taxicab hostile

The atlas restriction is load-bearing:

```text
86^3+344^3=197^3+323^3=41,343,640.                     (21)
```

The first pair is `86(1,4)` and its primitive ratio lies in `R`.  The competing
pair is primitive and outside `R`:

```text
197+323=520=2^3*5*13.                                  (22)
```

Moreover, the actual first sum is

```text
86+344=430=2*5*43,                                     (23)
```

so THM-3793's all-inert physical-sum singleton hypothesis does not apply.

The two THM-778 schedules are different:

```text
(86,344):  gcd 86, reduced 1:4, CF (4),
           zero ties, 430 wall blocks;
(197,323): gcd 1, both odd, CF (1,1,1,1,3,2,3,2),
           one tie at 1/2, 519 wall blocks.             (24)
```

Therefore the raw integer cube sum is not an unrestricted scale, word, or
phase carrier.  The decoder returns the unique **atlas member**, not every
positive two-cube representation.

## 7. A same-packet global-arrival hostile

Put

```text
S_+={1,4,21,21^2,...,21^11},
S_-={1,4,21,21^2,...,21^10,5*21^11}.                   (25)
```

Both are primitive sets of thirteen distinct positive speeds.  In each one,
the unique nonzero relation of `l1` norm at most five, up to sign, is

```text
4(1)-1(4)=0.                                           (26)
```

To prove uniqueness, suppose a relation of norm at most five uses a largest
term `21^k`.  Its remaining coefficient mass is at most four, so for `k>=2`
its absolute contribution is at most `4*21^(k-1)<21^k`; for `k=1` it is at
most `4*4<21`.  If the largest term is `5*21^11`, the same bound is stronger.
Thus only `1` and `4` can occur.  The integer solutions of `c_0+4c_1=0`
have norm at least five, with equality only at (26).

Nevertheless, at the common phase `t=1/5`,

```text
min_(v in S_+) ||v/5||=1/5,
min_(v in S_-) ||v/5||=0.                              (27)
```

Indeed, `21^k=1 mod 5`, while the last speed of `S_-` is zero modulo five.
The two rows have the same unique minimum cube address `65`, scale, labelled
pair, covector, facet pair, THM-778 word, internal phase, and tie packet, but
opposite safety at the same global phase.  Hence no sidecar depending only on
the selected pair can supply semantic arrival.

These are interface hostiles, not LRC counterexamples.

## 8. The lawful ambient-residue sidecar

For an ambient labelled row containing (14), put

```text
D=n_i+n_j=g(p+q),
d_D(x)=min(x mod D, D-(x mod D)),
A_D=(n_l mod D)_(l notin {i,j}).                       (28)
```

For `m in Z/DZ`, define the selected-pair safe phases

```text
G_(i,j)={m:14 d_D(m n_i)>=D}.                          (29)
```

Because `n_j=-n_i mod D`, the same set is obtained from owner `j`.  For every
other runner define

```text
B_l={m:14 d_D(m n_l)<D}.                               (30)
```

Since

```text
||m n_l/D||=d_D(m n_l)/D,                              (31)
```

the complete safe schedule on this lattice is

```text
K_D(n)=G_(i,j) minus union_(l notin {i,j}) B_l.        (31a)
```

there is a lonely witness on the pair-sum lattice `D^(-1)Z` if and only if

```text
G_(i,j) is not contained in union_(l notin {i,j}) B_l. (32)
```

Consequently every hypothetical counterexample in this branch must satisfy
the exact finite covering obligation

```text
G_(i,j) subset union_(l notin {i,j}) B_l.              (33)
```

The eleven labelled residues in (28) determine every set in (30).  Under a
common global dilation by `c`, `(D,A_D)` becomes `(cD,cA_D)` and the schedule
repeats `c` times, so the cover predicate is equivariant.  In (25), the final
residue modulo five is `1` for `S_+` and `0` for `S_-`; thus (28) is exactly
the first coordinate that distinguishes the hostile.

Equation (33) is a genuine reduction, but only on one discrete pair-sum
lattice.  A row satisfying it can still have an off-lattice witness.  It also
does not form THM-3731's owner/root/deep tensor or repair THM-3718's excluded
zero-charge target.

The grid restriction is sharp even when (31a) is known exactly.  For

```text
n=(1,2,...,13),             selected pair=(1,4),       D=5,
```

runner five makes `K_5(n)` empty, whereas the off-grid phase `t=1/14` is a
lonely witness.  Conversely, the independent positive control
`(1,2,...,10,12,13,14)` with selected pair `(1,10)` has all ten nonzero
classes modulo eleven in `K_11`.  Thus the residue packet is a lossless
description of its declared lattice and no more.

## 9. Typed connection and boundary

```text
source:      scaled support-two THM-3743 branch on the 5,855 inert ratios
target:      cube class, integer address, labelled covector/facets,
             and pair-sum residue schedule
map:         {gp,gq} -> [p^3+q^3], g^3(p^3+q^3),
             (n_l mod g(p+q))_l
preserved:   primitive ratio, scale, selected flatness width/facets,
             complete two-clock word, midpoint cocycle, ties,
             and safety on the selected pair-sum lattice
destroyed:   off-lattice phase, later flatness slices, semantic root,
             ancestry, complete atom, and positive arrival
next test:   apply (33) to actual rank-eleven/star or finite-terminal rows;
             only then form a same-base owner/root/deep tensor.             (34)
```

The finite atlas has not been intersected with every rank-eleven star space,
and (33) has not removed a row.  LRC(14) remains **OPEN**.

## 10. Reproduction and audit status

Run

```bash
python3 -B 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.py
python3 -B -O 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.py
python3 -B 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_independent_audit_thm3818.py
python3 -B -O 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_independent_audit_thm3818.py
```

Both streams reproduce

```text
05-knowledge/results/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.out
```

byte for byte, with the independent commands reproducing the correspondingly
named independent transcript.  The companions check `1,775,955` and
`2,764,096` active gates without Python assertions.  The proof of Sections
3--9 is algebraic conditional on the exact finite separation in Sections 1--2
and on the cited THM-778 decoder.

An independent implementation generated the 94 admissible sums
multiplicatively, enumerated the 5,855 ratios, and tested all `17,137,585`
unordered row pairs using the criterion that, after division by the gcd, both
remainders must be integer cubes.  It found zero cube-class collisions and
also recovered every row and scale in `585,500` controls with `1<=g<=100`.
The same audit checked the interval, THM-778 packet, both hostiles, the exact
residue-cover iff, replay, hashes, and absence of executable assertions.  It
found and repaired the minimal-selection quantifier above.  The independently
written preserved companion separately checks the sum-divisor decoder,
labelled placements, full grid identity, and the split/exponent/off-grid
boundaries.  No LRC(14) conclusion follows.
**QED.**
