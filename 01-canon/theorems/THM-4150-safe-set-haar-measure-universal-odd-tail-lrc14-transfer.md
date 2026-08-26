---
id: THM-4150
title: "Safe-set Haar-measure transfer for universal odd-tail LRC(14) completion"
status: >
  PROVED FULL-SAFE-SET HAAR/CROSS-COMB TRANSFER + VERIFIED-EXACT +
  INDEPENDENT IMPLEMENTATION AUDIT; LRC(14) OPEN. Every finite positive body
  H whose complete 1/14-safe set has Haar measure at least 4/63 accepts every
  pair of distinct positive odd tails after doubling. The cross-comb bound
  4/63 is sharp and has unique primitive equality ratio (1,9). The explicit
  33-speed pool has exactly 193,536,720 eleven-body subfamilies, of which
  193,328,720 fail THM-4148's stated min/max width gate, but THM-4154 proves
  that this concrete pool-family safety was already inherited from the
  small-denominator/divisor sieve. Arbitrary bodies, parity-class entry, and
  LRC(14) remain open.
source: codex-lrc-multiwindow-probe-20260825
depends_on:
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
related:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4142-common-safe-arc-clock-pool-universal-odd-tail-lrc14-completion
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  - THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction
script: 04-computation/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150.py
output: 05-knowledge/results/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150.out
independent_audit_script: 04-computation/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150_independent_audit.out
script_sha256: 4e81f65385ea81137bb3a8465ca6b7bb635592873960df57001f7d10e5c67283
output_sha256: e35db11b0e380bf9391f29d8d4d7221548f227fdde058625c6716d1401d3e16d
semantic_sha256: 9b08f96a0b6771ddadc366232f0cd7e3d07b208d3328668890b6a6cc24a5390f
independent_audit_script_sha256: fed4a4ffbc750b66b88a2161319812e843d876c98cacb5db6773091c739aecca
independent_audit_output_sha256: f8aa43bb153b637e04f5b7653e532dc08fcedefca4b07744cc38bc2bcd0a5cba
independent_semantic_sha256: 1686df22845483de1378e8c91e9fb1260bfd826fb56e3652e607e8eb67ccf81e
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic reconstructs the Bernoulli formula and a
  separate strict-wall cross-comb census for all 1,053 primitive odd ratios
  through q=101. It rebuilds all 2,472 pool walls, 46 positive safe
  components, the exact Haar measure, the 193,536,720-family count, the
  208,000-family THM-4148 subcount, and both positive and hostile controls.
  Normal, optimized, and hash-seeded replays byte-match the frozen output.
independent_audit: >
  ACCEPT. A no-import implementation uses direct clipped danger-tooth and
  repeated safe-tooth interval intersections, checks 2,350 primitive ratios
  through q=151, independently reproduces the pool and fragmented-body
  geometry, and recomputes the family split. Normal, optimized, and
  hash-seeded replays byte-match its frozen output.
---

# THM-4150 -- full-safe-set Haar transfer

**PROVED FULL-SAFE-SET HAAR/CROSS-COMB TRANSFER + VERIFIED-EXACT +
INDEPENDENT IMPLEMENTATION AUDIT; LRC(14) REMAINS OPEN.**

The closest proved measure mechanism is THM-2061, Sections 1 and 4: its
folded-seam obstruction already gives the sharp nonstrict cap `4/63`, with
primitive equality ratio `(1,9)`. The new abstract content here is the strict
upgrade

```text
failed dyadic seam  =>  mu(G_H)<4/63,
```

together with a Fourier--Bernoulli dual proof. The exact pool geometry and
census below remain valid, but THM-4154 corrects their significance: every
label in that concrete pool is nonzero modulo `6`, so its family-safety
corollary was already covered by THM-366 and THM-2061.

THM-4148 keeps one connected first window and forgets every other component
of the body-safe set. The useful dual operation is to retain all components
but forget their addresses: a failed completion would have to pack the whole
closed safe set into the pullback of one open two-sheet cross-comb. Haar
measure makes that packing impossible above a sharp threshold.

## 1. The transfer theorem

Let `H` be a nonempty finite set of positive integers and put

```text
G_H={y in R/Z: min_(h in H)||hy||>=1/14}.                (1)
```

Write `mu` for normalized Haar measure on `R/Z`.

> **Theorem.** If
>
> ```text
> mu(G_H)>=4/63,                                         (2)
> ```
>
> then for every two distinct positive odd integers `a,b` there is an
> `x in R/Z` such that
>
> ```text
> min_(v in 2H union {a,b})||vx||>=1/14.                 (3)
> ```

The body may have any cardinality. The thirteen-speed corollary below takes
eleven-element subbodies.

## 2. Exact Haar measure of the cross-comb

Put `delta=1/14` and

```text
A(z)=1_(||z||<delta),              D_r={z:A(rz)=1}.       (4)
```

For coprime odd integers `0<p<q`, THM-4136 reduces simultaneous failure of
the two physical half-lifts to the two cross terms. If

```text
B_(p,q)=D_p intersect (D_q-1/2),                         (5)
```

then the other term is `B_(p,q)+1/2`. These terms are disjoint because
`D_p intersect (D_p-1/2)` is empty. More precisely, the inverse image of the
open bad quotient under doubling is their disjoint union,

```text
m_2^(-1)(C_(p,q))=B_(p,q) disjoint_union (B_(p,q)+1/2).
```

Haar invariance under the surjective doubling map therefore gives

```text
mu(C_(p,q))=2mu(B_(p,q)).                                (6)
```

The Fourier coefficients of `A` are

```text
Ahat(0)=1/7,
Ahat(n)=sin(pi n/7)/(pi n),                  n!=0.        (7)
```

Orthogonality in the integral
`int A(pz)A(q(z+1/2)) dz` forces the frequency pair `(qk,-pk)`.
Because `p` is odd, the half-shift contributes `(-1)^k`. Thus

```text
mu(B_(p,q))
 =1/49 + 2/(pi^2 p q)
    sum_(k>=1) (-1)^k sin(pi qk/7)sin(pi pk/7)/k^2.      (8)
```

Let

```text
B_2(u)=u^2-u+1/6,                         0<=u<1,
u_-={1/2+(q-p)/14},          u_+={1/2+(q+p)/14},         (9)
```

where braces mean fractional part. Product-to-sum in `(8)` and

```text
sum_(k>=1) cos(2pi k u)/k^2=pi^2 B_2({u})               (10)
```

give the exact rational formula

```text
mu(B_(p,q))=1/49+[B_2(u_-)-B_2(u_+)]/(pq).              (11)
```

Since `-1/12<=B_2<=1/6`, the bracket in `(11)` is at most
`1/4`. Hence for `pq>=23`,

```text
mu(C_(p,q))
 <=2/49+1/(2pq)
 <=2/49+1/46
 <4/63.                                                 (12)
```

The complete coprime odd universe with `p<q` and `pq<=22` is:

| `(p,q)` | `mu(C_(p,q))` | `(p,q)` | `mu(C_(p,q))` |
|---:|---:|---:|---:|
| `(1,3)` | `0` | `(1,5)` | `0` |
| `(1,7)` | `2/49` | `(1,9)` | `4/63` |
| `(1,11)` | `4/77` | `(1,13)` | `4/91` |
| `(1,15)` | `4/105` | `(1,17)` | `4/119` |
| `(1,19)` | `4/133` | `(1,21)` | `2/49` |
| `(3,5)` | `2/105` | `(3,7)` | `2/49` |

Substitution in `(11)` proves

```text
mu(C_(p,q))<=4/63,                                      (13)
```

with equality only at `(p,q)=(1,9)`. The constant is therefore sharp for
the cross-comb itself, not merely an artifact of an inequality.

## 3. Compact closed set versus open pullback

Fix distinct positive odd tails, interchange them if necessary, and write

```text
a=pt,          b=qt,          t=gcd(a,b),
0<p<q,         gcd(p,q)=1.                               (14)
```

All of `p,q,t` are odd. For every `y in G_H`, both physical lifts

```text
x_0=y/2,                         x_1=(y+1)/2             (15)
```

preserve every doubled-body gap because `||2h x_i||=||hy||`. If `(3)`
were false, both lifts would be tail-bad for every `y in G_H`. In the
THM-4136 quotient this says

```text
G_H subset m_t^(-1)(C_(p,q)),             m_t(y)=ty.     (16)
```

The map `m_t` is a surjective compact-group endomorphism. The pushforward
of normalized Haar measure is again normalized Haar measure, and therefore

```text
mu(m_t^(-1)(C_(p,q)))=mu(C_(p,q)).                       (17)
```

Now `G_H` is compact and `m_t^(-1)(C_(p,q))` is open. Moreover it is proper:
`0 notin C_(p,q)`, because at quotient zero the zero half-lift is bad but the
half-turn lift is safe for both odd primitive tails. From `(2)`, `(13)`,
`(16)`, and `(17)`, inclusion would force equality of all measures.

That equality is impossible. If an open proper set `U` contains a nonempty
compact set `K` and `mu(U)=mu(K)`, then `U\K` must be empty because every
nonempty open subset of the circle has positive Haar measure. But then
`U=K` is a nonempty proper clopen subset of the connected circle. Taking
`K=G_H` and `U=m_t^(-1)(C_(p,q))` contradicts `(16)`. This proves `(3)`.
**QED.**

## 4. Exact geometry of a 33-label pool; safety inherited

Let

```text
P={1,2,4,5,8,10,16,17,19,20,23,25,29,31,32,34,38,
   40,41,43,47,50,51,53,58,62,64,69,71,73,75,76,80}.    (18)
```

Literal reduction modulo `6` gives `P intersect 6Z=empty`. Therefore
THM-4154's fixed phase `x=1/12` already gives body clearance at least `1/6`
and every odd tail clearance at least `1/12`. The Haar argument below is a
second, logically independent certificate for this pool and supplies the
fragmented-safe-set geometry used to test the abstract transfer theorem; it
is not the first proof of the pool-family safety.

The complete exact wall arrangement of `G_P` has `2,472` distinct walls
and `46` positive-length components. Direct rational summation gives

```text
mu(G_P)
 =110551382435042260737/1702610555154297252800
 =4/63
  +3148874954907620719/2189070713769810753600
 >4/63.                                                 (19)
```

Its largest component has length only

```text
3/700 < 2/189.                                          (20)
```

For every eleven-element `H subset P`, inclusion `G_P subset G_H` and
`(19)` invoke the theorem. Hence every row

```text
2H union {a,b},                    a,b distinct odd,     (21)
```

is `1/14`-safe. The pool supplies exactly

```text
binom(33,11)=193,536,720.                               (22)
```

This is an exact hereditary census, but not a new coverage count relative to
the proved divisor sieve: all `193,536,720` bodies were already safe by the
common `x=1/12` certificate of THM-4154.

The comparison with THM-4148 is made **inside precisely the universe of all
eleven-subsets of `P`**. For fixed minimum and maximum positions `i<j`,
there are `binom(j-i-1,9)` such bodies. Summing this number over exactly the
pairs satisfying THM-4148's stated conditions

```text
27(13 min H-max H)>=4 min H max H.                       (23)
```

gives `208,000`. Therefore the complementary

```text
193,536,720-208,000=193,328,720                         (24)
```

families fail that exact min/max width gate. This is not a claim that they
evade every other theorem or every body-specific certificate. Also
`|P intersect P_THM-4142|=9`, so no eleven-subset in `(22)` is one of
THM-4142's common-pool bodies.

## 5. Hostile controls and equality boundary

The eleven-body

```text
H_*={1,17,31,32,41,47,50,51,58,62,71}                  (25)
```

has `100` positive-length safe components and

```text
mu(G_(H_*))=1142176622583/5854727790800,
max_component_width=21/4100 < 2/189.                    (26)
```

Its first min/max window is empty because
`13/(14*71)<1/14`. Thus `(25)` is a literal fragmented-body hostile to
reading this theorem as a disguised first-window argument.

The threshold is sufficient, not necessary. For `H_0={1,...,11}`,

```text
mu(G_(H_0))=10931/194040 < 4/63,                         (27)
```

yet the particular row `2H_0 union {1,3}` has clearance `1/13` at
`x=1/13`. Falling below `(2)` is therefore uncertified, not unsafe.

At the other boundary, `(p,q)=(1,9)` attains cross-comb measure `4/63`.
The non-strict hypothesis in `(2)` remains valid only because a compact
closed safe set cannot fill an equal-measure proper open pullback. Dropping
the closed/open distinction would lose the equality case.

## 6. Information ledger and scope

```text
source:       the complete closed body-safe set G_H
target:       every doubled-body row 2H union {a,b}, a,b distinct odd
map:          y -> {y/2,(y+1)/2}, then y -> gcd(a,b)y
preserved:    every body clearance, both physical sheets, strict tail danger
destroyed:    safe-component addresses, owners, and internal body labels
sidecar:      exact open cross-comb Haar measure and compact/open equality
hostile:      H_* has 100 components and no component of width 2/189
decisive test: mu(G_H) versus sharp sup_(p,q)mu(C_(p,q))=4/63. (28)
```

This theorem proves a full-safe-set transfer criterion and includes
minimum-one bodies in its exact test pool. THM-4154 supersedes only the claim
that the explicit pool contributes previously uncovered families. This
theorem does not prove a uniform lower bound `mu(G_H)>=4/63`, classify
bodies below the threshold, handle mixed/even tail parities, or provide
physical entry into the `11+2` odd-tail branch. It does not prove LRC(14).

## 7. Exact replay

The primary implementation compares `(11)` against a separately rebuilt
strict-wall cross-comb on every primitive odd ratio through `q=101`; the
independent implementation instead intersects clipped periodic teeth through
`q=151`. They use disjoint body-safe-set representations and agree on every
stated value. Replay with

```text
python3 -B 04-computation/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150.py
python3 -B -O 04-computation/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150.py
PYTHONHASHSEED=271828 python3 -B 04-computation/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150.py

python3 -B 04-computation/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150_independent_audit.py
python3 -B -O 04-computation/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150_independent_audit.py
PYTHONHASHSEED=8675309 python3 -B 04-computation/lrc14_safe_set_haar_measure_odd_tail_transfer_thm4150_independent_audit.py
```

All streams byte-match their frozen outputs. **QED.**
