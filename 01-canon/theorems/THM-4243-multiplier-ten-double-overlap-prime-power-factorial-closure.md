---
id: THM-4243
title: "Multiplier-ten double-overlap and prime-power factorial closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. The exact THM-3474
  candidate barcode for multiplier ten is classified by two binary carry
  overlaps. Empty barcode is equivalent to simultaneous overlap of {2H,8H}
  and {4H,6H}, proving the stated factorial gcd and exact quadratic
  three-moment conclusion. The twenty modulus-64 classes are exactly the
  uniformly suffix-forcing classes, not all closing heights.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-4237-multiplier-six-binary-adjacency-prime-power-factorial-closure
  - THM-3475-factorial-all-divisor-digit-polygon-and-pair-ledger-compiler
mistake_firewall:
  - MISTAKE-350
  - MISTAKE-363
scripts:
  - 04-computation/factorial_multiplier_ten_double_overlap_thm4243.py
  - 04-computation/factorial_multiplier_ten_double_overlap_independent_audit_thm4243.py
outputs:
  - 05-knowledge/results/factorial_multiplier_ten_double_overlap_thm4243.out
  - 05-knowledge/results/factorial_multiplier_ten_double_overlap_independent_audit_thm4243.out
script_sha256:
  - 18ada8fca613264c7b04b5e29d22289a20112672ec91c39d224924766ae263cd
  - eaf7fee670800f0ce5df344a38921153fbc638cd57df89147ec42e874a813055
output_sha256:
  - 5ade804f3895eb7c2d354ceb8d9a92181612a11eb998ce27e38ad2b9dd31d9ef
  - 58d315e0d524fe729c06bc785cafbe65445c18578ac40370c2ee60fe1bc2f526
hash_basis: raw LF bytes
audit: >
  Hostile audit PASS. Carry-free splitting gives the two complementary even
  blocks, and the retained lower representative in every reset range proves
  the iff emptiness criterion even for p=11,13. Independent bit-support and
  moment-polynomial engines reproduce the classification, all twenty uniform
  suffix classes, their twelve hostile representatives, and ten good-prime
  gcd rows. Normal and optimized outputs byte-match; no surviving barcode is
  promoted to a factor.
---

# THM-4243 -- multiplier-ten double-overlap and prime-power factorial closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

Let `p>=11` be prime, `k>=1`, and put

```text
H=p^k,                       N=10H,
T=min(9,(p-1)/2).                                             (1)
```

For nonnegative integers, write `x subset_2 y` when every binary one-bit of
`x` is also a one-bit of `y`.  The positive multiplier form of THM-3474's
cross-place candidate barcode is

```text
C(p,k)={t:1<=t<=T and tH subset_2 10H}.                      (2)
```

Define the two carry-overlap predicates

```text
P_28(H): (2H)&(8H) != 0,
P_46(H): (4H)&(6H) != 0.                                    (3)
```

They also admit the shifted forms

```text
(2H)&(8H) = (H&(4H))<<1,
(4H)&(6H) = ((2H)&(3H))<<1.                                 (4)
```

Then the complete exact classification of (2) is

```text
C(p,k)
 = ({2,8} intersect [1,T]  if not P_28(H), else empty)
   union
   ({4,6} intersect [1,T]  if not P_46(H), else empty).      (5)
```

In particular,

```text
C(p,k)=empty  iff  P_28(H) and P_46(H).                     (6)
```

The exact reset-range cases are

```text
p=11:     T=5,  C=({2} if not P_28) union ({4} if not P_46),
p=13:     T=6,  C=({2} if not P_28) union ({4,6} if not P_46),
p=17:     T=8,  C=({2,8} if not P_28) union ({4,6} if not P_46),
p>=19:    T=9,  the same four even candidates as p=17.       (7)
```

Thus the complementary carry pairs are always `{2H,8H}` and `{4H,6H}`.
For `p=11,13`, THM-3474's reset range omits one or both upper complementary
multipliers, but retains `t=2,4`, one representative of each pair.  This is
why the same two overlap tests still decide emptiness in every admissible
prime case.

Whenever (6) holds,

```text
gcd_Q(A_(10p^k-1)^(10p^k+1), A_(10p^k)^(10p^k+1))=1,       (8)
```

where

```text
L(x^r)=r!,            A_n^d(v)=L((d-x+v x^2)^n).           (9)
```

Consequently, for every exact quadratic

```text
q(x)=alpha+beta*x+gamma*x^2,       alpha*beta*gamma != 0,   (10)
```

at least one of

```text
L(q^(10p^k-1)),       L(q^(10p^k)),       L(q^(10p^k+1))   (11)
```

is nonzero whenever (6) holds.

## 2. Binary submasks are carry-free splittings

For `0<=x<=y`,

```text
x subset_2 y  iff  x & (y-x) = 0.                          (12)
```

If `x subset_2 y`, subtracting `x` simply removes its bits from `y`, so the
supports of `x` and `y-x` are disjoint.  Conversely, if those supports are
disjoint, their sum is carry-free and equals `y`; hence every bit of `x`
occurs in `y`.

Apply (12) to `x=tH`, `y=10H`.  Then

```text
tH subset_2 10H  iff  (tH)&((10-t)H)=0.                    (13)
```

This equivalence is exact and retains the candidate degree `tH`; it is not a
heuristic about digit density.

## 3. Proof of the exact classification

The height `H=p^k` is odd.  If `t` is odd, then both `tH` and `(10-t)H` are
odd, so their binary supports overlap at bit zero.  Equation (13) fails.
Thus every surviving multiplier is even.

The even integers between one and nine are `2,4,6,8`.  Complementation
`t -> 10-t` partitions them into

```text
{2,8},                         {4,6}.                       (14)
```

For `t=2` or `8`, equation (13) is precisely the disjointness of `2H` and
`8H`; for `t=4` or `6`, it is precisely the disjointness of `4H` and `6H`.
Intersecting these two blocks with the exact reset interval `[1,T]` proves
(5).

Since `p>=11`, one has `T>=5`.  Therefore `2` and `4` always lie in the reset
range.  If either pair is carry-free, at least its lower representative
survives.  Conversely, if both pairs overlap, every possible even multiplier
fails.  This proves the iff statement (6), including the truncated cases in
(7).

The shifted identities (4) follow because bitwise intersection commutes with
a common left shift.  Thus `P_28` detects two one-bits of `H` separated by two
positions.  The second test remains carry-sensitive because `3H=H+2H` may
itself involve carries.  Both coordinates are load-bearing.

## 4. Transfer through THM-3474

THM-3474's hypotheses are exactly satisfied:

```text
p is prime, p>=11,
a=10 is even,
2<=a<p,
k>=1,
N=a p^k=10H is even.                                      (15)
```

For

```text
d=N+1,
F=A_(N-1)^d,
G=A_N^d,                                                   (16)
```

that theorem gives the complete prime-power reset barcode

```text
{0,H,2H,...,TH}                                            (17)
```

and the complete binary necessary degree barcode for `G`, namely the binary
submasks of `10H`.  Every positive rational common-factor degree must
therefore be `tH` for some `t` in (2).  The coordinate-root capacity is zero.

If (6) holds, there is neither a positive degree address nor a coordinate
root contribution.  THM-3474's exact degree compiler therefore gives (8).

This logical direction is only

```text
empty exact necessary-degree intersection  =>  rational coprimality.       (18)
```

A nonempty set (2) is merely an address left by this compiler.  It does not
construct a common factor and is not evidence that (8) fails.

## 5. Scaling to the exact quadratic three-moment window

Suppose all three moments in (11) vanish and write `N=10p^k`.  THM-3124's
symbolic recurrence, applied at the middle-to-last step, forces the unique
resonance

```text
beta/alpha=-1/(N+1).                                      (19)
```

Put

```text
v=(N+1)gamma/alpha.                                       (20)
```

Then

```text
q(x)=(alpha/(N+1))*((N+1)-x+v*x^2),                       (21)
```

so the first two vanishing moments say

```text
A_(N-1)^(N+1)(v)=A_N^(N+1)(v)=0.                          (22)
```

The two polynomials have rational coefficients.  Any common complex root is
algebraic, and its minimal polynomial over `Q` would be a positive-degree
common rational factor, contradicting (8).  Thus the three simultaneous
zeros are impossible.  This proves (11).

The scaling uses `v=(N+1)gamma/alpha`, not `gamma/alpha`; the factor `N+1`
is load-bearing.

## 6. Congruence families

### 6.1 The root family

If `H=5 mod 8`, reducing the relevant products modulo 16 gives

```text
2H=10, 8H=8, intersection=8;
4H=4,  6H=14, intersection=4.                              (23)
```

If `H=7 mod 8`, the corresponding row is

```text
2H=14, 8H=8, intersection=8;
4H=12, 6H=10, intersection=8.                              (24)
```

Both overlap predicates hold.  Hence `H=5 or 7 mod 8` closes the compiler.
For odd `k`,

```text
p=5 mod 8 => p^k=5 mod 8,
p=7 mod 8 => p^k=7 mod 8.                                 (25)
```

Therefore every prime `p>=11` with `p=5 or 7 mod 8` and every odd `k`
satisfies (8)--(11).

### 6.2 Stronger suffix families

Two further low-bit certificates are

```text
H=27 mod 32,                                               (26)
H=41 or 57 mod 64.                                        (27)
```

For (26), reduce the products modulo 64.  At the representative `27`,

```text
54 & 24 = 16,                 44 & 34 = 32.                (28)
```

For (27), reduce modulo 128:

```text
H=41:   82 & 72 = 64,        36 & 118 = 36;
H=57:  114 & 72 = 64,       100 &  86 = 68.                (29)
```

Thus the following union is a sufficient congruence family:

```text
H=5 or 7 mod 8,
or H=27 mod 32,
or H=41 or 57 mod 64.                                     (30)
```

The extra `27 mod 32` family can be written directly in prime/exponent
coordinates.  For odd `k`, the exact unit calculation modulo 32 is

```text
(k mod 8,p mod 32)=(1,27),(3,3),(5,11),(7,19)
  iff p^k=27 mod 32.                                      (31)
```

A simple even-exponent corollary is

```text
p=11,13,19,21 mod 32 and k=2
  => p^2=57,41,41,57 mod 64, respectively.                 (32)
```

These families are sufficient shadows of the exact bit test, not replacements
for it.

## 7. Exact maximality of the modulus-64 suffix list

This section has a deliberately narrower quantifier than (6).

Call an odd residue `r mod 64` **uniformly suffix-forcing** if every positive
odd `H=r mod 64` satisfies both exact overlap predicates in (3).  Then the
complete list is

```text
R_64={5,7,13,15,21,23,27,29,31,37,
      39,41,45,47,53,55,57,59,61,63}.                      (33)
```

Equivalently,

```text
R_64={r:r=5 or 7 mod 8}
     union {r:r=27 mod 32}
     union {41,57 mod 64}.                                 (34)
```

### Proof of sufficiency

If `H=r+64q` and `m` is one of `2,4,6,8`, then

```text
mH=mr (mod 128),                                           (35)
```

because every such `m` is even.  Hence the lower seven bits of all four
products, and any overlaps visible there, depend only on `r mod 64`.
Equations (23)--(24), lifted through all congruent residues, cover the sixteen
classes `5,7 mod 8`; equation (28) covers `27,59`; and equation (29) covers
`41,57`.  Therefore every class in (33) uniformly forces both overlaps.

### Proof of maximality

There are twelve remaining odd residues.  Their least positive
representatives have the following exact multiplier-nine candidate sets:

```text
r= 1: {2,8}          r= 3: {2,4,6,8}
r= 9: {2,8}          r=11: {4,6}
r=17: {2,8}          r=19: {2,8}
r=25: {2,8}          r=33: {2,8}
r=35: {2,8}          r=43: {4,6}
r=49: {2,8}          r=51: {2,4,6,8}.                     (36)
```

Every set in (36) is nonempty, so the representative `H=r` fails at least one
overlap predicate.  Consequently none of these twelve residue classes can
have the universal quantifier in the definition of suffix-forcing.  Together
with sufficiency, this proves that (33) is maximal.

The two companions reconstruct (33)--(36) independently: the primary uses
integer bit operations modulo 128; the second uses explicit sets of occupied
binary positions.

### Why this is not the unbounded classification

The exact iff criterion (6) inspects every bit of `H`.  Membership in `R_64`
inspects only the suffix and is merely a uniform sufficient condition.

The admissible prime-power height

```text
H=11^3=1331=51 mod 64                                     (37)
```

lies outside `R_64`, yet

```text
(2H)&(8H)=2048,
(4H)&(6H)=5120,                                            (38)
```

so its exact barcode is empty.  By contrast, its least residue representative
`H=51` has all four even multipliers surviving.  Higher bits therefore turn a
nonforcing suffix class into a closing individual height.  This witness
strictly separates the exact unbounded criterion from the maximal uniform
modulus-64 certificate.

## 8. Hostiles and equality boundary

The two overlap coordinates are independent:

```text
H=11: P_28 true,  P_46 false; multiplier block {4,6} survives.
H=17: P_28 false, P_46 true;  multiplier block {2,8} survives.
H=51: both false; all four even multipliers survive.
H=5:  both true; the barcode is empty.                     (39)
```

Reset truncation is visible at actual prime rows:

```text
(p,k)=(11,1): H=11, T=5, C={4};
(p,k)=(13,1): H=13, T=6, C=empty;
(p,k)=(17,1): H=17, T=8, C={2,8};
(p,k)=(43,1): H=43, T=9, C={4,6}.                          (40)
```

The independent companion constructs the actual moment-polynomial pairs for
the rows in (40), together with `(23,1)`, modulo each of `1000003` and
`1000033`.  Both primes exceed twice the largest tested order and preserve the
leading factorial coefficients.  Every gcd has degree zero.  In particular,
the nonempty barcodes at `p=11,17,43` do not correspond to actual common
factors in these rows.  This is a direct hostile to reversing (18).

## 9. Verification

The primary companion checks:

```text
odd heights:  1<=H<2^18, H odd, T in {5,6,8,9}, 524288 cells;
prime powers: 11<=p<1000 prime, 1<=k<=30,          4920 cells;
suffix table: all 32 odd residues modulo 64;
congruences:  the complete p^k=27 mod 32 odd-exponent table and square family.
```

It verifies the direct submask definition against (5), odd-bit-zero death,
the iff emptiness criterion, the root family, the exact maximal suffix list,
and the strict witness (37).

The independent companion checks:

```text
explicit support sets: odd H<2^15 and four T values,       65536 cells;
prime powers:          primes 11<=p<350, k<=19,             1254 cells;
suffix table:          independently reconstructed from occupied positions;
moment gcds:           five rows at each of two large primes.
```

It also compares its moment coefficient engine with an explicit multinomial
enumeration through order eight and checks the positive polynomial-gcd control

```text
gcd((x+1)^2,(x+1)(x+2))=x+1.                              (41)
```

Normal and optimized runs byte-match the frozen outputs.

Reproduce with

```bash
python3 -B 04-computation/factorial_multiplier_ten_double_overlap_thm4243.py
python3 -B -O 04-computation/factorial_multiplier_ten_double_overlap_thm4243.py
python3 -B 04-computation/factorial_multiplier_ten_double_overlap_independent_audit_thm4243.py
python3 -B -O 04-computation/factorial_multiplier_ten_double_overlap_independent_audit_thm4243.py
```

Byte-compare commands are

```bash
python3 -B 04-computation/factorial_multiplier_ten_double_overlap_thm4243.py | cmp - 05-knowledge/results/factorial_multiplier_ten_double_overlap_thm4243.out
python3 -B -O 04-computation/factorial_multiplier_ten_double_overlap_thm4243.py | cmp - 05-knowledge/results/factorial_multiplier_ten_double_overlap_thm4243.out
python3 -B 04-computation/factorial_multiplier_ten_double_overlap_independent_audit_thm4243.py | cmp - 05-knowledge/results/factorial_multiplier_ten_double_overlap_independent_audit_thm4243.out
python3 -B -O 04-computation/factorial_multiplier_ten_double_overlap_independent_audit_thm4243.py | cmp - 05-knowledge/results/factorial_multiplier_ten_double_overlap_independent_audit_thm4243.out
```

Raw-LF SHA-256 hashes are

```text
primary script      18ada8fca613264c7b04b5e29d22289a20112672ec91c39d224924766ae263cd
primary output      5ade804f3895eb7c2d354ceb8d9a92181612a11eb998ce27e38ad2b9dd31d9ef
independent script  eaf7fee670800f0ce5df344a38921153fbc638cd57df89147ec42e874a813055
independent output  58d315e0d524fe729c06bc785cafbe65445c18578ac40370c2ee60fe1bc2f526
```

The finite computations audit the elementary classification, congruence
tables, and hostile examples.  The unbounded result follows symbolically from
(12)--(18) and the cited proved theorems.

## 10. Scope firewalls and nonclaims

1. **Candidate barcode versus factors.**  Formula (5) exactly classifies the
   THM-3474 intersection of necessary degree barcodes.  A surviving degree is
   not an actual factor.  Only empty intersection is used for coprimality.
2. **Exact criterion versus congruence shadow.**  Equation (6) is the exact
   all-bit criterion for every odd height.  Equation (33) is only the maximal
   list of modulo-64 classes that force closure uniformly.  Heights outside
   that list may close, as (37) proves.
3. **Admissible primes.**  The THM-3474 transfer requires `a=10<p`; hence the
   theorem begins at prime `p=11`.  Statements about arbitrary odd `H` in the
   bit lemma do not silently extend the factorial conclusion to `p<=7`.
4. **Exact support and nonzero coefficients.**  The moment conclusion concerns
   the univariate exact support `{0,1,2}` with all three coefficients nonzero.
   It is a three-slot slice of `SFC(1)`, not `SFC(3)`; see MISTAKE-350.
5. **No arbitrary-support conclusion.**  No claim is made for translated or
   gapped supports, missing slots, arbitrary multipliers, `SFC(1)` in full,
   ambient `SFC(3)`, or the multivariable Factorial Conjecture.
6. **No converse.**  Nonempty `C(p,k)` does not imply a common factor, a common
   complex root, or a bad three-moment window.  The modular hostiles explicitly
   demonstrate this failure of converse.
7. **No use of a bounded audit as an unbounded proof.**  The exhaustive runs
   are verification and regression controls.  The proof of (5)--(6) is the
   carry-free splitting argument, and the factorial transfer is THM-3474 plus
   THM-3124.
8. **Modular validity.**  The hostile gcd primes exceed every tested factorial
   index, preserve both polynomial degrees, and avoid the modular-factorial
   failure genus recorded in MISTAKE-363.
9. **No dependency on barcode completeness elsewhere.**  THM-3483's warning
   that an all-divisor barcode can be incomplete at `d=9996` does not affect
   this theorem.  Here only the proved exact local barcodes and one-way empty
   intersection compiler of THM-3474 are used.
10. **Canon scope.**  The proved statement is the exact multiplier-ten
    compiler closure and its one-way factorial consequence. It does not
    promote a nonempty barcode, a bounded audit, or a suffix residue outside
    its stated quantifier to a factor theorem.

## 11. Duplication and method ledger

Exact searches over canon, current navigation, hypotheses, reflections, and
the mistakes ledger found no existing multiplier-ten factorial theorem or
candidate.  The closest proved mechanism is THM-3474; THM-4237 is the direct
multiplier-six model.  The least-used relevant sidecar is THM-3475's complete
pair-ledger perspective, which reinforces that surviving degrees are addresses
rather than factors.

The selected reusable method cards are: search the statement before the
method; test structured adversaries; use redundant paths as detectors; and
canonicalize mechanics while retaining the predicate.  No new meta-pattern is
proposed from this single thread.
