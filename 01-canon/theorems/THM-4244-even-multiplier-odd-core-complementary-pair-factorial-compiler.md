---
id: THM-4244
title: "Even-multiplier odd-core complementary-pair factorial compiler"
status: >
  PROVED ANALYTIC COMPILER + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every even multiplier a<p at prime-power height H=p^k, the exact
  THM-3474 candidate barcode decomposes into complementary pairs controlled
  solely by the odd core of a. Empty barcode gives rational coprimality and
  an exact three-moment nonvanishing window. A universal binary suffix and
  the multiplier-fourteen mod-eight family give infinite exponent families.
  Candidate survival is not a factor certificate.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3475-factorial-all-divisor-digit-polygon-and-pair-ledger-compiler
  - THM-4237-multiplier-six-binary-adjacency-prime-power-factorial-closure
  - THM-4243-multiplier-ten-double-overlap-prime-power-factorial-closure
primary_script: 04-computation/factorial_even_multiplier_odd_core_compiler_thm4244.py
primary_output: 05-knowledge/results/factorial_even_multiplier_odd_core_compiler_thm4244.out
independent_audit_script: 04-computation/factorial_even_multiplier_odd_core_compiler_independent_audit_thm4244.py
independent_audit_output: 05-knowledge/results/factorial_even_multiplier_odd_core_compiler_independent_audit_thm4244.out
primary_script_sha256: a1249897a322d80d1263b4123cb9e41c059ee9d9fe9781c0413c4347d71429e1
primary_output_sha256: 09e9c3e900d55cb2808159dd622ff30d4096776bf778b97a0b6b4b728238a9bd
independent_audit_script_sha256: 4b6659a93b4b28ae7fb6e77200f0dd560c0820a8308bc7986cab2118e167e9fc
independent_audit_output_sha256: a2f8d99a87d921f6227f1cc7ba01b8a6255d072ad20b225a34f8371b16efd623
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT AFTER ONE SCOPE REPAIR. Hostile review corrected the pure-bit
  cutoff to a/2<=T<=a-1; the factorial cutoff already satisfied it. Four
  normal/-O replays byte-match, an extra all-cutoff sweep checked 266,240
  cells and 5,523,456 pair tests, and an independent third polynomial engine
  reproduced the modular nonconverse. No duplicate, correction, or MISTAKE
  conflict was found.
---

# THM-4244 -- even-multiplier odd-core complementary-pair factorial compiler

**PROVED ANALYTIC COMPILER + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Dependencies and scope

The only proved mathematical dependencies needed for the main statement and
its factorial consequence are:

- `THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families.md`
  (**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED**), for the exact local
  candidate barcode and the one-way empty-barcode-to-coprimality compiler;
- `THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census.md`
  (**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**),
  for the unbounded symbolic resonance reduction.

The closest proved model is
`THM-4237-multiplier-six-binary-adjacency-prime-power-factorial-closure.md`.
The least-used relevant sidecar is
`THM-3475-factorial-all-divisor-digit-polygon-and-pair-ledger-compiler.md`,
which reinforces that a surviving degree is an address, not a factor.
`THM-4243-multiplier-ten-double-overlap-prime-power-factorial-closure.md`
is now a proved and independently audited special case. It is **not** needed
as a dependency of the general proof.

The factorial conclusion below retains the hypotheses of THM-3474:

```text
p >= 5 prime,       k >= 1,       2 <= a < p with a even.       (1)
```

The pure bit lemma is valid more generally for odd `H` and any cutoff `T`
with `a/2<=T<=a-1`, but this theorem makes no factorial claim at `p=3`,
outside the proved scope of THM-3474.

## 2. Exact complementary-block theorem

Put

```text
H=p^k,                N=aH,                T=min(a-1,(p-1)/2).  (2)
```

For nonnegative integers write `x subset_2 y` when every one-bit of `x` is a
one-bit of `y`. Define the exact positive THM-3474 candidate multiplier set

```text
C_a(p,k)={t:1<=t<=T and tH subset_2 aH}.                    (3)
```

The corresponding positive candidate **degree** barcode is `H*C_a(p,k)`.
It is the exact intersection of the two exact local necessary-degree
observers in THM-3474; it is not the set of actual factor degrees.

Choose one representative from each nonfixed complementary even orbit:

```text
R_a={t even: 2<=t<a/2}.                                    (4)
```

For `t in R_a`, define

```text
O_t(H)=(tH)&((a-t)H),
B_t(T)={t,a-t} intersect [1,T].                             (5)
```

Then the exact decomposition is the disjoint union

```text
C_a(p,k) = disjoint_union_(t in R_a, O_t(H)=0) B_t(T).      (6)
```

In particular,

```text
C_a(p,k)=empty
  iff O_t(H)!=0 for every t in R_a.                         (7)
```

Thus the candidate barcode is empty exactly when every nontrivial
complementary even multiplier pair overlaps in binary support. Equation (7)
is vacuous, and correctly true, for `a=2` and `a=4`.

The sharpened normal form writes

```text
a=2^q b,              q=nu_2(a)>=1,              b odd.   (7a)
```

Then the exact candidate set is

```text
C_a(p,k)
 = disjoint_union_(1<=s<=(b-1)/2, (sH)&((b-s)H)=0)
     ({2^q s,2^q(b-s)} intersect [1,T]).                   (7b)
```

Consequently its emptiness depends only on the odd core `b` and the height
`H`, not on the power `2^q`:

```text
C_a(p,k)=empty
  iff (sH)&((b-s)H)!=0 for every 1<=s<=(b-1)/2.            (7c)
```

The scale `2^q` changes visible multiplier labels and reset truncation, but
only left-shifts every reduced overlap and cannot change closure.

### Proof

For `0<=x<=y`,

```text
x subset_2 y    iff    x & (y-x)=0.                         (8)
```

Indeed, a submask is removed without borrowing, while disjoint supports add
without carrying. Apply (8) with `x=tH` and `y=aH`:

```text
tH subset_2 aH    iff    (tH)&((a-t)H)=0.                  (9)
```

Because `a` is even and `H` is odd, an odd `t` makes both `tH` and
`(a-t)H` odd. They overlap at bit zero, so no odd multiplier survives.

Complementation `t -> a-t` is an involution on the even multipliers in
`[1,a-1]`. Its nonfixed orbits are the pairs represented uniquely by (4),
and (9) is constant on each orbit. Its sole fixed point in `[1,a-1]` is

```text
s=a/2.                                                      (10)
```

The fixed point never survives, because

```text
(sH)&((a-s)H)=(sH)&(sH)=sH!=0.                             (11)
```

If `a=2 mod 4`, this fixed point is odd and was already killed at bit zero;
if `a=0 mod 4`, it is an even self-complement and (11) is the necessary
extra check.

It remains to audit reset truncation. Since `p>a` and `p` is odd,

```text
(p-1)/2 >= a/2,       hence       T>=a/2.                  (12)
```

Therefore the fixed point and the lower representative `t<a/2` of every
nontrivial orbit are always in the reset range. The upper representative has
the exact visibility criterion

```text
a-t <= T    iff    a-t <= (p-1)/2
            iff    p >= 2(a-t)+1.                          (13)
```

Thus truncation can change a surviving block from two visible multipliers to
one, but can never erase the whole block. Combining (9)--(13) proves (6) and
both directions of (7). QED.

### Proof of the odd-core normal form

Write `q=nu_2(a)`. If `t in R_a` has `nu_2(t)<q`, then

```text
nu_2(a-t)=nu_2(t).                                         (14)
```

Since `H` is odd, both products in (5) have the same least occupied bit, so
`O_t(H)!=0` automatically. Hence only representatives divisible by `2^q`
need a genuine all-bit test. After writing

```text
a=2^q b,       t=2^q u,       b odd,                       (15)
```

the remaining overlap is a common shift of

```text
(uH)&((b-u)H),                                              (16)
```

where the complementary reduced multipliers sum to the odd core `b`. More
precisely,

```text
(2^q uH)&(2^q(b-u)H)=2^q*((uH)&((b-u)H)).                 (16a)
```

The representatives `2^q u<a/2` are exactly
`1<=u<=(b-1)/2`, and the exact reset block is the one displayed in (7b).
This proves the normal form and invariance under adjoining powers of two.
This recovers THM-3474's power-of-two multiplier family immediately: if
`a=2^q`, no positive `t<a/2` is divisible by `2^q`, so every nontrivial pair
is automatically overlapping.

## 3. Coprimality and exact quadratic window

Use

```text
L(x^r)=r!,
A_n^d(v)=L((d-x+v*x^2)^n),
d=N+1,
F=A_(N-1)^d,       G=A_N^d.                                (17)
```

THM-3474 proves that every positive rational common-factor degree of `(F,G)`
must lie in `H*C_a(p,k)`, and separately proves that the common coordinate-root
capacity is zero. Consequently (7) gives

```text
if O_t(H)!=0 for every t in R_a, then gcd_Q(F,G)=1.         (18)
```

Now let

```text
f(x)=alpha+beta*x+gamma*x^2,       alpha*beta*gamma!=0.     (19)
```

If its moments at `N-1,N,N+1` all vanished, THM-3124's symbolic recurrence at
window start `r=N-1` would force

```text
beta/alpha=-1/(N+1).                                       (20)
```

Put `v=(N+1)gamma/alpha`. Then

```text
f(x)=(alpha/(N+1))*((N+1)-x+v*x^2),                        (21)
```

so the first two vanishing moments make `v` a common complex root of `F` and
`G`. A common complex root of rational polynomials has a positive-degree
minimal polynomial common to both, contradicting (18). Therefore

```text
at least one of L(f^(ap^k-1)), L(f^(ap^k)), L(f^(ap^k+1))
is nonzero whenever every pair in (7) overlaps.             (22)
```

The factor `N+1` in the definition of `v` is load-bearing.

## 4. A universal suffix certificate

Let `ell>=1` satisfy

```text
2^(ell-1) >= b,          where a=2^q b and b is odd,        (23)
```

and let any odd height `H` satisfy

```text
H=-1 mod 2^ell.                                             (24)
```

Then every nontrivial complementary even pair overlaps, so (7), (18), and
(22) all hold in the factorial setting.

### Proof

Work modulo `2^ell` on the reduced odd-core pairs. For every
`1<=s<=(b-1)/2`, equation (24) gives

```text
sH=-s mod 2^ell,        (b-s)H=-(b-s) mod 2^ell.           (25)
```

```text
0<s<b<=2^(ell-1),       0<b-s<b<=2^(ell-1).                (26)
```

Thus the least nonnegative residues in (25) are `2^ell-s` and
`2^ell-(b-s)`, both strictly between `2^(ell-1)` and `2^ell`. Their bit
`ell-1` is one, so the reduced products overlap there. Equation (16a) shifts
that common bit to `q+ell-1` in the original multiplier pair. This holds for
every reduced pair. QED.

In prime-power coordinates, a sufficient family is

```text
p=-1 mod 2^ell,       p>a prime,       k odd.              (27)
```

For each fixed admissible `p`, the infinitely many odd `k` give infinitely
many exact quadratic windows. For example, at `a=14` one has `q=1,b=7`;
choose `ell=4` and `p=31`; every odd `k` closes.

More generally, a residue `r mod 2^ell` is a sufficient suffix certificate
whenever the low `ell`-bit residues of `sr` and `(b-s)r` overlap for every
`1<=s<=(b-1)/2`. The original overlaps are their common `q`-bit shifts. This
is sufficient, not an all-bit classification.

## 5. The multiplier-fourteen family

For `a=14`, the nontrivial blocks are

```text
{2,12},              {4,10},              {6,8}.           (28)
```

The exact reset cases are

```text
p=17:     T=8,   visible blocks {2}, {4}, {6,8};
p=19:     T=9,   visible blocks {2}, {4}, {6,8};
p=23:     T=11,  visible blocks {2}, {4,10}, {6,8};
p>=29:    T=13,  all three complete blocks visible.         (29)
```

If `H=7 mod 8`, every pair overlaps. Since the multipliers are even, their
products modulo `16` depend only on `H mod 8`; at the residue `7`,

```text
(2H mod16)&(12H mod16)=14&4 =4,
(4H mod16)&(10H mod16)=12&6 =4,
(6H mod16)& (8H mod16)=10&8 =8.                            (30)
```

Therefore every prime

```text
p>14,       p=7 mod 8,       k odd                        (31)
```

gives coprimality (18) and the nonzero-window conclusion (22) at multiplier
fourteen. The explicit fixed-prime family `p=23`, all odd `k`, is already
infinite in the exponent.

## 6. Equality boundaries and hostiles

1. **Vacuous small multipliers.** At `(a,p,k)=(2,5,1)` and `(4,5,1)`, there
   are no nontrivial complementary even pairs; the candidate set is empty.
   For `a=4`, the visible even fixed point `t=2` dies by self-overlap.

2. **Self-complement is never a singleton candidate.** At
   `(a,p,k)=(12,17,1)`, `T=8`. The block `{2,10}` overlaps, `{4,8}` is
   disjoint, and the exact candidate set is `{4,8}`. The visible self-point
   `t=6` is absent.

3. **Reset truncation changes block size, not emptiness logic.** At
   `(a,p,k)=(14,17,1)`, `H=17,T=8`; all three overlaps are zero, but the exact
   set is `{2,4,6,8}`, not all six even multipliers. The upper members `12,10`
   are outside the reset range while every lower representative remains.

4. **Pair coordinates are independent.** At `(a,p,k)=(14,29,2)`, `H=841`,
   the three overlaps are `(1536,0,4608)`, leaving exactly the middle block
   `{4,10}`.

5. **Odd exponent is load-bearing in (31).** At `(a,p,k)=(14,23,2)`,
   `H=529,T=11`; all three overlaps vanish and the truncated candidate set is
   `{2,4,6,8,10}`. No all-exponent extension follows from `p=7 mod 8`.

6. **Suffix certificates are not necessary.** At `(a,p,k)=(14,19,1)`, all
   three overlaps are nonzero `(36,12,16)` and the barcode closes although
   `H=3 mod 8`, outside (31). At `H=23`, closure also lies outside the
   universal `-1 mod 16` certificate.

7. **Nonempty barcode does not construct a factor.** The exact hostile row
   `(a,p,k)=(14,17,1)` has `N=238,d=239` and candidate set `{2,4,6,8}`. Two
   independent exact coefficient engines computed the actual pair
   `A_237^239,A_238^239` over each of `F_1000003` and `F_1000033`; in both
   fields the degrees remain `237,238` and the gcd degree is zero. Both primes
   exceed `2N`, so all leading factorial coefficients survive. Any common
   rational factor would retain positive degree modulo either prime, a
   contradiction. Hence this row has `gcd_Q=1` despite its nonempty candidate
   barcode. This is a `FINITE-EXACT` hostile to reversing the compiler, not an
   ingredient in the unbounded proof.

## 7. Exact verification

The primary checker uses integer bit operations and a hand-written modular
polynomial Euclidean algorithm. It verifies:

```text
43,104 prime-power barcode cells,
14,172 reset-truncated cells,
43,104 fixed-point checks,
247,152 automatic valuation-pair checks,
16,384 sharpened odd-core suffix heights / 341,376 reduced-pair checks,
26,112 power-of-two invariance cells / 208,896 exact overlap shifts,
32,768 multiplier-fourteen suffix heights,
770 prime/odd-exponent multiplier-fourteen cells,
the nine named boundary/hostile rows,
and the two-prime modular non-converse row.                 (32)
```

The independent checker uses explicit sets of occupied bit positions,
constructs complement orbits directly, expands the moment polynomial by
three-slot multinomial counts, checks that expansion against a word-by-word
state engine through order eight, and delegates polynomial gcd to SymPy. It
verifies:

```text
23,715 independently shaped prime-power cells,
74,196 one-visible and 248,202 two-visible block checks,
9,216 sharpened odd-core suffix heights / 142,944 reduced-orbit checks,
12,480 power-of-two invariance cells / 74,880 explicit support shifts,
16,384 multiplier-fourteen suffix heights,
504 prime/odd-exponent cells,
the same nine hostile rows,
and the same two modular gcd certificates.                 (33)
```

Neither script contains a Python optimization-removable truth gate. Normal
and `-O` runs are byte-identical. The two canonical transcripts and all raw
file hashes are frozen in the repository.

Reproduce from the repository root with

```bash
python3 -B \
  04-computation/factorial_even_multiplier_odd_core_compiler_thm4244.py
python3 -B -O \
  04-computation/factorial_even_multiplier_odd_core_compiler_thm4244.py
python3 -B \
  04-computation/factorial_even_multiplier_odd_core_compiler_independent_audit_thm4244.py
python3 -B -O \
  04-computation/factorial_even_multiplier_odd_core_compiler_independent_audit_thm4244.py
```

Each primary stream must byte-match
`05-knowledge/results/factorial_even_multiplier_odd_core_compiler_thm4244.out`;
each independent stream must byte-match
`05-knowledge/results/factorial_even_multiplier_odd_core_compiler_independent_audit_thm4244.out`.

## 8. Duplication audit

Exact searches covered canon theorem files, current navigation, active
hypotheses, reflections, and targeted mistake entries, using the theorem IDs,
`binary submask`, `carry-free splitting`, `complementary multiplier`,
`self-complement`, `a-t`, reset truncation, and the suffix congruences.

- The proved general input is THM-3474, which states the direct bit test but
  does not orbit-decompose it or handle the all-even empty-iff criterion.
- THM-4237 is exactly the `a=6` one-pair specialization.
- Proved THM-4243 is exactly the `a=10` two-pair specialization.
- No earlier proved/candidate file contains (6)--(7), the self-point and exact
  visibility audit (10)--(13), the odd-core reduction (14)--(16), the universal
  suffix theorem (23)--(27), or the multiplier-fourteen family (28)--(31).


## 9. Nonclaims and firewalls

- Formula (6) classifies the exact THM-3474 **candidate intersection**, not
  actual factor degrees. Only the empty direction is used.
- A nonempty candidate barcode does not imply a common factor, common root, or
  bad moment window; hostile 7 proves the converse false.
- The general suffix criterion and the `a=14` congruence are sufficient
  shadows. They are not necessary and do not replace the exact all-bit test.
- The factorial result requires `p>=5`, prime-power height, even `a<p`, and the
  exact support `{0,1,2}` with all coefficients nonzero.
- The moment conclusion is a three-slot slice of ambient `SFC(1)`, not
  `SFC(3)`, arbitrary-support `SFC(1)`, or the multivariable Factorial
  Conjecture (MISTAKE-350 firewall).
- No claim is made for translated/gapped supports, missing slots, arbitrary
  odd multipliers, or a converse to THM-3474.
- The modular hostile uses primes larger than every factorial index and
  preserves degrees, respecting the MISTAKE-363 modular-factorial firewall.
- Proved THM-4243 is a corroborating special case, not a dependency.
