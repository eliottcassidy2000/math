---
id: THM-3343
title: "Shifted-donor rotation bisects exactly the dyadic critical-run annuli"
status: >
  PROVED + VERIFIED-EXACT.  For 1<=d<M, the length-M words whose first d
  bits are equal but which are nonconstant admit an exact bisection in every
  Hamming layer if and only if d and M are powers of two.  When they are,
  label every critical value n>d heads for n-d odd and tails for n-d even;
  left rotation injects the resulting layer defect into the critical-d donor
  branch, which has exactly enough capacity for a deterministic repair.  For
  consecutive dyadic boundaries this gives one uniform fair extractor with
  T(n)=n+1 off the powers of two and T(n)=2n on the powers of two.  Under the
  prescribed alternating labels, the dyadic donor deadline M is sharp.
source: codex-kps-2026-08-12-shifted-donor-handoff
depends_on: []
related:
  - THM-2160-dyadic-checksum-extracts-a-fair-bit-under-the-critical-run-deadline
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors
  - THM-3342-sublinear-deadline-excess-is-impossible-for-fair-critical-run-extractors
  - THM-3344-orientation-splitting-saves-exactly-one-dyadic-donor-bit
script: 04-computation/amm12592_shifted_donor_dyadic_annuli_thm3343.py
output: 05-knowledge/results/amm12592_shifted_donor_dyadic_annuli_thm3343.out
script_sha256: 6e27307311e093f96b6a935392d6fbb13b72ddfb2344124e6c996169f1ed9a40
output_sha256: fa44f18891fdfc59caa27d82d5cb0217f30feae773d5bef9973c5c68cd944dba
hash_basis: working-tree bytes (LF)
---

# THM-3343 -- shifted donors and dyadic annuli

Let `1<=d<M`.  For a nonconstant binary word `x` of length `M`, write
`n(x)` for its initial constant-run length and put

```text
S_(d,M)={x:n(x)>=d}.                                                (1)
```

Thus the first `d` bits are equal and the first `M` bits are not constant.
The `S_(2^r,2^(r+1))` partition every nonconstant infinite stream.

## 1. Exact parity characterization

The Hamming-layer enumerator of `S_(d,M)` is

```text
B_(d,M)(z)=(1+z)^(M-d)(1+z^d)-1-z^M.                              (2)
```

Indeed, the two initial orientations contribute `(1+z)^(M-d)` and
`z^d(1+z)^(M-d)`, and the two constant words are deleted.  A Hamming layer
can be bisected only if its size is even.  All interior coefficients in (2)
are even exactly when

```text
(1+z)^(M-d)(1+z^d)=1+z^M                    in F_2[z].              (3)
```

**Lemma.** Equation (3) holds if and only if `d` and `M` are powers of two.

*Proof.* Sufficiency is Frobenius.  For necessity, write

```text
d=2^a e,       M=2^b f,       e,f odd,
U_j(z)=(1+z^j)/(1+z)          for odd j.                            (4)
```

The polynomials `U_e,U_f` are squarefree and do not vanish at `z=1`.  Unique
factorization in (3) gives

```text
(1+z)^(M-d+2^a) U_e(z)^(2^a)
  =(1+z)^(2^b) U_f(z)^(2^b).                                      (5)
```

If either `U` is nonconstant, equality of irreducible-factor multiplicities
forces `a=b`; it then forces `U_e=U_f`, hence `e=f` and `d=M`, a
contradiction.  Therefore `e=f=1`, so both integers are powers of two. QED

This proves necessity not merely for our construction but for every
composition-exact coloring of the annulus.

## 2. The shifted rotation defect

Assume henceforth that `d` and `M` are powers of two.  For every critical
value strictly above the lower boundary, prescribe

```text
n-d odd:       heads,
n-d even:      tails,                           d<n<M,              (6)
```

independently of the initial-bit orientation.  Fix a Hamming weight `w` and
let

```text
A_w = number of critical-d donor words,
E_w = number prescribed heads in (6),
O_w = number prescribed tails in (6),
D_w = E_w-O_w.                                                       (7)
```

**Rotation lemma.** For every layer,

```text
0<=D_w<=A_w.                                                        (8)
```

*Proof.* Left rotation by one bit is a weight-preserving bijection

```text
{n-d positive and even}
  <--> {n-d positive and odd, and x_M=x_1}.                         (9)
```

Indeed it lowers the critical value by one and appends the old initial bit;
right rotation is the inverse.  Hence `D_w` counts precisely the prescribed
heads ending in the bit opposite their initial bit.

For such an unmatched word of critical value `n`, rotate left by `n-d`
places.  The image has critical value `d`.  It ends in exactly `n-d` copies
of its initial bit, preceded by the opposite bit `x_M`; this terminal run
recovers `n-d`, and rotating back recovers the source.  The map is therefore
an injection into the donor layer, proving (8). QED

## 3. Exact donor repair

Let `B_w=A_w+E_w+O_w` be the full layer size.  Section 1 makes `B_w` even.
Define

```text
r_w=B_w/2-E_w=(A_w-D_w)/2.                                        (10)
```

It is an integer, and (8) gives

```text
0<=r_w<=A_w.                                                       (11)
```

Order the critical-`d` words of weight `w` lexicographically and declare the
first `r_w` of them heads.  Together with (6), the number of heads in layer
`w` is exactly

```text
E_w+r_w=B_w/2.                                                     (12)
```

This is a literal deterministic coloring and proves sufficiency in Section
1.  For calculation without word enumeration, the number of weight-`w`
words of critical value `n` is

```text
C(M-n-1,w-1)+C(M-n-1,w-n),                                        (13)
```

with out-of-range binomial coefficients zero.

## 4. One uniform extractor

Use the construction on every consecutive annulus

```text
S_(2^r,2^(r+1)),                r>=0.                              (14)
```

Each Hamming layer is bisected, so heads and tails have equal Bernoulli mass
inside every annulus for every unknown `p`.  The annuli partition the
nonconstant streams, while the two constant rays have probability zero.
The resulting extractor is exactly fair and pathwise total.  Its deadline is

```text
T(n)=2n,          n a power of two,
T(n)=n+1,         otherwise.                                      (15)
```

At an interior value, (6) decides as soon as the first disagreement reveals
`n`.  At the lower boundary, the lexicographic donor repair reads through the
upper endpoint.

This profile is Pareto-incomparable with THM-2160: it meets every non-dyadic
pointwise floor, but pays the full factor two at the sparse dyadic donors.
It makes the uniform obstruction concrete: the cost is concentrated on a
zero-density handoff spine rather than spread over every critical value.
THM-3344 subsequently splits the two upper-boundary orientations and saves
one donor bit while retaining every non-dyadic floor.

## 5. Sharpness of the donor inside this architecture

Let

```text
C_n(z)=(z+z^n)(1+z)^(M-n-1),
Delta(z)=sum_(n=d+1)^(M-1) (-1)^(n-d+1) C_n(z).                    (16)
```

The required donor-head enumerator is

```text
R(z)=(C_d(z)-Delta(z))/2=sum_w r_w z^w.                            (17)
```

At `z=-1`, only the last critical row can survive.  Directly,

```text
R(-1)=-1,        d=1,
R(-1)= 1,        d>=2.                                             (18)
```

In particular `1+z` does not divide `R`.  If every donor word were decided
before bit `M`, refining its terminal prefix to length `M` would multiply
each contribution by at least one factor `1+z`; their sum `R` would be
divisible by `1+z`, contradicting (18).  Thus some critical-`d` stream must
read bit `M`.  The donor deadline in (15) is sharp once (6) is fixed.

## 6. Exact audit

```bash
python 04-computation/amm12592_shifted_donor_dyadic_annuli_thm3343.py
python -O 04-computation/amm12592_shifted_donor_dyadic_annuli_thm3343.py
```

Both runs pass.  The companion checks all binomial repair inequalities for
45 dyadic pairs through `M=512`; the iff parity characterization on all 2016
pairs `1<=d<M<=64`; both rotation maps and the literal coloring on 107,432
words for every dyadic pair through `M=16`; and the nonzero value (18) for all
45 large pairs.  The proof above is uniform; the computation is a hostile
referee. QED.
