# Independent audit: golden period fifteen, the 341 quotient, and source mass

**Verdict: ACCEPTED, with the explicit typing guards below.** The ring and map identities are PROVED by the elementary calculations here. The complete 1364-cell and 1364-word checks are FINITE-EXACT independent controls. This is a bounded supplement to the already frozen arithmetic packet; that packet was not changed.

## 1. Exact ring factorization, not an order coincidence

Let `O=Z[phi]`, `phi^2=phi+1`, and `R_15=O/(phi^15-1)`. Since

```text
phi^15-1 = 376+610phi,
M = [[376,610],[610,986]],    det(M)=-1364,
```

this quotient has `1364=4*11*31` elements. The determinant is negative; `1364` is its absolute value.

There are three ring evaluations:

```text
O -> F_4,   phi -> alpha,  alpha^2+alpha+1=0;
O -> F_11,  phi -> 4;
O -> F_31,  phi -> 19.
```

Their kernels are `(2)`, `(11,phi-4)`, `(31,phi-19)` and are pairwise comaximal. The evaluations of `phi` have respective orders `3,5,15`, so each kills `phi^15-1`. CRT makes the product evaluation onto a ring of exactly `1364` elements. Since its kernel contains `(phi^15-1)` and has the same index, the two ideals coincide. Thus

```text
R_15 ~= F_4 x F_11 x F_31,     phi -> (alpha,4,19).       (A)
```

The evaluation `phi -> 81 modulo341` is a ring map: `81^2-81-1=19*341` and `81^15=1 mod341`. Under (A) and ordinary integer CRT it forgets the `F_4` coordinate. It is onto, every fibre has four elements, and its kernel is exactly `F_4 x {0} x {0}`.

**Clock typing guard.** Scalar `4` has coordinates `(0,4,4)`. It is **not a unit or an automorphism of the full ring R_15**: its image has only `341` cells and its kernel has four. It becomes a unit clock on the `Z/341` quotient. Golden time `phi` is a unit already on the full source ring. The joint group `<4,81>` belongs to the quotient clock system, not to two invertible full-ring clocks.

## 2. The actual eleven-source lift lies on an axis

The period-five quotient is `R_5=O/(2+5phi) ~= F_11`, with `phi -> 4`. Its repetition into period fifteen uses

```text
Q=1+phi^5+phi^10=38+60phi,
(2+5phi)Q=376+610phi,
Q -> (0,3,0) under (A).
```

Consequently the source-preserving additive lift is

```text
J(a)=Qa -> (0,3a,0),    a in F_11,
projection(J(a)) = 124a mod341.                         (B)
```

Here `a` is the original **golden label**, not the literal rational Collatz source value. The eleven projected labels are exactly the multiples of `31`, including zero. They are not the 75 unit square residues and not the full 341-cell quotient.

The golden quotient clock is multiplication by `81`. On these selected labels, it is identical to multiplication by `4`, because both reduce to `4` modulo11 and the modulo31 coordinate is zero. The full group

```text
K=<4,81>=QR_11 x QR_31,  |K|=75
```

therefore acts on the eleven-source image only through its order-five `QR_11` quotient. Its action kernel is the order-fifteen `QR_31` factor. Each nonzero source orbit has size five; zero is fixed. The original eleven cells do not see the additional 31-clock coordinate.

**Multiplication guard.** `J` is additive and O-linear, not a ring embedding:

```text
Q^2=3Q,
J(a)J(b)=J(3ab),
Q^2 != Q.
```

The image is an ideal whose local identity is `e_11=4Q`, with coordinates `(0,1,0)`. Rescaling by four can give a nonunital embedding of the field into that ideal, but it changes the chosen source-label map and must not be performed silently when weights are attached.

## 3. The exact weighted bridge, and an explicit hostile to uniformization

For any source mass `mu` on the 341 labels, the arithmetic packet's square-residue ceiling observation is the linear functional

```text
L(mu)=sum_{x mod341} mu(x)*[-x^2]_341.
```

The bracket formula uses **counting mass on all 341 residues**. It does not use arbitrary source mass. Exact independent values are:

| Original source and its transported mass | Total mass | L(mu) | Mean |
|---|---:|---:|---:|
| all 341 quotient labels, counting mass | 341 | 59334 | 174 |
| all 1364 ring cells, then quotient pushforward | 1364 | 237336 | 174 |
| original eleven cells transported by (B) | 11 | 2046 | 186 |

Thus even normalizing total mass does not turn the original eleven sources into the ambient residue distribution. The means differ by twelve. The full ring contributes four copies of every quotient label, whereas the eleven-source lift contributes one copy of every label on one conductor axis.

For arbitrary original weights `w_a`, the lawful replacement is

```text
L(J_*w)=sum_{a in F_11} w_a*[-(124a)^2]_341.
```

It is not a count over `K`, and it is not ambient averaging. The checker also tests a nonuniform weight vector to keep this interface literal.

This supplies a source-preserving connection between the golden ring, the old 341 clock, and the new weighted bracket observation. It supplies no equality of their original measures, no original-doubling explanation for `{2,3,11}`, and no new prime-count or Collatz-convergence result.

## 4. The four-cell kernel is the period-three image

The repetition multiplier from period three to fifteen is

```text
Q_(3,5)=1+phi^3+phi^6+phi^9+phi^12=117+188phi.
```

Under (A), it has coordinates `(1,0,0)`: `phi^3=1` on the F4 factor, and the other two geometric sums vanish. Also `R_3=O/(phi^3-1)=O/(2phi) ~= F_4`, because `phi` is a unit. Therefore `J_(3,5)` has image exactly the four-cell kernel of the projection to 341. Here its multiplier is idempotent and gives the local identity of that ideal.

An independent exhaustive census of all `1364` legal rooted binary words of length15 confirms that the only words whose golden labels project to zero are

| Word | Rational standard-map fixed point |
|---|---:|
| `000000000000000` | 0 |
| `001001001001001` | 4 |
| `010010010010010` | 2 |
| `100100100100100` | 1 |

The rational decoder independently composes the affine branches `x/2` and `3x+1`, then solves their exact fixed-point equation. The label reader uses polynomial Horner evaluation modulo `phi^15-1`; it does not infer labels from the rational values.

This identifies the kernel quartet with the standard-map sources `{0,1,2,4}`. It does **not** make another fibre into a literal additive coset of those four rational numbers. Ring addition is addition of golden labels, and the rational source decoder has not been asserted additive.

## 5. Independent implementation and reproduction

The exact quotient-ring computation does not use the producer's CRT model as its representation. It derives the integer ideal-lattice basis

```text
(682,0), (-162,2)
```

from the multiplication matrix using the unimodular coefficient matrix with determinant `-1`. Its `1364` canonical representatives are `(a,b)` with `0<=a<682`, `b in {0,1}`. For arbitrary `(a,b)`, put `b=2k+r`; reduction is

```text
(a+b*phi) -> ((a+162k) mod682) + r*phi.
```

Multiplication is first performed in `Z[phi]` and then reduced in this lattice. The CRT coordinates, quotient projection, both clock actions, all source cells, all eleven lift cells, all 121 source products, all 75 joint multipliers, and all 1364 legal rooted words are checked against this independent ring representation.

```bash
python3 check_crt.py > receipt.replayed.json
python3 -O check_crt.py > receipt.replayed.optimized.json
cmp receipt.replayed.json receipt.normal.json
cmp receipt.replayed.optimized.json receipt.optimized.json
sha256sum -c SHA256SUMS
```

There are **8387 explicit checks**, using exceptions rather than removable assertions. The normal and optimized receipts are byte-identical. The code uses only the standard library and has no external file or absolute-path dependencies.

## 6. Producer review receipt

The producer supplement `session74-clocks-supplement/RESULT.md` was reviewed in full and accepted, including its alternative principal-ideal CRT proof, exact factorization, rooted-source interpretation, phase-sensitive weight hostile, and all four typing guards. Its `verify.py` was independently rerun and byte-matched the retained receipt. `PRODUCER-REVIEW.json` records the exact candidate bytes reviewed; a later status-only promotion does not change the mathematical content audited here. No substantive error was found.
