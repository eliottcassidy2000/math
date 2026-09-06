# Signed (1,1,2) sharp one-ray classification

Status: **PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT +
INDEPENDENTLY AUDITED (THM-4444); LRC(14) OPEN.**

For every primitive sorted distinct positive ternary-unit triple with a signed
coefficient-magnitude \((1,1,2)\) relation, the coefficient-two coordinate is
pointwise the physical envelope. Consequently
\[
 \min_iE_i=\mu(F)=E_{\text{coefficient-two}}\le {11\over140},
\]
sharply at \((2,11,20)\).

The complete \(6/77\) locus is:
\[
 \min E>6/77\iff w=(2,11,20),\qquad
 \min E=6/77\iff w=(1,5,11).
\]

## Three cones

| cone | ray | sharp value | equality |
|---|---|---:|---|
| \(c=2a+b\) | \((2,1,-1)\) | \(58/833\) | \((5,7,17)\) |
| \(c=2b-a\) | \((1,-2,1)\) | \(11/140\) | \((2,11,20)\) |
| \(c=a+2b\) | \((1,2,-1)\) | \(6/77\) | \((1,5,11)\) |

The carrier-zonotope covector is especially rigid here:
\[
 |\phi(C)|<(3/14)\|u\|_1=6/7<1,
\]
so every integral carrier lies on \(\mathbb Zu\), without the mod-three
repair needed at coefficient norm five. Exact profile integration gives the
same positive-half selected and physical integral \(9/196\) in every cone.
Deleted-third sampling then yields
\[
 {3\over49}-{4\over7c}\le\min E=\mu(F)
 <{3\over49}+{4\over7c}.
\]
Thus all heights \(c\ge35\) are strictly below \(6/77\).

The exact head through height 68 has 361 rows. A height-611 extension has
28,438 rows, split \(9,482/14,214/4,742\), and reproduces the same sharp
leaders.

## Reproduction

```powershell
python -B 04-computation/lrc14_signed_112_sharp_classification_thm4444.py
python -B -O 04-computation/lrc14_signed_112_sharp_classification_thm4444.py
python -B 04-computation/lrc14_signed_112_sharp_classification_thm4444_independent.py
python -B -O 04-computation/lrc14_signed_112_sharp_classification_thm4444_independent.py
```

See [THM-4444](../../01-canon/theorems/THM-4444-lrc14-signed-112-sharp-one-ray-classification.md)
for the proof and scope.
