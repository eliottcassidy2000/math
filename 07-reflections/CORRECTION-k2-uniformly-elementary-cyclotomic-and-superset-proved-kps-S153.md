---
source: kind-pasteur-2026-07-24-S153 (Opus 4.8)
status: MAJOR CORRECTION + RESULT. (1) RETRACTS kps-S151 (the "signature-specific, C_{1/5}={1}" claim was a
  FALSE NEGATIVE: S_{1/5}(2) IS elementary). (2) RETRACTS kps-S152 section 4 (the "no uniform k=2 formula"
  obstruction). (3) PROVES the correct result: S_a(2) is UNIFORMLY elementary for all a != 1/2 via a cyclotomic
  (Gauss-digamma) closed form, derived cleanly from Mehler-Dirichlet. This PROVES the superset direction
  {1,...,P(n)-1} subset C_{1/n} for all four signatures. (4) Consequently the prime-power law's NON-elementary
  side (k>=P(n)) rests on PSLQ nulls of the same unreliable type and is NOT established.
tags: [hypergeometric, closed-forms, cyclotomic, digamma, correction, false-negative, superset, series]
related: [kps-S146, kps-S150, kps-S151, kps-S152, verify-inherited-residues]
corrects: [kps-S151, kps-S152]
retracts: [kps-S151 signature-specificity, kps-S152 section 4]
---

# CORRECTION: k=2 is uniformly elementary (cyclotomic); superset proved; prime-power hard-half unsupported

## 1. The error (a textbook false negative)
kps-S151 claimed `S_{1/5}(2)` non-elementary (hence `C_{1/5}={1}`, "the law is signature-specific") from a PSLQ
**null against a golden-field basis**. **That was a false negative.** The true closed form lives at level 20
(10th-roots-of-unity constants), which the golden basis does not span. Concretely
`S_{1/5}(2)=(4 sin(pi/5)/pi)*(25/6) int_0^1 (t-t^7)/(1+t^{10})dt` -- an integral of a rational function, hence a
`Q`-combination of `log`s and `arctan`s. **`S_{1/5}(2)` is elementary.** (Verified 80+ digits.)

This is exactly the failure my own memory warns about (`verify-inherited-residues`: "report the detection floor
next to any negative"): a PSLQ null certifies non-elementarity **only against a complete basis for the relevant
field**. I did not establish basis-completeness. Retract kps-S151 in full.

## 2. The correct result: S_a(2) is uniformly elementary (proof)
> **THEOREM.** For `0<a<1`, `a != 1/2`, with `c=1-2a`:
> `S_a(2)=int_0^1 2F1(a,1-a;1;x^2)dx = (4 sin(pi a)/pi) sum_{n>=0} (-1)^n/((2n+1)^2-c^2)`
> `= (sin(pi a)/(2 pi c)) [psi((3-c)/4)-psi((1-c)/4)-psi((3+c)/4)+psi((1+c)/4)]`,
> which by **Gauss's digamma theorem** is a `Q`-combination of `pi`, `log`(algebraic), `arctan`(algebraic) for
> every rational `a`. (Verified to 115 digits at `a=1/3,1/4,1/5,1/6`.)

**Proof (clean, Mehler-Dirichlet).**
1. `2F1(a,1-a;1;t)=P_{-a}(1-2t)` (Legendre), so `S_a(2)=int_0^1 P_{-a}(1-2x^2)dx = int_0^{pi/2}P_{-a}(cos 2th)cos th dth` (`x=sin th`).
2. Mehler-Dirichlet `P_{-a}(cos psi)=(sqrt2/pi)int_0^psi cos((1/2-a)u)/sqrt(cos u-cos psi)du`; substitute, swap
   order, and the inner `th`-integral evaluates to `sqrt2 log(cot(t/4))`. Result (verified to 0.0):
   > `S_a(2)=(1/pi) int_0^pi cos((1/2-a)t) log(cot(t/4)) dt`.
3. Fourier `log(cot(t/4))=2 sum_{n>=0} cos((2n+1)t/2)/(2n+1)`; integrate term-by-term -> the sum. QED.

**Degeneration.** At `a=1/2` (`c=0`) the sum is `sum(-1)^n/(2n+1)^2 = Catalan G` -- the one (conjecturally)
non-elementary point at `k=2`. So `k=2` is elementary for **all `a` except `a=1/2`**.

**It specializes to the three signature values** (via Gauss's digamma theorem, verified `d=0.0`):
`a=1/4 -> 4log(1+sqrt2)/pi`; `a=1/3 -> 3sqrt3 log2/pi`; `a=1/6 -> 3sqrt3 log(2+sqrt3)/(2pi)`.
So the "three hard separate derivations" of kps-S152 are ONE uniform theorem. (kps-S152 section 4, which argued
no uniform formula could exist, is retracted: it inferred that from the false premise `S_{1/5}(2)` non-elem.)

## 3. Consequence A -- the SUPERSET direction is now fully proved
`{1,...,P(n)-1} subset C_{1/n}` for all four signatures, from three proved facts:
- `k=1`, all `a`: `S_a(1)=sin(pi a)/(pi a(1-a))` (Gauss). [covers n=2,3,4,6]
- `k=2`, all `a != 1/2`: section 2 theorem. [covers the k=2 need of n=3,4,6]
- `k=3`, `lambda=1/4`: `S_{1/4}(3)=-(sqrt3 log(5-2sqrt6)+2 arctan(sqrt2/5))/pi` (verified positive). [covers n=4]

n=2 needs only k=1; n=3,6 need k=1,2; n=4 needs k=1,2,3 -- **all covered. Superset: DONE.** (The user's request.)

## 4. Consequence B -- the prime-power law's HARD half is unsupported
The law's non-elementary claims -- `S_{1/3}(3), S_{1/6}(3)` non-elem (kps-S149), `S_{1/4}(4)=S(4)` non-elem
(kps-S146/S148) -- ALL rest on PSLQ nulls of exactly the type that just failed. Re-tests against level-9/16
cyclotomic bases still find no relation, but a **null is now uninformative** (it was for `S_{1/5}(2)` too). So:
> **The prime-power law is NOT established.** Its easy half (k=1,2 elementary, and `k=3` for `lambda=1/4`) is
> proved and even universal; its substantive claim (non-elementary at `k>=P(n)`) has exactly ONE unverified-
> direction of evidence and one conjectural anchor (Catalan at `a=1/2`). I over-claimed it. It may be true,
> false, or have its boundary elsewhere -- open, and honestly so.

Note the k=2 mechanism is **cyclotomic, not CM/elliptic**: it is pure Gauss-digamma / roots of unity, universal
in `a`, with no signature or elliptic input. So the "CM signature" framing (kps-S147/S150) was premature at least
for `k<=2`; genuine elliptic/CM structure, if it enters at all, can only be at the unverified `k>=3` layer.

## 5. Lesson + a Jacobian-thread echo
The meta-lesson (MISTAKES-worthy): **a coarse degeneracy detector reports the wrong boundary.** A PSLQ null is
the "zero-locus" of elementarity; it is too coarse unless the basis is complete for the field. This is precisely
the lesson the repo's Jacobian thread reached independently (HYP-9027 / THM-2455, surfaced this session): the
resolvent *discriminant's zero-locus* "can never separate anything -- only odd-order **valuations** can"; the
genuine degeneracy is caught by a finer (parity) invariant, not the coarse locus. Same shape of mistake, same
fix: replace the coarse null-test by a complete/finer certificate. (The JC degree-monoid HYP-9030,
`{3^k} subset KDeg`, remains a tempting prime-power analogy for a *future* correct form of the law -- but only
once the non-elementary side has a real certificate, which it does not yet.)

Files: `/tmp/{CRISIS,k2thm,fast,retest2}.py`. Corrects/retracts kps-S151 (whole) and kps-S152 (section 4).
