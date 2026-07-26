---
id: HYP-1994
status: OPEN
source: oracle-2026-06-01-S521
related:
  - HYP-1965
  - HYP-1966
  - HYP-1987
  - THM-2413-prime-index-affine-drift-and-twin-center-weld
---

# HYP-1994: twin-prime Goldbach reduces to binary Goldbach on the twin-center necklace; last exception is 4208

**Setup.** Twin prime = prime in a twin pair; T={3,5,7,11,13,17,19,29,31,...}.
Even numbers NOT of the form t1+t2 (t1,t2 in T) = OEIS A007534.

**Verified (exhaustive to 2,000,000):** exactly **35** such even numbers, largest
**4208**. Apart from {2,4} they are **11 triples** of consecutive evens
{6m-2,6m,6m+2}.

**Structural claim (VERIFIED exactly on even n in (8,6000], 0 mismatches).**
Twin pairs are centered at multiples of 6 (pair {6k-1,6k+1} has center 6k; only
{3,5}, center 4, is anomalous). THM-2413 identifies the full center set
`C_all` exactly with OEIS A014574 and with the plateau edges of the slope-two
prime-index drift. Put `C_6=C_all minus {4}` and let K = C_6/6 =
{k : 6k±1 both prime} = {1,2,3,5,7,10,12,17,18,...} (the "twin-center necklace").
A sum chosen from two nonexceptional twin pairs lands in
`{c1+c2-2,c1+c2,c1+c2+2}`; since `C_6 subset 6Z`, the covering center sum is a
multiple of 6. Thus each mod-6 triple `{6m-2,6m,6m+2}` is represented by this
packet iff `6m in C_6+C_6` iff `m in K+K`. On the stated finite range, the
separate exhaustive check confirms that the exceptional pair `{3,5}` creates
no additional representability status. Hence, within that verified scope:

  **even n>8 is a sum of two twin primes  ⇔  round(n/6) ∈ K+K.**

**Internal-necklace ancestry (THM-2422).** The tempting local law
`k_i-k_(i-1) in K` is false: its first post-startup failure is
`58=52+6` with `6 notin K`, i.e. A014574 term `348=312+36`. The stronger
distinct-parent fibre nevertheless satisfies

```text
every k in K, 3<=k<=16,666,598, has k=a+b with a<b in K
```

by a FINITE-EXACT, independently reproduced census through center
`100,000,000` (`440,309` targets). The all-`k` extension is OPEN and is
strictly narrower than the present whole-integer `K+K` conjecture.

There is a real local explanation for the difference. Apart from the
prime-equals-five startup `k=1`, one has

```text
K mod 5 subset {0,2,3}.
```

The ordered parent-channel counts from `{0,2,3}^2` into residues
`0,1,2,3,4` are `3,1,2,2,1`. Thus an internal `K` target avoids the two
one-channel residue classes, whereas whole-integer `K+K` must cover them;
ten of the eleven recorded holes below are `1 mod 5`. This is an exact local
thinning mechanism, not a proof that no further holes exist.

The 11 exception-triples ⇔ the 11 integers not in K+K:
m = 16, 67, 86, 131, 151, 186, 191, 211, 226, 541, **701** (×6 → 96,402,516,
786,906,1116,1146,1266,1356,3246,4206). The triples come in threes purely as the
±2 thickening of one center-sum 6m; the genuine content is the 11-hole comb of
K+K.

**Conjecture (OPEN):** m=701 is the LAST hole of K+K — i.e. every m>701 is a sum
of two twin-center necklace elements, so twin-Goldbach has no exception above
4208. K has density ~ 2·C2·n/log^2 n (sparser than primes); the holes are a
finite startup deficit before the representation count r_K(m) stays positive.

**To do:** (a) verify r_K(m)>0 for m well past 701 (push the bound to 10^7+);
(b) lower-bound r_K(m) via Hardy–Littlewood twin-prime density to make the
"no holes past 701" claim conditional-provable; (c) explain the specific 11
m-values (do they sit in local thinnings of K?).

**Evidence:** `04-computation/twin_goldbach_exceptions_s521.py`,
`04-computation/twin_goldbach_necklace_reduction_s521.py` (+ `.out` in
`05-knowledge/results/`). **Reflection:**
`07-reflections/twin-goldbach-necklace-triples-s521.md`.
**Lens link:** the center coordinate 6k is exactly the pair-lens coordinate of
HYP-1965/1966 (`pair-first-twin-prime-lens-s502`); twin-Goldbach is Goldbach
"on pairs."
