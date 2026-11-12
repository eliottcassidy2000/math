# A twin clock can be arbitrarily near a unit seed without being one

**PROVED, pending independent review.** The literal seed and its defect
polynomial are inherited; the precise valuation comparison below is this
session's arithmetic consequence. This is not a Chan–Pak obstruction outside
this one-port family.

## 1. Source and exact physical length set

The two seven-vertex examples in
[the five-stress atlas, SEEDS §6](../../five_stress_connection_atlas_20261112/SEEDS.md)
have outgoing rows, selected masks and unselected ports

```
(90,28,9,16,36,79,30), O=72, v=2;
(54,28,24,17,32,78,31), O=68, v=0.
```

Replace v by an unselected transitive module of positive length m. The
inherited complete skeleton calculation gives

    δ(m)−1 = −4 [3^m−3·2^m+3] = −4·3!S(m,3).       (1)

Therefore the **entire physical unit-length set is {1,2}**. At m≥3, the
surjection interpretation gives S(m,3)>0, so δ(m)<1. The verifier supplies a
separate vertex-level dynamic program for both cores through m=5; this tests
the inherited polynomial but is not the proof of its all-length identity.

## 2. Opposite prime views of the same failed length

For every m≥3,

    v_2(δ(m)−1) = 3 if m is odd, 4 if m is even;     (2)
    v_3(δ(m)−1) = 1 if m is odd,
                     2+v_3(m) if m is even.         (3)

**Proof of (2).** Write δ−1=−4·3^m+12·2^m−12. For m≥3 the middle term is
zero modulo32. Modulo8, 3^m is3 for odd m and1 for even m. Consequently
δ−1 is8 modulo32 for odd m and16 modulo32 for even m. This proves both
exact valuations.

**Proof of (3).** Write δ−1=12(2^m−1)−4·3^m. For odd m≥3,
2^m−1 is1 modulo3, giving valuation1 versus m. For even m≥4,
LTE applied to4^(m/2)−1 gives v_3(2^m−1)=1+v_3(m). Its first term has
valuation2+v_3(m), strictly below m, the second term's valuation. This proves
(3). The exceptional physical roots m=1,2 have infinite valuation and are
not included in these formulas. ∎

In particular, the physical lengths m=2·3^a (a≥1) have

    v_3(δ(m)−1)=a+2,

but never δ(m)=1. They give arbitrarily accurate *ternary near-seeds*, not
new seeds. By contrast, one dyadic modulus32 detects every failure at once:
for every m≥3, δ(m) is not1 modulo32.

## 3. Good-prime approximation is equally insufficient

For any prime p≥5 and a≥1, take

    m=2+(p−1)p^(a−1).

Euler's theorem gives2^m=4 and3^m=9 modulo p^a, hence δ(m)=1 modulo p^a.
Yet m≥6 and (1) proves δ(m)<1. Thus every such prime admits arbitrary local
precision from nonunit physical lengths. The precise dyadic exclusion and
finite true length set are retained rather than erased by local lifting.

This is an elementary explicit counterpart of Session287's warning that
p-adic closure is not absolute integer realization. Here the failure is
stronger: an entire unbounded tower of p-adic near-solutions can lie outside
a two-point exact admissibility set.

## 4. Connection contract and construction consequence

Source: actual unselected transitive-module inflations of a marked seed.
Target: the signed length response and its prime-adic valuation.
Map: count H and the selected two-path tax, then subtract1.
Preserved: actual tournament existence, selection, literal defect, module
length and congruence accuracy. Lost if only valuations are retained: exact
unit equality and the ordinary magnitude/sign of the defect. Needed sidecar:
the complete exponential response, or an independent modulus such as32.

The constructive lesson is not to abandon long blocks. It is to distinguish
free length axes from coupled directions that cancel the complete response.
The synchronized physical family in `../wildcard/PROOF.md` implements that
second possibility rather than approximating δ=1 prime by prime.

## Reproduction

`python3 verify.py` and `python3 -O verify.py` print identical JSON. All checks
use explicit exceptions, not assertions disabled under optimization. Literal
control universe: both listed cores, each length1–5. Arithmetic controls:
(2)–(3) for m3–500, ternary near-seeds through a7, and good-prime residues for
p in{5,7,11,13}, precisions1–6. The displayed proofs, not these finite ranges,
establish the all-order claims.
