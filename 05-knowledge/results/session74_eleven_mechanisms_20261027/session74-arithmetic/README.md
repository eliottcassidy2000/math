# Session74: the bracket-prime eleven, a weighted residue bridge, and its prime-power lift

**Status:** PROVED (elementary identities and proofs below); FINITE-EXACT (the explicitly bounded checks and the 341 stratum totals). No literature-priority claim and no canon promotion. Source pin: `4959b4933cba7e6561b9f19b1d68278ecf6cb2f5`. This packet does not prove a prime-distribution theorem, Collatz, or a Hamiltonian-count theorem.

## Result first

The old `{2,3,11}` is genuinely special for one **scale predicate**: these are exactly the primes whose double remains in the same odd-square bracket. It is not a multiquadratic or multiplicative-order classification.

There is, however, a new exact way to leverage its parent construction. Relax doubling along the natural integer-reciprocal family

```text
r_q = 1 + 1/(q-1) = q/(q-1),     q >= 2,
S_q = {n >= 2 : n and r_q*n lie in the same odd-square bracket}.
```

Doubling is the first member, `q=2`. The complete integer cardinality is a quadratic polynomial **minus a weighted local square-residue correction modulo q**. For an odd prime `p`, put

```text
S(p) = sum_{a=1}^{p-1} a*chi_p(a),
```

where `a` is the least positive representative and `chi_p` is the quadratic character. Then

```text
|S_p| = (16p^2 - 27p + 11)/3 + 2*S(p)/p.                 (1)
```

Here `2*S(p)/p` is an ordinary rational expression, **not** a Legendre symbol. At `p=11`, `QR_11={1,3,4,5,9}`, so `S(11)=22-33=-11` and `|S_11|=548`. The same weighted defect has an exact prime-power lift:

```text
q=p^e,  Z=p^floor(e/2),  h=ceil(e/2),
|S_q| = (16q^2-27q+8)/3 + Z
        + (2*S(p)/p)*(p^h-1)/(p-1).                    (2)
```

For `p=11`, this gives `548, 77008, 9436336, 1143115696` for exponents `1,2,3,4`.

At `q=341`, the unit square-image is literally `QR_11 x QR_31`, the same subgroup `<4,81>` recovered in Session24. The bracket count is `617092`. It needs the source root multiplicities, least-residue weights, and both zero-coordinate axes. Keeping the unit correction alone gives `617128`, an exact error of `36`.

These are real operation-preserving links. They do **not** say that the original doubling exception `11` occurs because of a Paley group, or that any of the displayed integer counts is a prime count.

## 1. Inheritance and concept board

Pinned source snapshots and full source-path/hash mapping are in `sources/` and `SOURCE-MANIFEST.json`.

- **Closest proved mechanism:** `procgen_numerology_20260923_odd_square_brackets.md` §2 and `procgen_brackets_20260924_pairings_transitions.md` §1: exact odd-square brackets, multiplier cutoff, and threshold function.
- **Corrected near miss:** the real interval cutoff is not the integer cutoff; right-closed upper endpoints are load-bearing. `decoder_prime_square_20260925.md` §2 independently restates the all-height classification.
- **Canonical hostile:** `11` and `44` have the same squareclass, but only `11` survives doubling in its own bracket. Replacing the basis triple by `{2,3,5}` leaves its abstract Fano geometry unchanged but makes the third prime fail the source predicate.
- **Least-used sidecar:** the least residue of the terminal odd square modulo the multiplier denominator; it retains the floor/ceiling error that an unweighted squareclass support forgets.
- **Related recovered carrier:** `arithmetic_clock_tournaments_20261021/session24-arithmetic/README.md` §§2–3: complete conductor strata at 341, and the warning that arbitrary independent prime triples have the same Fano/Hamming labeling.

| Object | Operation / question | Preserved coordinate | Hostile or boundary |
|---|---|---|---|
| Prime bracket escape | multiply within a square interval | actual magnitude and endpoints | `11` vs `44` |
| Multiplier family | `r_q=1+1/(q-1)` | exact floor counts | no prime count from integer count |
| Square residues | push terminal odd roots modulo q | least-residue weight and root multiplicity | support alone loses weights |
| Prime-power lift | separate unit roots from p-divisible roots | zero and valuation strata | omit `Z=p^floor(e/2)` |
| Joint 341 clocks | CRT and square map | all four conductor strata | units-only error `36` |

## 2. Exact recovery of the original eleven

Write

```text
B_m=((2m-1)^2,(2m+1)^2],   m>=1,
```

and keep `{1}` separate. A positive integer `n in B_m` has a nontrivial integer multiple in `B_m` iff its double is there. For `m>=3`,

```text
2*((2m-1)^2+1) - (2m+1)^2 = 4m^2-12m+3 > 0.
```

Thus only the first two brackets contribute:

```text
all integers: {2,3,4,10,11,12};
prime integers: {2,3,11}.
```

The complete prime-multiplier pairs `(k,p)`, `k>=2`, are

```text
(2,2),(3,2),(4,2),(2,3),(3,3),(2,11).
```

Under a strict upper endpoint only `(3,3)` disappears; the prime set is unchanged. The source threshold `rho(n)=(2m(n)+1)^2/n` gives precisely the prime plateau `{2,3,11}` for `r in (25/13,25/11]`. At the lower endpoint `13` enters; immediately above the upper endpoint `11` leaves.

This is an inherited complete theorem, not an extrapolation from a prime search. Its mechanism is the mismatch between the linear square-gap width and multiplicative growth. The prime filters are applied only after that scale cutoff.

## 3. Universal denominator formula, including even and composite moduli

For every integer `q>=2`, define `S_q` as above. Its contribution from bracket `m` is the integer interval

```text
(2m-1)^2 + 1 <= n <= floor((q-1)*(2m+1)^2/q).           (3)
```

### Theorem A: exact bracket range and weighted floor formula

Exactly the brackets `1 <= m <= 2q-2` contribute, and each is nonempty. Put

```text
T_q = sum_{m=0}^{q-1} [-(2m+1)^2]_q,
```

where `[z]_q` is the unique residue in `0,...,q-1`. Then

```text
|S_q| = (16q^2-24q+8)/3 - 2*T_q/q.                    (4)
```

**Proof of the range.** Nonemptiness of (3) is equivalent to

```text
D_q(m) = -4m^2+(8q-4)m-q-1 >= 0.
```

At both `m=1` and `m=2q-2`, this is `7q-9>0`; concavity gives all intervening integers. At `m=2q-1`, it is `-q-1<0`, and it decreases thereafter. Thus the range is complete.

**Proof of the count.** The number of source integers in (3) is

```text
W_q(m)=8m-ceil((2m+1)^2/q).
```

At `m=0` and `m=2q-1`, the formal value is `-1`. Therefore add these two formal terms and then add `2` back. The polynomial part sums to

```text
sum_{m=0}^{2q-1} (8m-(2m+1)^2/q)
  = (16q^2-24q+2)/3.
```

The ceiling correction is `[ -(2m+1)^2 ]_q/q`, periodic in `m` with period `q`. Its sum over the two periods is `2*T_q/q`. Adding the endpoint correction gives (4). This proof works for **every** `q>=2`; primality and invertibility of `2` are not required. QED.

For odd `q`, multiplication by `2` is invertible. Define

```text
R_q = sum_{x=0}^{q-1} [x^2]_q,
Z_q = #{x mod q : x^2=0 mod q}.
```

Then `T_q=q*(q-Z_q)-R_q`. This is the form in which root multiplicities and zero strata are visible.

## 4. Prime-specific weighted character defect

For an odd prime `p`, the nonzero squares each have exactly two roots, and only zero has a zero square. Hence

```text
R_p = 2*sum_{a in QR_p} a
    = p*(p-1)/2 + S(p),
T_p = p*(p-1)-R_p.
```

Substituting in (4) proves (1). If `p=1 mod4`, pairing `a` with `p-a` gives `S(p)=0`. For `p=3 mod4`, it need not vanish and is not uniformly `-p`:

| p | weighted character sum S(p) | correction 2*S(p)/p | count |
|---:|---:|---:|---:|
| 3 | -1 | -2/3 | 24 |
| 5 | 0 | 0 | 92 |
| 7 | -7 | -2 | 200 |
| 11 | -11 | -2 | 548 |
| 19 | -19 | -2 | 1756 |
| 23 | -69 | -6 | 2612 |
| 31 | -93 | -6 | 4844 |

No class-number formula is used. In particular, `11` is not uniquely selected by a correction of `-2`. The actual prime-specific ingredient is the **weighted** residue split, with the local modulus named.

At `q=11`, a separate exhaustive primality filter on the exact source intervals finds `100` primes, the largest `1523`. That is a FINITE-EXACT side result only; (1) counts `548` integers and does not by itself count or predict those primes.

## 5. Prime powers: an exact source-stratified lift

Fix an odd prime `p` and an exponent `e>=1`. Let `R_e=R_{p^e}` and let `U_e` be the same residue sum restricted to roots not divisible by `p`. Set `R_0=0` for the one-element ring modulo `1`.

### Theorem B: weighted square-residue lift

Writing `q=p^e`, `Z=p^floor(e/2)`, `h=ceil(e/2)`, one has

```text
U_e = (p-1)*p^(2e-1)/2 + p^(e-1)*S(p),                (5)
R_e = U_e + p^3*R_(e-2)                 (e>=2),       (6)
R_e/q = (q-Z)/2 + (S(p)/p)*(p^h-1)/(p-1).            (7)
```

Consequently formula (2) holds.

**Proof.** For a unit `x modulo p^(e-1)`, its `p` lifts are `x+k*p^(e-1)`. Their square residues modulo `p^e` share the same low part, while their high digits run through all `0,...,p-1`: the coefficient `2x` is nonzero modulo `p`. Thus, for `e>=2`,

```text
U_e=p*U_(e-1)+(p-1)^2*p^(2e-2)/2.
```

Starting with `U_1=R_p=p*(p-1)/2+S(p)` proves (5).

A nonunit root is `x=p*y`, with `y modulo p^(e-1)`. Its square residue is `p^2*[y^2]_(p^(e-2))`, and every lower root repeats `p` times. This proves (6). The zero roots are precisely those divisible by `p^ceil(e/2)`, so their number is `Z`.

After dividing (6) by `p^e`, the recurrence is

```text
R_e/p^e=(p-1)*p^(e-1)/2+S(p)/p+p*R_(e-2)/p^(e-2).
```

Summing its `h` terms gives (7). Finally use `T_q=q*(q-Z)-R_e` in (4), which yields (2). QED.

This lift preserves the actual weights and the p-divisible source strata. It is not an identification of the mod-11 field with the mod-121 ring or the split 121-element ring; those remain different objects, as the inherited Session24 warning requires.

## 6. CRT, 341, and the exact price of deleting the axes

For coprime odd moduli `u,v`, let `nu_u(a)` count roots of `x^2=a modulo u`, and similarly for `v`. Put `kappa_q(0)=0` and `kappa_q(a)=q-a` for `1<=a<q`. CRT gives the general weighted transport law

```text
T_(uv) = sum_{a mod u,b mod v}
         nu_u(a)*nu_v(b)*kappa_(uv)(CRT(a,b)).           (8)
```

The `CRT(a,b)` on the right is the least residue in `0,...,uv-1`; replacing it by the pair without this weight loses information. The proof is simply the CRT bijection on source roots, followed by squaring.

At `u=11,v=31`, the nonzero unit image is `QR_11 x QR_31`. This is the inherited joint-clock subgroup `K=<4,81>` in `Z/341`:

- `4` has order `5` modulo both primes;
- `81/4` has coordinates `(1,28)`;
- `28^3=4 mod31`, `28^5=5 mod31`, and `28^15=1 mod31`, so it has order `15`;
- both generators are squares, and these orders generate all `QR_11 x QR_31`.

Thus the same subgroup is obtained by an intrinsic operation on the bracket source, **squaring terminal odd roots**, not by matching a number or an abstract group order.

The complete weighted CRT audit is:

| Source roots modulo341 | root count | square-image size | roots per image | sum of square residues | ceiling tax |
|---|---:|---:|---:|---:|---:|
| zero | 1 | 1 | 1 | 0 | 0 |
| nonzero multiples31 (additive order11) | 10 | 5 | 2 | 1364 | 2046 |
| nonzero multiples11 (additive order31) | 30 | 15 | 2 | 6138 | 4092 |
| units | 300 | 75 | 4 | 49104 | 53196 |
| total | 341 | 96 | varying | 56606 | 59334 |

The integer count from (4) is therefore `617092`. Deleting both axes changes the tax by `6138=18*341` and changes the count by `2*6138/341=36`. The false units-only count is `617128`.

This is precisely the kind of source-stratification issue seen in Session24's invariant-skew-section test, but it is a **different preserved predicate**: a weighted bracket cardinality, rather than the existence of an invariant tournament. No Hamiltonian count or Collatz stopping-time consequence is imported.

## 7. Connection contract and stopping boundary

**Source:** labelled integer intervals (3), retaining the primality tag only when it is separately tested.

**Target:** square-residue distributions modulo the chosen denominator `q`; for odd prime `q`, their weighted quadratic-character defect; for `q=341`, the full CRT-stratified square image.

**Map:** terminal odd root `2m+1` to its square modulo `q`, pushed forward with source root multiplicity and least-residue ceiling weight.

**Preserved:** the exact count of all integers surviving the multiplier `q/(q-1)`, and the prime-power/CRT transport formulas for that count.

**Lost if unadorned:** bracket height, the prime tag, individual source integers, least-residue magnitude, and root multiplicity. The residue support alone or an abstract Paley/Fano isomorphism does not preserve the count.

**Restoration sidecars:** exact finite bracket range, original endpoint formula, local modulus, least-residue weight, root multiplicity, and zero/nonunit strata. A prime count additionally needs the actual source integers and a primality test.

**Cheapest decisive tests:** `q=2` recovers `{2,3,4,10,11,12}`; `p=3` tests the rational correction; `p=23` defeats a universal `-2`; `11` versus `44` defeats squareclass-only escape; `q=341` defeats units-only transport by exactly `36`.

The new map creates a rigorous way to use the older 11/341 arithmetic. It does not retroactively explain the original doubling set by a quadratic-character theorem. No claim of infinitely many primes, a prime in every bracket, global Collatz convergence, an LRC certificate, or a Hamiltonian gap closure is made.

## 8. Reproduction and audit scope

From this packet directory:

```bash
python3 check_brackets.py > receipt.replayed.json
python3 -O check_brackets.py > receipt.replayed.optimized.json
cmp receipt.replayed.json receipt.normal.json
cmp receipt.replayed.optimized.json receipt.optimized.json
sha256sum -c SHA256SUMS
```

The standard-library checker has no repository or absolute-path dependency and uses exceptions, not removable assertions. Normal and optimized outputs are byte-identical. It performs **2023 explicit checks**:

- every denominator `2<=q<=512`: exact bracket intervals versus the weighted residue formula, including all cutoff endpoints;
- every `2<=q<=64`: independent literal integer membership using integer square roots versus the interval lists;
- **67 prime-power cases**, prime bases `3..101`, exponents at most `6`, modulus at most `20000`: unit lifting, nonunit recurrence, zero-root count, direct residue totals, interval totals, and the closed formula;
- exact recovery of the doubling prime set, full integer set, integer-multiplier pairs, and both plateau boundaries;
- all `341` source roots, all CRT strata and square-map fiber multiplicities, and direct comparison with `<4,81>`;
- the finite `q=11` prime list and the stated hostile controls.

The algebra proves the all-height/all-modulus statements. The checks audit their conventions and consequences; they are not substituted for those proofs. An independently authored root checker is maintained outside this packet and may be attached by the publishing coordinator.
