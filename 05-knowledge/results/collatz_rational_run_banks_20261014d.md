# Rational run schedules enlarge a native Collatz descent bank

**Status:** PROVED elementary finite-certificate theorem and the three scoped banks; FINITE-EXACT declared replay controls. Universal coverage and arbitrary child grounding remain OPEN. This is a specialization and strict enlargement of one declared sufficient bank, not a new general theory of descent cylinders or an external-priority claim.

Companions: [independent checker](../../04-computation/experiments/collatz_rational_run_banks_20261014d.py), [saved output](collatz_rational_run_banks_20261014d.json).

## Inheritance and connection contract

The closest mechanism is the [native successor bank](collatz_successor_holonomy_20261014b.md), pinned initially at `bd2d0d6660b3889179e226172a9572351bd122bc`. Its exact native cells, affine carry and source-preserving receipt semantics are inherited. The [source-paid cylinder compiler](collatz_paid_exit_127_20261012.md) already handles general guarded affine words. The contribution here is a finite rational-schedule certificate and explicit enlargement of the `1^k 2^(2k)` bank.

The canonical hostile is the nontrivial Mersenne family. The concurrent [joint-guard analysis](mersenne_joint_guards_20261014c.md), inspected at `200076ca71ea3abd4e6e6fb6bc8f1c2fd55092d2`, explains why its actual first two runs cannot supply descent. The smaller prescribed schedule here still cannot be substituted for that source's actual word.

Board: exact run index / native residue / affine fixed point / finite schedule seed / natural-density tail / unchanged child. The source-to-certificate map preserves the source, its exact valuations and the strict inequality at the actual endpoint. Its target is a smaller-child dependency. It does not preserve or supply a proof that the child reaches 1; that must remain a separate supplied terminal sidecar.

## 1. Independently proved rational-run certificate

Let `U(n)=(3n+1)/2^v2(3n+1)` on positive odd integers, stopped when it reaches 1. Letters of a word are exact positive valuations, and superscripts on letters mean repetition.

For `k,m>=1`, the word `w(k,m)=1^k 2^m` has affine carrier

`F(x)=(P x+B)/Q`,

`P=3^(k+m)`, `Q=2^(k+2m)`,

`B=Q+P−2^(k+1)3^m`.

Its native positive sources are exactly

`n=2^(k+1)u−1`, where `u>0` and `3^k u=1 (mod 4^m)`.

In particular, `u` is odd, so `v2(n+1)=k+1`; the index k is recoverable from the unchanged source.

### Native-word proof

The first k valuation-one steps are the expansion about −1:

`x=−1+(3/2)^k(n+1)=2*3^k*u−1`.

A valuation-two step acts as `x -> 1+(3/4)(x−1)`. Thus m successive exact twos require and suffice for `x−1` to be divisible by `2*4^m`, which is the stated congruence. All intermediate values are positive odd integers. The endpoint is

`h=1+2*3^m*(3^k*u−1)/4^m >1`.

The quotient is a positive integer, so no earlier state can be 1. The full affine expression follows by composing these two maps.

If `P<Q`, then `F(n)<n` if and only if `n>B/(Q−P)`. Write `z=B/(Q−P)` and let `n0` be the least positive native source. Because every other native source is `n0+2Q t`, `t>=0`, the complete native cell pays if and only if `n0>z`. This exact condition is stronger than merely observing that the slope is contracting.

### Finite all-height certification of a rational schedule

Fix positive integers a,b satisfying

`3^(a+b)<2^(b+2a)`,

and set `m(k)=ceil(a k/b)`. For `r=1,...,b`, put `m_r=ceil(a r/b)` and

`lambda_r=3^(r+m_r)/2^(r+2m_r)`.

The slope condition implies every `lambda_r<1`. Set

`L=max_r lambda_r`, `Z=(1+L)/(1−L)`.

For `k=r+j b`, one has `m(k)=m_r+j a` and

`lambda_k=lambda_r*(3^(a+b)/2^(b+2a))^j <=L`.

Since `B/Q=1+lambda_k−2*(3/4)^m`,

`z_k < (1+lambda_k)/(1−lambda_k) <= Z`.

Choose any K with `2^(K+1)−1>=Z`. Every native source with `k>=K` is at least `2^(k+1)−1`, hence strictly exceeds `z_k`. It therefore suffices to check the finite exact inequalities `n0(k)>z_k` for `1<=k<K`. If any fails, that full bank has not been certified; this theorem does not silently discard the failure or replace the source.

This is an independently proved certificate theorem. No claim of external novelty is made for the elementary affine method.

## 2. The concrete enlarged banks

For `(a,b)=(3,2)`, the maximum multiplier is `243/256`. Thus `Z=499/13<63`, and `K=5` suffices. The four required exact rows are:

| k | m | Least native n0 | Fixed point z |
|---:|---:|---:|---:|
| 1 | 2 | 43 | 23/5 |
| 2 | 3 | 455 | 283/13 |
| 3 | 5 | 8495 | 10865/1631 |
| 4 | 6 | 71199 | 7789/499 |

Each n0 exceeds its z. Therefore all positive native sources of every word `1^k 2^ceil(3k/2)` have `1<h<n`.

The old bank is contained in this bank: `ceil(3k/2)<=2k`, so its exact word is a prefix of the old word. Containment is strict. The source 967 follows

`967 -> 1451 -> 2177 -> 1633 -> 1225 -> 919`,

with exact word `(1,1,2,2,2)`; the next valuation is 1, so it misses the old word `(1,1,2,2,2,2)`. The original source is unchanged throughout the test.

The schedule `(a,b)=(10,7)` is also certified. Here

`L=129140163/134217728`, `Z=263357891/5077565<63`,

and the same four native threshold rows suffice. It contains the `3/2` bank because `ceil(10k/7)<=ceil(3k/2)`. This yields a further small enlargement. It is not needed for the simpler main result.

### Exact density and tail justification

For each k, the native cell is one odd residue modulo `2^(k+2m(k)+1)` and has relative density `2^−(k+2m(k))` among positive odd integers. Distinct k cells are disjoint because their `v2(n+1)` values differ.

All cells with `k>=K` lie in `n=−1 mod 2^(K+1)`, which has relative odd density `2^−K`. This uniform tail bound proves existence of the natural density of the countable union and justifies summing it. For any certified rational schedule its density is exactly

`[sum_(r=1)^b 2^−(r+2 ceil(a r/b))] / [1−2^−(b+2a)]`.

The three densities are:

- old schedule `2k`: `1/31`;
- schedule `ceil(3k/2)`: `3/85`;
- schedule `ceil(10k/7)`: `4737097/134217727`.

The first enlargement adds exactly `8/2635` relative to all odd integers. Its density is `3/170` relative to all positive integers. Combined with the disjoint immediate-descent class `n=1 mod4, n>1`, its relative odd paid-dependency density is `1/2+3/85=91/170`. The root 1 is handled separately, with no nonempty receipt appended at 1. This is only a comparison with the declared old bank and immediate-descent baseline. The repository's larger adaptive banks may already cover these sources. These are not ROOT densities.

### Mersenne hostile and stopping boundary

For a Mersenne source `n=2^E−1`, `E>=2`, the only possible index is `k=E−1` and `u=1`. Membership would require `v2(3^k−1)>=2m(k)`. But this valuation equals 1 for odd k and `2+v2(k)` for even k. Each schedule here has `m(k)>=k+1`, so the required inequality is impossible. Thus every nontrivial Mersenne still misses all three banks.

The rational schedule improves a previously declared sufficient domain. It gives no source-preserving selector for all odd inputs, no bound on the unresolved complement's stopping times, and no grounding of arbitrary emitted children. Pursuing ever closer rational slopes has sharply diminishing density returns and does not repair the Mersenne obstruction; the next substantive task would need another native word shape or a different child interface.

## 3. Exact verification

The independent script imports no other repository module. It checks three rational schedules; eight literal native sources for each k=1..128; every odd residue in the complete native moduli for k=1..4; Mersenne exponents2..128 as finite controls in addition to the all-exponent proof; the strict source967 witness; and the exceptional root1 without padding. All checks stay active under optimized Python.

Run from the repository root:

```
python3 -B 04-computation/experiments/collatz_rational_run_banks_20261014d.py
python3 -B -O 04-computation/experiments/collatz_rational_run_banks_20261014d.py
```

Both runs pass **7,322 checks** and produce byte-identical output with SHA256
`0daf49e1275efda67db90af967e163213ded662736ab958c0056da34c21835aa`.
