# Loss-aware modular observers for scaled counts

**Status: PROVED** elementary statements below; finite checks are in `verify_modular.py`. No factoring-time or ETH claim is proved here. The inspiration is Rudich–Rousseau, *Deterministic Factoring in Time N^{o(1)}*, hexagon:2610.00600v1, §4, steps 3–4: the engine returns a known multiple `V=2^E D` of a nonnegative count and decides whether it vanishes from small modular evaluations. We use only this elementary interface, not the claimed counting engine.

## 1. The exact information retained by any collection of moduli

Let `s>0`, let `m_1,...,m_t>0`, put `L=lcm(m_i)`, and observe `sD mod m_i` for an unknown integer `D`. The observations of `D` and `D'` agree if and only if

    D ≡ D' (mod Q),       Q = L/gcd(L,s).

**Proof.** Agreement means every `m_i` divides `s(D−D')`, equivalently `L` divides it. Divide by `gcd(L,s)`; the remaining multiplier is coprime to `Q`. This proves both implications. Equivalently,

    Q = lcm_i [m_i/gcd(m_i,s)].

The last identity follows prime by prime: both sides have exponent `max(0,max_i v_p(m_i)−v_p(s))`.

This is an exact quotient, not an estimate. On any interval of `B+1` consecutive integers, all observations determine `D` uniquely iff `Q>B`. For nonnegative `0≤D≤B`, zero is separated from every positive value iff `Q>B`. Necessity is witnessed by `D=0,Q` when `Q≤B`. For signed `−B≤D≤B`, zero detection still needs only `Q>B`, whereas reconstruction of arbitrary `D` requires `Q>2B`.

For constructive recovery, let `r_i` be the observed residue and `g_i=gcd(s,m_i)`. Necessarily `g_i|r_i`, and

    D ≡ (r_i/g_i) (s/g_i)^(-1) (mod m_i/g_i).

These compatible generalized CRT equations recover `D mod Q`. The inverse exists. Reducing modulo `m_i` and attempting to invert `s` there is invalid when they are not coprime.

**Hostile:** `s=100`, moduli `25,50,100,200` have `L=200` but `Q=2`. All even counts, including 0 and 2, are indistinguishable. Four apparently large modular readings provide only one parity bit about the count. With the same four moduli, scales `25,50,100,200,5000` retain effective moduli `8,4,2,1,1`, respectively. In the final two cases every count is invisible to these four observations.

## 2. Remove the scale's prime support, not its enormous exponent

Let `P` be the set of primes dividing `s`. Retain only moduli coprime to `s`, up to `K`, and put

    L_P(K) = lcm{m : 1≤m≤K, gcd(m,s)=1}.

Then the observations of `sD` retain **exactly** `D mod L_P(K)`. Zero detection needs no modular inverse: since `s` is a unit modulo each chosen modulus, `sD=0` there iff `D=0` there. Consequently the zero test remains valid even if the exponent of each prime in `s` is unknown. The set `P`, rather than the size of `s`, controls the observation loss.

For `K=2M`, `M≥1`,

    L_P(2M) ≥ 4^M / [(2M+1) product_{p in P, p≤2M} p^floor(log_p(2M))]
             ≥ 4^M / [(2M+1)(2M)^|P|].

**Proof.** For any prime `p`, Legendre's formula expresses `v_p(binomial(2M,M))` as a sum of `floor(2M/p^j)−2floor(M/p^j)`, each 0 or 1, and zero for `p^j>2M`. Thus the central binomial coefficient divides `lcm(1,...,2M)`. Its size is at least `4^M/(2M+1)`, since it is the largest of the `2M+1` terms summing to `4^M`. Deleting the prime-power factors supported on `P` yields `L_P(2M)`, proving the bound.

For a **fixed** number `r=|P|` of scale primes, a sufficient cutoff for `D≤2^n` is therefore

    2M − log2(2M+1) − r log2(2M) > n.

This permits `2M=n+O((r+1) log(n+2))` as `n` grows with fixed `r`. This is a modulus-size statement, not a claim that the underlying modular evaluations are cheap.

### Explicit dyadic corollary

For every integer `n≥0`, every `E≥0`, and every `0≤D≤2^n`, evaluating `2^E D` modulo the odd positive integers at most `2n+8` determines whether `D=0`. If `E` is known, these readings recover `D` itself.

Take `M=n+4` in the displayed bound with `P={2}`. It suffices that `2^(n+8)>(2n+8)(2n+9)`. This holds at `n=0`; increasing `n` doubles the left side while the ratio of successive right sides is less than 2. Hence the effective modulus is strictly greater than `2^n`.

The cutoff is deliberately simple, not optimal. It depends on the size of `D`, **not on `E`**. In contrast, observing only powers of two up to `2^E` gives identically zero readings for every `D`.

## 3. Burnside's numerator: an actual bridge to symmetry counts

For a finite group `G` acting on a finite set `X`, let

    V = sum_{g in G} |Fix_X(g)| = |G| · |X/G|.

The equality follows by counting pairs `(g,x)` with `gx=x`: each orbit contributes `|G|` because its size times the stabilizer size is `|G|`. Thus Sections 1–2 apply to the orbit count `D=|X/G|` with scale `s=|G|`.

For the order-5,000 normalizer from Session 111, `s=2^3·5^4`. Modular orbit counting can use moduli coprime to 10, and the sufficient bound loses only `(2M)^2`, independent of the exponents 3 and 4. For group orders 25, 50, 100 and 200, the bad prime supports are respectively `{5}`, `{2,5}`, `{2,5}`, `{2,5}`. Equal prime support gives the same coprime-modulus interface even though the groups and their actions can be entirely different.

This does **not** identify the signed and unsigned actions, assert a free group action, replace the fixed-point sums by `|X|`, or identify bank orbits with tournament isomorphism classes. Stabilizers are retained inside Burnside's numerator. There is no claim that an order-200 group is determined by its order.

A second exact strategy, useful when primes dividing `|G|` cannot be avoided, is to compute `V mod |G|m`. The residue is divisible by `|G|`, and dividing it gives `D mod m`. Thus one may pay a modulus-amplification cost instead of discarding prime support. Computing only `V mod m` and dividing its representative by `|G|` is generally wrong.

## 4. Robust modular observations and the exact deletion criterion

The scale-aware observations also form an error-detecting code. Let `q_i=m_i/gcd(m_i,s)` and restrict `D` to an interval of `B+1` consecutive integers. For two distinct admissible values with difference `delta`, their observation vectors differ in precisely those positions `i` for which `q_i` does not divide `delta`.

They can therefore be uniquely reconstructed in the presence of at most `e` arbitrarily corrupted readings **iff**

    lcm{q_i : i not in S} > B    for every S with |S|≤2e.

**Proof.** Unique correction of `e` errors is equivalent to minimum Hamming distance exceeding `2e`: otherwise two radius-`e` balls overlap by splitting the differing coordinates; if distance exceeds `2e`, the triangle inequality prevents overlap. Distance at most `2e` means some nonzero difference of magnitude at most `B` is divisible by all `q_i` outside a set `S` of at most `2e` positions. Such a difference exists exactly when the lcm in the condition is at most `B`; use the lcm itself as the difference. The proof also covers deletion of all coordinates, whose lcm is 1.

For **erasures** in known positions, replace `2e` by the number of erased positions. Repeated or highly overlapping moduli are not automatically useful redundancy. This sharp criterion uses the effective moduli after scaling, not their unreduced sizes.

## 5. Prefix-sensitive divisor observers

Let `N` have binary length `n`. A prefix of length `ell≤n` specifies exactly `2^(n−ell)` nonnegative `n`-bit integers (leading zeros allowed). Let `D_x` count the divisors `d>1` of `N` in that prefix interval. Then

    0≤D_x≤2^(n−ell).

If a correct modular engine returns `2^E D_x mod m`, the odd moduli at most `2(n−ell)+8` suffice to decide whether that prefix contains a divisor, regardless of `E`. Thus the observation budget shrinks as the prefix is refined. At a full-length prefix, `D_x∈{0,1}`, and the **single modulus 3** suffices: `2^E` is never zero modulo 3. No inversion and no knowledge of `E` are needed for this Boolean decision.

If the input has `n` bits and the oracle is queried after selecting `ell=1,...,n` bits, the simple cutoff policy uses at most

    sum_(ell=1)^n [n−ell+4] = n(n−1)/2+4n

odd-modulus evaluations, counting modulus 1 as a harmless redundant evaluation. Removing it saves `n` calls. This is only the cost of the modular observation schedule, not a bound on the counting engine.

The same width argument applies to any binary prefix search on an actual finite set: once a correct nonnegative extension count exists, its upper bound is the number of remaining binary strings. This candidate-count bound fails if each object carries an uncontrolled number of auxiliary witnesses, even when positivity still preserves the zero/nonzero predicate. The paper explicitly states unique auxiliary completion in its multiplication network; that is the relevant guard. Orbit representatives likewise require either a genuine selector or a stabilizer-correct Burnside numerator; raw labelled completions are not automatically orbit counts.

**Hostile to an unbounded zero test:** for any finite family of odd moduli, the positive count equal to their lcm has the same all-zero observations as zero under every dyadic scale. Removing the scale exponent does not remove the required bound on the count.

## Scope and provenance

The effective-modulus and robust-code arguments are elementary CRT/Hamming-distance deductions, not claims of literature priority. The new synthesis is to make the paper's scale factor explicit as a quotient, remove only its prime support, and apply the same interface to stabilizer-sensitive orbit counts. The paper's subexponential counting engine and Lean build are not dependencies. No runtime speedup follows until the cost of evaluating the required residues is separately established.

**Formal-engine boundary:** the incoming [source audit](../../scaled_count_fft_carriers_20261030/SOURCE_AUDIT.md)
reports a generic Lean multiplier `(t+1)^A * 2^(t*R)`. The dyadic corollary
requires its stated purely dyadic scale. Applying the general theorem to
the full multiplier must retain the prime support and effective-modulus
loss of the additional factor; no engine speedup is inferred here.
