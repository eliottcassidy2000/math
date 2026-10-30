# Prime-power-plus-one symmetry: full moments and signed Frobenius families

2026-10-30. **PROVED scoped deductions / FINITE-EXACT controls / CITED
classical inputs / OPEN general composition.** Same-agent proof review,
with alternative computational formulations; no independent referee or
literature-priority claim.

The main advance is a bound on **every fixed positive moment of the
residual automorphism order**, including nonuniform IID color palettes.
It closes the weighted-tail gap for the specified recursive tournament
actions. A second all-depth law converts factors of p^k+1 into complete
banks of circulant orientations with a prescribed converse Frobenius.
These are different symmetry mechanisms at related integer sizes.

## The exact family inventory, with concurrent credit

Incoming [T1--T2 in GROWTH.md](../prime_power_plus_one_20261029/GROWTH.md),
commit51cd593a3e, independently completed the scalar and rigid-child
recurrences while this session was in progress. Those counts are shared
inherited results here. We retain the full stabilizer-order polynomial
and compare the hierarchical and flat parents at the same size.

For the iterated carousel T_(p,k)=R_p^[k], write
g_k=(p^k-1)/(p-1), so its full group has order p^g_k. The Q-color orbit
polynomial and unmarked one-vertex-child inventory satisfy

    P_0(Q,z)=Q,
    P_(k+1)(Q,z)=[P_k(Q,z)^p-P_k(Q,z^p)]/p+z P_k(Q,z^p),
    E_k(z)=P_k(2,z)-z^(g_k-k).

The final subtraction is the regular-parent clone correction, not an
approximate division by a presumed group order. Complete sector counts:

| child size | parent | all child classes | rigid child classes |
|---|---|---:|---:|
|10|hierarchical C3[C3]|23|2|
|10|flat R9|59|55|
|26|hierarchical R5[R5]|6559|1554|
|26|flat R25|1342183|1342175|
|50|hierarchical R7[R7]|182857159|87460002|

For the recursive prime bases, regular-color counts obey
A_0=2, A_(k+1)=(A_k^p-A_k)/p. They stay at2 only for p=3, and grow
for every prime p>=5. Nevertheless, the probability that a uniformly
chosen binary cut is rigid tends to zero for EVERY fixed p. Counting
rigid classes and sampling labelled cuts are different measures.

## The weighted law and the actual collision statistics

For a fixed faithful base action G on b points, let h=|G| and let d be
the minimum of b minus the number of cycles of a nonidentity element.
Color b^k leaves independently with palette probabilities pi, and set
a=max pi. G3 in [PROOFS.md](PROOFS.md) proves an explicit finite bound

    M_(1,s)^(b^(k-1)) <= E|Stab|^s
       <= M_(1,s)^(b^(k-1)) exp(epsilon_s),
    epsilon_s=O_(G,s)(b^k a^(bd)) as a->0,

for each fixed s>=1. The proof tilts the child-orbit distribution by
the child stabilizer weight before bounding collisions. This keeps
rare, large stabilizers in the estimate.

Put c_j=sum pi_a^j. The resulting exact critical laws are:

| full base action | condition as n=b^k grows | limiting residual group | limiting s-th moment, s>0 |
|---|---|---|---|
|C_p on prime p vertices|n c_p -> lambda|C_p^Z, Z~Poisson(lambda/p)|exp((p^s-1)lambda/p)|
|C_p semidirect C_l, l odd prime|n c_l^((p-1)/l) -> lambda|C_l^Z, Z~Poisson(lambda)|exp((l^s-1)lambda)|

Paley7, Paley11 and the recovered local13 have critical statistics
c3^2, c5^2 and c3^4. The carousel R7 instead uses c7. For uniform
Q colors the corresponding exponents are4,8,8 and6. One heavy color
cannot be hidden by adding many negligible colors: palette cardinality
alone does not govern symmetry breaking.

G3 also gives exact recursive moment evaluators. All parameters of the
base action and the moment order are fixed in the asymptotic theorem.
Multiple ports are individually named. An arithmetic application must
establish the joint leaf law and a weighted error guard; it cannot
replace a deterministic source by fresh independent colors.

## An infinite arithmetic bank, including its reversal phase

Suppose odd ell divides p^k+1. Then ord_ell(p) has exactly the same
2-part as2k. Set ord_ell(p)=2^(a+1)u_1 with u_1 odd, and retain
delta=v_ell(p^ord_ell(p)-1). At EVERY depth e,

    u_e=u_1 ell^max(0,e-delta),
    r=p^(2^a) mod ell^e,  H_e=<r^2>.

On Z/ell^e, H_e-invariant skew connection sets are precisely the
circulant tournaments for which multiplication by r reverses all arcs.
Their complete labelled bank has2^b_e orientations, where

    b_e=sum_(j=1..e) phi(ell^j)/(2u_j)
       =(ell^min(e,delta)-1)/(2u_1)
        +max(0,e-delta) ell^(delta-1)(ell-1)/(2u_1).

Thus the number of new orientation bits per level becomes constant
after the exact lifting delay. Examples:

| p | ell | number of independent bits b_e |
|---|---|---|
|3 or5|7|e|
|3|5|2e|
|7|5|2 at e=1; 10e-8 thereafter|
|5|13|6e|
|7|43|7e|

The delay at7/5 comes from7^2+1=50: modulus25 allows all4096 skew
circulant orientations, and modulus125 has22 independent bits.
These are labelled banks, with a certified translation/multiplier
subgroup; extra graph symmetries remain possible.

The map x -> zeta^x to primitive ell^e-th roots in characteristic p
makes exponent multiplication literally Frobenius. In this carrier
Frob^(2^a) reverses orientation, and its square preserves it. The
F729 implementation checks every arc on its seven root vertices.
This complements the incoming norm-one conic and paired-lattice
theorems; neither equal order nor the word Frobenius identifies the
three carriers.

## Inheritance, evidence and next targets

The Anchor/Niche/Wildcard board is: **recursive child inventories /
weighted collision laws / prime-adic lifting delay / converse phase /
actual joint-port realization**. The closest mechanisms are regular-parent
augmentation, the ternary weighted necklace recurrence, and the old
Fermat--Wieferich split/grow tower. The hostile control remains the
nonregular81/27 pair, whose binary CORE stabilizers differ1620/1728.
The corrected near miss is transferring a bounded-observable limit to
unbounded symmetry weights without uniform integrability. The least-used
sidecar recovered here is the first nonzero prime-adic lifting defect.

The verifier passes132374 explicit gates, including1208 actual child
cuts,162 exact moment queries, complete binary color laws on7/11/13,
50 factorizations, six four-level arithmetic banks, and the literal
F729 phase. Normal and optimized runs must give identical receipts.
See [AUDIT.md](AUDIT.md) for precise finite universes and source pins.

From the repository root:

    python -X utf8 -B 05-knowledge/results/prime_power_plus_one_symmetry_20261030/check.py
    python -X utf8 -O -B 05-knowledge/results/prime_power_plus_one_symmetry_20261030/check.py --output /path/to/separate-receipt.json
    python -X utf8 -B 05-knowledge/results/prime_power_plus_one_symmetry_20261030/seal.py

The strongest next targets are weighted arithmetic transfer under an
actual growing whole-state bound; joint action inventories beyond group
order; and integration of regular-parent augmentation into the current
enumerator with a common-universe benchmark. General native complement
closure and a new full A000568 term remain OPEN. No further binary-color
necklace sweep is needed to establish the stated all-depth recurrences.
