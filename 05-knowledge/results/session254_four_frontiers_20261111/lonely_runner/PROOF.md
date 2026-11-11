# Session254: phase-uniform critical packs and affine seven-runner attachment

**Status: PROVED IN NOTE (elementary); the arbitrary six-speed core corollary is PROVED RELATIVE TO the cited lower-dimensional Lonely Runner theorem. FINITE-EXACT controls supplied. LRC(14) remains OPEN. No general literature-priority claim.**

## 1. Native question and inheritance

Can a critical seven-runner packet retain a **uniform positive safe set under arbitrary phase shifts**, and can that robustness be localized inside an actual safe interval of six other runners? This is a native LRC question, rather than a safe-word/Mahler encoding.

Closest proved mechanisms are THM-4112, `01-canon/theorems/THM-4112-antipodal-component-ancestry-chain-and-scale-separated-lrc-families.md`, and THM-4448, `01-canon/theorems/THM-4448-lrc14-general-shore-attachment-and-decoder-pair-cones.md`. The former keeps component ancestry under ratio gates; the latter attaches scaled shores inside protected arcs. THM-4066, `01-canon/theorems/THM-4066-lrc14-diagonal-intercept-pullback-and-exact-affine-ray-closure.md`, preserves a moving intercept in two selected exceptional shapes.

The least-used relevant sidecar is the **phase of an isolated divisor frequency**, combined with the fact that danger occupancy is an integer. The historical THM-594, `01-canon/theorems/THM-594-pair-overlap-law-mirsky-newman-floor.md`, points toward a divisor-minimal Fourier coefficient and an unshifted quadratic floor. Its file is CLAIMED in startup routing and contains a corrected/refuted pair-overlap branch. We use it only as idea provenance: every load-bearing formula below is proved from scratch, with no pair-overlap formula or status promotion.

The canonical hostile for this lane is seven equal slopes with phases j/7. Their open danger arcs at radius 1/14 tile the circle except for seven boundary points; every radius greater than 1/14 covers it. This is why arbitrary phases and slope multiplicity cannot be discarded. The corrected near miss is replacing a whole-period mass statement by a statement on a prescribed physical core component; Theorem 3 explicitly retains the component and bounds the within-block drift.

Concept board: (i) integer occupancy and L1 defect; (ii) divisor-isolated Fourier frequency; (iii) phase-retaining affine blow-up; (iv) protected physical interval; (v) repeated-slope tiling hostile; (vi) positive target fibres and cost in Chan–Pak. The resulting theorem connects (i)–(v); the exact but limited finite-count transport to (vi) is stated in Section 7.

## 2. An integer-occupancy certificate

Let C be any nonnegative integer-valued integrable function on a probability space. Set mu=E C and p=Pr(C=0). Then

    E |C-1| = 2p + mu - 1.                                      (1)

Indeed C-1 is -1 on C=0 and nonnegative elsewhere. Thus its negative part has mean p, its positive part has mean p+mu-1, and their sum gives (1).

If a complex observable chi has E chi=0 and |chi|<=1, then

    p >= (|E(C chi)| - mu + 1)/2.                               (2)

This follows by applying the triangle inequality to E((C-1)chi). The right side may be negative, in which case the trivial p>=0 is stronger. No independence, equidistribution, pairwise-overlap estimate, or probability model beyond the specified measure is assumed.

## 3. A shifted-packet theorem

Let a_1,...,a_m be positive integers, with repetitions permitted, and let phi_i be arbitrary real phases. For 0<d<1/2 define

    C(u)=sum_i 1_{||a_i u+phi_i||<d},  u in R/Z,
    S_d={u: ||a_i u+phi_i||>=d for every i}.

For every nonzero integer k, direct substitution on the a_i equal subintervals gives

    c_k = integral C(u) exp(-2 pi i k u) du
        = sum_{i:a_i divides k}
          [sin(2 pi d k/a_i)/(pi k/a_i)]
          exp(2 pi i (k/a_i) phi_i).                           (3)

Terms with a_i not dividing k vanish by the geometric sum of a_i roots of unity. Each danger comb has measure 2d, so mu=2md. By (2),

    measure(S_d) >= (|c_k| - (2md-1))/2.                       (4)

**Theorem 1 (isolated-divisor phase-uniform floor).** Suppose one index * has the property

    for every j != *, a_j does not divide a_* .                (5)

Then, uniformly in all phases,

    measure(S_d) >= [sin(2 pi d)/pi - (2md-1)]/2.              (6)

Proof: take k=a_* in (3). There is exactly one contribution, of modulus sin(2 pi d)/pi. Apply (4). QED.

A unique smallest slope suffices, but is not necessary. The slopes (2,2,2,2,2,2,3) also satisfy (5), with a_*=3. Seven pairwise distinct slopes always satisfy (5), using their minimum. Repeated slopes other than the isolated one cause no difficulty.

**Corollary 2 (robust critical seven-pack).** For m=7 under (5):

    measure(S_(1/14)) >= sin(pi/7)/(2 pi),                     (7)
    measure(S_(2/25)) > 311/18750.                            (8)

For (8), use sin x >= x-x^3/6 at x=4pi/25 and pi^2<10:

    sin(4pi/25)/pi > 4/25 - 64/9375.

Here 14*(2/25)-1=3/25, so (6) is strictly greater than

    (4/25 - 64/9375 - 3/25)/2 = 311/18750.

The weak threshold 1/14 has therefore been replaced by a strictly larger safe radius 2/25, uniformly over phases and over the sizes, ratios, and number of repeats of all nonisolated slopes. The margin is exactly

    2/25 - 1/14 = 3/350.                                    (9)

This is a deliberately convenient rational constant, not a sharp claim. The one-frequency method works whenever sin(2pi d)/pi > 14d-1.

## 4. Localizing inside a physical safe interval

**Theorem 3 (affine critical-pack attachment).** Let I be a closed real interval of length L, and let a slow core be 1/14-safe throughout I. Let

    v_i=h a_i+b_i,  i=1,...,7,
    h>0, a_i positive integers satisfying (5),
    B=max_i |b_i|,
    h L>=1,  3h>=175 B.                                    (10)

The b_i may be real for this interval theorem. Assume the v_i are positive. Then there is a time in I safe for both the core and all seven v_i. More quantitatively, in every closed subinterval of I of length 1/h the simultaneous safe set has measure strictly greater than

    311/(18750 h).                                          (11)

Proof: write the subinterval as t=t_0+u/h, -1/2<=u<=1/2. Retain the exact phases

    phi_i=(h a_i+b_i)t_0.

Then

    v_i t = phi_i+a_i u+(b_i/h)u.                            (12)

By Corollary 2, a set of u of measure >311/18750 satisfies ||phi_i+a_i u||>=2/25 for all i. Its last perturbation in (12) has absolute value at most B/(2h)<=3/350. Since distance to Z is one-Lipschitz, all those u give ||v_i t||>=1/14. The core is safe because t remains in I. Scaling du to dt gives (11). Equality in the drift gate is allowed because the target inequalities are closed. QED.

Disjoint blocks show that the total safe measure in I is greater than

    floor(h L) * 311/(18750 h)

when floor(h L)>=1. This addresses local positive mass on this explicit affine-scale family; it is not an arbitrary-core recurrence theorem.

**Corollary 4 (an explicit thirteen-speed family).** Let C be any set of six distinct positive integer speeds, M=max C. Let a_i be seven positive integer slopes satisfying (5), let b_i be integers, B=max |b_i|, and let h be a positive integer satisfying

    h>=7M,  3h>=175B.                                       (13)

If C together with the seven h a_i+b_i consists of thirteen distinct positive speeds, it satisfies LRC(14).

Proof: cited lower-dimensional LRC supplies t_0 with ||c t_0||>=1/7 for every c in C. The closed interval centered at t_0 with radius 1/(14M) stays 1/14-safe for C, by the Lipschitz inequality. Its length is 1/(7M). Apply Theorem 3. The core theorem is the sole non-elementary imported input. Alternatively, supplying any explicit safe interval I makes Theorem 3 wholly elementary and removes that dependency. QED.

For C={1,...,6}, t_0=1/7 is explicit. Examples in the exact output include

    {1,2,3,4,5,6,172,348,524,700,876,1052,1228}

(h=175, a_i=i, b_i=i-4), and the much more clustered carrier

    C union {175,347,348,349,351,352,353}

(h=175, slopes (1,2,2,2,2,2,2), offsets (0,-3,-2,-1,1,2,3)). In each case the physical safe mass in [1/7-1/(2h),1/7+1/(2h)] is printed, rather than just a Fourier proxy.

## 5. What has changed relative to the inherited mechanism

THM-4448's protected-arc mechanism is retained, not replaced. Its common-dilation shore has one fixed danger-component profile. Here offsets alter the profile and phases vary with the block center. The new ingredient is a **phase-uniform, strict-radius packet supplier**, which absorbs that change quantitatively through (12). No selected component is silently changed, and no adjacent-slope ratio or upper bound on max a_i is required.

THM-4112's ancestry ratios remain useful for other packets. The present hypothesis instead allows six equal carrier slopes beside one isolated slope, with distinct physical speeds supplied by the offsets. It therefore targets a different operational cell: a bounded-width moving cluster plus an isolated carrier.

The old quadratic unshifted Fourier floor is not being reverified. The new proof uses L1 and integer occupancy, permits arbitrary phases, retains a positive safety margin, and proves an actual localized physical conclusion.

## 6. Failure boundaries and frontier

1. **Repeated-slope hostile.** With all a_i=1 and phi_i=i/7, i=0,...,6, the 1/14-safe set is exactly seven points, hence has zero measure. For any d>1/14 the safe set is empty. In (3) the first-harmonic phasors cancel. Thus erasing slope multiplicity or phase destroys the theorem's supplier.
2. **Not an eight-pack claim.** At m=8,d=1/14, the coefficient sin(pi/7)/pi is smaller than 2md-1=1/7, so (6) gives no positive mass. This is failure of this certificate, not an eight-runner counterexample. Additional modes, their joint phase, or another integer-occupancy test are needed.
3. **Not arbitrary affine slopes.** Condition (5) is sufficient, not necessary. Without it, (3) remains exact, but the phasor sum can vanish. Investigate phase paths phi_i(t_0), rather than substituting phasor counts.
4. **Not arbitrary LRC(14).** A general speed set has not been shown to admit the scale, bounded offset, isolated carrier, and actual safe core interval in (10)/(13). This representation/entry problem remains open.
5. **Not primitive-line recurrence.** We construct witnesses in named finite blocks; no Haar-mass-to-recurrence inference is used.
6. **No universal denominator cutoff.** All witnesses are physical interval points. No bounded-denominator claim is made.

Next native question: for repeated primitive carrier slopes with cancelling first phasor, can the actual offset path phi_i(t_0) force a nonzero higher/divisor coefficient somewhere on a protected core interval? The exact formula (3) is the retained sidecar. The equispaced equal-slope hostile determines what cannot be forced without offset or interval hypotheses.

## 7. Honest Chan–Pak transport

There is one exact transport, and an important limitation. On a finite address set Omega, replace the integral by the uniform average. Let C(x) count the number of **actually legal** constructions landing at x, and let a character/test vector chi have zero sum and |chi(x)|<=1. Then (1)–(2) remain exact:

    # {x:C(x)=0} / #Omega
      >= (|sum_x C(x)chi(x)| - sum_x C(x) + #Omega)/(2 #Omega).

Source: integer incidence of forbidden combs. Target: integer incidence of legal output constructions. Map: replace Lebesgue integration and Fourier characters by normalized finite summation and a specified bounded zero-mean test vector. Preserved predicate: a signed occupancy imbalance forces an empty fibre when total load is near one. Lost data: legality of the constructions, integer output values, order cost, and the relation between a quotient address and its full target fibre. Required sidecar: an actual legal-construction incidence table retaining those data. Cheapest decisive test: evaluate both sides on such a table, not on a residue histogram alone.

This is an obstruction/certificate transport, **not** a positive-fibre construction for Chan–Pak and not a claim that its presently missing sidecars have been supplied. The useful design lesson is precise: an isolated observable must be retained **before** phase or address summation; a quotient with the same counts need not retain it. Native LRC progress above does not depend on a forced Chan–Pak identification.

## 8. Reproduction and source pins

From repository root:

    python3 05-knowledge/results/session254_four_frontiers_20261111/lonely_runner/reproduce.py

`output.json` is the matching output. The code uses exact Fraction endpoints and interval union, tests 156 phase rows across six carrier multisets, four physical thirteen-speed attachments, and the repeated-slope hostile. A distinct wall-midpoint implementation rechecks the least-safe tested row of each carrier family. The normal and optimized Python outputs agree byte for byte. These are finite controls of the new statements, not empirical substitutes for their proofs and not a census of all phases or speed sets.

`SOURCE_PINS.json` records inherited source paths and scope. `SHA256SUMS` pins this packet's proof, code, output, and source list. No historical theorem has been promoted, and no canon theorem ID is reserved by this packet.
