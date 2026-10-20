# Independent audit of the neighboring-disjoint all-length exclusion

20 October 2026. Status: **VERIFIED proof application and FINITE-EXACT independent certificates.** No mathematical or scope defect was found. The proposed all-length exclusion is justified by the separately audited full-content theorem, the cited Matveev estimate and frozen absolute-gap certificate, and the independent finite certificates below. The producer's candidate may be promoted to **PROVED** while retaining these dependencies.

This audit writes only into this new audit directory. It does not modify the math repository, the producer, or either inherited frozen proof packet.

## 1. Exact verdict and inherited mechanisms

For primitive lower mechanical sources with

    p>=2, 1<=r<p, gcd(p,r)=1, q=p+r, n=p-r,

the two specified classes of disjoint actual swaps of consecutive source ones are

    plus:  residues {u,u+n}, 0<=u<r-n;
    minus: residues {u,u+n}, n<=u<r.

Both orders are included. Every resulting branch word is nonintegral, for either clock sign and at every length.

The nearest proved mechanism is the neighboring-disjoint full-content estimate from `../../collatz_two_swap_neighbor_disjoint_content_20261020/NEIGHBOR_DISJOINT_CONTENT.md`, independently audited in that packet's `audit/NEIGHBOR_DISJOINT_AUDIT.md`. It bounds the entire gcd by `121^(4 p^(2/3))` or `65^(4 p^(2/3))`, not only the radical. Its proof and all stated scope restrictions were reviewed for this application; its separate content audit is retained as a dependency rather than replaced by the present residual computation.

The canonical hostile is the genuine negative cycle with shape `(p,q)=(7,11)`: it cannot be excluded by a new-cycle minimum estimate and is not a convergent in the residual list. The corrected near miss is transferring a positive-clock-only inequality to negative clocks. The least-used sidecars are the sign of the clock, known-cycle shapes, and full physical source labels at the original cut.

The concept board was actual content, the absolute multiplicative gap, signed minimum magnitude, ordinary continued fractions, physical edit paths, and branch parity. The bridge is

    full content bound -> absolute clock gap -> finite weight range
                       -> signed minimum bound -> finite rational shapes
                       -> complete actual paths -> nonintegral fixed points.

The continued-fraction quotient forgets known cycles, path legality and cut data; the signed cycle controls and physical-label enumeration restore them. The method card used is “Audit and close sections under their next native operation,” with labels, chronology, and branch integrality retained.

## 2. Large-weight transfer is valid for both signs

For every p>=8000=20^3, one has p^(2/3)<=p/20. The exact bounds

    log121<29/6, log3>13/12

therefore give

    log(g/3^p)
      <=4 p^(2/3)log121-p log3
       <(29/30-13/12)p
       =-7p/60
       <-p/9.

The audit independently certifies these rational comparisons. The first 21 nonnegative Taylor terms of exp(29/6) exceed 121; its independent log3 interval exceeds 13/12.

For an integral output, g=|D| and

    U=2^q 3^(-p)-1, L=q log2-p log3,
    0<|U|=g/3^p<exp(-p/9).

Unique factorization gives U!=0. Since U>-1 and |U|<1/2,

    |L|=|log(1+U)|<=|U|/(1-|U|)<2 exp(-p/9).

This is an absolute-value estimate, so no negative clock is omitted.

The frozen endpoint large-range argument has the following reusable arithmetic conclusion:

    There are no coprime p,q with p>=100, p<q<2p, and
    0<|2^q 3^(-p)-1|<exp(-p/9).

Inspection confirms that none of its arithmetic steps uses an endpoint edit after the displayed gap has been supplied. Its cited real Matveev bound is applied to U, not L; it forces p<10^12. Legendre's criterion then selects ordinary convergents, and exact signed intervals exclude all eighteen convergents with 100<=p<10^12. The next denominator exceeds 10^12. Both signs are tested.

Both frozen large-range implementations were rerun unchanged. The producer output matches its stored receipt byte for byte, as does the independent endpoint audit output. These replays are stored here as `inherited.large-range.replay.json` and `inherited.large-range.audit-replay.json`. The primary-source Matveev theorem remains an inherited cited dependency, already inspected in the frozen endpoint audit; this audit does not claim a new logarithmic-form theorem.

Thus every possible integer output in the present family has 2<=p<8000. No endpoint-specific content formula or resultant estimate was transferred.

## 3. Independent complete signed descent through 2^26

The independent source `audit_signed_descent.cpp` is an exact copy of the prior independently audited full-half-step uint64 program, with only its inclusive limit changed from 2^24 to 2^26. The one-line difference is stored in `descent.source-difference.diff`. The inherited source is unchanged.

Unlike the producer's accelerated odd-only uint128 computation, this program:

- visits every integer starting magnitude from 1 through N=67,108,864;
- applies every individual half-step, including the even steps;
- uses checked uint64 arithmetic;
- accepts a non-base start only after it reaches a strictly smaller positive magnitude;
- skips only the separately verified cycle minima 1 for the positive map and 1,5,17 for the negative-magnitude map.

The verified cycles are `(1,2)` on positive integers and

    (1),
    (5,7,10),
    (17,25,37,55,82,41,61,91,136,68,34)

on negative magnitudes. The program checks every cycle edge, distinctness of its entries, and that the chosen base is the minimum. Strong induction on the start proves that all 2N signed nonzero starts in the claimed bounded universe enter one of these four cycles. Trajectories are not assumed to stay below N.

Optimized and O0+UBSan executions are byte-identical:

| Map | Base starts | Certified descents | Half-steps | Longest descent | First start attaining it |
|---|---:|---:|---:|---:|---:|
| positive | 1 | 67,108,863 | 234,384,451 | 376 | 63,728,127 |
| negative magnitude | 3 | 67,108,861 | 234,370,991 | 276 | 30,430,361 |

The largest post-half-step magnitudes are 237,318,849,425,546 and 127,143,512,668,792, at starts 38,595,583 and 58,882,625, respectively. Doubling these values agrees exactly with the producer's largest prehalving numerators. Full-step and accelerated-step totals are different observables and are not falsely compared.

Every multiplication is preceded by an explicit overflow check. Self-tests reject zero magnitude and UINT64_MAX for both signs and verify safe odd-step values near the uint64 boundary. The step cap is fatal rather than a pruning rule. Recompiling with `AUDIT_STEP_LIMIT=1` exits with code 1 and `FAIL: step limit at start 3`; its stdout is empty. The ordinary runs hit neither overflow nor the million-step cap. UBSan reports no undefined behavior.

## 4. Small-weight completeness and known cycles

If an integer cycle is not one of the four verified cycles, its minimum absolute entry M exceeds N. Otherwise its orbit from an entry of magnitude at most N would enter a known cycle and, being periodic, would be that known cycle. Because M is integral, M>=N+1.

Multiplying the odd-step transitions with their intervening powers of two gives

    2^q=product_i (3+1/a_i) for positive cycles,
    2^q=product_i (3-1/a_i) for negative cycles,

where p is the number of odd entries and all a_i>=M. The q here is the total half-step period, not the number of odd entries.

For alpha=log3/log2, the positive logarithm bound log(1+t)<t and the negative bound -log(1-t)<t/(1-t) give

    |q/p-alpha| < 1/((3M-1)log2)
                < 3/[2(3N+2)]
                < 1/(2N).

The strict middle inequality uses log2>2/3. Since p<=7999<8192 and N=8192^2, the right side is below 1/(2p^2). As gcd(p,q)=gcd(p,r)=1, Legendre applies to the unreduced q/p directly.

The independent logarithm calculation uses 48 positive-series terms for log2 and log(3/2), with

    log3=log2+log(3/2).

It computes the ordinary rational continued fractions of both exact interval endpoints independently and checks their digit lists agree through the first denominator at least 8000. It does not import the producer or use its direct log3 series. This certifies the prefix

    [1;1,1,2,2,3,1,5,2,23]

and exactly the residual shapes

    (2,3),(5,8),(12,19),(41,65),(53,84),(306,485),(665,1054).

The next denominator is 15601, from 24727/15601.

The known cycles have shapes `(1,2),(1,1),(2,3),(7,11)`. Their shapes and branch transitions are independently checked as positive controls. Any multiple traversal would multiply both p and q by the same integer, impossible under gcd(p,q)=1. The p=1 shapes are outside the theorem; `(2,3)` is already listed and has no paths; `(7,11)` is explicitly added before enumeration. No known cycle is assigned the false minimum premise M>N.

## 5. Independent physical paths and all 740 record comparisons

`audit_small_neighbor_candidates.py` starts from source positions

    ceil(k*q/p)-1, 1<=k<=p,

rather than the producer's floor-difference source generation. It assigns source labels, takes consecutive physical source ones cyclically, and attempts both common directions. A path is admitted only if both destination bits are zero and the two affected physical bit sets are disjoint. Both move orders are actually executed with labels. Thus no residue interval selects the independent universe.

Residues are attached only after a physical path has been generated, to permit comparison with the producer's record convention. For each output, rational affine-map composition computes the exact fixed point and carry. Its reduced denominator is verified as |D|/g. This differs from both of the producer's integer carry implementations.

The physical audit obtains:

| p | q | Ordered paths | Maximum actual g |
|---:|---:|---:|---:|
| 2 | 3 | 0 | 0 |
| 5 | 8 | 4 | 1 |
| 7 | 11 | 4 | 1 |
| 12 | 19 | 8 | 311 |
| 41 | 65 | 28 | 19 |
| 53 | 84 | 36 | 11 |
| 306 | 485 | 208 | 1 |
| 665 | 1054 | 452 | 37 |
| total | | 740 | |

All 740 exact fixed points are nonintegral. There are 14 ordered paths crossing the original cut. The partial content 311 at `(12,19)` is not confused with the clock magnitude 7153. All four paths at the known negative-17 shape fail integrality.

All fields of every record, including p,q,r, direction, residue, ordered edges, word, B,D and g, match the producer's full trace digest:

    a089dddaa4a1c19eb6846ee3671b5938fe6ab8939f43899fc403dbf3272f4aaf

The complete independent records are saved in `independent_physical_records.jsonl`. Normal and optimized Python runs are byte-identical, with 8,575 explicit checks; no load-bearing check is a Python assert.

## 6. Integer branch validity, clock sign, and scope

For a length-q binary word of weight p, composition of its affine branches is

    x -> (3^p x+B)/2^q.

Its only rational fixed point is B/D. Since B>0 and D!=0, its sign is exactly the sign of D; it is never the zero cycle. The integer condition is precisely g=|D|.

This divisibility test does not forget branch parity. If x=B/D is integral, its fixed-point equation modulo 2 forces x congruent to the first bit: B has that parity and 3^p is odd. The first prescribed branch is therefore an actual integer half-step. Its integer image is the fixed point of the cyclically rotated affine composition, so the same argument supplies the next parity. Induction verifies all q branches and closure. Conversely, an actual integer cycle satisfies the same affine fixed-point equation. Thus exact fixed-point nonintegrality really excludes integer-cycle branch words. Integrality is also rotation-invariant.

The exact all-length theorem remains limited to primitive q=p+r with 0<r<p, disjoint physical edges, consecutive source ones, and one common direction. It does not classify mixed directions, nonneighboring ones, repeated sources, arbitrary mechanical clocks, or all two-edit words. It proves no global signed or positive Collatz convergence theorem.

## 7. Reproduction and audit status

From this audit directory:

    g++ -std=c++17 -O2 -Wall -Wextra -Wconversion -pedantic audit_signed_descent.cpp -o /tmp/neighbor-audit-O2
    /tmp/neighbor-audit-O2 > descent.replay.O2.txt
    g++ -std=c++17 -O0 -g -fsanitize=undefined -fno-sanitize-recover=all -Wall -Wextra -Wconversion -pedantic audit_signed_descent.cpp -o /tmp/neighbor-audit-O0
    /tmp/neighbor-audit-O0 > descent.replay.O0.txt
    cmp descent.O2.txt descent.replay.O2.txt
    cmp descent.O2.txt descent.replay.O0.txt
    python3 audit_small_neighbor_candidates.py > candidates.replay.json
    python3 -O audit_small_neighbor_candidates.py > candidates.replay.optimized.json
    cmp candidates.audit.json candidates.replay.json
    cmp candidates.audit.json candidates.replay.optimized.json

The deliberately failing cap control is:

    g++ -std=c++17 -O2 -DAUDIT_STEP_LIMIT=1 audit_signed_descent.cpp -o /tmp/neighbor-audit-hostile
    /tmp/neighbor-audit-hostile

It must return exit code 1, not success.

`AUDIT_MANIFEST.json` records exact source, receipt, and dependency hashes. The repository worktree was clean on startup and at the audit close; no publication was attempted. The coordination status command could not run because this shared checkout has no `.machine-id`; no registration file was invented or written. This administrative limitation did not affect the proof or computation.

Direct progress is a complete independent audit of the claimed exclusion. The reusable arithmetic conclusion is the family-independent absolute-gap gate. The remaining mathematical frontier is the other disjoint-edit families, not an untested range inside this theorem.
