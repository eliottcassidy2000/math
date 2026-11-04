# Independent odd TT2-quotient trace audit

**PASS.** The all-order modulo-eight proof, the physical query battery and the two-level reversal-cube corollary are valid at their stated scope. The final producer sources and retained data are pinned by `retained.json`; no even-quotient separation or general ambient transfer is inferred by this audit of this packet.

## Primary input and analytic checks

The primary PDF of Grinberg and Stanley, [The Redei–Berge symmetric function of a directed graph](https://arxiv.org/pdf/2307.05569), arXiv v1, was opened and checked directly. Its Theorem1.39 on PDF page14 has the stated odd-cycle permutation expansion. Setting one color to1 and the others to0 counts paths in the complementary tournament, whose scalar Hamilton count equals the original by path reversal. Thus the displayed odd-cycle collection formula has the correct factor2 for each nontrivial directed odd cycle, counted up to rotation. This classical identity is an input, not a new theorem of the packet.

For a single odd cycle containing both ports, the number of old physical names is odd, so a singleton pair-label occurs. Free pair-label assignments multiply its count by independent factors2. Modulo8 only exactly one free block survives after the cycle weight2. That block must occur singly; all others are forced consecutive pairs. The port supplies a distinguished cyclic origin, preventing a rotational-count ambiguity. Every collapsed quotient cycle with k old labels has exactly k possible singleton blocks and two assignments each, hence contribution4k. This leaves exactly odd quotient cycles containing both ports modulo8.

With two odd cycles and both ports on one cycle, an old block occurs singly on that distinguished cycle. Swapping its physical names, even when the other name lies in the other cycle, is a valid fixed-point-free involution. That contribution is therefore zero modulo8.

For one port on each of two cycles, a free block similarly gives a factor2 against cycle weight4. Surviving blocks occur as whole consecutive pairs. The collapsed cycles cannot have only one old label, because uniform attachments forbid a two-cycle through a port and one block. Thus the removed old paths each have at least two vertices and all four endpoints are distinct. Direct expansion of their two contact factors gives the stated pair contrast: two mixed quadratic terms, no shared-label cubic, and cancelling quartics over F2. Its total mixed-quadratic trace is zero. The proof does not claim that arbitrary physical quartics vanish.

The quotient contrast divided by2 is integral before reduction. Combining the two collection types gives the coefficientwise mixed identity with factor4, and incoming E1 yields total trace4q modulo8. Every individual L_j is divisible by4. Odd q>=3 is used precisely to make the total trace nonzero. The q=1 identity comparison and the zero even-q trace are not mistaken for separators or collisions.

The battery contains all lower-subset queries needed for the stated finite differences. Its count is2q²+q+1. Coarse queries are exactly insensitive to pair-interior reversals. Once quotient orientation is known, isolated pair probes read the remaining bits without interference. The adaptive and fixed upper bounds2q²+2q+1 and2q²+3q+1 correctly count actual queries and all2^(q+1) distinct physical states.

## Independent exact computation

`independent.cpp` constructs both compared physical tournaments explicitly and uses literal integer Hamilton-path subset DP, not the producer's complement-query shortcut or the odd-cycle expansion. Its comparison cycle-pair polynomial is generated from the full contact products for disjoint directed paths in Q and Q-op separately, not by importing the simplified cancellation formula.

The independent run checks:

- all complete coarse coefficients for every quotient through order4 and six specified order5 quotients:81 cases,8,057 mixed coefficients in the modulo-eight identity;
- direct cancellation of the cycle-pair quartics, all shared-index cubics, and the total mixed-quadratic trace;
- additional exact trace controls at one order6 and two order7 quotients, retaining the even-order zero-trace boundary;
- 48,182 literal integer queries and46,576 interpolated coefficients in total;
- all484 physical named coefficients for all16 C3/pair-interior states, retained in `full_named_cube.json`, with exact vector distinctness and exact coarse-query invariance under all child bits.

The C++ `-O2` and `-O0` receipts and complete cube files are byte-identical. At the largest tested extension there are16 vertices, so exact path counts and finite differences fit signed64-bit arithmetic; no floating-point or hash-equality test is used.

`compare_retained.py` compares the independently reconstructed complete physical vectors against every one of the producer's8,192 source/target named truth values. All968 source/target physical coefficients and all300 nonzero contrast coefficients agree, including the36 quadratic,132 cubic and132 quartic terms. It adds44 literal backward-DP values on the producer's coarse controls, checks every final source/data integrity pin, and verifies the paired producer receipts. Its own normal and optimized receipts agree.

The producer's broader bounded trace census and sample results were inspected, but the independent census above is stated separately rather than presented as a second reproduction of every producer quotient.

## Resolved implementation hardening

The two producer trace programs originally counted an individual-divisibility failure without making that count affect the exit status. All saved counts were already zero, so no numerical result was wrong. Before this audit froze, both exit gates were changed to fail on either the total-trace failure count or the individual-divisibility failure count. Final normal/optimized source-bound receipts match, and the complete retained mathematical data are unchanged.

The final proof hash is `6d741cbef4d88c31a0d39f8c0df28ca50bebfd4c3195d6b33591bd5f39fc18c4`. Exact final checker and data pins are in `retained.json`.
