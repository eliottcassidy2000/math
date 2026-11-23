# Seven Kaprekar walks in one exact base

At **base 109394**, one can run all seven solved width-two orbit walks from [base 10](KAPREKAR-METAGRAPH.md), [base 12](BASE12-CYCLE.md), [base 14](BASE14-CYCLES.md), and [base 16](BASE16-CYCLE.md) (including [its second orbit](BASE16-SECOND-CYCLE.md)). Every original arc schedule, checkpoint, and detour is transported unchanged except for a scalar multiplier on the exact Hamiltonian-path count. The transported lengths are upper bounds, not minimum claims in the new base.

## Minimal arithmetic base, and the role of the detours

In an even source base \(b\), write \(m=b+1\); its positive periodic gaps use the odd representatives of the nonzero sign classes modulo \(m\). The update is \(g\mapsto |2g-m|\), and the exact output count is \((b-1)g\). In target base \(B\), let \(M=B+1\). If \(M=dm\), the embedding \(g\mapsto dg\) commutes with the update, and the corresponding count multiplier is
\[
K_b=\frac{(M-2)d}{b-1}=\frac{M(M-2)}{b^2-1}.
\]
The ordinal join of a tournament of count \(K_b\) and the source walk realizes this multiplication on every literal stage.

For any positive uniform scaling of a complete first orbit, the image of gap 1 is some positive integer gap \(d\), and its successor must be \(d(b-1)\). Therefore \(|2d-M|=d(b-1)\). For \(b\ge4\) this forces \(M=(b+1)d\); for \(b=2\) the other algebraic possibility \(M=d\) is outside the legal gap range. An odd target base has no nonzero odd outputs, so cannot host these walks.

For a finite set of source bases containing their first orbits, set \(L_0=\operatorname{lcm}(b+1)\) and \(A=\operatorname{lcm}(b^2-1)\). The exact **arithmetic** conditions for simultaneous integral multipliers are
\[
L_0\mid M,\qquad A\mid M(M-2),\qquad M \text{ odd}.
\]
Since \(\gcd(M,M-2)=1\), each odd prime power \(p^a\Vert A\) must divide exactly the relevant side: \(M\equiv0\) or \(2\pmod{p^a}\), with \(0\) forced when \(p\mid L_0\). These local alternatives combine independently by the Chinese remainder theorem. Attaining the integer multipliers as Hamiltonian-path counts is the additional construction problem. In particular an individual lift \(M=(b+1)((b-1)x+1)\) is on the \(2\)-branch at every prime dividing \(b-1\), so it cannot be shared with a base whose \(b'+1\) forces such a prime onto the \(0\)-branch.

For \(b=10,12,14,16\), every prime in \(A\) already divides \(L_0\), and thus \(M\) must be an odd multiple of
\[
L=\operatorname{lcm}(99,143,195,255)=109395.
\]
The lower bound survives even if a multiplier is initially allowed to be a positive real number, provided **all stages of these same packets** must scale to integer counts. The first three selected packets have coprime readings already at their starting checkpoint and one detour: \((9,83),(11,35),(13,85)\). Any real multiplier sending those two integers to integers is itself integral. The resulting factors \(9,11,13\) in \(M\), together with \(L_0\mid M\), already force \(L\mid M\), regardless of the multiplier for base 16. Arithmetic gap paths alone would embed already at \(M=L_0=36465\); it is the absolute decimal detour, for example \(83\), which detects the missing factor three. This is not a lower bound on other methods that change the packets or do not use uniform scaling.

## Exact count factors at the minimum

At \(M=L\), write \(K_b=C_b h\) with \(h=M-2=109393=139\cdot787\) and \(C_b=L/(b^2-1)\). The common factor is realized by two ordinal blocks:

| Height | Outgoing row masks |
|---:|---|
| 139 | \((10,108,73,112,39,69,17)\) |
| 787 | \((478,284,296,304,36,323,158,318,80)\) |

| Source base | \(C_b\) | \(K_b\) |
|---:|---:|---:|
| 10 | \(1105=5\cdot13\cdot17\) | 120879265 |
| 12 | \(765=5\cdot9\cdot17\) | 83685645 |
| 14 | \(561=3\cdot11\cdot17\) | 61369473 |
| 16 | \(429=3\cdot11\cdot13\) | 46929597 |

The seven cycles in the large base are distinct. For a gap \(g\), its reduced denominator \(q=M/\gcd(g,M)\) is invariant under doubling up to sign. The seven source cycles have \(q=11,13,15,5,3,17,17\), and the two cycles for 17 are disjoint. In fact these are the **complete cycles at these six denominators**: in the denominator-\(q\) class the number of cycles is \(\varphi(q)/(2\operatorname{ord}^{\pm}_q(2))\), where \(\operatorname{ord}^{\pm}_q(2)\) is the first positive \(t\) with \(2^t\equiv\pm1\bmod q\). It follows by identifying the class with \((\mathbb Z/q)^\times/\{\pm1\}\). Altogether base 109394 has 701 nonzero gap cycles; only the seven specified here have the supplied tournament walks. The base-2 and base-4 imports coincide with the base-14 fixed and two-phase cycles, respectively, rather than adding two new cycles.

## One host, two ways to program the multiplier

The ordinal layout has two fixed blocks of total order 16, a selector of height \(C_b\), and a moving block on ten vertices. Append ordinal singleton padding inside the moving block for the shorter packets. To switch packets, reverse the disagreeing arcs of this block; its cross arcs to all other blocks never change. Then perform the chosen packet with the selector held fixed. Program edits are not counted as part of the orbit replay.

**Five program flips (41 vertices total).** In the selector ordinally join three slots
\[
(6,12,8,1),\qquad(10,20,25,18,1),\qquad(40,53,57,50,1,16)
\]
of initial heights \(5,13,17\). Reverse slot-local arcs \((0,1)\) in the second slot; then \((0,1)\) in the first and \((2,4)\) in the second; then \((4,5)\), \((0,1)\) in the third. This gives the exact programming trajectory
\[
1105\to765\to459\to561\to297\to429.
\]
A later program can undo any prefix. The intermediate states 459 and 297 are programming states, not source phases.

**Five program flips on a smaller host (35 vertices total).** A single nine-vertex selector starts from
\[
(444,401,138,322,76,414,423,24,148).
\]
Reverse \((1,2),(1,4),(3,6),(0,5),(7,8)\) in order. Its exact counts are
\[
1105\to765\to955\to561\to405\to429.
\]
This retains the same program length while removing six vertices from the factored selector, so the host uses \(16+9+10=35\) vertices. Nine is minimal for the **standalone selector** since its count 1105 exceeds the [exact order-eight maximum 661](DECIMAL-SHELF.md). The five flips are an upper bound for the program, not a proved minimum, and there is no minimum claim for the entire host.

## The next admissible base, and a sufficient condition beyond it

For any odd \(y>0\) with attained ordinary counts \(y\) and \(Ly-2\), take \(M=Ly\) and \(B=Ly-1\). Use the same selector and packets with the fixed factor \(y(Ly-2)\). This sufficient condition separates the problem into two exact counts; no necessity or infinitude is claimed.

The next admissible \(y\) is 3. At base \(328184\), \(Ly-2=328183=307\cdot1069\); join a directed triangle and the following two blocks to obtain the fixed factor.

| Height | Outgoing row masks |
|---:|---|
| 307 | \((222,124,40,240,36,65,148,54)\) |
| 1069 | \((478,316,408,48,352,325,142,314,72)\) |

Three more exact common bases follow from the same condition. The row masks for all their factors, including the \(y\)-blocks, are in the replay.

| \(y\) | \(B\) | Factorization of \(B-1\) | Vertices for fixed factor \(y(B-1)\) |
|---:|---:|---|---:|
| 9 | 984554 | \(809\cdot1217\) | 23 |
| 15 | 1640924 | \(31\cdot43\cdot1231\) | 26 |
| 27 | 2953664 | \(577\cdot5119\) | 25 |

Thus all **five verified bases** fit in one 45-vertex layout: reserve 26 vertices for a programmable common-factor block, nine for the selector, and ten for the moving packet. Ordinal padding preserves the count of any common-factor program. Programming flips only internal arcs of that block and carries no stated cost bound. These are examples, not the first five feasible bases: \(y=5\) and \(y=7\), for instance, are not settled here.

Run [`node kaprekar-common-base.js`](kaprekar-common-base.js) to verify the fixed factors, selector, all literal source walks and exact BigInt heights in both common bases.
