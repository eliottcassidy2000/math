# Run receipts in the inherited credit potential

2026-10-11. **PROVED:** sharp one-credit run interface, immutable-source payment and finite episodes in the stated H/G/R library, exact native-language separation, and the pre-funded entry ray. **FINITE-EXACT:** the declared controls. **OPEN:** arbitrary source coverage and ROOT for unspecified terminal ports.

Artifacts: [program](../../04-computation/experiments/collatz_run_credit_interface_20261011.py), [output](collatz_run_credit_interface_20261011.out).

## 1. Inheritance and the missing native interface

The [adaptive potential theorem](adaptive_credit_potential_20261004.md), §§2–4, already supplies
\[
H(x)=\frac{729x+669}{1024},\quad x\equiv155\pmod{2048};
\qquad G(x)=\frac{9x+5}{8},\quad x\equiv11\pmod{16}.
\]
H earns two credits, G spends one, and the immutable-source account is
\[
\mathcal E(x,k)=(x+5)(9/8)^k.
\]
H is an authenticated common-future dependency; G is the actual valuation word $(1,2)$. The prior theorem gives exact neutrality for G and sharp factor $\kappa_H=2349/2560$ for H.

The new [projective run family](collatz_projective_receipts_20261009.md) supplies
\[
R_r(n)=a_r(n+1/2)-1,\qquad a_r=\frac{2^{r+2}}{3^r},
\quad n\equiv3\pmod4,\quad4\le r\le v_3(2n+1).                \tag{1}
\]
Its actual receipt is $n\xrightarrow{(1)}J\xleftarrow{(1^r,3)}R_r(n)$, with $1<R_r(n)<n<J$. Its source domain is exactly $n\equiv283\pmod{324}$. This is already contained in the old $G_5$ guard: this integration adds no source-guard mass.

The interface preserves source, native residue, potential and actual common-future words. Its scalar inequality forgets which next action is legal. The needed sidecar is the native source state and labelled action history. The cheap hostile is a proposed second run credit; the stronger stopping boundary is an exactly neutral receipt cycle. Our live board is **immutable source / credit potential / native language / actual receipt / terminal ROOT obligation**.

## 2. Sharp one-credit import

Every native run input is at least 283, and $a_r\le64/81$. For fixed $n>0$, the quotient
\[
\frac{R_r(n)+5}{n+5}=\frac{a_r(n+1/2)+4}{n+5}
\]
increases with $a_r$. Consequently
\[
\frac{R_r(n)+5}{n+5}
\le\frac{64}{81}+\frac{4}{9(n+5)}
\le\frac{19}{24}.
\]
Equality is attained at $r=4,n=283,R_4(n)=223$. Awarding one credit therefore gives the sharp uniform bound
\[
\boxed{\frac{\mathcal E(R_r(n),k+1)}{\mathcal E(n,k)}
\le\frac{57}{64}=:\kappa_R<1.}                              \tag{2}
\]

Two credits fail this individual potential test throughout the native $r=4$ family:
\[
\frac{\mathcal E(R_4(n),k+2)}{\mathcal E(n,k)}
=1+\frac9{16(n+5)}>1.                                      \tag{3}
\]
At 283 this is $513/512$, decreasing to 1 along the infinite native ray. This is **not** a legal RGG size-overspending witness: the first G after any R is already nonnative. Nor does (3) prove a maximal budget for all possible guarded continuations or other potentials.

## 3. Unified payment, with an exact language boundary

Start at immutable source $N$, current value $x=N$, and credit zero. Permit only native H, G and $R_r$ actions, with credit changes $+2,-1,+1$, respectively, and no negative prefix credit. Let $h,g,s$ count these actions. At every nonempty completed-action prefix,
\[
k=2h+s-g\ge0,
\qquad116\le\mathcal E(x,k)
\le(N+5)\kappa_H^h\kappa_R^s<N+5,
\qquad x<N.                                                \tag{4}
\]
The first action must be H or R. Their outputs are at least 111 and 223, respectively; G increases its input, giving the lower bound. Induction on the three exact potential contracts proves the upper bound. Since $k\ge0$, $x+5\le\mathcal E(x,k)$, proving payment against the unchanged original source. These assertions concern action endpoints, not every internal Collatz state of a receipt.

Because $\kappa_R<\kappa_H<1$, the number $h+s$ is bounded by the largest integer $b$ satisfying
\[
116\le(N+5)\kappa_H^b.
\]
The program computes this bound using exact fractions. Also $g\le2h+s\le2(h+s)$, so there are at most $3b$ executed actions. No infinite legal episode exists in this library. This bounds executed actions, not the running time of an unspecified search for a next action, and does not assert ROOT at a stopped state.

**The native grammar is more restrictive than (4).** Every R output is
\[
R_r(n)=2^{r+1}t-1\equiv15\pmod{16},
\]
whereas H and G require $11\pmod{16}$. They cannot follow R. Conversely H outputs are $3\pmod9$, and G outputs are $4\pmod9$. For native G,
\[
2G(x)+1=9(x+1)/4.                                         \tag{5}
\]
After an H/G action, $x\equiv3$ or $4\pmod9$, hence a further G output has exactly two ternary factors in $2G(x)+1$, below the four required for R. An H output itself is not $4\pmod9$ and cannot enter R. Thus every nonempty zero-credit native episode is either an H/G episode starting with H, or a pure R episode. The algebraically unified account does **not** create new mixed switching paths; its R credits are stranded in this specific library.

The pure R family does genuinely repeat: for every $t\ge0$,
\[
13891+26244t\ \xrightarrow{R_4}\ 10975+20736t
\ \xrightarrow{R_4}\ 8671+16384t.
\]
The two credits remain attached to the same original source. The inherited run packet carries the exact composed actual words and raw clock differences $(-12,-8)$.

## 4. The five-letter boundary and a pre-funded incoming ray

The [native carry interfaces](collatz_carry_interfaces_20261004.md), controller table, give the larger old library:

| Letter | Affine action | Source guard | Output modulo 9 |
|---|---|---|---:|
| H | $(729x+669)/1024$ | $155\pmod{2048}$ | 3 |
| G | $(9x+5)/8$ | $11\pmod{16}$ | 4 |
| A | $(81x+85)/128$ | $187\pmod{256}$ | 2 |
| B | $(81x+73)/128$ | $7\pmod{256}$ | 5 |
| L | $(9x-3)/16$ | $219\pmod{256}$ | 6 |

R outputs $15\pmod{16}$, blocking all five native inputs. For an old letter immediately before R, only G has the required output residue $4\pmod9$. If at least one other old letter precedes that G, its input lies in $\{3,4,2,5,6\}\pmod9$, never 8. Equation (5) then gives ternary fuel at most three. Therefore **no native old-word prefix of length at least two can enter R**. This guard statement does not need a credit policy and is not an assertion about all other known rules.

A lone G can enter R when it starts from an independently funded state. For every $t\ge0$,
\[
251+288t\ \xrightarrow{G}\ 283+324t
\ \xrightarrow{R_4}\ 223+256t.                            \tag{6}
\]
All native guards are exact. With entry credit one, the account evolves $1\to0\to1$, and its potential contracts by at most $57/64$. To use this mid-episode under an original source $N$, one still needs an authenticated path to that entry state and the entry bound $\mathcal E(251+288t,1)\le N+5$. Formula (6) alone supplies neither. Starting the zero-credit episode at its first state would reject G.

The distinction is visible in a neutral recycling loop. On positive odd $n\equiv4\pmod9$, the inherited inverse-G dependency
\[
C(n)=(8n-5)/9,
\qquad C(n)+5=(8/9)(n+5)
\]
could earn one credit with **equality** in the potential. Its child is $11\pmod{16}$ and G returns it exactly to $n$. Thus
\[
(283,0)\xrightarrow{C}(251,1)\xrightarrow{G}(283,0)
\]
can repeat forever with authenticated actual common-future steps and constant potential. C is deliberately not part of (4). Adding it without a separate cycle budget or a proved cycle-removal procedure would invalidate the finite-episode conclusion and strict payment at every nonempty prefix. This is receipt recycling, not progress towards a new terminal proof. It also explains why the [common-rank policy certificate](collatz_receipt_flow_geometry_20261009.md) requires strict expected decrease rather than only nonincreasing energy.

## 5. Receipt-aware adapter and verification

`Action` retains H/G type or the chosen run $r\ge4$. `transition` checks the exact source guard and returns both credit change and an actual common-future receipt. It does not by itself authorize spending credit. `execute` starts with zero credit, keeps $N$ immutable, and checks every credit/potential prefix. `continue_episode` authenticates the existing packet and extends the same source history. Caches cannot replace their provenance: `audit_episode` replays the entire labelled packet and rejects boolean/float aliases.

`discharge` first verifies a supplied strict ROOT word for the final state. It then applies the inherited clock-bounded receipt-path exporter. There is no orbit search for an open terminal. Controls supply the frozen ROOT(223) word and an explicit ROOT(111) word extracted from the already certified ROOT(27) suffix through 167; production only replays them. Outputs ROOT(283) and ROOT(155) have respective rank/cost $(20,40)$ and $(30,55)$. Empty ROOT(1) is allowed; padding, wrong suffixes and forged states are not.

The declared universe is all funded H/G words of lengths 1–6 at their native lifts 0, 1 and 7 (168 episodes); 845 native R instances, $r=4,\ldots,16$ with parameters 0–64; 100 two-return rays; exact five-letter symbolic output residues with 64 literal lifts each; and 128 instances of (6) and the C/G cycle. Positive receipts, incompatible native transitions, invalid types, insufficient credits and ROOT padding are separate controls. This finite universe tests the APIs; the unbounded claims are proved above.

    python -B -X utf8 04-computation/experiments/collatz_run_credit_interface_20261011.py
    python -B -O -X utf8 04-computation/experiments/collatz_run_credit_interface_20261011.py

The result imports a sound potential contract and identifies why the current native languages still do not mix. It does not enlarge the proved source domain, ground a giant child, or prove general Collatz convergence.

Normal and optimized runs agree exactly with the saved output: **10,450 checks**.
LF-normalized SHA256: `0060693288a6cb6f82c9457067a066445126d23ed01ad606d6194debc319b96e`.
