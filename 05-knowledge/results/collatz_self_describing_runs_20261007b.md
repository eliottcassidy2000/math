# Self-describing run clocks and exact Collatz prefix compilation

**Status:** CITED sequence definitions and supplied-paper mechanisms; PROVED elementary clock/compiler/height statements; FINITE-EXACT prefix counts and small exhaustive proof certificate. No limiting Kolakoski density, global Collatz convergence, or new paid-guard coverage is asserted.

## 1. Two sequences, and the count convention

[OEIS A000002](https://oeis.org/A000002) is the self-run-length sequence on \(\{1,2\}\) starting \(122112122\ldots\). [A071820](https://oeis.org/A071820) is the different sequence on \(\{2,3\}\) starting \(2233222333\ldots\). In each, the \(j\)-th term is the full length of the \(j\)-th maximal constant run. The starting symbol fixes the alternating run labels. The conjectured one-half frequency for A000002 is still recorded as open by the current OEIS entry; a finite count is not a proof of this limit.

The cumulative count of ones in A000002 is [A156077](https://oeis.org/A156077), not A071820. Its signed discrepancy \(E(n)=2O_n-n\) is [A088568](https://oeis.org/A088568).

The exact reproduction gives:

| Sequence | Prefix length | Runs intersecting prefix | Symbol transitions | Fully completed runs | Count of smaller symbol | Final observed/full run length |
|---|---:|---:|---:|---:|---:|---|
| A000002 | 200,000 | **133,321** | **133,320** | 133,321 | 100,010 | 1/1 |
| A000002 | 10,000,000 | **6,666,660** | **6,666,659** | 6,666,660 | 5,000,046 | 2/2 |
| A071820 | 200,000 | 80,001 | 80,000 | 80,001 | 100,007 | 3/3 |

Thus the stated 200,000-term number counts **transitions**, while the stated ten-million-term number counts **runs**. Every displayed large prefix ends at a complete run, so truncation does not explain the first off-by-one. At the genuinely truncated prefix \(12\) of A000002 there are two intersecting runs, one transition, and only one completed run: the second full run has length two but only its first symbol has been included.

The generator is checked against an independent two-block substitution for its first10,000 A000002 terms, and against the published A071820 initial terms. Every completed run in each reported prefix is checked against its purported self-describing length; a lookahead symbol verifies the final boundary. These are finite exact statements.

## 2. The two clocks and their boundary marker

For a self-run-length sequence \(k_1,k_2,\ldots\) on positive integers \(\{a,b\}\), \(a<b\), put

\[
S_j=\sum_{i=1}^j k_i,\qquad
R_N=\min\{j:S_j\ge N\},\qquad
\delta_N=S_{R_N}-N.
\tag{1}
\]

Then \(R_N\) counts runs intersecting the first \(N\) symbols, \(R_N-1\) counts transitions, and \(R_N-\mathbf1_{\delta_N>0}\) counts completed runs. The deficit obeys \(0\le\delta_N\le b-1\).

Let \(O_j\) count the smaller symbol among the **first \(j\) sequence terms**. The exact inverse-clock identity is

\[
N=bR_N-(b-a)O_{R_N}-\delta_N.
\tag{2}
\]

For A000002 this becomes

\[
3R_N-2N=2O_{R_N}-R_N+2\delta_N.
\tag{3}
\]

The frequency on the right is evaluated at the run index \(R_N\), not at \(N\). This keeps the self-description from becoming a spurious closed balance equation. If a limiting smaller-symbol frequency \(f\) exists, monotone inversion of \(S_j\) gives

\[
\frac{R_N}{N}\longrightarrow\frac1{b-(b-a)f}.
\tag{4}
\]

Conversely such a run-rate limit implies the corresponding letter-frequency limit, by evaluating it at \(N=S_j\). Thus the expected rates \(2/3\) for A000002 and \(2/5\) for A071820 both correspond to a one-half letter frequency, but they concern distinct clocks. The A000002 clock and inverse-clock viewpoint is also treated in [Bordellès–Cloitre, *Bounds for the Kolakoski Sequence*, definitions and Theorems3–4](https://cs.uwaterloo.ca/journals/JIS/VOL14/Bordelles/bordelles7r.pdf); the exact finite boundary convention is retained here.

There is a useful second observer, rather than a closed scalar contraction. For A000002 define

\[
E(r)=2O_r-r,\qquad D(r)=\sum_{i=1}^r(-1)^{i+1}k_i.
\]

At complete run endpoints, and then at arbitrary truncated endpoints,

\[
S_r=\frac{3r-E(r)}2,\qquad E(S_r)=D(r),\qquad
E(N)=D(R_N)-(-1)^{R_N+1}\delta_N.
\tag{4a}
\]

The proof counts the full odd-numbered runs as ones and the even-numbered runs as twos; truncating the last run removes its signed contribution. Self-description therefore changes ordinary discrepancy into an **alternating-index** observable. The run phase and deficit remain necessary, and (4a) gives no contraction of \(E\) in terms of itself. Exact controls cover run indices1 through2,048 and all endpoints1 through4,096.

## 3. A lossless ordered-run compiler

Closest proved mechanisms are [ordered affine routing](eleven_squares_routing_bridges_20261007.md), [valuation information/carry separation](valuation_information_spectrum_20261005.md), and [the exact child-prefix decoder](paper_child_closure_transfers_20261007.md). Their shared requirement is to retain chronology, the source guard, and the endpoint boundary. The present compiler specializes this to alternating runs of valuation symbols.

A run packet stores the alphabet, its first run label, the ordered full run lengths, and the final deficit \(\delta\). Only the last length is shortened by \(\delta\). For a valuation symbol \(c\) repeated \(d\) times, the formal Collatz carrier is

\[
F_c^d(x)=\frac{3^d x+B_{c,d}}{2^{cd}},\qquad
B_{c,d}=\frac{3^d-2^{cd}}{3-2^c}\in\mathbb Z.
\tag{5}
\]

It is obtained by powering \(\begin{pmatrix}3&1\\0&2^c\end{pmatrix}\); the geometric sum proves integrality. Chronological composition of carriers \((P,Q,B)\) and \((p,q,b)\) is

\[
(P,Q,B);(p,q,b)\longmapsto(pP,qQ,pB+bQ).
\tag{6}
\]

Consequently run powers compile exactly the same \((P,Q,B)\) as expanding every valuation. A nonempty positive valuation word with \(P=3^r,Q=2^A\) has the exact odd-source cylinder

\[
n\equiv(Q-B)P^{-1}\pmod {2Q}.
\tag{7}
\]

Final oddness enforces all prefix integrality and actual valuations. The standard proof peels the last affine step backwards, or inducts on the nested oddness cylinders. Source positivity and strict first-hit ROOT conventions are checked separately; the formal root loop must not be exported as an extra step.

Two minimal losses are explicit:

* Words \((1,2)\) and \((2,1)\) have the same valuation sum and symbol counts, but carries5 and7, and native cells \(11\pmod {16}\) and \(9\pmod {16}\). An unmarked starting run label loses the actual source guard.
* Full run lengths \((1,2)\) and starting label1 expand to \(122\). With terminal deficit1 they instead encode \(12\). Omitting that boundary marker changes the carrier and source cylinder.

This is an exact compiler, not a claim of bounded-memory validation or of new convergence coverage. A self-describing prefix is already algorithmically determined by its length; checking its arithmetic effect still requires the ordered carrier and an actual compatible source. The implementation exposes a distinct `authentic_self_prefix` test because a generic alternating-run packet need not be a prefix of either named sequence.

## 4. A finite forbidden-factor proof of the nine-letter bound

The following elementary recovery avoids assuming any limiting frequency. In A000002, \(111\) and \(222\) are impossible because all run lengths are1 or2. Any finite factor has a consecutive list of certainly complete runs: discard a length-one run at either exposed end, but retain a length-two end run since it cannot extend farther. Their lengths form a factor of the original self-describing sequence.

This gives a short deduction chain:

| Forbidden factor(s) | Certainly complete run-length factor forcing exclusion |
|---|---|
| \(12121,21212\) | \(111\) |
| \(112211,221122\) | \(222\) |
| Each of the six cubes \(u^3\), where \(u\in\{1,2\}^3\) is nonconstant | \(12121\) or \(21212\) |

Together with \(111,222\), these are twelve forbidden factors. The script exhausts **all512** length-nine binary words and verifies that the42 surviving necessary candidates have four or five ones. This is a finite proof certificate of the implication; it does not claim that every survivor occurs in A000002. Similar forbidden factors appear in [Bordellès–Cloitre, Lemma8](https://cs.uwaterloo.ca/journals/JIS/VOL14/Bordelles/bordelles7r.pdf). Broader factor languages of iterated run differentiation are studied in the current primary paper [Cassaigne–Henry, *The complexity of smooth words over binary alphabets*, version4,2026](https://arxiv.org/abs/2603.10733); no factor-language theorem there is needed for this twelve-word certificate.

Reading an A000002 factor as nine Collatz valuations therefore gives total valuation at most14 and coefficient

\[
\frac{3^9}{2^A}\ge\frac{19683}{16384}>1.
\tag{8}
\]

Every actual nine-step block of this prescribed infinite valuation stream would increase its positive input by at least this multiplicative factor, since its affine carry is positive. Thus a positive integer realizing the entire stream would have an exponentially growing subsequence. No such positive realization is asserted. Early shorter prefixes may still pay: the first three values \((1,2,2)\) take43 to37.

Necessary local tests are not a self-description certificate. The word \(12212\) has the correct initial seed \(122\) and avoids all twelve forbidden factors, but it is not the length-five A000002 prefix \(12211\). Its third run closes too early. This is the cheap hostile against replacing an exact run-clock check by a local forbidden-factor screen.

## 5. Opposite drift, and a genuine finite/infinite boundary

For any actual positive odd step whose valuation is at least two,

\[
U(n)-1\le\frac34(n-1).
\tag{9}
\]

A strict ROOT hit from a positive odd input greater than1 has last valuation even and at least4: its predecessor is \((2^a-1)/3\), which is integral only for even \(a\), and \(a=2\) gives the already-root input1.

It follows that an actual nonroot prefix with every valuation in \(\{2,3\}\) cannot hit ROOT during the prefix. At length \(r\), its endpoint is at least3, so (9) gives

\[
n\ge1+2(4/3)^r.
\tag{10}
\]

**Therefore no positive integer realizes the entire A071820 stream as actual Collatz valuations.** This statement is proved by an elementary height rank; it uses no conjectural letter density and no general Collatz convergence premise. Nevertheless every finite prefix has an infinite native arithmetic progression (7). For the initial all-two prefixes, the root-padding point1 must be removed; the other positive members remain.

These coherent native cylinders have a unique intersection in the odd2-adic integers, because their binary precision tends to infinity. Equation (10) shows that the intersection contains no positive integer. Finite guards at every depth therefore do not furnish one immutable integer satisfying every depth. The A000002 stream also cannot itself be a strict finite ROOT word using only1 and2; its global positive realization is left open here, and (8) describes the consequence if it existed.

The preserved object is the exact finite valuation word and its native cylinder. The destroyed information in a mere run rate is the ordered carry and the supplied source. The required sidecars are the source, start label, and final boundary deficit. This gives a reusable prefix compiler and a sharp hostile to a purported closure argument, not a universal terminal rule.

## 6. The two supplied papers: what can actually transfer

The user accepts the following headline theorems. We inspect their constructions rather than re-audit those theorems.

| Paper and inspected interface | Mechanism | Exact transfer boundary |
|---|---|---|
| [*Polynomial removal fails for ordered binary matrices*, paper-22.pdf](C:/Users/Eliott/Downloads/paper-22.pdf),12pp; pp.1–11, Figure1 p.4 visually inspected | Both row/column orders and zeros/ones remain part of an occurrence. A set meeting every old copy need not repair the matrix: editing those cells creates new copies (Remark4.7,p.10). | Our chronological run carrier keeps order and boundary data. The \(12212\) control proves that eliminating the listed local obstructions does not authenticate the global self-description. The matrix theorem itself is not promoted to a complexity lower bound for this compiler or for Collatz. |
| [*Sharp binary-information contraction on the discrete cube*, main-4.pdf](C:/Users/Eliott/Downloads/main-4.pdf),29pp; pp.1–6,10–11,15–17, Figure1 p.11 visually inspected | The local comparison retains the output mean and the original child entropies. A common-output mixture transports information along a stated probability model; it does not silently identify the mixture with the original channel. | Exact source guards remain authenticated data. In an explicitly introduced independent bit-flip observation model with error \(0<\epsilon<1\), every noisy output has positive probability from every input. Hence no observation alone can give a zero-error positive certificate for a nonconstant residue guard. The script checks this full-support obstruction exactly at \(\epsilon=1/4\) on three bits. This is a separate elementary statement, not a stochastic model for an actual Collatz orbit. |

For the last obstruction, if a decoder accepts some observed word with positive probability, the same observed word has positive probability from an invalid input. Thus zero false positives forces it never to accept. Additional retained exact data or a separately authenticated receipt is necessary. Average information bounds do not alter this support argument.

## 7. Reproduction and scope

[Exact program](../../04-computation/experiments/collatz_self_describing_runs_20261007b.py) and [saved output](collatz_self_describing_runs_20261007b.out).

```text
python -B -X utf8 04-computation/experiments/collatz_self_describing_runs_20261007b.py
python -B -O -X utf8 04-computation/experiments/collatz_self_describing_runs_20261007b.py
```

The explicit universes are the reported10-million/200,000-digit prefixes; the complete512 length-nine words and twelve-factor deduction chain; both named prefix families at every length0 through256; two positive native source lifts per compiled packet; and all64 input/output pairs of a three-bit noise channel. The514 ordered packets are compared against independent symbol-by-symbol affine composition and literal Collatz replay. Types, ROOT padding, missing start marking and terminal truncation have explicit hostile controls. Optimization does not remove any check.

Both modes perform **11,905 exact checks** and reproduce the saved output. Normalized-LF SHA256: `9482db77f55a7bb89db0c1dc1d0441c94378ad6273195083d0c88be2143dcc1b`.

The output is a finite exact measurement and a proof of stated compiler/rank properties. It supplies no new actual paid-guard coverage and no empirical-to-universal density inference.
