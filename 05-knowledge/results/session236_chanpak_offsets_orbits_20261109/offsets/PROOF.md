# Genuine tournament offset transport and the periodic-word barrier

Status: PROVED inherited identities, rederivations and scoped corollaries; FINITE-EXACT independent replay. The near-join law and mixed-word cyclic drift are recovered prior work, not new theorems. Chan–Pak remains OPEN.

This recovery preserves the mathematical content and provenance corrections of the lost Session236 proof. Its single-final-bridge obstruction does not apply to multiple bridges; in particular it does not contradict concurrent S19 multibridge negative-offset constructions.

## 0. Inheritance and scope

The bridge identity is inherited from THM-4875, `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`. The nonnegative port state is THM-5070, `01-canon/theorems/THM-5070-port-chain-transfer-for-odd-cycle-packing-counts-and-ternary-coverage-to-a-million.md`.

THM-4865(d), `01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci-h-k-equals-h-k-1-plus-h-k-2-plus-h-k-3-for-every-rooted-base-the-lineage-c2-law-at-any-vertex-and-the-three-root-response-numerations-attain-every-odd-number-in-69-702189.md`, explicitly records HYP-9453's near-join formula h+2^k B(T,u). Its B(T,u) is exactly R=E+O below. The singleton-final-block case and its shift transport are therefore inherited.

THM-4882(1),(3), `01-canon/theorems/THM-4882-signed-one-apex-responses-slope-and-parity-classification-session198-clock-merged-positive-count-barrier.md`, and its source `05-knowledge/results/signed_responses_slope_classification_clock_barrier_20261105.md`, section4, paragraph beginning “On a strong input” (lines166–170 in the audited source), already prove positive cyclic drift by run decomposition, crediting Session198. Section4 below rederives that mechanism in word-integer coordinates and states its fixed-offset consequence. It does not claim a new all-length positivity theorem.

THM-5055/5056's graph offset moves are a comparison, not a dependency or a transfer to tournaments. Session233's moving-core windows already show that offset variation does not imply arbitrary target coverage. The hostile examples are forbidden complete counts7 and21 and equal complete counts with unequal port responses. The retained sidecar is the marked port, its response, and the padding location.

## 1. Recovered fixed-port offset move and arbitrary-B composition

Let H(T) count directed Hamilton paths, with H(empty)=1. For a port u put

E(T,u) = sum H(T-V(P)) over even-order paths P starting at u,
O(T,u) = sum H(T-V(P)) over odd-order paths P starting at u.

The singleton path is included. These are THM-5070's x=2 coordinates. Let R=E+O. Let T_m be the ordinal sum of T with a transitive tournament of m sinks, keeping the old port u. For every m>=1,

(H,E,O)(T_m,u) = (H(T), 2^(m-1)R, 2^(m-1)R).                 (1)

Proof. Adding a sink fixes H of each remaining induced tournament. A path from u either stays in T or extends an old path by the new sink. The alternatives exchange parity, giving E'=O'=E+O. Iterate. This is fixed-port padding, not the moving-port tower.

Fix a nonempty tournament B and entry w. Form J_m from T_m=>B by reversing exactly u->w. Define q_e,q_o as the sums H(B-V(Q)) over even-/odd-order paths Q ending at w. The inherited seam law gives

H(J_m)=H(T_m)H(B)+2(E(T_m,u)q_o+O(T_m,u)q_e).

Consequently, for every m>=1,

H(J_m)=C+2^m A,
C=H(T)H(B),  A=R(q_e+q_o).                               (2)

Appending r more sinks before the same final bridge transports C from A2^m to A2^(m+r). This gives actual physical shifts2 and3, and any positive shift, on this decorated family. The singleton-B scalar law is already THM-4865(d)/HYP-9453. Formula(2) is its arbitrary-B composition with the inherited seam identity, not an independent completeness engine.

## 2. Exact obstruction for this single-final-bridge family

Both A and C are positive integers. Moreover, C=H(T=>B) is itself a tournament count. Therefore this particular parameterization cannot carry a negative offset or C=7 or21; the exclusions use THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`. This neither excludes targets represented at another scale nor restricts multibridge constructions.

The unit-coefficient case A=1 forces C=1. Every term in R is positive. The singleton path contributes H(T-u)>=1; an out-neighbor of u supplies another positive term. Thus R=1 iff u is a sink and H(T-u)=1. A tournament has H=1 iff transitive: the transitive case has one path, whereas any nontransitive tournament contains a directed triangle and has H>=3 by the odd-cycle formula. Hence R=1 iff T is transitive and u its sink. By converse, q_e+q_o=1 iff B is transitive and w its source. Both have H=1.

Clock shifts A=2^r allow more rays but retain the C barrier. Two explicit five-vertex seeds (adjacency bitmasks, vertices0..4) are:

- T rows [12,29,24,16,1], port1: H=5, R=32. With B singleton, H(J_m)=2^(m+5)+5 for m>=1, on m+6 vertices.
- T rows [14,28,24,16,1], port2: H=9, R=16. With B singleton, H(J_m)=2^(m+4)+9 for m>=1, on m+6 vertices.

These are constructions, not new scalar-value claims. This note supplies the arbitrary-B composition, explicit family/clock corollaries and independent replay; the underlying near-join transport is inherited. Equal H does not determine A: the source and sink ports of a transitive two-vertex tournament have R=2 and1.

## 3. Inherited one-apex zero tails and coefficient loss

Take a transitive chain v_1->...->v_k and apex u. The i-th bit of W, read most significant first, is1 iff u beats v_i. Appending a bit on the right changes the integer to2W+b. THM-5070's recurrence is H(W)=1+2a(W), where a(0)=0, a(2W+1)=a(W), a(2W)=a(W)+ceil(W/2). Thus

H(2W+1)=H(W),
H(2W)=H(W)+W+(W mod2).                                  (3)

For even W and m>=0 this gives

H(2^m W)=2^m W+[H(W)-W].                                 (4)

This is inherited algebra. Negative offsets genuinely occur: W=6 has H=5, so H(6*2^m)=6*2^m-1. The leading coefficient6 is not a power of two. If positive even W is a power of two, H(W)=W+1, so the normalized pure-dyadic ray has offset1 only. A scalar offset analogy loses this leading-coefficient sidecar.

## 4. Recovered cyclic drift and fixed-offset corollary

This is a coordinate-level rederivation of Session198/THM-4882's cyclic positive drift. Their coordinates satisfy U=W+1 and epsilon_i=1-b_i: bit1 is O and bit0 is Z. Their H increment is epsilon_i*(U-epsilon_(i-1)), where epsilon_(i-1) is the cyclic ending indicator. This is exactly the increment in(3). The drift D below is their Delta_cyc. The complemented-digit run proof and all-period quantifier are inherited, not new.

Fix a binary word b=(b_1,...,b_k), q=2^k, and its integer value d. Repeatedly append b to W_0. First assume W_0 mod2=b_k. Otherwise discard the first block; the condition then holds. Fixed nonnegative integers alpha,beta satisfy

W_(m+1)=qW_m+d,
H_(m+1)=H_m+alpha W_m+beta.                              (5)

To compute them initialize c=1,e=0,p=b_k,alpha=beta=0. At each digit b_i=0 add c to alpha and e+p to beta. Then set (c,e,p)=(2c,2e+b_i,b_i). This is direct substitution in(3).

Define

A=alpha*(W_0+d/(q-1))/(q-1),
D=beta-alpha*d/(q-1).

Summing(5) yields, for every m>=0,

H_m=A q^m + Dm + (H_0-A).                                (6)

Inherited drift and corollary: if b contains both digits, D>0. Thus no mixed periodic suffix of any length transports a constant offset from its dyadic exponential clock. This also holds on every infinite arithmetic progression of repetition indices. An all-zero suffix has D=0 and gives(4); an all-one suffix keeps H constant.

Proof of strict positivity. Evaluate the block increment at the formal fixed point W*=-d/(q-1). Its value is D. For a mixed word this fixed point and all cyclically rotated fixed points lie strictly in(-1,0), being negatives of mixed periodic binary fractions. Group zeros cyclically into maximal runs, each preceded by1. For a run of length l write the value before the run as -r. The first zero contributes1-r, and successive zeros contribute -2r,-4r,...,-2^(l-1)r. Their sum is1-(2^l-1)r. The value after the run is -2^l r in(-1,0), so0<r<2^(-l) and the contribution is greater than2^(-l)>0. Ones contribute zero. At least one run exists, proving D>0.

The fixed-offset conclusion follows from(6): the recurrence f_(m+2)-(q+1)f_(m+1)+qf_m vanishes for an exponential plus constant, but equals(1-q)D here. A subsequence with positive integer step retains a positive linear drift equal to that step times D.

For the smallest mixed examples01 and10, respectively (q,d,alpha,beta)=(4,1,1,1) and(4,2,2,2), both give D=2/3. The initial parity condition matters: W_0=0 with suffix01 has actual first increment0 rather than beta=1. Discarding the first block restores the condition.

This restricted corollary does not obstruct the full15-matrix orbit, nonperiodic digit choices, varying blocks, multibridge constructions, or other tournament operations. It applies the already proved cyclic drift to interpret finite suffix searches. No new all-length drift theorem or proof of Chan–Pak is claimed.

## 5. Reproduction and status

Run `python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets/replay.py`.

The replay and results survived on main906c04d and were not rewritten during recovery. They use no repository imports or external packages. Hamilton paths are counted by subset endpoint DP independently of the identities; path-response sums are directly enumerated. The complete rooted labelled universe of orders1..3 checks87 padding identities and1682 bridges, with50 deterministic random larger bridges. Both concrete rays are checked at m=1..4. The replay also checks2046 periodic words through length10 by rational arithmetic and254 one-apex tournaments through word length7 by independent DP. Hostile controls reject m=0 extrapolation, loss of response and initial parity, and confusion of negative offsets with a unit dyadic coefficient. Normal and optimized-Python runs agreed. The finite checks corroborate the separately proved all-order identities; they do not imply completeness.
