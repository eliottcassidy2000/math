# Genuine tournament offset transport and an all-length periodic-word barrier

Status: PROVED identities and scoped obstructions; FINITE-EXACT independent replay.
Chan–Pak remains OPEN. No claim here transfers graph dominating-set results to tournaments.

## Inheritance and concept board

The inherited bridge is THM-4875, `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`.
The nonnegative state is THM-5070, `01-canon/theorems/THM-5070-port-chain-transfer-for-odd-cycle-packing-counts-and-ternary-coverage-to-a-million.md`.
The graph comparison is THM-5055/5056, not a dependency. Session233's moving-core windows already show that genuine offset variation does not imply target coverage. The hostile examples are the forbidden *complete* counts 7 and21, and equal complete counts with unequal port responses. The least-used sidecar here is fixed-port padding, which is outside the moving-port singleton tower.

Live concepts: graph dyadic offsets; rooted packing response; one-arc joins; fixed-port ordinal padding; one-apex digit words. The maps below preserve actual complete Hamilton-path counts, but scalar projection loses the response coefficient. No artificial tournament relation is introduced: vertices and arcs are the actual constructed tournament.

## 1. Fixed-port padding is an actual offset move

Write H(T) for the number of directed Hamilton paths, with H(empty)=1. For a port u define

E(T,u) = sum H(T-V(P)) over paths P starting at u of even order,
O(T,u) = the analogous sum over odd-order paths (including the singleton).

These are THM-5070's x=2 coordinates, by its odd-cycle formula. The proofs below use only the displayed Hamilton-path sums and the inherited seam identity.

Let T_m be the ordinal sum of T with a transitive tournament of m sinks, and **keep the old port u**. Then for every m>=1,

(H,E,O)(T_m,u) = (H(T), 2^(m-1)R, 2^(m-1)R),  R=E(T,u)+O(T,u).                 (1)

Proof: one appended sink fixes every remainder's H. A path from u either stays in T or extends a path from u by the new sink. The two alternatives exchange parity, so E'=O'=E+O. Iteration proves (1). This is not the moving-port tower matrix.

Take any fixed nonempty tournament B, an entry w, and form J_m from T_m => B by reversing exactly u->w. Let q_e,q_o be the sums H(B-V(Q)) over even/odd paths Q ending at w. Then

H(J_m) = C + 2^m A,  C=H(T)H(B),  A=R(q_e+q_o),  m>=1.                 (2)

Indeed, the seam identity is H(T_m)H(B)+2(E(T_m,u)q_o+O(T_m,u)q_e). Substitute (1). The actual operation of appending r further sinks before the same final bridge transports C from A2^m to A2^(m+r). Thus shifts2 and3 are genuine physical tournament moves on this decorated family, just as arbitrary positive shifts are. The decoration T,u,B,w and padding location must be retained.

## 2. Exact barrier for this offset mechanism

Both A and C are positive integers. Moreover, C itself is a tournament Hamilton-path count: it is H(T=>B). Consequently this mechanism cannot supply a negative offset or an offset7 or21, using the all-order exclusion THM-1370 (`h-spectrum-omits-7-21-all-n` theorem). This does not exclude any such numerical target from being reached by a different choice of scale or another construction.

The unit-coefficient case A=1 is completely rigid: C=1. To see this, every term in R is positive. The singleton path contributes H(T-u)>=1; if u has an out-neighbor, a second path contributes positively. Thus R=1 iff u is a sink and H(T-u)=1. A tournament has H=1 iff transitive: a directed triangle supplies at least3 through the odd-cycle formula, whereas a transitive tournament has one path. Hence T is transitive and u its sink. The converse argument gives q_e+q_o=1 iff B is transitive and w its source. Both have H=1.

Allowing a clock shift A=2^r gives more rays, but does not remove the C barrier. Concrete five-vertex seeds, with adjacency encoded by bitmasks and vertices numbered0..4, are:

- T rows [12,29,24,16,1], u=1: H=5 and R=32. With B singleton, (2) gives 2^(m+5)+5 for m>=1, on m+6 vertices.
- T rows [14,28,24,16,1], u=2: H=9 and R=16. It gives 2^(m+4)+9 for m>=1, on m+6 vertices.

These are constructions, not assertions of new scalar values. The new point is the exact transport and its family-level obstruction. Equal H does not determine A: the source and sink ports of the same transitive two-vertex tournament have R=2 and1.

## 3. One-apex zero tails: inherited transport and coefficient loss

Use THM-5070's one-apex word convention: take a chain v_1->...->v_k and an apex u; u beats v_i exactly when the i-th bit of W, read most significant first, is1. Appending a bit means adding that bit on the right of the written word, so the integer changes to2W+b. In this convention H(W)=1+2a(W), where a(0)=0,
a(2W+1)=a(W), a(2W)=a(W)+ceil(W/2). Thus

H(2W+1)=H(W), H(2W)=H(W)+W+(W mod2).                  (3)

For even W and m>=0, iteration gives

H(2^m W)=2^m W+[H(W)-W].                              (4)

This is inherited algebra, not a new theorem. Negative offsets genuinely occur: W=6 gives H=5 and hence H(6*2^m)=6*2^m-1. They are not pure dyadic leading scales. If positive even W is a power of two, then H(W)=W+1, so the normalized pure-dyadic ray has offset1 only. This explains exactly which sidecar a scalar graph analogy discards.

## 4. New all-length obstruction: a mixed periodic suffix always drifts

Fix a binary word b=(b_1,...,b_k), let q=2^k and d its integer value, and repeatedly append b to a one-apex word W_0. Assume the parity of W_0 is b_k; without this assumption discard the first block, after which it holds. Then fixed nonnegative integers alpha,beta satisfy

W_{m+1}=qW_m+d,
H_{m+1}=H_m+alpha W_m+beta.                            (5)

Specifically, start c=1,e=0,p=b_k, alpha=beta=0. At digit b_i=0 add c to alpha and e+p to beta; then replace (c,e,p) by (2c,2e+b_i,b_i). This is just (3), with the previous digit p carrying parity.

Put

A=alpha*(W_0+d/(q-1))/(q-1),
D=beta-alpha*d/(q-1).

Summing (5) proves the exact formula

H_m=A q^m + Dm + (H_0-A).                              (6)

**Theorem.** If b contains both digits, D>0. Therefore no mixed periodic suffix, of any length, transports a constant offset from a dyadic exponential clock. This holds even on every infinite arithmetic progression of repetitions. If b consists only of zeros, D=0 and (4) applies; if it consists only of ones, H is constant.

For the smallest mixed examples01 and10, respectively (q,d,alpha,beta)=(4,1,1,1) and(4,2,2,2), and both have D=2/3. The initial parity sidecar is essential: for W_0=0 and suffix01 the actual increment is0, not beta=1. Dropping the first block resolves this mismatch.

Proof of strict positivity: evaluate the affine block increment at its formal fixed point W*=-d/(q-1). Its value is D. For a mixed word, W* and every cyclically rotated fixed point lie strictly between -1 and0: each is minus a periodic binary fraction with both digits. Traverse the digits cyclically and partition the zeros into maximal runs. A run of length l is preceded by1. Write the formal value immediately before the run as -r. By (3) the first zero contributes1-r, the following zeros contribute -2r,-4r,...,-2^(l-1)r. The total is 1-(2^l-1)r. Immediately after the run the formal value is -2^l r in (-1,0), so 0<r<2^(-l), and the run contribution is strictly greater than2^(-l)>0. Ones contribute zero. At least one zero run exists, proving D>0. Formula(6) now excludes a fixed-offset exponential: its first differences divided by q^m have a nonzero term D/q^m, or directly the second-order recurrence for an exponential plus constant fails by (q-1)D. Passing to a repetition subsequence only multiplies the linear drift by its positive step.

This is a restricted all-order obstruction, not an obstruction to the full15-matrix orbit, arbitrary nonperiodic words, varying blocks, or other tournament moves. It strengthens a finite search over suffix lengths to a proof for every length within the one-apex digit grammar.

## 5. Reproduction and status

Run `python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets/replay.py`.
No repository imports or external packages. Hamilton paths are counted independently by subset endpoint DP. Path-response sums are enumerated directly. The complete rooted labelled universe of orders1..3 checks87 padding identities and1682 bridges, plus50 deterministic random larger bridges. Both concrete rays are checked at m=1..4. Hostile controls reject extrapolation to m=0 and scalar-H-only response transport. Periodic-word tests check every word through length10 with exact rational arithmetic, plus independent one-apex Hamilton-path DP on short words. All arithmetic proofs above are all-order; enumerations are only checks.
