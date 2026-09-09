---
id: THM-4464
title: "Checksum projective tournament and output sheet cocycle"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED. For a checksum response of
  even order 2q, oriented sections of its q projective directions define
  intrinsic tie-free determinant tournaments. For every r dividing q the
  labelled section action fixes no tournament if q/r is even and exactly
  2^(r-1) if q/r is odd. Its minimum period is the two-part of q, and every
  dyadic section has period q. At step q all valued pair determinants and
  projective addresses return while selected checksum verdicts reverse.
  A root sign restores a tracked input's output when that input lies in
  the transported section; an arbitrary fixed section additionally needs
  the input's relative sheet. Common GL2 frame covariance is retained.
  Separately, rank-two tie-free determinant amplitudes on at least four
  vertices cannot be absorbed into positive vertex weights. No coin
  deadline is improved and no physical LRC owner action is constructed.
source: tournament-continuation-20260908
depends_on:
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
related:
  - LEM-004-tournaments-are-odd-functions
  - THM-1475-the-pfaffian-is-the-odd-function-and-its-spectrum-has-gaps
  - THM-2235-response-antipode-barriers-for-lrc-sheets-tournament-cycles-and-knot-kernels
  - THM-2253-online-dyadic-contrast-tournament-extractor
  - THM-4461-checksum-orientation-transcendence-and-nonlinear-dyadic-carry
script: 04-computation/tournament_orientation_transport_20260908.py
output: 05-knowledge/results/tournament_orientation_transport_20260908.out
script_sha256: 451bb4c89d87aaa01094a4dcc87c25497374b05f9bc9025b6390440f8f940318
output_sha256: 66f5a9e2f6ffcf6e81b6385ab97fef35237ecef7fb71856a924598f112f29495
hash_basis: working-tree bytes; output LF
independent_audit:
  - tournament-continuation-20260908-metric-agent
  - tournament-continuation-20260908-root
---

# THM-4464 -- Checksum projective tournament and output sheet cocycle

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** A checksum response admits a real determinant tournament on its projective directions after an explicit orientation section is chosen. Its continuation is exact: if the response has `2q` phases, every dyadic `q`-vertex section tournament has period exactly `q`, while step `q` reverses every selected verdict and the physical checksum response has period `2q`. For general `q`, the minimum section period is exactly the two-part of `q`. For a tracked input lying in the transported section, a single retained root orientation restores the lost output. Separately, from four vertices onward no positive vertex weights can restore a rank-two determinant field from its tournament signs.

These are statements about actual cyclic response and its faithful or lossy representations. They do not improve a coin deadline, construct an LRC owner exchange, prove an LRC profile exclusion, or transfer the checksum to a thirteen-power sheet group.

## Inheritance and the live concept board

The inherited baseline is `c53b407b7`. The closest proved mechanisms are:

- **THM-2225**, [dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection](../../01-canon/theorems/THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection.md): rotation at tail weight `j` changes the checksum by `j`; its response coset contains the output antipode precisely when `gcd(j,m)` divides `m/2`.
- **THM-2294**, [anchored-plucker-tournament-and-kakeya-address-bank](../../01-canon/theorems/THM-2294-anchored-plucker-tournament-and-kakeya-address-bank.md): an alternating rank-two determinant field gives a locally transitive tournament when tie-free. It retains more than symmetric quadratic-character colours at thirteen, and its magnitudes and affine offsets remain sidecars.
- **LEM-004**, [tournaments-are-odd-functions](../../01-canon/theorems/LEM-004-tournaments-are-odd-functions.md), entry (a): odd circulant kernels are tournaments, whereas even cyclic order has an unavoidable antipodal tie. Only this elementary dictionary is inherited, not its historical asymptotic discussions.
- **THM-2253**, [online-dyadic-contrast-tournament-extractor](../../01-canon/theorems/THM-2253-online-dyadic-contrast-tournament-extractor.md): at the first homogeneous sibling contrast, a genuine alphabet tournament reverses under the exact exchange involution. Its exchangeability and stopping-node sidecars are essential. This is a positive example of a tournament acting on an actual extractor, different from the checksum quotient below.
- **THM-2235**, [response-antipode-barriers-for-lrc-sheets-tournament-cycles-and-knot-kernels](../../01-canon/theorems/THM-2235-response-antipode-barriers-for-lrc-sheets-tournament-cycles-and-knot-kernels.md): odd-order fixed-response groups cannot supply a Boolean antipode. A signed double cover cannot be renamed an action of its odd quotient.

The canonical hostile is the non-dyadic checksum `m=6,j=2`, with head/tail counts `7/8`. The corrected near miss is treating quadratic characters at thirteen as antisymmetric, or treating a fixed-XOR odd-group response as a local exchange. The least-used sidecar in this lane is the common orientation of all vectors after a closed projective return; pair determinants erase it even before signs are taken.

The Anchor is the requested tournament continuation of the coin operation; the Niche is a sharp amplitude-reconstruction obstruction; the Wildcard is the complete period spectrum of the orientation-section action. The five concepts are:

| Concept | Predicate and current operation | Coordinate retained or lost |
|---|---|---|
| Checksum response | Literal rotation preserves weight and changes checksum by `j` | Antipodal output sheet |
| Plucker tournament | Pair determinant determines an honest arc after choosing a section | Magnitudes and common vector orientation |
| Circulant parity | Seek a fixed tournament under the projected rotation | Two-part of the response order |
| First dyadic contrast | Exchange terminal cylinders at a selected node | Exchangeability, node address, stopping time |
| Finite continuation | Compare tournament return with output return | The signed lift of the return |

Targeted searches found the cited odd-circulant, antipode and switching mechanisms, but not the checksum section-period count below. No global novelty or literature priority claim is made. The elementary proofs, rather than a historical analogy, supply the claimed connections.

## 1. The physical response is a signed double cover

Let `m` be even and let a binary tail `y` have length `m` and weight `j`, with `1<=j<m`. Number positions `1,...,m` and set

```
s(y)=sum_i i*y_i mod m,     0<=s(y)<m.
```

Let `rho` move each position forward by one cyclically. Then

```
s(rho y)=s(y)+j mod m.                                  (1)
```

Set `g=gcd(j,m)`. Assume `g` divides `m/2`, and define

```
2q=m/g,             k=j/g,             gcd(k,2q)=1.      (2)
```

In particular `k` is odd. On one response coset choose `0<=s_0<g` and write

```
s=s_0+g*r,      r in Z/(2q).
```

The actual checksum rule is heads iff `r<q`; this uses the same half-open boundary convention as THM-2225. Its response action is `r->r+k`. The antipodal quotient is `x=r mod q`, and its extra sheet is `epsilon=floor(r/q) in {0,1}`. Explicitly,

```
x'=(x+k) mod q,
epsilon'=epsilon+floor((x+k)/q) mod 2.                 (3)
```

After `q` applications, `r` increases by `q*k`, which is `q mod 2q`. Thus the quotient address returns and the output flips. After `2q`, the response itself returns. A literal word orbit can be longer than `2q`, but its length is a multiple of `2q`, because (1) has exact response order `2q`.

This recovers the existing checksum bisection mechanism without asserting that `rho^q` is an involution on literal words. It is an antipode on their response. If `g` does not divide `m/2`, the response orbit has odd order and misses its antipode; the excluded `m=6,j=2` control lies exactly there. For dyadic `m`, every nonconstant weight satisfies (2).

## 2. The vertices and determinant observable are explicit

Represent the `q` projective directions by

```
u_x=(cos(pi*x/q), sin(pi*x/q)),        0<=x<q.
```

The offset `s_0` amounts to one common rotation and does not change pair determinants. An **orientation section** is a sign vector `sigma in {+1,-1}^q`, selecting one of the two antipodal response phases over each direction. Its chosen vectors are `v_x=sigma_x*u_x`.

The vertices are these labelled projective directions with their selected oriented lifts. Their intrinsic pair observable is

```
h_sigma(x,y)=det(v_x,v_y)
             =sigma_x*sigma_y*sin(pi*(y-x)/q).          (4)
```

For distinct `x,y`, it is nonzero. Orient `x->y` iff (4) is positive. This is the real rank-two determinant construction of THM-2294, and every resulting tournament is locally transitive. No quadratic character or label ordering is substituted for the observable.

On the full `2q` response phases, antipodal pairs have determinant zero. They remain actual ties. Passing to projective directions and then choosing a section is an explicitly declared quotient and lift; it is not a tie-breaking rule claimed to be canonical. Global replacement `sigma->-sigma` leaves the complete valued field (4), and hence the tournament, unchanged.

For `x<y`, the sine in (4) is positive, so

```
S_sigma(x,y)=sigma_x*sigma_y.                            (5)
```

Consequently two sections give the same labelled tournament exactly when they differ by one common sign. There are `2^(q-1)` section tournaments. The statement remains valid for `q=1`, with its one trivial tournament.

## 3. Exact signed rotation and the complete fixed-period count

Write

```
pi(x)=(x+k) mod q,
tau_x=(-1)^floor((x+k)/q).
```

Rotation through angle `pi*k/q` sends `u_x` to `tau_x*u_pi(x)`. Thus the physical action on sections is the signed permutation

```
(T sigma)_pi(x)=tau_x*sigma_x.                          (6)
```

It obeys

```
T^q sigma=-sigma,                                      (7)
h_(T sigma)(pi(x),pi(y))=h_sigma(x,y).                  (8)
```

Equation (8) preserves the full determinant amplitudes, not just the tournament signs. Relative to a fixed canonical section, the same formula appears as a cyclic relabelling together with vertex switching by `tau`. On sections modulo their common sign, `T` has order dividing `q`.

This is not a diagonal-gauge assumption. Under any common real invertible frame change `A`, determinants multiply by `det(A)` and the rotation conjugates to `A*Rot(pi*k/q)*A^(-1)`; the signed section action (6) is unchanged. A negative determinant globally reverses the tournament. The physical heads half-plane and its half-open boundary mark must also be transported if the frame is changed. The exact `q=2` shear `A=[[1,1],[0,1]]` gives conjugated quarter rotation `[[1,-2],[1,-1]]`, whose square is `-I`; it sends `(1,0)` to `(1,1)` and then to `(-1,0)`. Thus the negative return persists in a nonorthogonal, nondiagonal frame.

**Theorem A (fixed sections).** For every positive divisor `r` of `q`, the number of section tournaments fixed by `T^r` is

```
Fix(T^r)=0                 if q/r is even,
         2^(r-1)           if q/r is odd.              (9)
```

**Proof.** By (5), being fixed as a labelled tournament is equivalent to `T^r sigma=lambda*sigma` for one common `lambda in {+1,-1}`. The permutation `pi^r` has `r` cycles, each of length `q/r`. By (7), the sign product around each such cycle is `-1`. Therefore consistency requires `lambda^(q/r)=-1`. This is impossible when `q/r` is even. When `q/r` is odd it forces `lambda=-1`; one initial sign may then be chosen freely on each of the `r` cycles. There are `2^r` sections, and quotienting by their common sign leaves `2^(r-1)` tournaments. QED.

**Corollaries.**

1. Write `q=2^a b`, `b` odd. The minimum possible section period is exactly `2^a`. It is attained because (9) has `2^(2^a-1)` fixed sections at that divisor; no smaller divisor can support a return.
2. If `q` is a power of two, **every section tournament has period exactly `q`**. There are exactly `2^(q-1)/q` labelled section orbits.
3. A section tournament fixed by one rotation exists iff `q` is odd, and then it is unique. It is the regular cyclic tournament, in the alternating section `sigma_x=(-1)^x`.

For the last claim, direct substitution gives `T sigma=-sigma` when `q` is odd. Its arcs from zero go to the positive even residues `2,4,...,q-1`. Every vertex has outdegree `(q-1)/2`; multiplication by `2` identifies it with the usual cyclic half-interval tournament. This is a specific instance of LEM-004's odd-function dictionary.

The dyadic failure is stronger than the absence of a fixed tournament: every section is forced to run through the full `q`-step cycle. For example, `q=4,8,16` give respectively `2,16,2048` section orbits, each of lengths `4,8,16`. The statement is about labelled response sections under the specified rotation, not about isomorphism classes of tournaments or all finite-state extractors.

At odd `q`, the physical signed action still has the missing two-sheet information. The regular tournament does not create a Boolean antipode in an odd-order group: the output flip lives in the physical `2q`-cycle. Thus the positive odd-carousel case respects THM-2235.

## 4. A closed tournament continuation can reverse the actual fair bit

**Theorem B (one necessary and sufficient output bit).** The actual checksum verdict does not factor through the pair determinant field and the projective checksum address. One root orientation bit, together with the labelled section tournament and that address, restores it exactly.

**Proof.** By (7), after `q` rotations every selected vector is negated. All pair determinants and all projective addresses are unchanged, while each chosen phase crosses to its antipode and its checksum verdict reverses. Hence even the complete valued edge field cannot determine the verdict.

Conversely, retain `eta=sigma_0`. Equation (5) recovers

```
sigma_0=eta,
sigma_x=eta*S_sigma(0,x),                 1<=x<q.      (10)
```

For an input on the selected lift at address `x`, its checksum verdict is heads iff `sigma_x=+1`. More generally the actual input's sheet relative to the section must also be specified; changing that relative sheet flips the decoded verdict. Equivalently, choose the transported section to contain the tracked input, and (10) gives the required single global bit. QED.

This typing matters: a tournament on projective directions does not identify which of the two physical inputs over a direction was observed. The one-bit claim is exactly recovery of a section from its edge field, with the tracked input placed in that section. It is not a claim that an arbitrary fixed section and no input-sheet information determine every response.

There is a nontransitive literal control. Take `m=8`, weight `j=1`, so `q=4`, and the section

```
sigma=(+1,-1,+1,-1).
```

Its determinant tournament is strongly connected, with two directed triangles. The tail `00000001` has checksum `0` and gives heads. Four forward rotations give `00010000`, checksum `4`, and tails. The transported section is now `-sigma`; its entire determinant field and the input's projective address are identical. The boundary convention is the literal half-open checksum convention, not an arbitrarily assigned sign at a tied anchor.

Thus a closed tournament continuation is not automatically a closed output continuation. The missing sign is exposed in a strong tournament as well as in the transitive section.

## 5. Four vertices already force genuinely edgewise amplitudes

The same determinant carrier has a second exact loss boundary relevant to weighted tournament transport.

**Theorem C (sharp vertex-amplitude obstruction).** Let `h_ij=det(v_i,v_j)` be a tie-free rank-two real edge field on at least four vertices, and `S_ij=sign(h_ij)`. There are no positive vertex weights `lambda_i` with

```
h_ij=lambda_i*lambda_j*S_ij.                            (11)
```

The boundary is sharp: every positive edge-amplitude assignment on at most three vertices has such a vertex factorization.

**Proof.** On any four vertices, the Plucker identity gives

```
h_12*h_34-h_13*h_24+h_14*h_23=0.                       (12)
```

Under (11), the left side equals the positive factor `lambda_1*lambda_2*lambda_3*lambda_4` times

```
S_12*S_34-S_13*S_24+S_14*S_23.
```

This is a sum of three numbers in `{+1,-1}` and is odd, hence nonzero. Contradiction. This combines THM-2294's decomposable two-form with the four-vertex instance of **THM-1475**, [the-pfaffian-is-the-odd-function-and-its-spectrum-has-gaps](../../01-canon/theorems/THM-1475-the-pfaffian-is-the-odd-function-and-its-spectrum-has-gaps.md). For three vertices with positive amplitudes `a_12,a_13,a_23`, set `lambda_1=sqrt(a_12*a_13/a_23)` and recover the other two by division; one and two vertices are immediate. QED.

The integer control `v_1=(1,0)`, `v_2=(1,1)`, `v_3=(0,1)`, `v_4=(-1,1)` gives six positive determinants

```
(h_12,h_13,h_14,h_23,h_24,h_34)=(1,1,1,1,2,1).
```

Their Plucker expression is `1-2+1=0`, while their all-positive sign tournament gives `1-1+1=1`. Thus a node-weighted tournament cannot reconstruct even this smallest valued determinant packet. Ties, zero weights, or non-rank-two fields lie outside the theorem.

No quantitative LRC capacity follows from unweighted tournament degree or strength alone. In particular, if an actual positive relation `sum_i w_i v_i=0` is retained, the weighted edge flow `w_i*w_j*h_ij` has zero divergence, since its row sum is `w_i det(v_i,sum_j w_jv_j)=0`. The signs orient this flow, but (11) proves that its determinant amplitudes cannot in general be absorbed into vertex masses. An edge-amplitude sidecar is mathematically necessary from four vertices onward.

## 6. Connection contracts and stopping boundaries

**Coin response -> determinant tournament continuation.** The source is the literal checksum response (1) on a fixed weight and coset. The map is the antipodal quotient, an explicitly chosen orientation section, and pair determinant (4). It preserves rotation as signed permutation (6), full edge amplitudes through (8), and the exact period census (9). It loses the common orientation of all lifts and the input sheet. A root sign and the tracked input's relation to the section restore the checksum verdict. The cheapest decisive test is the eight-bit hostile in Section 4; exhaustive controls also replay every response-compatible word through length 12.

**Odd circulant dictionary -> the dyadic obstruction.** The source is the signed rotation of the projective frame; the target is a tournament invariant under its projected cyclic action. The map takes fixed sections modulo global sign. It preserves the actual angular pair orientation, not just an assigned odd function. Its existence is exactly the oddness condition in (9) at `r=1`. The source has a physical double cover; discarding it loses the output antipode. The decisive test is the fixed-section count, including mixed `q=6,10,12` controls, rather than only powers of two.

**Plucker field -> weighted tournament.** The map takes signs. It preserves an intrinsic pair orientation and local transitivity but loses amplitudes. Equation (12) and the odd sign Pfaffian prove that vertex masses cannot repair this loss. The needed sidecar is genuinely edgewise magnitude, together with the affine offsets required by THM-2294's original rank-three packet. The explicit four integer vectors are the minimal size control.

**First contrast and LRC remain separate.** THM-2253 already gives a lawful tournament on actual categorical source symbols because the selected sibling exchange preserves the terminal cylinder's probability. The checksum's projective section does not replace that exchangeability mechanism. No demonstrated LRC moving-density or owner-word action intertwines with (6); its thirteen-power fixed-response obstruction remains live. No physical LRC predicate is asserted to descend to these section tournaments.

After these tests the board changes as follows: the coin lane gains exact tournament continuation with a stated output sidecar; the Plucker lane gains the sharp four-vertex amplitude obstruction; the odd-circulant lane gains the full two-part period law; the contrast lane remains the positive control for an actual exchangeable-source tournament; and finite continuation now distinguishes a closed pair observer from an open output sheet. The stopping boundary is precise rather than a claim that tournaments are universally unsuitable.

## 7. Exact evidence

Run from the repository root:

```
python 04-computation/tournament_orientation_transport_20260908.py
python -O 04-computation/tournament_orientation_transport_20260908.py
```

The standard-library script imports no repository producer. It enumerates all **5,448** nonconstant tails at every even length `2..12`; **3,918** satisfy the exact antipode condition. It enumerates **36,863** labelled orientation sections modulo their common sign for every `q=1..12` and `q=16`, including all section periods and fixed counts. For `q<=8` every unit response step modulo `2q` is also checked. Determinant signs use exact residue order, not floating trigonometry. Independent literal tail rotations check the physical response, and direct sign-field iteration checks the section action. All **64** four-vertex sign fields and the integer Plucker hostile are retained. The source prints the actual head/tail return witness, not only a graph invariant.

All **223,007 always-active exact gates** pass. The all-parameter fixed-period, one-bit and amplitude statements are proved above; the finite computation is an adversarial check. A next useful question is whether a particular physical continuation consumer admits an explicit edge-amplitude and sheet-preserving map. The present result supplies the exact test such a proposed map must pass.

Both ordinary and optimized Python runs pass with the same transcript. Working-tree SHA256: source `451bb4c89d87aaa01094a4dcc87c25497374b05f9bc9025b6390440f8f940318`; LF output `66f5a9e2f6ffcf6e81b6385ab97fef35237ecef7fb71856a924598f112f29495`.

The metric-continuation agent independently accepted the fixed-section proof, `q=1,2` boundaries, minimum two-part period, alternating odd-order section for every admissible step, the precisely scoped tracked-input decoder, common `GL2` frame covariance and the explicit shear, and the four-vertex amplitude obstruction with its three-vertex sharp boundary. No mathematical correction was required.

The root tournament-continuation agent independently read and accepted the full draft as well. Promotion preserves the labelled section-action scope, the distinction between response return and Boolean verdict repetition, the tracked-input condition on the one-bit decoder, and the edge-amplitude boundary. No deadline improvement or physical LRC map is asserted.
