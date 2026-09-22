# The greedy 3-adic inverse map G: symbolic conjugacy, exact moments, G_b census, hostiles, and the Collatz reframe

**Status: PROVED (scoped) + FINITE-EXACT + CITED classical ingredients; no novelty claim for the
symbolic-dynamics and Markov-chain machinery, nor for the level-2 residue chain, its stationary law and
`E_pi[k]=1`, which are wave-one Theorem 4.2. Audited adversarially on 2026-09-21 (audit script and output
named in the reproduction block); four draft statements were refuted and corrected, see the near-miss list.**
PROVED: the level-J transfer law and Markov structure of
the greedy map `G(m)=(2^k m-1)/3` on the 3-adic units (level 2 inherited, levels `J>=3` new); topological conjugacy of `(Z_3^x, G)` to a
6-state subshift of finite type with entropy `log 3`, the k-word being an injective 2-block recoding;
`G`-invariance of `mu=(2/3)Haar(1+3Z_3)+(1/3)Haar(2+3Z_3)`, primitivity at every level, stationary drift
exactly `log(2/3)`; the exact exponential moment `E_Haar[u^{K_J}]=(u+1)(u^2+2)/6*((u^2+u+1)/3)^(J-1)` for
EVERY tilt `u` (the wave-one tilted matrix has determinant zero for every `u`), the exact
large-deviation rate function, and the fact that the wave-one `(7/9)^(J-1)` is the `u=2` Chernoff bound
and not the sharp exponent (`0.758751`); `sigma=sigma_res` for every `m>=2` whose greedy stopping time
is `<=12` (exact per-level threshold) and, from the same enumeration, no positive `G`-cycle of length
`<=12` other than `{1}` (S3.9); `E[j]=1` for every signed `G_b` (`3` not dividing `b`), the
`G_b` cycle gate `m_0(2^J-3^L)=bB'(w)` with `B'(w)=B(reversed w)`, the scaling lemma `G_b(bm)=bG_1(m)`;
the hostile families `3^j+1` (lands on `4^(j-1)`), the sharp peak bound `K_n<=2n+[k_1=3]`, and the
impossibility of two consecutive `8/3`-letters. FINITE-EXACT: all `m<=10^7` reach 1 (max 93 steps, attained at
`m=8751065` and `m=9454814`; worst `peak/m=133.033` at `m=4847486`, full word printed); `d_J` exact for `J<=12`;
`G_b`/`T_b` cycle census for all odd `|b|<=49`, `3` not dividing `b`, `|m|<=10^5`, no escape to `10^18`.
REFUTED: the reversal correspondence as a bijection between `G_b`-cycles and `T_b`- or `T_{-b}`-cycles
(only the two fixed points are common; census counts `#G_1=3` vs `#T_1=4`); the draft's "exactly one class
mod `3^(n+3)`" for the `(8,5,1^(n-1))` family (it is two classes, witness `8,35 mod 81`); the draft's
"primitive `= q=|b|`" reading of the census column (it counts content `d<|b|`; `b=25,35,49` differ); and the
draft's stopping-boundary conjecture `m*(w) <` least member of the class (witness word `(1,2,2)`, class
`2 mod 81`, `m*=37/5`). SCOPE: `Collatz => Q1` is trivial, `Q1 => Collatz` is not established here; no map
found from the census to Bott/octonion structure or from `6/pi^2` to any `G`-quantity; an abstract bijection
between the full (unknown) cycle sets is not a claim. HEURISTIC: the observed prefix-tail ratio `0.6358`
below `exp(-I(c))=0.7588` is a finite-`J`/lattice effect. OPEN: no positive `G`-cycle of length `>12` and no
divergent positive `G`-orbit (equivalently Q2), `sigma=sigma_res` for all `m`, and whether `exp(-I(c))` is
the exact decay rate of `1-d_J`.

Inheritance and concept board. Closest proved mechanism: the wave-one E-graph lane
([extended_collatz_scc](collatz_mod6_20260917_extended_collatz_scc.md), S2/S4/S6/S7): the first `J`
greedy letters depend only on `m mod 3^(J+1)` (sharp), the exact residue chain mod 9 with stationary law
`2/9` on `{1,4,7}`, `1/9` on `{2,5,8}` and `E_pi[k]=1` (its Theorem 4.2; Corollary 1.2 and Theorem 2.2 below
are its level-`J` generalization), `E[2^{K_J}]=3(7/3)^(J-1)` via the rank-one
tilted matrix `[[5/3,1/3],[10/3,2/3]]`, the density-1 stopping theorem, the hostile word of `3^j+1`,
the leaf identity of the `even -> 3n+1` arrows, and the `G_-` two-cycle `{4,11}`; none of these is
re-proved here. The cycle-gate polynomial `n_0=bB/(2^K-3^L)` and the minimal parameter
`q=|Delta|/gcd(B,|Delta|)` are inherited from
[arithmetic_braids_20260917_collatz](arithmetic_braids_20260917_collatz.md) (B7) and
[arithmetic_braids2_20260917_signed_cycles](arithmetic_braids2_20260917_signed_cycles.md) (8)-(10);
the affine word/carry model `W(n)=(2^r n-B)/3^m` from
[collatz_blueprint_20260921_affine](collatz_blueprint_20260921_affine.md) Sec. 1 and 3. Canonical
hostile: the worst `m=4847486`, whose word prefix `3 2^5 0 3 2^9 0` stacks two `(5,1^n)` runs and
reaches `2^34/3^17=133.03` of `m` at step 17 (the peak is a 3-adic accident: `v_3(G^8(m)-1)=10`).
Corrected near miss: the recovered draft asserted the bijective-lift law at `J=1` and a sharpness
witness of the form `a, a+3` -- both false at `J=1` (the three lifts of a class mod 3 have `k`
`(2,0,0)` resp. `(1,3,1)`, images `(1,1,2)` mod 3; the witness is `G(1)=1`, `G(7)=2` mod 3); it also
claimed the k-word map is not injective (it is, S1.3), and quoted `d_12=0.9997` (true value
`530885/531441=0.998954`). The 2026-09-21 audit refuted four more draft statements, each now corrected
in the script and asserted by an explicit check: the pair-resolving clause of Theorem 1.3(ii) (draft:
`k_(n+1) in {0,2} iff X_n in {1,2,4,5}`; truth: `iff X_n in {1,4,7}`, witness `X_n=7`); the class count of
the `(8,5,1^(n-1))` family (draft: one class mod `3^(n+3)`; truth: two); the "primitive `= q=|b|`" label
of the census column (truth: content `d<|b|`, which splits for `b=25,35,49`); and the stopping-boundary
conjecture (refuted by the note's own candidate list). It also found the tie `9454814` for the maximal
93 steps and downgraded the Lind-Marcus theorem numbers to UNCITED-RECOLLECTION. Least-used sidecar: the exact threshold `m*(w)=B_i/(2^{K_i}-3^i)` of S3.7, the
rational fixed point of the word's affine map, which is simultaneously the cycle gate and the boundary
between "descent" and "residue says growth".

## 1. The map, the transfer law, and the symbolic conjugacy

`G(m)=(2^k m-1)/3` with `k>=0` minimal such that `2^k m` is `4` or `7` mod 9; `k` depends only on
`m mod 9`: `k(1,2,4,5,7,8)=(2,1,0,3,0,1)`. Since `2^k m` is never `1` mod 9 the image is a 3-adic unit,
and `3G(m)+1=2^k m` says that a `G`-orbit reversed is an E-path (one `3n+1` arrow and `k` halvings per
step; the letter `k=0` is the E-graph's extra `even -> 3n+1` arrow, wave-one S1.1). `G(1)=1`.

**Theorem 1.1 (transfer law; PROVED, S1.1 of the output).** For `J>=2` and a unit class `a mod 3^J`,
its three lifts `a+t3^J` (`t=0,1,2`) mod `3^(J+1)` share the same `k` and
`G(a+t3^J) = G(a)+t2^k 3^(J-1)` mod `3^J`; since `2^k` is a unit mod 3, the three lifts map bijectively
onto the three lifts (mod `3^J`) of the class `G(a) mod 3^(J-1)`. Hence `G(m) mod 3^J` is a function of
`m mod 3^(J+1)`, and `G` maps each cylinder `a+3^(J+1)Z_3` affinely and bijectively onto the cylinder
`G(a)+3^J Z_3` (3-adic scale factor 3). At `J=1` the statement is different and is the truth asserted by
the script: the three lifts mod 9 of either class mod 3 have `k=(2,0,0)` resp. `(1,3,1)` and images
`(1,1,2)` mod 3, i.e. `P_1=[[2/3,1/3],[2/3,1/3]]`. Verified for every unit class at `J=1..7`; the modulus
`3^(J+1)` is sharp with witnesses `(a, a', G mod 3^J) = (1,7 -> 1,2)`, `(1,10 -> 1,4)`, `(1,28 -> 1,10)`,
`(1,82 -> 1,28)` at `J=1..4`.

**Corollary 1.2 (Markov structure; PROVED, S1.1c).** For `m` Haar-distributed on `Z_3^x` and `J>=2`, the
residue sequence `X_n=G^n(m) mod 3^J` is a Markov chain with kernel `P_J(a,b)=1/3` iff
`b = G(a) mod 3^(J-1)`; at `J=1` the kernel is `P_1=[[2/3,1/3],[2/3,1/3]]` (Theorem 1.1). Level 2 is
wave-one Theorem 4.2; levels `J>=3` are new here.
Proof: the history event `{X_0=a_0,...,X_n=a_n}` is a disjoint union of classes mod `3^(J+n)`, each
mapped by `G^n` affinely onto the full class `a_n+3^J Z_3` (Theorem 1.1 iterated), so conditionally on
the history the next digit is uniform. At `J=2` this is the session lead's law: from `{1,2,4,5}` the next
residue is uniform on `{1,4,7}`, from `{7,8}` uniform on `{2,5,8}`.

**Theorem 1.3 (conjugacy and entropy; PROVED, S1.3-S1.4).** Let `A` be the `6x6` 0/1 matrix with
`A(a,b)=1` iff `b = G(a) mod 3`, i.e.
`rows {1,2,4,5} -> {1,4,7}`, `rows {7,8} -> {2,5,8}`, and `Sigma_A` the one-sided SFT.
(i) `Phi(m)=(G^n(m) mod 9)_n` is a homeomorphism `Z_3^x -> Sigma_A` with `Phi o G = shift o Phi`.
Injectivity: by induction the block `(X_0..X_(J-1))` determines `m mod 3^(J+1)` (`J=1` trivial; knowing
`m mod 3^(J+1)` and `G(m) mod 3^(J+1)`, the level-`(J+1)` law pins `m mod 3^(J+2)` because the three
lifts have three distinct images). Surjectivity: admissible blocks of length `J` number
`6*3^(J-1)=2*3^J`, the number of unit classes mod `3^(J+1)`, so the injective block map is bijective at
every length. Verified for `J<=7`: blocks `6,18,54,162,486,1458,4374`.
(ii) The k-word is a 2-block recoding of `Phi`: `k_n` determines `X_(n-1)` except the pairs `{4,7}`
(`k=0`) and `{2,8}` (`k=1`), and `k_(n+1) in {0,2}` iff `X_n in {1,4,7}` iff `X_(n-1) in {1,2,4,5}` resolves
the pair (`4`, `2` are the heavy members, `7`, `8` the light ones; verified on all units mod 81, S1.3 of
the output; the draft's clause "iff `X_n in {1,2,4,5}`" is REFUTED by `X_n=7`, which has `k=0`). So the
k-word map is injective (a homeomorphism onto a proper SFT of `{0,1,2,3}^N`); admissible k-words of
length `J` number `4*3^(J-1)` (`4,12,36,108,324,972,2916` for `J<=7`), a `(J+1)`-letter word pins
`m mod 3^(J+1)` exactly. This corrects the draft's claim that the k-word map is not injective; what is
lost against Lagarias' parity map (CITED: Lagarias 1985, `Q: Z_2 -> Z_2` bijective, Haar-preserving) is
surjectivity onto the full shift, not injectivity.
(iii) `A` has constant row sum 3, so `h_top(Sigma_A)=log 3=1.098612`; its Parry measure has kernel
`A/3=P_2` and stationary law the left Perron vector `(2,1,2,1,2,1)/9=pi_2` (checked exactly). The same
holds at every level (`P_J=A_J/3`), so `mu` (Section 2) is the measure of maximal entropy and
`h_mu(G)=h_top(G)=log 3` (CITED for the Parry-measure facts: W. Parry, Intrinsic Markov chains, Trans.
Amer. Math. Soc. 112 (1964); the Lind-Marcus theorem and section numbers quoted by the draft, Thm 4.4.4
and Sec. 13.3, are UNCITED-RECOLLECTION and are not relied on).

## 2. Invariant measure, primitivity, drift

**Theorem 2.1 (PROVED, S2.1).** With `w(1)=4/3`, `w(2)=2/3` and `mu = w(m mod 3)*Haar` on `Z_3^x`
(each unit class mod `3^(J+1)` of Haar mass `1/(2*3^J)`), `mu` is `G`-invariant: for a cylinder
`B=b+3^J Z_3`, the preimage classes `a mod 3^(J+1)` are `a=2^(-k)(3b+1)` for the `k` with `k(a)=k`; for
`b=1 mod 3` all four candidates `4,2,1,5 mod 9` (`k=0,1,2,3`) are valid, `mu(G^-1 B)=4/(2*3^J)=mu(B)`;
for `b=2 mod 3` only `7,8` (`k=0,1`) are valid, `mu(G^-1 B)=2/(2*3^J)=mu(B)`. Equivalently `w` is the
left Perron vector of the branch matrix `[[2,1],[2,1]]`. Verified on all cylinders of level `J<=7`.
`G` is 4-to-1 over `1+3Z_3` and 2-to-1 over `2+3Z_3`. This is the session lead's
`mu=(2/3)Haar(1+3Z_3)+(1/3)Haar(2+3Z_3)`.

**Theorem 2.2 (PROVED, S2.3).** `pi_J(b)=w(b mod 3)/(2*3^(J-1))` is `P_J`-stationary; at `J=2` it is
`2/9` on `1,4,7` and `1/9` on `2,5,8`. For `n<J` the support of `P_J^n(a,.)` is exactly the `3^n` classes
over `G^n(a) mod 3^(J-n)`, and for `n>=J` every entry is positive: `P_J` is primitive (irreducible and
aperiodic) at every level; verified `J<=7`. Drift: `E_pi[k]=2(2/9)+1(1/9)+3(1/9)+1(1/9)=1`, so the
stationary log-growth per step is exactly `log(2/3)=-0.405465`; the uniform-residue average `7/6` would
give the wrong `-0.2899`. On E-arrows: at stationarity one tripling is undone per halving (forward Collatz
uses two halvings per tripling).

## 3. Stopping time, exact densities, exact moments, sharp exponent, threshold

Notation: `G^i(m)=(2^{K_i} m-B_i)/3^i` with `B_i=sum_{s<=i} 3^(s-1) 2^(K_i-K_s)>0`;
`sigma(m)=min{i: G^i(m)<m}`, `sigma_res(m)=min{i: 2^{K_i}<3^i}`.

**Theorem 3.1 (PROVED, S3.1; the density theorem is inherited).** `sigma<=sigma_res`; a mismatch at level
`i` forces `m<B_i`; `{sigma_res<=J}` is a union of unit classes mod `3^(J+1)`, so
`d_J := dens{sigma<=J} = dens{sigma_res<=J} = #good classes/(2*3^J)`, and `d_J -> 1` by the wave-one
Markov bound `1-d_J <= (7/9)^(J-1)` (S4.5 there) or by the ergodic theorem for the primitive chain
(CITED: Norris, Markov Chains, Thm 1.10.2).

Exact values (S3.2), `J<=10` equal to the wave-one `f_J`, `J=11,12` new:

| J | d_J | 1-d_J | g_J = P(2^{K_J}<3^J) |
|---|---|---|---|
| 1 | 2/3 | 1/3 | 2/3 |
| 2 | 8/9 | 1/9 | 5/6 |
| 3 | 25/27 | 2/27 | 43/54 |
| 4 | 26/27 | 1/27 | 73/81 |
| 5 | 236/243 | 7/243 | 214/243 |
| 6 | 239/243 | 4/243 | 686/729 |
| 7 | 241/243 | 2/243 | 2126/2187 |
| 8 | 2173/2187 | 14/2187 | 12647/13122 |
| 9 | 19609/19683 | 74/19683 | 19339/19683 |
| 10 | 58868/59049 | 181/59049 | 6413/6561 |
| 11 | 176809/177147 | 338/177147 | 58396/59049 |
| 12 | 530885/531441 | 556/531441 | 1057277/1062882 |

FINITE-EXACT (S3.3): for `m<=10^6`, `sigma=sigma_res` for every `m>1` (no mismatch), `max sigma = 31`
at `m=128669`. Terras densities (CITED: Terras 1976; values computed mod `2^J`, S3.4) reach only
`0.9448` at `J=12` and `0.9739` at `J=20`, against `d_7=0.991770`, `d_12=0.998954`: the per-step drifts
are `log(2/3)=-0.405` versus `(1/2)log(3/4)=-0.144`, and the 3-adic carry helps descent while the
2-adic carry hinders it.

**Theorem 3.2 (exact exponential moments for every tilt; PROVED, S3.5-S3.6).** Lump the level-2 chain to
the class `c_n=G^n(m) mod 3`. Given `c_n`, `X_n` is uniform on its class, and the pairs
`(k(X_n), c_(n+1))` are `(2,1),(0,1),(0,2)` for class 1 and `(1,1),(3,1),(1,2)` for class 2 -- the same
multiset shifted by one halving. Hence for `u=e^theta`
`M_u=(1/3)[[u^2+1,1],[u+u^3,u]] = v w^T`, `v=(1,u)`, `w=((u^2+1)/3,1/3)`, of rank one for every `u`,
eigenvalue `rho(u)=(u^2+u+1)/3`, and with `X_0` Haar
`E[u^{K_J}] = (u+1)(u^2+2)/6 * ((u^2+u+1)/3)^(J-1)` exactly for all `J>=1`.
Verified by enumeration over units mod `3^(J+1)` at `J in {1,2,3,5,8}`: `u=2: 3(7/3)^(J-1)` (wave one),
`u=3: (22/3)(13/3)^(J-1)`, `u=5: 27(31/3)^(J-1)`, `u=7: 68*19^(J-1)`. Consequences: (a) the scaled cumulant
generating function is `Lambda(theta)=log((u^2+u+1)/3)` exactly, so `K_J` satisfies the LDP with rate
`I(c)=sup_theta(theta c-Lambda(theta))` (CITED: Dembo-Zeitouni, Large Deviations Techniques and
Applications, Sec. 3.1 and Thm 2.3.6). (b) `1-g_J=P(K_J>=cJ)` with `c=log_2 3=1.584963`; the wave-one
`(7/9)^(J-1)` is Chernoff at `u=2` (per-step `rho(2)/3=7/9`), valid but not sharp: `Lambda'(log 2)=10/7=1.4286<c`.
The optimizer is `u*=2.782079` (root of `(2-c)u^2+(1-c)u-c=0`), `rho(u*)=3.840680`, sharp per-step
factor `exp(-I(c))=rho(u*)/u*^c=0.758751` versus `7/9=0.777778`. (c) `1-d_J<=1-g_J`, so the prefix tail
decays at least at this rate; the observed ratios `(1-d_(J+1))/(1-d_J)` for `J=1..11` are
`0.333,0.667,0.500,0.778,0.571,0.500,0.778,0.587,0.815,0.622,0.548` (geometric mean `0.6358` over
`J=8..12`). PROVED: `limsup (1/J) log(1-d_J) <= -I(c)`; OPEN: whether `exp(-I(c))=0.7588` is the exact
rate of `1-d_J`; HEURISTIC: the observed `0.6358` sits below it through finite-`J` and lattice effects
(the ratios of `1-g_J` alternate between about `0.48` and about `1.25` with the fractional part of
`J log_2 3`). So the precise
statement of the task's item: the tilted matrix gives the EXACT moment of the word weight for every
tilt, and the exact LD rate; the specific bound `(7/9)^(J-1)` is not the exact exponent.

**Theorem 3.3 (exact threshold; PROVED per level, S3.7-S3.8).** On a class `a mod 3^(i+1)` with word `w`
and `2^{K_i}>3^i`, `G^i` is affine and `G^i(m)<m` iff `m<m*(w):=B_i/(2^{K_i}-3^i)`. A mismatch
`sigma(m)=i<sigma_res(m)` forces `m<m*(w)` and no earlier descent (`G^l(m)>=m`, `l<i`); conversely a
candidate with both properties is a mismatch, because no earlier descent excludes `2^{K_l}<3^l` for `l<i`
by Theorem 3.1(a). Both are finite checks. Enumerating all words of length `i<=12`:

| i | classes with 2^K>3^i | max m*(w) | at (a,K,B) | candidates m>=2 below m* | true mismatches |
|---|---|---|---|---|---|
| 1 | 2 | 1 | (1,2,1) | 0 | 0 |
| 2 | 3 | 11/7 | (8,4,11) | 0 | 0 |
| 3 | 11 | 53/5 | (43,5,53) | 1 | 0 |
| 4 | 16 | 239/47 | (124,7,239) | 1 | 0 |
| 5 | 58 | 973/13 | (187,8,973) | 3 | 0 |
| 6 | 86 | 827/59 | (1645,10,4135) | 2 | 0 |
| 7 | 122 | 17269/1909 | (1645,12,17269) | 2 | 0 |
| 8 | 475 | 63071/1631 | (9850,13,63071) | 4 | 0 |
| 9 | 688 | 51769/2617 | (29533,15,258845) | 4 | 0 |
| 10 | 2664 | 940375/6487 | (44302,16,940375) | 8 | 0 |
| 11 | 3918 | 3820549/84997 | (398596,18,3820549) | 8 | 0 |
| 12 | 5605 | 15459343/517135 | (398596,20,15459343) | 8 | 0 |

Every candidate (e.g. `(i,m)=(3,2),(5,4),(5,14),(8,5),(8,7),(10,13)`) descended at an earlier step. Hence
`sigma(m)=sigma_res(m)` for every `m>=2` with `sigma(m)<=12`, and `d_J` is the exact density of
`{sigma<=J}` with no exceptional set for `J<=12`; the class of `m=1` (word `2^i`) has `m*=1` exactly.
OPEN: all levels (FINITE-EXACT to `10^6`, where `sigma_res<=31`).

**Corollary 3.4 (cycle boundary; PROVED + FINITE-EXACT, S3.9).** Every element `m_0` of a positive
`G`-cycle of length `L` is the rational fixed point `m*(w)` of its own word `w` (a growth word, since
`m_0(2^{K_L}-3^L)=B_L>0`), so positive cycles of length `L` are exactly the integers `m*(w)` lying in
their own class `a mod 3^(L+1)`. The enumeration finds, for every `L<=12`, only `(L, m*)=(L,1)`: no
positive `G`-cycle of length `<=12` other than `{1}`. (This is stronger than the `10^7` sweep for these
lengths: `m*(w)<=145` bounds every candidate cycle element.) Note the two regimes: candidates `m<m*(w)`
exist at every level `i>=3` (they descend earlier); integer boundary points `m=m*(w)` do not.

## 4. Signed maps G_b: drift, cycle gate, scaling lemma, census, duality

`G_b(m)=(2^j m-b)/3`, `j>=0` minimal with `2^j m in {b+3,b+6} mod 9`, on all integers coprime to 3.

**Theorem 4.1 (PROVED, S4.1).** For every `b` with `3` not dividing `b` the transfer law holds verbatim,
the invariant density is `4/3` on `m = b mod 3` and `2/3` on `m = -b mod 3`, and `E_pi[j]=1`, so the drift
is `log(2/3)` for every `b`. Proof: with `T_1=b+3`, `T_2=b+6` and `u` the discrete log of `T_1/m` base 2
mod 9, `T_2/T_1` is `2^2` (`b=1 mod 3`) or `2^4` (`b=2 mod 3`); the case `b=1 mod 3` gives exactly the
`b=1` table under the relabelling `u=0..5 <-> 4,2,1,5,7,8`, and `G_{-b}(-m)=-G_b(m)` conjugates the other
case. All six `b mod 9` verified: `j`-multiset `{0,0,1,1,2,3}`, `E_pi[j]=1`.

**Theorem 4.2 (cycle gate; PROVED, S4.3).** For a `G_b`-cycle `m_0 -> ... -> m_L=m_0` with word
`(j_1..j_L)`, `J=sum j_i`: `m_0(2^J-3^L) = b B'`, `B'=sum_{i=1}^L 3^(i-1) 2^(J-J_i) = B(j_L,...,j_1)`,
where `B` is the inherited `T_b` gate polynomial (`n_0(2^K-3^L)=bB(k_1..k_L)`, braids (B7)). The two gates
coincide up to word reversal. Consequences: (i) `sign(m_0)=sign(b)sign(2^J-3^L)`, so drift-side words
(`2^J<3^L`) give cycles of sign `-sign(b)`; (ii) a `G_b`-cycle is a reversed `T_b`-cycle iff all its
elements are odd, and a `T_b`-cycle is a reversed `G_b`-cycle iff every exponent is the greedy minimum
at its node (both verified on every cycle found); (iii) `G_{-b}(-m)=-G_b(m)`. With the 2-adic gate: the
roles of 2 and 3 are exchanged in the carry (`B'` accumulates `3^(i-1)` weighted by later powers of 2;
`B` accumulates `3^(L-1-i)` weighted by earlier powers of 2); what is NOT dual is branch selection:
`T_b` takes the forced 2-adic valuation `k=v_2(3n+b)`, `G_b` takes the deterministic 3-adic minimum
`j` (the E-graph offers all `j` in two classes mod 6, the greedy rule picks the least). Example
`b=1`: `{-4,-11}` has word `(3,0)`, `-4(2^3-3^2)=4=B'(3,0)`.

**Theorem 4.3 (scaling lemma; PROVED, S4.5).** `G_b(bm)=bG_1(m)` with the same `j` (the target set
`{1+3b^{-1},1+6b^{-1}} mod 9` is `{4,7}` for every unit `b`). Hence every `G_b` has the three universal
cycles `b{1}`, `b{-1}`, `b{-4,-11}`, exactly as `T_b(bn)=bT_1(n)` gives `T_b` the four universal cycles
`b{1}`, `b{-1}`, `b{-5,...}`, `b{-17,...}` (braids2: content `d=gcd(cycle,b)=|b|/q`, universal iff `q=1`).
Cycles not inside `bZ` have content `d<|b|` (`q>1`; "non-universal" below); the braids2 primitive cycles are
those with `d=1` (`q=|b|`). The two notions coincide for prime `|b|` and differ for `b=25,35,49`, where a
non-universal cycle can be a dilation of a primitive cycle of a proper divisor (the draft's "primitive
`= q=|b|`" label for the census column is REFUTED: `G_25` has `5*{-4,...}`, `5*{16,...}` from `G_5`).
Verified for all `|b|<=49`, `|m|<=200`.

FINITE-EXACT census (S4.4, S4.6; `|m|<=10^5` for `G_b`, odd `|n|<=2*10^5` for `T_b`, escape bound
`10^18`, no escapes, negation symmetry verified, F3 counts `1:4, 5:9, 7:5, 11:7, 13:13` reproduced, gate
identity verified on every cycle):

| b | #G_b | non-universal G_b (d<b) | primitive G_b (d=1) | #T_b | non-universal T_b | primitive T_b (d=1) | non-universal G_b cycles (minimum, content d) |
|---|---|---|---|---|---|---|---|
| 1 | 3 | 0 | 0 | 4 | 0 | 0 | -- |
| 5 | 5 | 2 | 2 | 9 | 5 | 5 | (-4,1), (16,1) |
| 7 | 4 | 1 | 1 | 5 | 1 | 1 | (-4,1) |
| 11 | 5 | 2 | 2 | 7 | 3 | 3 | (-16,1), (-17,1) |
| 13 | 7 | 4 | 4 | 13 | 9 | 9 | (-7,1), (209,1), (236,1), (256,1) |
| 17 | 6 | 3 | 3 | 8 | 4 | 4 | (-10,1), (-59,1), (-64,1) |
| 19 | 6 | 3 | 3 | 6 | 2 | 2 | (-13,1), (-16,1), (-85,1) |
| 23 | 8 | 5 | 5 | 15 | 11 | 11 | (-13,1), (-14,1), (-29,1), (-2621,1), (-3620,1) |
| 25 | 7 | 4 | 2 | 12 | 8 | 3 | (-13,1), (-19,1), (-20,5), (80,5) |
| 29 | 4 | 1 | 1 | 9 | 5 | 5 | (-20,1) |
| 31 | 4 | 1 | 1 | 6 | 2 | 2 | (-16,1) |
| 35 | 7 | 4 | 1 | 12 | 8 | 2 | (-19,1), (-20,5), (-28,7), (112,7) |
| 37 | 4 | 1 | 1 | 7 | 3 | 3 | (-28,1) |
| 41 | 6 | 3 | 3 | 5 | 1 | 1 | (-23,1), (-26,1), (-101,1) |
| 43 | 5 | 2 | 2 | 5 | 1 | 1 | (-22,1), (-64,1) |
| 47 | 6 | 3 | 3 | 11 | 7 | 7 | (-25,1), (-29,1), (64,1) |
| 49 | 8 | 5 | 4 | 7 | 3 | 2 | (-28,7), (-43,1), (-52,1), (-59,1), (-61,1) |

`b<0` is the negative of `b>0`. Within `|m|<=10^5`, `G_1` on negatives has exactly the 2-cycle `{-11,-4}`
and the fixed point `-1`, and the only positive `G_1`-cycle is `{1}` (Corollary 3.4 proves this for all
lengths `<=12`). The common cycles of `G_b` and `T_b` are exactly the two fixed points `{b},{-b}` for every
`b` (criterion (ii)). REFUTED: the reversal correspondence (ii) as a bijection `G_b`-cycles `<->` `T_b`- or
`T_{-b}`-cycles (intersection = two fixed points; census counts `#G_1=3` vs `#T_1=4`, `b=13`: 7 vs 13).
SCOPE: an abstract bijection between the full, census-bounded cycle sets is not a mathematical claim and
is not refuted. Sample words:
`b=13` positive cycles `[337,445,589,781,256]` word `(2,2,2,0,2)`, `[209,553,733,973,320]` word
`(3,2,2,0,1)`, `[236,625,829,272,721]` word `(3,2,0,3,0)`, all with `2^8-3^5=13` and `m_0=B'`.

## 5. Hostiles

**Theorem 5.1 (PROVED, S5.1).** For `m=3^j+1`: `G^i(m)=4^i 3^(j-i)+1` for `i<=j-1` (each is `1 mod 9`
while `j-i>=2`, letter 2), `G^(j-1)(m)=3*4^(j-1)+1 = 4 mod 9` (letter 0), and `G^j(m)=4^(j-1)`.
So `peak/m >= (3*4^(j-1)+1)/(3^j+1) -> infinity` like `(4/3)^(j-1)`: `j=15`, `m=14348908`,
`G^14(m)=805306369`, ratio `56.123`; all reach 1 (e.g. 60 steps at `j=15`). This is the 3-adic mirror of
Collatz's `2^L-1` (`T^L(2^L-1)=3^L-1`).

**Theorem 5.2 (sharp peak bound; PROVED, S5.2).** Along every word `K_n <= 2n+[k_1=3]`, hence
`G^n(m) < 2*(4/3)^n m` (and `<(4/3)^n m` unless `m=5 mod 9`). Proof: letter 3 is residue 5 mod 9, which
is never followed by 5 (`G(5+9t)=13+24t = 4+6t mod 9`) and is entered only from `{7,8}` (letters 0,1);
pair each non-initial 3 with its predecessor (sum `<=4`), the pairs are disjoint, unpaired letters are
`<=2`. The family `3^j+1` attains `K_(j-1)=2(j-1)`, so the rate `4/3` is optimal (Collatz: `3/2`).
Consequently "J consecutive `k=3` letters" is impossible for `J>=2`: the only `8/3`-patterns are
`(7 or 8),5,1^n`. The residue pattern `8,5,1^(n-1),{4 or 7}` (word `1,3,2^(n-1),0`, `n+2` letters, factor
`(16/9)(4/3)^(n-1)` after `n+1` steps) is realized by exactly TWO classes mod `3^(n+3)`: `m=8 mod 27`,
`G^2(m)=13+48s`, and `v_3(G^2(m)-1)=n` exactly iff `v_3(1+4s)=n-1` exactly, one class of `s mod 3^(n-1)`
and two of its three lifts mod `3^n` (the third lift continues the run of 2's); e.g. `8, 35 mod 81` at
`n=1`, `143, 224 mod 243` at `n=2` (S5.3; the draft's "exactly one class" and "`1^n`" were off by one).
Smallest members `8,143,62,1277,548,11483,4922,103337,44288,930023,398582` for `n=1..11`, growth ratio e.g.
`31.5692` at `n=11` against `(16/9)(4/3)^10=31.5693` (S5.3-S5.4).

FINITE-EXACT (S5): all `6666667` non-multiples of 3 in `[1,10^7]` reach 1; max steps 93, attained
exactly twice, at `m=8751065` (the session lead's F1 names this first attainer) and `m=9454814`; max
`peak/m=133.033` at `m=4847486`, peak `644874241` at step 17; `K_n<=2n+[k_1=3]` held
along every orbit. Full greedy word of `m=4847486` (76 letters):
`3,2,2,2,2,2,0,3,2,2,2,2,2,2,2,2,2,0,2,0,0,3,0,0,2,2,0,3,2,2,2,0,0,1,2,2,2,2,0,0,0,2,0,0,3,2,2,0,3,0,0,1,3,0,3,2,2,0,3,0,0,1,0,3,2,0,0,0,2,0,0,0,1,3,0,0`;
the 3-adic reason for the peak: `m=5 mod 9`, `v_3(G(m)-1)=6` (run `2^5 0`), `G^7(m)=5 mod 9`,
`v_3(G^8(m)-1)=10` (run `2^9 0` ending at the peak), `K_17=34`, `2^34/3^17=133.03`.

## 6. The reframe and the typed dualities

**Theorem 6.1 (PROVED, S6.1).** Let `U={n>=1: 3` not dividing `n}`; Q1: every `n>=1` reaches 1 in E; Q2:
1 reaches every `m in U` in E. Then `(Q1 and Q2) <=> (U is one strongly connected component of E
containing 1, and every multiple of 3 is a transient singleton feeding it)`. (`=>`: `m -> 1 -> m'`;
`3Z` is never entered and leaves via `3t -> 9t+1`, wave-one S1.3. `<=`: Q2 is reachability inside the
SCC; Q1 for `n=3t` goes through `9t+1`.) SCOPE: `Collatz => Q1` (the deterministic path is an E-path);
`Q1 => Collatz` is not established here because E has the extra `even -> 3n+1` arrows. Q2 carries the
`G`-certificate: a `G`-orbit `m -> ... -> 1` reversed is an explicit E-path `1 -> m`, so Q2 follows from
"every `m in U` reaches 1 under `G`" (FINITE-EXACT to `10^7`; stopping-time density 1 by Section 3;
OPEN in general, as for Collatz: no positive `G`-cycle but `{1}` and no divergent positive `G`-orbit).

Typed analogies (S6 of the output).
A1 Source: Collatz `T` on `Z_2` (CITED Lagarias 1985, Terras 1976). Target: `G` on `Z_3^x`. Map: 2-adic
digits `<->` 3-adic digits, parity word `<->` k-word, `n mod 2^J <-> m mod 3^(J+1)`, gate
`n_0(2^K-3^L)=bB(w) <-> m_0(2^J-3^L)=bB(rev w)`. Preserved: word determined by the residue at depth `J(+1)`,
gate shape, Terras-type density theorem. Lost: Haar invariance and surjectivity of the word map onto the
full shift (`G` is 4:1/2:1, invariant density `4/3, 2/3`; the k-word map is injective onto a proper SFT
of entropy `log 3 < log 4`); the carry direction. Sidecar: the rank-one `M_u` for every `u`. Test:
`E[k]=1` versus `E[parity]=1/2`.
A2 Source: `6/pi^2=prod_p(1-p^-2)`. Target: the law `(4/3,2/3)`. Preserved: density of a residue set
equals its Haar mass. Lost: multiplicativity over primes -- `G` lives at the single prime 3. No map found
from `6/pi^2` to any `G`-quantity (SCOPE).
A3 The three pieces mod 3 of E: `0 mod 3` = transient leaves; `1 mod 3` = heavy class (`4/3`); `2 mod 3` =
light class (`2/3`), glued by the transfer law; a statement about the greedy inverse walk, not about
forward Collatz (whose image residues mod 3 are `(-1)^k`).
A4 `3n+b` versions: same chain and drift for every `b` (Theorem 4.1); typical cycles small with sign
`-sign(b)`. No map found from the census to Bott periodicity or octonions (SCOPE; the only periodicity
in the data is the order 6 of 2 mod 9).
A5 Source: the blueprint audit's affine word model (affine note, Sec. 1 and 3): chronological D/E
words, `W(n)=(2^r n-B)/3^m`, legality iff `2^r n = B mod 3^m` (one legal source class mod `3^m`), fixed
point `B/(2^r-3^m)`, and the ternary/dyadic progression bijection `s+3^m j <-> W(s)+2^r j`. Target: the
greedy word of `G`. Map: letter `k>=1 -> E o D^(k-1)` (`3G+1=2^k m`); letter `k=0 ->` the E-graph's extra
arrow, outside the D/E semigroup (wave-one leaf identity, fibre index `j=-1`). Preserved: the carry
formula (`r=K_J, m=J, B=B_J`), the legality congruence, the fixed-point gate, freeness of the word
semigroup. Lost: the dyadic endpoint address mod `2^r` and the uniform word mass `2^-r` -- `G` chooses
the word from the 3-adic side (`m mod 3^(J+1)`, one extra ternary digit beyond the legality class
mod `3^J`), so word masses are the Markov products `2*3^J P_chain(w)` (wave-one S4). Sidecar: the
blueprint's spectral hostile (DEDE vs DDEE, equal trace and determinant, different carry) is precisely
the `B'(w)` dependence of the gate in Theorem 4.2. Test: "same words, opposite adic completion" -- the
blueprint counts words by dyadic endpoints (E-count binomial), `G` by ternary sources (`E[k]=1`).

## Reproduction

```bash
python3 04-computation/experiments/collatz_mod6_20260917_three_adic_g_map.py > 05-knowledge/results/collatz_mod6_20260917_three_adic_g_map.out
```

Runs in about 17 s with peak RSS about 0.55 GB (the `10^7` sweep is chunked); every check is an explicit
`raise` and the `python3 -O` replay is identical modulo timing lines. Output ends with `ALL CHECKS PASSED`
and the source sha256 (`ffcd4f63571787101d0cd723c17e417f2a5880cf70206c7320a8e45bd27edf58`). Sections S0-S6
of the `.out` correspond to Sections 1-6 above (S1.3 conjugacy, S3.5-S3.9 moments/threshold/cycle
boundary, S4.5-S4.6 scaling/non-universal table, S6.1 reframe).

Independent audit (2026-09-21, separately written code, 66 s, 0.57 GB):

```bash
python3 04-computation/experiments/collatz_mod6_20260917_three_adic_g_map_audit.py > 05-knowledge/results/collatz_mod6_20260917_three_adic_g_map_audit.out
```

It recomputes the transfer law, the resolving clause of Theorem 1.3(ii) (and the draft's false variant), the
Markov property at level 3 from conditional laws, `mu`-invariance, `d_J`, `g_J` and Terras `F(J)` for
`J<=12`, `sigma=sigma_res` to `10^6` (pure Python), the moment identity also at `u=1/2`, the threshold and
cycle-boundary enumeration, the census with contents `d`, the fixed points `{b,-b}`, the `10^7` sweep with
ties, the two-class count of the `(8,5,1^(n-1))` family for `n<=7`, and the E-path reversal; it ends with
`ALL AUDIT CHECKS PASSED`.

## Stopping boundary / next question

The 3-adic side is now completely symbolic: `(Z_3^x,G)` is a 6-state SFT with Parry measure `mu`,
exact moment generating function `log((u^2+u+1)/3)`, and exact per-level thresholds. What is not
symbolic is the integer side: positive `G`-cycles other than `{1}` and divergent positive orbits
(Q2), exactly as for Collatz. The one new lever is Theorem 3.3 with Corollary 3.4: a positive `G`-cycle
of length `L` with word `w` is an integer `m*(w)=B_L/(2^{K_L}-3^L)` lying in its own class (none for
`L<=12` but `m=1`), and a mismatch `sigma<sigma_res` is an integer of its class strictly below `m*(w)`
that has not descended earlier. The draft's proposed next question, "prove `m*(w) <` least positive
member of the class of `w` for every growth word", is REFUTED by the enumeration itself: the word
`(1,2,2)` at level 3 is the class `2 mod 81` with `m*=37/5=7.4>2`; integer members below `m*(w)` exist
at every level `i>=3` (the candidate column of Theorem 3.3 sums to 41 for `i<=12`, all descending earlier). The correct next question is
the two-regime statement: (a) every integer `m>=2` of a growth class below `m*(w)` descends at an
earlier step (this is `sigma=sigma_res` for all `m`), and (b) no growth word other than `2^L` has an
integer fixed point in its own class (this is "no positive `G`-cycle but `{1}`"). For (b) the first
place to test is `b=13`, where every word with `(J,L)=(8,5)` has `2^J-3^L=13` and the gate
`m_0=13B'/13=B'` is automatically an integer, so a positive cycle only needs `B'(w)` to land in the class
of `w` mod `3^6` -- which happens for the three census cycles `337,209,236` and their rotations. For `b=1`
the denominator `2^J-3^L` equals 1 only at `(J,L)=(2,1)` for `L>=1` (PROVED: `2^J-3^L=1` with `L>=1`
forces `L` odd by mod 4, then `2^J = 3^L+1 = 4 mod 8`, so `J=2`), so every other positive `G_1`-cycle needs the divisibility
`2^J-3^L | B'(w)` with `2^J-3^L>1`, which the gate does not control; the enumeration of Corollary 3.4
shows it fails for every growth word of length `<=12`.
