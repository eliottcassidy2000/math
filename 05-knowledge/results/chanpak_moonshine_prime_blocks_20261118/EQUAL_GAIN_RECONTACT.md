# Destroying the terminal cut: equal gains and unequal arc usage

Status: PROVED universal accounting identity; FINITE-EXACT physical 131/133 pair and scoped inflation tests. No all-order repetition, minimum order, or new discovery of these scalar counts is claimed.

The companion `ORDINAL_RECONTACT_DETERMINANT.md` classifies the simplest selective recontact: its output c(P_2(A)+H(A)) has a terminal factor c, and its neutral-swap ±2 case is only 3↔5. Here we change the operation, not merely increase its bank: reverse one cross-seam arc from the terminal C block back into A. This removes the terminal ordinal cut.

## All-order accounting

For distinct u,v in any tournament T define the directed seam cover

J_T(u,v) = sum_{u in S, v notin S} E_u(T[S]) S_v(T[V\S]),

where E_u counts paths ending at u and S_v paths starting at v. The seam itself is not required. If u→v is an arc, J_T(u,v) is its actual Hamilton-path usage. Reversing only that arc gives

H(T^{uv})=H(T)−J_T(u,v)+J_T(v,u).                 (1)

Proof: unaffected paths cancel. A path using the changed arc is uniquely a prefix ending at its tail and a complementary suffix starting at its head. The changed arc is never inside either induced part, so their path counts are the same before and after reversal.

Consequently, if two original arcs e=(u,v),f=(s,t) have equal reverse gain J_T(v,u)=J_T(t,s)=B and forward usages differing by two, reversing e or f gives counts differing by exactly two. More generally a bank with common B and forward usages forming a length-r step-two interval gives a physical step-two output interval, based at H(T)+B minus that interval.

This is a recovered exact accounting identity, not by itself a selector. To hit prescribed N the bank must contain usage H(T)+B−N; equal gains alone do not position the bank absolutely.

## Literal escape from the factorized family

Let A be the two-vertex transitive tournament, outgoing masks [0,1]. Let C have outgoing masks [16,1,3,7,14]. Its profile is (H,P_2)=(9,84). Put A on vertices0,1, singleton b=2, C on3..7, and x=8. Orient A→b→C and x→A,C, b→x.

The common base has H=45. Two choices of a single return contact give:

| original arc reversed | removed paths | added paths | new H |
|---|---:|---:|---:|
| 1→6 (A1→C3) | 4 | 90 | 131 |
| 0→3 (A0→C0) | 2 | 90 | 133=7·19 |

Thus both modules are nontrivial, the new backward gain is exactly equal, and the two-path correction is carried solely by a difference in old forward usage. The two final tournaments are strong; this is checked by exact reachability. Their outgoing masks are

131: [252,189,504,128,8,24,58,112,251]
133: [244,253,504,129,8,24,56,112,251].

These are two different arc reversals from a common host, not a claim that the final tournaments are adjacent by a single flip, and not a claim about a merged-metagraph line. Their P_2 and endpoint vectors differ; the full continuation is not preserved automatically.

## What does and does not repeat

The transitive-twin inflation probe replaced each of C's five vertices in turn by a transitive module of length1..4 and retained the same two indicated port choices (the first copy of an inflated port). Every nontrivial such inflation destroys the +2 gap. For vertex1 or2, the pairs are

(131,133), (239,269), (455,541), (887,1085).

The reverse gains become (90,90),(162,186),(306,378),(594,762): the failure is loss of equal gain, not merely the center moving. For the fixed C and transitive A of sizes1..8, an exhaustive scan of all A–C return ports finds a common-gain/usage-gap-two bank only at size2. These statements are FINITE-EXACT only.

The useful repeatability target is therefore concrete: **change C or the selected ports so that the two backward seam covers remain equal while the old arc usages retain a gap of two, and simultaneously select H(T)+B from the target**. This is a three-coordinate marked-profile requirement. It allows the high counts to move together rather than freezing them. The finite example shows the requirement is physically possible after removing the terminal cut; it does not yet show a family at arbitrary heights.

## Reproduction

`python3 contact_replay.py` writes the exact profiles, usage/gain receipts, reachability checks, and all five twin-inflation tables. `python3 -O contact_replay.py` gives byte-identical JSON. The larger transitive-A scan is a separate exploratory computation. The standalone contact replay does not import the exploratory scripts.
