# Nonunit loops become unit windows under coprime composition

Status: PROVED, independently audited ([supplement audit](../graft_audit/SUPPLEMENT_AUDIT.md)). This strengthens the stable two-unit-core theorem in PROOF.md. The new advance is to reuse individually nonunit moving loops rather than search only for unit loops. Chan-Pak completeness remains OPEN.

## 1. General weighted two-core theorem

Suppose physical cut-labelled families P_i, 0<=i<=K, and Q_j, 0<=j<=L, have fixed full interfaces except T, with deletion path counts h and k and T increments respectively 2u and 2v per index, where u,v are nonzero integers. Ordinal composition has T increments 2ku and 2hv. Reverse either finite index if needed, so that closure counts have positive offsets

    2(a i+b j), a=k|u|, b=h|v|.

If a,b>=2, gcd(a,b)=1, K>=b-1 and L>=a-1, the bounded coin-box lemma in PROOF.md gives every odd count in the translated interval with half-count endpoints C and aK+bL-C, C=(a-1)(b-1). Its size is aK+bL-2C+1. This is valid even if neither original loop is unit-spaced. The whole output remains a common full physical interface, with h-output=hk, because only T changes. Strongness requires the chosen closing construction; the protected construction below provides it.

This supplies a constructive sufficient escape from the repeated-identical-core divisor obstruction. Coprime weighted increments plus a large enough physical rectangle, rather than gcd alone, are the sufficient condition.

## 2. Two explicit nonunit physical loops

Core C13 has outgoing masks

    [2,12,17,5,11,31,63,127], cut23,
    state=(13,832,73,584,46,0,11).

Use twelve-bit words

    X13=000101101001,
    Y13=000001101101.

The states of s1=>C13=>Y13=>s0 and s1=>X13=>C13=>s0 agree in every coordinate except T, and the latter has T smaller by 10. X13 ends1, Y13 starts0, matching the flanks. Hence

    P_(K,i)=s1=>X13^i=>C13=>Y13^(K-i)=>s0

has fixed full interface and T=T0-10i. Its deletion order is 12K+10 and h=13.

Core C29 has outgoing masks

    [8,1,3,22,39,79,31,127], cut14,
    state=(29,1160,133,668,60,0,21).

Use twenty-eight-bit words

    X29=0010001001111010010110000100,
    Y29=1100110100001100001000101001.

The states of s0=>C29=>Y29=>s1 and s0=>X29=>C29=>s1 differ only by -6 in T. X29 ends0,Y29 starts1. Hence

    R_(K,i)=s1=>s0=>X29^i=>C29=>Y29^(K-i)=>s1

has fixed full interface, T=T0-6i, deletion order 28K+11, and h=29.

Both finite identities are verified by independent literal core subset endpoint DP followed by the exact seven-coordinate ordinal law. The endpoint-sensitive context law proves every repetition, not a finite-K extrapolation. Neither core need itself be strong. In each protected deletion family its first vertex is selected and beats all subsequent vertices. After ordinal-composing any number of these families and the inherited h45 family, add v,r with deletion=>v=>r and r beating exactly its selected vertices. The closure is strong: r reaches the first vertex and hence every deletion vertex, while every deletion vertex reaches r via v.

## 3. Efficient pair corollaries

Let Q_(L,j) be the inherited twelve-bit h45 family of THM-4960, with T=T0+2j, deletion order 12L+11. Reverse j and close the ordinal composition.

* C13/Q45 has coefficients 225 and 13. Fix K=12; for every L>=224 it contains 13L-2675 consecutive odd counts on 12L+167 vertices.
* C29/Q45 has coefficients 135 and 29. Fix K=28; for every L>=134 it contains 29L-3723 consecutive odd counts on 12L+808 vertices.

The baselines are exactly the closure counts of P_(12,0)=>Q_(L,L), respectively R_(28,0)=>Q_(L,L), shifted by twice the relevant conductor. The conductors are 2688 and 3752. These are direct physical recurrence definitions, not estimated scalar positions. For example L=300 in the C13 pair yields 1225 consecutive odds on 3767 vertices. L=200 in the C29 pair yields 2077 on 3208 vertices.

## 4. CRT rectangle: a triple-core improvement

Compose P_(12,i)=>R_(28,j)=>Q_(L,k), then close. Relative to the baseline B_L at i=j=0,k=L, replace k by L-t. The exact offsets are

    2(6525i+1755j+377t),
    0<=i<=12, 0<=j<=28, 0<=t<=L.

The 377 representatives 6525i+1755j cover every residue modulo 377=13*29 exactly once: modulo 13 only the i term survives and its coefficient is invertible, and modulo 29 only the j term survives and its coefficient is invertible. Their minimum is 0 and maximum is

    R=6525*12+1755*28=127440.

A general residue-box lemma says that, if a finite nonnegative representative set contains every residue modulo c and has maximum R, then adding {0,c,...,Lc} covers every integer in[R-c+1,cL+c-1]. Proof: for the representative r of the target residue, the previous value r-c is at most R-c, while the first value beyond r+Lc is r+(L+1)c>=cL+c. Thus any target within the displayed interval lies between r and r+Lc in that residue. Duplicate representatives would not invalidate sufficiency; here all are unique by CRT.

Consequently, for every L>=337, one physical full interface supports

    377L-126687 consecutive odd counts

on

    12L+962 vertices.

The exact scalar interval is

    [B_L+254128, B_L+2(377L+376)] in odd steps.

L=337 yields 362 consecutive odds on 5006 vertices. L=1000 yields 250313 consecutive odd counts on 12962 vertices. The asymptotic width/order ratio is 377/12. This is 377 times the inherited h45 one-core family's asymptotic ratio, with no claim of global optimality. A target inside this interval is selected explicitly: solve i=(z*6525^-1) mod13, j=(z*1755^-1) mod29, set r=6525i+1755j, and t=(z-r)/377. The lemma guarantees0<=t<=L; the actual h45 family index is L-t.

## 5. Optional nested-box comparison

Combining h13 and h29 with K=38,L=149 gives 378 consecutive values in one h377 interface, using 4649 deletion vertices. Keep its first 377 members and compose with Q45. The bounded coin lemma gives 377L-16167 consecutive odd values on 12L+4662 vertices for L>=44. This has a smaller allowable L but larger fixed graph overhead. The CRT rectangle above is the simpler final headline. Both are retained in the verifier for independent cross-checks of the composition mechanism.

## 6. Verification and scope

Run python3 verify_weighted.py and python3 -O verify_weighted.py. The checker uses literal core DP, the two finite flanked identities, selected family positions through K=100, every arithmetic selector target in the reported pair/triple windows, and representative selected physical state compositions. It does not enumerate full graphs on thousands of vertices. Both modes agree. nonunit13.json and nonunit29.json are standalone physical certificates. The inherited h45 certificate remains core45.json.

The discovery scan regenerated all 6880 order-eight tournament classes and 26956 saturated core/cut states, then inspected weighted reachable SCCs. That bank was used only to FIND the small certificates; neither the positive theorem nor its verifier depends on census completeness or the discovery search. In particular, h341 was not promoted from a finite search failure to an all-tournament impossibility.

The three families are saturated, not unsaturated. The theorem covers real translated intervals with explicit physical selectors inside those intervals, not arbitrary odd targets. Their baselines still grow exponentially with L. It neither proves eventual coverage of all admissible targets nor Chan-Pak completeness, and it does not use path-orbit omissions as existence statements.

### Exact endpoint receipts

Hashes below are SHA-256 of the lower endpoint as unsigned big-endian bytes, with no leading zero byte. The recurrence and selector in verify_weighted.py reproduce them.

* L=337: lower endpoint bit length 5010; hash 4b5c21120fc4667f69cd62cb31363576d25b48b92726cf3b85fd8e62a0cb8e1e.
* L=400: lower endpoint bit length 5766; hash 4ee3250a06d4c51e6e8f01e78da5b1f98419d00537dbd7b635f2804dded61084.
* L=1000: lower endpoint bit length 12966; hash fe78f2c31a532384fc58a3329368bcdf8bcb3cc22d382b7f8ba3bb9f1fc22ac1.
