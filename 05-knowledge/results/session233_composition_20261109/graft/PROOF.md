# Coprime moving-core composition restores consecutive windows

Status: PROVED, independently audited (see ../graft_audit/AUDIT.md in the integrated packet). No completeness or arbitrary-target placement claim.

## Inheritance

Session231 `fibers/UNIT_INTERVALS.md` proves the h=37 full-interface family using 1188-bit words. THM-4960, `01-canon/theorems/THM-4960-twelve-bit-moving-core-consecutive-odd-intervals.md`, and Session231 `fibers/SHORT_INTERVALS.md` prove its h=45 counterpart using twelve-bit words. Their different deletion counts are the key unused sidecar: 37 and 45 are coprime. Session229 `05-knowledge/results/session229_combined_constructions_20261109/dynamic/PROOF.md` gives the seven-coordinate ordinal law and warns that repetitions of the same h create a growing divisor. MISTAKE-736 distinguishes scalar gap two from saturation; both cores here are saturated and their composition remains saturated.

## 1. A bounded coin-box lemma

Let a,b>=2 be coprime. If K>=b-1 and L>=a-1, set C=(a-1)(b-1) and M=aK+bL. Then

    [C,M-C] intersect Z is contained in {ai+bj: 0<=i<=K, 0<=j<=L}.

This is an actual constructive interval, not a gcd-only assertion. For each residue r=0,...,a-1 of j modulo a, write j=r+a ell. The attainable quotients in that residue are

    t=br+a(i+b ell), 0<=i<=K, 0<=ell<=floor((L-r)/a).

Since K>=b-1, the integer intervals [b ell,b ell+K] meet or are adjacent. They therefore fill every integer quotient from 0 through K+b floor((L-r)/a). The smallest attainable t in this residue is br; the integer immediately before it is br-a <= b(a-1)-a=C-1. The largest attainable t is br+aK+ab floor((L-r)/a). Write L=r+a ell+s, 0<=s<a. This upper endpoint is M-bs, so its next integer in the same residue is M-bs+a >= M-b(a-1)+a=M-C+1. Thus every integer of [C,M-C] in every residue is covered.

An explicit selector is r=t*b^(-1) mod a, i0=(t-br)/a, ell=max(0,ceil((i0-K)/b)), i=i0-b ell, j=r+a ell. The proof gives the required bounds. The lower conductor is sharp for the unrestricted numerical semigroup: C-1=ab-a-b is absent. For a=45,b=37, C=1584 and 1583 is absent.

## 2. Physical state composition

Use the seven coordinates (h,q,A,G,T,S,E) and the inherited exact law

    (h,q,A,G,T,S,E)*(j,r,B,F,U,V,W)
      =(hj,qr,Ar+hB,Fq+jG,AF+Tj+hU+EV,Sj,hW).

Both inputs below are actual cut-labelled deletion tournaments. Let P_(K,i) denote the h=37 deletion family in UNIT_INTERVALS, and Q_(L,j) the h=45 family in SHORT_INTERVALS. In each fixed K or L, every coordinate except T is independent of the position, and

    T(P_(K,i))=T(P_(K,0))-2i,
    T(Q_(L,j))=T(Q_(L,0))+2j.

Their deletion orders are 1188K+10 and 12L+11. Ordinal-compose P before Q, then perform the star closure by adding v,r with all deletion vertices pointing to v, v->r, and r beating precisely the selected cut. This has path count 2A+h-T. It is strong: the first deletion vertex is selected and reaches all subsequent deletion vertices; r reaches it, and every deletion vertex reaches r via v. The closing vertices also lie in this same strong component.

The common interface after composition is fixed; h=37*45=1665. The only position dependence of T is

    -90i+74j.

Let B_(K,L) be the closure count of P_(K,0)*Q_(L,L). Reindex j=L-u. The complete rectangle has exactly the scalar count set

    {B_(K,L)+2(45i+37u): 0<=i<=K, 0<=u<=L}.

This statement includes all values, even collisions; it does not assert distinctness of all graph classes.

## 3. Consecutive-window theorem

For every K>=36 and L>=44 the single common physical interface above contains every odd count from

    B_(K,L)+3168

to

    B_(K,L)+2(45K+37L-1584),

inclusive. Thus it supports 45K+37L-3167 consecutive odd values on exactly

    1188K+12L+23

vertices. Arbitrary common ordinal singleton contexts preserve the same interface difference and unit spacing; strongness must still be checked in the chosen context (the displayed protected construction supplies it directly).

In particular, fix K=36. For every L>=44 there are

    37L-1547 consecutive odd values on 12L+42791 vertices.

The width/order ratio tends to 37/12 rather than the inherited 1/12. This is a 37-fold asymptotic improvement within this explicit family, with a substantial finite overhead, not a global optimum claim. The first instance L=44 gives 81 consecutive odd counts on 43,319 vertices; L=100 gives 2,153 counts on 43,991 vertices. The exact baselines are defined by the physical recurrence and computed by verify.py; they have approximately 43,000 bits already. They are not close to small unsolved scalar targets.

## 4. Exact input certificates and reproduction

core37.json and core45.json are copied verbatim from the inherited Session231 word certificates. The standalone verify.py independently computes both seven-coordinate core states using subset endpoint DP and ordered split sums. It checks each finite flanked word identity, the actual endpoint bits needed for repeated physical movement, family states at selected positions through k=100, and the composed count identity at each tested parameter pair. It exhaustively checks the general coin lemma on 474 small coprime rectangles, and tests the explicit residue selector at every target in each of five displayed large-core windows. Representative target counts are independently recomposed from the selected actual family states. No full-graph subset DP on 43,000 vertices is claimed.

Run `python3 verify.py > CHECK.json` and `python3 -O verify.py > CHECK.optimized.json`; outputs agree. Exact baseline endpoints can be recovered using the exported combined function, avoiding huge printed decimal certificates. CHECK.json gives binary lengths and SHA-256 hashes of the unsigned big-endian lower endpoint. The code uses explicit failure checks, not optimizable assertions.

## 5. Limits and next operation

Repeated copies with the same h still have the old divisor obstruction. Here coprimality removes it only after the bounded rectangle is wide enough, and the residue selector is the physical restoration sidecar missing from a bare gcd argument. Both cores remain saturated. The interval baselines grow exponentially, while their widths grow linearly; there is no all-target baseline selector, interval-overlap theorem, or proof of Chan-Pak completeness. No use is made of automorphism-orbit omissions as a coverage theorem.

A useful next optimization is a shorter unit loop with deletion h coprime to45, or another coprime physical moving core. The theorem needs only two audited full-interface families and their scalar h values, so any such improvement can be substituted without changing the arithmetic proof.
