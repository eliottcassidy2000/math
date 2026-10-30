# Owner promotion and harmonic completion for second-band cores

Session 116, 2026-10-30. **PROVED for the explicitly defined family below.**
Here “second-band” means precisely the explicitly defined family below; no claim about unverified historical terminology is needed. No lower-chamber or general LRC closure is claimed.

## Concurrent provenance

The incoming commit `98f961a7aa142a739403de1ba343720371d6125c` publishes the same second-band owner-promotion mechanism and harmonic transfer in [the concurrent owner-family packet](../../odd_eighth_owner_families_20261030/README.md). It arrived during this session and is credited here. The present proof and independent audit were developed separately; they are not an exclusive discovery claim. The explicit atlas/minimum census, surgery formulas and boundary controls below give this session's verified presentation. The incoming finite q=8 classifications and native signed-cycle results are not claimed as results of this packet.

## Family and complete atlas

Let q>=3, U=(Z/qZ)^* represented in {1,...,q-1}, and

 C(q,A)={1,...,q-1} union {q+b:b in A},
 1 in A subset {1,...,q-1}.

Define e_A(u)=u+q 1_A(u), E_A={e_A(u):u in U}, s_A=min E_A.
For every 1/(q+2)<r<=1/q, put epsilon=1/q-r. The full good set is the disjoint union over primitive a/q of

 [a/q-epsilon/e_A(p), a/q+epsilon/e_A(q-p)], p=a^-1 mod q.

In particular M(C)=1/q. Every owner e in E_A occurs once as a left owner and once as a right owner. Nonunit second-band attachments are completely invisible in this open chamber, not necessarily below it.

### Proof

The Session114 Farey-tent argument confines the AP good set above the stated floor to denominator-q and denominator-(q+1) tents. Since q+1 is present, all denominator-(q+1) tents are excluded. Its argument also replaces owner1 by q+1 in the two relevant directions. Thus we need only prove the additional walls do exactly the stated promotions inside those known complete windows.

Consider a left direction at a/q and let p=a^-1 mod q. On the candidate window write t=a/q-d, 0<=d<=epsilon/e_A(p). An added speed v=q+b has fractional phase c/q at the center, c=ab mod q in {1,...,q-1}. (Nonzero follows because a is a unit and b is nonzero.) A residue c=1 is equivalent to b=p. For p>=2, throughout the entire unpromoted AP side, (q+p)d<2(q+p)/[q(q+2)p]<=1/q; its phase therefore cannot cross zero and re-enter safety after wrapping. Its lower wall is exactly d<=epsilon/(q+p), proving necessity of that promotion. For p=1 the already-proved q+1 window supplies this same wall (and its earlier wrap exclusion). The remaining lower walls must satisfy c/q-vd>=r. Upper walls are automatically safe since c/q<=1-1/q<=1-r and the phase decreases.

If p>=2 and c>=2, v=q+(cp mod q)<=q+cp. Consequently

 2(v-p)<=2q+2p(c-1)<=p(c-1)(q+2),

because p(c-1)>=2. Since epsilon<2/[q(q+2)], this implies

 (v-p)epsilon <= p(c-1)/q,
 c/q-v epsilon/p >= 1/q-epsilon=r.

As e_A(p)>=p, the actual candidate window is no wider. If p=1, its owner is already q+1; for c>=2, v<2q and

 v epsilon/(q+1)<4/[(q+1)(q+2)] < 1/q <= (c-1)/q+epsilon.

This proves the needed lower wall too. No hidden wrap is possible because the raw phase has just been bounded between r and 1-r. Reflection proves the right directions. All promoted owner constraints are necessary and all other constraints safe, proving exact equality of the complete good sets.

## Complete finite classification of owner atlases

Only A intersect U matters. There are exactly 2^(phi(q)-1) distinct atlases, each represented by 2^(q-1-phi(q)) actual subsets A. All distinctness statements follow because e_A(u) retains residue u modulo q.

The minimum owner is the least unpromoted unit if any exists, otherwise q+1. Thus its complete possible list is

 (U minus {1}) union {q+1}.

For s=u in U minus {1}, exactly 2^#{v in U:v>u} atlases have minimum s; for s=q+1, exactly one does. Multiply these numbers by 2^(q-1-phi(q)) for actual cores. In particular the least-missing-prime rule is special to A={1}: q=7,A={1,2,3} has minimum owner4, not a prime. The fully promoted family has s=q+1 regardless of optional nonunit attachments.

## Uniform harmonic extension

Let V=C(q,A) union {qm,2qm,...,kqm}, m,k>=1.

1. For 1<=k<=q and 2m>s_A,

 M(V)=min{m/(qm+s_A),1/(k+1)}.

If the first term is strictly smaller, exactly two maximizers occur:
 t=(am+h)/(qm+s_A), 1-t,
 a=s_A^-1 mod q in {1,...,q-1}, h=(a s_A-1)/q.
Their terminal mass coefficient is B=2(1/s_A+1/(qm)). This formula does not cover equality of the two bounds.

2. For k=q, without any condition on m,

 M(V)=1/(q+1) if and only if m>=s_A.

For that preserving range all maximizers, with no duplicates, are

 t_(e,j)=a_e/q-j/[q(q+1)m], 1-t_(e,j),
 e in E_A, a_e=e^-1 mod q,
 1<=j<=floor(m/e), gcd(j,q+1)=1.

With n=q+1, their number and terminal mass coefficient are

 N_A(q,m)=2 sum_(e in E_A) sum_(d|n) mu(d) floor(m/(ed)),
 B_A(q,m)=2/(qm) sum_(e in E_A) sum_(1<=j<=floor(m/e), gcd(j,n)=1) [1/ell_j+1/(n-ell_j)],
 ell_j=j^-1 mod n in {1,...,n-1}.

These are the Session115 formulas with the actual promoted-owner set, not the original speed minimum. To verify: set y=q|t-a/q|. Core safety is r<=(1-e y)/q, first-tail safety requires r<=||my||<=my; equality uses the unique least owner. Setting y=r/m makes all k tails safe when r<=1/(k+1). For k=q, the q-tail AP can reach its maximum 1/n only at primitive phases my=j/n. Substitution in the complete windows gives exactly j<=m/e. At any such point the fastest tail owners are qm ell_j and qm(n-ell_j). If the core endpoint ties, ej=m, so e<=m<qm and the tail still wins. Summing their reciprocal physical speeds gives B. Completeness of the windows excludes other maxima.

Set H_E=sum_e 1/e, H_U=sum_(1<=u<n,gcd(u,n)=1)1/u and L=n lcm(E_A). Exact periods and limits are

 N_A(m+L)-N_A(m)=2L phi(n) H_E/n,
 (m+L)B_A(m+L)-m B_A(m)=4L H_E H_U/(qn),
 N_A(m)/m -> 2 phi(n) H_E/n,
 B_A(m) -> 4 H_E H_U/(qn).

At m=s_A there are exactly two maxima and B=2(q+1)/(q^2 s_A).

## Creative operations with exact consequences

**Owner-promotion surgery.** Adding q+u for an unpromoted unit u replaces only u by q+u. At fixed preserving m, the exact lost count is

 2 #{j: floor(m/(q+u))<j<=floor(m/u), gcd(j,q+1)=1}.

The lost mass is the same sum with summand 2 w_(q+1)(j)/(qm). No other phase pair is altered. In the common preserving range this is a monotone loss law, while the clearance remains exactly 1/(q+1). The limiting count density drops by 2q phi(q+1)/[(q+1)u(q+u)], and limiting B drops by 4 H_U/[(q+1)u(q+u)]. This separates constant clearance from tunable surviving measure.

**Arbitrary joint banks.** For a finite multiplier bank T define tau_T(r)=inf{x>=0: ||tx||>=r for every t in T}, infinity if empty. In the same open chamber, C(q,A) union qT has nonempty good set if and only if qr+s_A tau_T(r)<=1. This is the concurrent first-arrival theorem applied to a now-proved complete family of atlases. Joint arrival is essential; separate minima cannot be substituted.

**No false lower floor.** If 2m<=s_A we infer only M<=1/(q+2), not equality. For q=3,A={1,2},m=1,k=3, the source includes {1,2,3,4,5,6,9}; its maximum is 1/7<1/5 (exact checker). Additional core walls can destroy the old single-core floor. The chamber endpoint is genuinely excluded: for q=3,A={1}, the point t=1/5 is safe at r=1/5, but is absent from the displayed promoted-owner windows (their first interval starts at3/10). Higher-band additions require a new wall audit; the second-band proof does not authorize arbitrary positive attachments.

## Prior credit and verification

The Farey confinement and A={1} windows are Session114 `lrc/GENERAL_CORE.md`; the harmonic classification, fastest-tie rule, census and periods are Session115 `lrc/HARMONIC_COMPLETION.md`; joint arrival is `odd_eighths_uniform_families_20261030/LRC.md`, L1. The new work is the complete owner-promotion family proof and its classified transfer and surgery laws. No literature priority is claimed. `check_second_band.py` uses exact rational interval intersection independent of the window prediction; normal and optimized receipts must agree. Finite checks supplement the proof, not replace its quantifiers.

## Odd-eighth specialization survives every owner promotion

For q=7 and k=7, every member of the defined family has harmonic clearance 1/8 exactly when m>=s_A. The physical tail phases are then only 1/8,3/8,5/8,7/8. All four occur if and only if m>=3s_A: phase indices ±1 occur from j=1, while the first representative of ±3 is j=3, and no owner is smaller than s_A.

The phase weights are independent of which core was chosen: each point with phase ±1/8 contributes 8/(49m), and each with phase ±3/8 contributes 8/(105m). Core choice changes only which points survive and their multiplicities. With chi_2(j)=(2/j) on odd j, define

 A_A(m)=sum_(e in E_A) sum_(odd j<=floor(m/e)) chi_2(j).

Then the exact character-resolved formula is

 B_A(7,m)=88 N_A(7,m)/(735m)+64 A_A(m)/(735m).

Here A_A is bounded and periodic with period 8 lcm(E_A). Also N_A/m -> H_E and B_A ->88 H_E/735. Thus the eighth-character signature is universal across the whole second-band family, while owner surgery tunes its multiplicity and limiting mass. For q=7,A={1,2,3}, E_A={8,9,10,4,5,6}, s_A=4: preservation starts at m=4, all four phases start at m=12, and the threshold mass is B=4/49. Full promotion gives E_A={8,9,10,11,12,13}, s_A=8: preservation starts at m=8 and all four phases at m=24. These are boundary-family identities, not an identification of LRC physical weights with unrelated eighths in other subjects.

## Additive family coordinates versus nonlinear threshold

For fixed q and m>=q+1, every allowed core is in the preserving range. Regard optional unit promotions u in U minus {1} as Boolean coordinates. Both N_A and B_A are affine functions of these coordinates: start with A={1}, and subtract the independent single-owner losses displayed above. Hence every mixed second finite difference vanishes exactly. The limiting reciprocal-owner statistic has the particularly uniform expression

 H_E= sum_(u in U)1/u - q sum_(u in A intersect U)1/[u(q+u)].

The optimum-preserving threshold s_A is instead the first unpromoted unit (or q+1); it is an order statistic, not an affine function. Thus this family separates two kinds of data: a nonlinear bottleneck controls when the harmonic optimum appears, while an additive owner statistic controls its count and mass. Nonunit promotions are zero coordinates for all these chamber observables. This is an exact finite-family decomposition, not a claim that arbitrary LRC sources admit the same compression.
