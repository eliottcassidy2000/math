# A target-selected smaller source for a physical adjacent-count repair

**PROVED inverse characterization for a specified two-word construction.** This connects the new parent interface to the inherited near-floor theme; it does not prove that every requested profile exists.

## Coprime deletions alone do not give a common-word +2 repair

For the equal-(H,T) rerooted pair in `COPRIME_DESCENDING_PARENTS.md`, the difference after a common word W is -Delta J(W), where J(W)=K(W)+(W mod2) and Delta=2^s. Therefore Delta>=4 forbids a common-word difference ±2 at every length. If Delta=2, a common-word ±2 requires J(W)=1, which occurs exactly at numerical address W=1: K is nonnegative and even; every odd W>=3 has K(W)>=2. Appending trailing zeros does not create a different count. Thus the small +2 opportunity is not automatically repeatable merely because its parents are coprime.

There is a literal small control. The roots of `(20,1,2,7,14)` at1 and `(16,1,3,7,14)` at0 have profiles (9,3,9,3) and (9,1,9,1). Common1 gives15 and17, respectively. This proves the positive boundary case without asserting that the much larger dyadic differences have the same freedom.

## A different two-word repair uses one source and an exact floor condition

Let a physical tournament have a root with exactly one outneighbor. Write h for its count, d for the deletion count, and eta for deletion paths ending at that unique outneighbor. Its full profile is

    (h,e,t,d)=(h,d-eta,h+eta,d).

Indeed ending at the root permits every deletion endpoint except the unique outneighbor, giving e=d-eta. The weighted-start identity is t=h+eta, as in the inherited sole-exit transfer.

Compare the two actual binary words with numerical addresses2 and3 (chronological01 and11). Their counts are exactly

    F(2)=3h+2eta,
    F(3)=F(2)+(h+2eta-3d).

Consequently they form an adjacent odd pair with signed gap g=±2 precisely when

    h+2eta-3d=g.

This is the exact meaning of the near-floor requirement for this particular repair. It uses a different source when necessary, rather than demanding the earlier coprime pair already have the right gap.

## Work backwards from the absolute requested baseline M

To produce M and M+g by these two words, the proposed smaller source MUST have

    eta=(M-3h)/2,
    d=(M-2h-g)/3,
    e=d-eta, t=h+eta.

The elementary nonnegativity constraints eta>=0, eta<=d<=h are equivalent to

    max((M+2g)/5,(M-g)/5) <= h <= M/3.

In addition h is positive odd and h=g-M modulo3. These conditions make eta a nonnegative integer and d a positive odd integer, with e>=0. (For very small M explicitly retain d>=1.) Every successful source is strictly smaller than M for M>0. Thus the target fixes a finite list of possible source heights in a concrete interval and fixes the ENTIRE rooted profile for each height. No free carry or approximate center remains.

**Exact reduction.** A prescribed pair (M,M+g) is obtainable by this sole-exit two-word repair if and only if at least one of these explicit profiles is physically realizable with a sole-exit root. The forward implication follows from the transfer; the reverse implication uses the supplied source and the two literal attachments. This is not a universal physical-supply theorem. It is a precise descending-parent obligation that can be tested or attacked structurally instead of searching arbitrary gadgets.

Examples:

* Baseline15, g=+2 selects h=5, eta=0, d=1. The actual source `(8,1,3,6)` rooted at0 has (5,1,5,1). The two words give15 and17.
* Baseline31, g=-2 selects h=9, eta=2, d=5. The actual source `(20,1,2,7,14)` rooted at2 has (9,3,11,5). The two words give31 and29.

The exclusion of7 or21 must still enter through physical source admissibility, not by ignoring a candidate's profile. For example a baseline5 with g=+2 would demand an impossible adjacent pair containing7. The numerical interval alone is not permission to assume the required source exists.

The Collatz-like feature is now exact and limited: the source height lies between roughly M/5 and M/3, so there is genuine descent IF its marked profile can be supplied. Proving physical supply at all admissible targets remains the missing theorem. The two-word identity itself is elementary and should not be promoted as a newly solved near-floor conjecture.

Replay: `python3 finite_edit_recontact/sole_exit_target_replay.py` independently checks the literal sources and their outputs, verifies the sole-exit profile identities for every labeled tournament through order5, and checks necessity/sufficiency of the arithmetic parameter constraints on finite integer controls.
