# Independent audit: two divisor gates and bounded source fibers

**PASS, without mathematical correction.** The packet proves the stated scoped elementary identities, complete divisor criteria, formal factorization classification, dyadic collapse, and bounded-both-g counting obstruction. It does not prove or refute uniform twelve-template occupancy. The hostile target is excluded only from two specified templates and is physically represented by J.

## Frozen inputs and independent path

`INPUTS.json` freezes the five original proof/code/receipt files, their source manifest, and the two already reviewed dependencies used in this audit. The principal original hashes are:

- `PROOF.md`: `e2ec1174853955ea959d0871285d0a28a3b548beddfe1c3e64b48ee56d2d8923`.
- `verify.py`: `557c217ba9d449702f98ee0deddeddc21893344fcf1a2029e45279a26b502f27`.
- Twelve-menu dependency: `f4c2fab8569aa682063f904a78b74026535f3a181a3e01bcc7ac23dc5a591275`.
- Physical core proof dependency: `49ad1123e80fbf751eac3e72ed0f5ad2ec39e7705f9ee89c77d4826ef5e1066c`.

Source originals were not edited. The independent checker imports no author code. It reconstructs every menu matrix from the O/Z/J operations, computes g through actual binary-source words rather than the author's recurrence, uses a different uniform scaled-factor inverse over four parity sheets, and certifies primes by a sieve of all possible prime divisors.

The independent and author checkers both pass normally and under `python -O`, with byte-identical paired receipts. The author replays also match its two frozen receipts exactly.

## 1. Actual g fibers and equality boundaries

The induction for odd, positive g and its upper bound is sound, including the h=2 endpoint in the even-h branch. Removing powers of two preserves g. For an odd core a=2h-1, the lower equality is possible precisely when h is a power of two; the odd-h branch has a strict extra increment. This gives a=2^k-1 with k>=2. Upper equality recurses through odd h until the exceptional endpoint a=3, giving exactly a=2^k+1 with k>=1.

Thus no hidden infinite collection of odd cores occurs in a fixed g fiber. For g<=G, every nontrivial odd core is at most 2G-3; allowing core one leaves at most G possibilities. This deliberately loose bound is valid also at G=1 and G=2. The restriction u>=2 removes only core one's zero-height occurrence.

The independent physical-word implementation checks the full profile tau(u), recurrence, parity, positivity, both equality classifications, and the unit fiber through u=65536. These are controls of the displayed induction, not its replacement.

## 2. J dyadic collapse

The literal J matrix gives g(a)g(b)+ab-1 when both heights are zero, and g(a)g(b)+2^k ab when their sum k is positive. The converse correctly reserves one positive height for each core equal to one. In particular, cores (1,1) require k>=2; a single singleton core requires k>=1; both nonsingleton cores permit k=0.

The statement is about an exact value-set description and legal index assignments, not uniqueness of core pairs or all-target coverage. The independent checker reconstructs every distribution of total height 0 through 12 for every pair of positive odd cores below 32, giving 3328 converse/boundary controls.

## 3. Exact divisor inverses

For c=3 and c=5, all sixteen matrix coefficients satisfy

    c F_c(u,v) = (c x(u)-p)(c x(v)-q) + (c-1)pq.

The three source-parity cases in the proof follow exactly. Source-sum membership is retained, so the criterion is stronger than factor congruences. The statements that X_e contains odd integers and X_o even integers are membership restrictions, not assertions that every integer of the stated parity occurs.

The source-sum inverse bounds have the correct orientations and rounding. Even u gives x in [u+1,3u/2], hence ceil(2x/3)<=u<=x-1. Odd u gives x in [(3u+3)/2,2u], hence ceil(x/2)<=u<=floor((2x-3)/3). Both domains require the actual equality x=u+g(u). The implementation enforces parity and u>=2; odd source candidates then start at three. Every fiber is finite, without a separate source-height cutoff.

The author inverse enumerates all ordered positive divisor pairs, including square factors only once, and handles both index orders of the mixed-parity case. Its final direct equality check is valid. In particular, it does not discard witnesses by imposing u<=v.

The independent inverse instead ranges over the four choices (p,q), factors c n-(c-1)pq, and tests the recovered source sums against an independently generated complete fiber bank. It agrees with complete direct physical searches for **every n from 1 through 5000**, including absent and even targets, for both forms. The direct domain is complete because these c-forms dominate uv: with x,y>=3 and c>=3,

    cxy-qx-py+pq >= 3xy-x-y >= xy >= uv.

The controls find 786 represented targets and 3611 ordered witnesses for JJZJ, and 474 targets and 1921 ordered witnesses for JJZZJ. A separate unrestricted source-bank comparison checks every source-sum fiber through x=4096.

## 4. Prime hostile and its physical survivor

The three integers 2135743, 6407227, and 10678711 have no prime divisor up to their respective square roots. The independent sieve explicitly constructs all 461 primes up to 3267, the largest required square-root bound. This is an exact primality certificate by exhaustive possible divisors, not a probabilistic test.

For a prime target, either even-source case is impossible by its displayed nontrivial factorization. In the odd/odd case, each factor c x-1 is at least 6c-1. Thus primality of c n-(c-1) excludes that case at every height. This proves the two-template absence for the hostile without relying on a finite source inventory.

The positive J witness is independently evaluated and compiled from actual binary-source words. Its precursor profile is

    (D,E,U,d) = (1051865,159305,1243183,159305),

so D-E+U=2135743. The literal chronological closing word is recorded in the independent receipt. Its physical meaning follows from the separately audited compiler; this audit does not claim a new literal Hamilton-path enumeration for that 23-vertex graph. The witness is above the occupancy threshold and shows that the two-template obstruction is not a full-menu gap.

## 5. Scope of the factorization classification

Every nonempty template is homogeneous in (A,p,u) and (B,q,v). The exact four-by-four matrices have zero constant row and column, as independently checked. For a separate product of linear forms after subtraction of lambda*pq, its three-by-three matrix must have rank one. Every template has a nonzero A*B pivot. That pivot uniquely determines the adjustable parity/parity entry; the independent checker then checks **every** two-by-two minor, not only the necessary core determinant.

The eleven listed determinants are correct. The even/even determinant-zero candidates are JJ, JJZJ, JJZZJ, and JJZJJ, but the other minors exclude JJ and JJZJJ. Exactly JJZJ and JJZZJ survive, with parity-product corrections 2/3 and 4/5 respectively. The norm-style identity has exactly the stated residual determinant coupling and correct cross-coefficient orientation.

This is a classification in the declared formal polynomial coordinates. It does not exclude identities using special facts about actual g values, coupled index transformations, nonlinear methods, or piecewise recipes outside that formal factorization class. No arithmetic norm theorem or surjectivity claim is used.

## 6. Counting barrier and quantifiers

The inherited product lower bound Psi_W>=uv/2 applies to every one of the twelve actual-source templates. As an independent coefficient control, the checker reconstructs the positive-coordinate ray matrices and verifies a minimum product coefficient of at least two for all twelve. Together with u=2^r a and v=2^s b, a value at most N forces r+s<=floor(log2(2N)).

There are at most G^2 core pairs and (L+1)(L+2)/2 height pairs, so the stated bound 12 G^2 (L+1)(L+2)/2 counts all possible presentations and therefore bounds distinct values. Removing illegal singleton occurrences or duplicate values can only improve it. J has only one potential value at each total height for each core pair, giving G^2(L+1).

For fixed G this is O_G(log^2 N), which cannot contain a tail of odd integers. At least floor(N/4) odd integers lie in (N/2,N], yielding the stated necessary growth scale G=Omega(sqrt(N)/log N) for interval coverage by witnesses with **both** g coordinates bounded by G. This is not a claim about every target, every representation, or a construction that bounds only one coordinate. Sixteen finite counting controls verify the inequalities independently but are not the proof of their all-order scope.

## Reproduction and release limitation

From this audit directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

The audit checker reads its copied, hash-frozen dependencies and is relocation-safe. The original author verifier instead reads a scratch-specific sibling path for MENUS.json. That is not a mathematical defect, but publication should preserve the original verifier and separately disclose/review a small dependency-path relocation change, or provide an explicitly compatible layout. No claim that the original author verifier already works in an arbitrary release layout is made here.

Universal twelve-template occupancy remains **OPEN**. No literature novelty, new coverage census, or universal representation theorem is certified.
