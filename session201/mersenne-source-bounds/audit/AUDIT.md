# Independent audit of the Mersenne bounded-source obstruction

Verdict: PASS for the corrected MERSENNE.md frozen in AUDITED-HASHES.json.
The all-height implication n=2^k-1=Psi_W(u,v), u,v>=2,
g(u),g(v)<=G, G>=1 => n<=18630 G^8 is proved by the stated argument.
This does not settle twelve-template occupancy or establish missing targets.

## Correction found and verified

The initial draft incorrectly asserted 0<=j<=1. Template J has j=-1.
The author corrected this to -1<=j<=1 before the audited freeze. Only j<=1
is used for an upper bound; no downstream bound changes. The independently
written checker enforces the corrected full range.

## All-height mathematical checks

The recurrence has g(1)=1, g(2h)=g(h),
g(2h-1)=g(h)+h-(h mod 2). Thus removing factors of two preserves g.
Strong induction gives 1<=g(n)<=n: for odd n=2h-1 with h odd,
g(h)+h-1<=2h-1; for h even, g(h)<=h/2, so g(n)<=3h/2<=2h-1.
For odd n=2h-1>=3, if h is even, g(h)>=1 gives g(n)>=h+1.
If h is odd then h>=3 and induction gives g(h)>=3, yielding
at least h+2. This proves g(n)>=(n+3)/2. The odd core therefore is
at most 2G, including core 1. The recurrence also proves
g(2^j+1)=2^j+1 by induction, starting at j=1.

All twelve matrices agree with the frozen release key 5, and both orientations
satisfy the coefficient bounds. Transposition is legitimate: swap sources and
transpose the form, without claiming every transpose is in the twelve-list.
The uniform estimates hold for the enlarged family. The fixed-source g
coefficient is positive for every actual source: q=0 gives lambda*B+c*v>0;
q=1 gives lambda*B+c*v-e>=lambda+3c-e>0, verified for every orientation.

The reviewed CORE_PROOF.md section 4 gives Psi>=uv/2, and the author's
finite cone coefficient certificate directly rechecks it for these matrices.
For the sorted valuations r>=s and r>=1, v>=2 implies n>=u>=2^r.
Consequently 2^k=n+1>2^r, giving k>r and justifying every modular reduction.
Positive C+1 in the nonempty branches permits divisor-to-magnitude bounds.
The even/even estimate totals 5+1260+84+17280=18629 after using G>=1.
The mixed estimate totals 939, safely below 940. The odd/odd bound is 66.

The empty-template formulas and constant terms are correct in both parity
branches. For s>=1 and Z=0, v=A+B+2>=4 is even, so v-1 is odd and at least
3 and cannot divide a power of two. For s=0 and Z=0, positive integer factors
(v-1)u=2^k force u to be a power of two and v=2^j+1 (j>=1), giving A=1,
B=v, and Z=2, a contradiction. This is not an illicit assumption Z!=0.
Once Z!=0, 2^r divides Z and 2^r<=abs(Z) is valid even for negative Z.
The empty-template constants 195 and 23 and the both-odd estimate are valid.

## Exact controls and reproduction

Run `python3 session201-mersenne-independent-audit/check.py` from the workspace.
Independent RESULT.json: PASS on 13,536,024 form/source checks using both
orientations, all indices 2..512, and additional a*2^r with a odd <=31,
1<=r<=20 (751 distinct sources). It checks every represented Mersenne target
in that universe against its branch bounds, all physical product bounds,
source lemmas through 65536, and g(2^j+1) for j<=100. There are 820 even/even
and 3168 even/odd zero-constant non-Mersenne controls (orientation duplicates
included). These controls confirm that Z=0 cannot simply be excluded for
arbitrary targets. Exact finite controls support, but do not replace, the proof.

The author's verifier was inspected and replayed with bytecode disabled.
The replay agrees byte-for-byte with MERSENNE-normal.json, which agrees with
MERSENNE-optimized.json. Its 189165 matrix/index tests and finite cone check
pass. The extracted TWELVE.json equals release MENUS.json key 5 exactly;
the release SHA256 matches dependency provenance. Its g-capped finite search
uses a height cap derived from the theorem and is therefore not independent
proof of the theorem; it is correctly presented as controls, not occupancy.
The independent dyadic-tower controls do not use that cap.

No author input, repository file, or external state was edited by this audit.
The dirty shared checkout was preserved. Publication remains with the parent.
