# Independent audit: quartic repair and congruent curves

## Verdict

The rational quartic repair is correct on ND != 0. The separate printed polynomial formula (5.1) is correctly transcribed and its sixteen polynomial square identities hold, including every integer boundary specialization. This does not imply that all values are distinct, nor does it provide a magic square of distinct integer squares.

The elliptic-center bridge and finite integer catalog are correct with d > 0, N=d/4, positive rational triangles, and nonzero y in the triangle map. An arithmetic progression in eligible centers, with positivity and nine-entry distinctness, is the relevant condition. Arbitrary curve x-coordinate arithmetic progressions are insufficient.

## Source verification

Visually inspected rendered PDF pages 11 and 12 from arXiv:2508.07806v1. All factors of (5.1) agree with polynomial_certificate.py. Page 12 genuinely states cross equality before (5.5); this is not a text extraction issue. The factor in (5.5) instead equals minus m² times the diagonal difference. Source URL: https://arxiv.org/pdf/2508.07806v1 . PDF digest supplied by the source worker is separately recorded in its repair.md. Source files and page renders are audit-only, not publication artifacts.

## Independent checks

check_independent.py directly verifies the cross harmonic and diagonal norm identities, the duplication identity h²=4x(2P) modulo y²=x³−N²x, the inverse triangle-to-curve map modulo h²=a²+b², and the center-to-triangle identities. Each uses explicit exceptions rather than assert.

All assertions in the three supplied certificates were mechanically converted into explicit conditional raises, preserving their calculations while preventing python -O from deleting verification. Each transformed certificate was run in normal and -O mode and the outputs match. These runs include 14 generic and 2 exceptional rational identities, 16 printed polynomial identities, the published numerical example, 24 rational curve multiples, and the integer catalog d≤20000 cross-checked by root enumeration through h=200.

The AST-transformed runs are robustness checks, not independent mathematical implementations. Independent symbolic derivations and rendered-source comparison supply the separate audit evidence.

## Domain and interpretation checks

For c−d=u², c=h², c+d=w² and d>0, take positive u,w,h. Then a=(w+u)/2 and b=(w−u)/2 are positive, a²+b²=c, and ab/2=d/4. Conversely c=h² for a positive rational triangle of area d/4 gives c±d=(a±b)². For rational nonzero d, c=d would require a rational square root of 2, so the boundary cannot occur.

For positive integer centers, u and w have the same parity since u²+w²=2h². Thus r=(w+u)/2, s=(w−u)/2 are positive integers with r>s, d=2rs, and c=r²+s² a square. This proves the finite factor catalog is complete for each fixed positive integer d; odd d gives no centers. It is not a global exclusion over all d.

For a nonconstant AP C−e,C,C+e in C_d, the nine numbers C+i d+j e with i,j in {−1,0,1} are squares. The supplied magic arrangement has line sum 3C. Positive entries require C>d+e after choosing d,e>0. Nine distinct entries require e not in {d,2d,d/2}. Common denominator scaling of roots multiplies all entries and gaps by L² and preserves distinctness. It scales every associated triangle area by L², rather than preserving its numerical value.

Points with y=0 cannot be used in the triangle formulas; the point at infinity is also excluded. Translating a point by rational 2-torsion leaves its double unchanged and therefore cannot provide distinct centers. Point multiples forming a group progression do not generally produce a progression of centers.

## Wording corrections

The source elliptic checker comment saying common scaling 'preserves the area' should say 'preserves a common area while multiplying it by L².' Center/triangle equivalence is understood modulo interchange/sign choices of legs, not as an unqualified one-to-one map from curve points.

## Final prose and citation review

Reviewed /tmp/session57-congruent-curve/PROOF.md after it arrived. Its positivity, rational inverse, sign convention, duplicated-center guard, complete fixed-integer-gap factor proof, and no-global-exclusion boundary are correct. Independently opened Stein's cited primary page and Cain's cited arXiv HTML: Stein contains the three-square and triangle/curve propositions; Cain Section 1 explicitly records the historical congruent-number/elliptic/right-triangle connection. The proof correctly does not claim inspection of Robertson's full text.

Recomputed source PDF SHA256 and obtained 964fdfc1046e0f34e652f97ce507e803bf50581e7a6042bb0c461d7a5bd3ad66, matching the source worker's pin.

## Final producer correction addendum

Verified final check.py comment now says common scaling preserves a common area across triangles while multiplying it by L². PROOF.md explicitly states unordered positive legs modulo root signs and leg order and a surjective, not one-to-one, point-to-center map. Final producer check.py also contains a new direct-half assertion at line 57. Regenerated the explicit-exception audit copy and reran it under both normal Python and -O: both PASS with byte-identical logs, including that new identity on all 24 samples. The original producer script still uses assert, so its own -O run is not a verification. Final independently recomputed producer hashes are recorded in final_elliptic_pins.json.

Additionally, direct_half_symbolic.py proves both coordinates of the new direct-half identity symbolically modulo h²=a²+b²: x(2P)=h²/4 and y(2P)=h(a²−b²)/8. For canonical a>b these equal (c/4,huw/8). Its explicit-exception checks pass identically in normal and -O modes. Thus the identity is supported by a symbolic proof as well as 24 finite samples.

## Third-checkpoint root-document review

Reviewed /tmp/session57-root/third/README.md and SYNTHESIS.md. Quartic text correctly distinguishes the erroneous intermediate reduction from the verified polynomial theorem and does not claim distinctness at boundary specializations. Elliptic text correctly retains the fixed curve, eligible doubled x-coordinates, common scaling, positivity and collision guards, and bounded-negative scope. Requested a fixed-d quantifier clarification in SYNTHESIS Section 6; confirmed it now says 'a rational magic square with this gap d exists precisely when'. Quartic/elliptic root-document scope is approved after that correction.
