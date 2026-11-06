# Independent review of exact prefix formula and dyadic gaps

Status: PASS. No mathematical correction required.

Reviewed `session204-multiport-inverse/RESULT.md` and `verify.py` against the independently audited physical four-register companion theorem.

## Exact algebra

Substitution into the two physical quotient transfers gives X'=2X-36b and Y'=13X/6-36-42b for X=D+u, Y=D+v. The initial X,Y values are 540,546. Consequently X_j=36K_j with K_j=15*2^j minus the ordinary binary prefix value.

Digit one changes u by X-Y-36. Initially this is -42. After a previous digit b, substituting the X,Y expressions gives -6K+6b. Digit zero leaves u unchanged. Therefore u=-23-6T and the stated formulas for D and B follow by induction. They are exact all-word identities, not extrapolated regressions.

## Bounds and endpoints

For every j>=0, 14*2^j+1<=K_j<=15*2^j. Thus each summand in T is nonnegative. For n>=1, the upper estimate is 7+15*(2^(n-1)-1)=15*2^(n-1)-8, including n=1, where it equals seven.

These inequalities yield exactly the advertised closed length-n bands:

    504*2^n+59 <= D <= 585*2^n-25,
    168*2^n+19 <= B <= 195*2^n-9.

The bands are valid but not claimed sharp. Their endpoints increase with n. All positive lengths at most n lie at or below the upper endpoint for n. All lengths at least n+1 lie at or above the lower endpoint for n+1. Hence the stated OPEN gaps exclude every positive length, including both adjacent lengths. The empty-word values D=563 and B=187 lie below the first gap and hence below all later gaps.

The B-gap length is 141*2^n+28 and its endpoints tend to infinity. Thus the gaps really contain arbitrarily long sequences of consecutive odd targets missing from the entire unchanged family. This is not merely failure of one particular binary interpretation of the words.

## Verifier and scope

The checker uses exact SymPy arithmetic and explicit exceptions, verifies five symbolic identities, and enumerates exactly 8191 words of lengths zero through twelve. Normal and optimized runs pass. The finite regression is appropriately separated from the induction proof. Physical realization is explicitly inherited from the companion certificate.

The conclusion excludes a tail only for this fixed seed, two fixed digit operations, and fixed closing context. Changing those ingredients or combining other families is outside its scope. Reinterpreting or selecting words differently cannot change their physical count values; this is the appropriate sense of the claimed reindexing obstruction. An arbitrary numerical relabeling of the outputs would not itself realize those relabeled numbers as Hamilton-path counts.

No conclusion about the full tournament spectrum, F(8,137), or all finite-control constructions follows.
