# The recontact menu is a rooted normalization, not a special scaffold

Let M be any nonempty tournament with y marked. Write h=H(M), g=H(M−y) (empty count 1), and w for the weighted starts at y: the number of two-slot covers with the second slot nonempty and starting at y. Add x beating only y in M, then z beating only x in M∪{x}. Call the result C, and append an ordinal sink s to get R. Now use R as the source in the same two-vertex return operation as in [LITERAL_LIFT](LITERAL_LIFT.md). Let J=H(R).

| Tags Q in R | New count |
|---|---|
| empty | J |
| s | 3J |
| x | 5J |
| y,s | 7J + 4(w−h−g) |
| x,y,s | 9J + 4(w−h−g) |

In particular, **one boundary equation w=h+g gives the whole menu**, and J becomes 3h+2g. Moreover the same signed defect shifts the last two rows, so the single recontact from the fourth to the fifth row always changes the count by exactly 2J, even off normalization. This makes the exact-target inverse explicit: pick a physical rooted seed with w=h+g and 3h+2g=p.

## Proof

H(C)=J=h+2w. If x is terminal, z must be just before it, giving h. Otherwise x→y occurs, with z just before x or at the very end, each giving w by splitting M at the rooted slot. For the weighted starts at y in C, the two tail vertices can be appended together as zx at either slot end, or separately at the ends in either assignment. The sole invalid case is x appended directly to the singleton marked slot y (g covers), giving W_C(y)=4w−g. The full ends at y are precisely full M−y followed by zxy, giving E_C(y)=g. Ends at x must have terminal zx and contribute E_C(x)=h.

Also W_C(x)=h+2w=J: if the rooted slot is singleton x, the other slot ends at the sink z after a full M. Otherwise the rooted slot continues with y, and z goes at either slot end, giving w each. For the joint tag set {x,y}, the U-count requires X to end at one tag and Z to start at the other. X cannot end at y while Z starts at x: in that case Z=x alone, making the sink z in X incompatible with ending at y. If X ends at x and Z starts at y, z can come just before x or finish Z, each contributing w. Thus U_C({x,y})=2w.

When the trailing sink s is untagged, W doubles and U is unchanged on passing from C to R. When it is tagged and the old tag set in C has data (W_C,U_C,E_C), one has W_R=2W_C+J and U_R=U_C+E_C+W_C. Substitute the above spellings into H(new)=J+2W_R−U_R. For singleton old tag y, U_C=0, yielding the fourth row. For both tags x,y, use W_C=J+4w−g and E_C=h+g, giving the fifth row. The other rows use the empty tag, the marked sink alone, or W_C(x)=J.

The second literal scaffold is a substitution into the marked five-vertex seed M=(8,5,17,6,11), y=1, retaining the x,z,s construction under that substitution. In contrast, choosing M with y the penultimate singleton in an ordinal decomposition gives h=g and w=2h, which supplies only J=5h. The first nontrivial inverse is not arithmetic representability as 3h+2g: the two counts must be realized by the *same* deletion pair with the required weighted start. Nothing here asserts such a normalized seed for every prime.

## Order-independent compression of the normalization

Set D=M−y and let A=N⁺(y) in D; use the cover quantities W,U,E on D with this mark set (E counts full paths ending in A). Then w=g+W, since the marked slot is either singleton y or continues into A. Likewise h=g−E+W−U, by separating paths ending y from paths with a tagged successor. Thus normalization means precisely

    U+E=g,    h=W.

In particular A is nonempty and lies in **one** strong component F of D. To prove this, suppose F and L are first and last tagged components and are distinct. Split full paths of F just before a tag, and full paths of L just after a tag, placing the F tag at the start of Z and the L tag at the end of X. In every intervening SCC independently split a full path into two slots; merge in ordinal order and attach outer full paths. These choices give U≥g times 2 to the power of the number of intervening components. If there is an intervening component this contradicts U≤g. If F,L are adjacent, there is an additional class of covers where X ends at a tag in F and Z starts at a tag in L. Such full paths in F and L exist: any vertex of a strong tournament with at least three vertices can be a start or end by cutting a directed Hamilton cycle (Camion's theorem); singleton cases are immediate. The additional class is disjoint and again makes U>g.

Write D=P→F→S in its unique SCC order (empty blocks allowed). With A inside F the test on (F,A) is U_F+E_F=H(F) if S is empty, and U_F=H(F) otherwise. Indeed U on D scales by H(P)H(S), while E is zero unless S is empty. The read for normalized seeds becomes

    J=H(P)[3W_F P₂(S)+2H(F)H(S)].

The outer prefix factor must be 1 for J prime, since the bracket is greater than 1. This is an exact single-component selection law, with a possible two-coordinate suffix, not a proof of the universal profile supply.
