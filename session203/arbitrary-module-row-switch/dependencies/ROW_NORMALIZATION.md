# A physical common-base family with two endpoint choices

**PROVED; independently audited with exact author controls.** This is a positive common-base family, not a proof of the staggered profile hypothesis or an eventual threshold. The O/Z operations and underlying rooted formulas are inherited; the contribution here is an explicit pair of different physical sources with the same padding state and different endpoint counts.

## 1. State and actual start-preserving padding

Let U have at least two vertices and root v. Write h=H(U), s=start_v(U), e=end_v(U), L for ordered root-ending two-path covers and R for ordered root-starting covers, with the other path allowed empty. Put T=U-v, d=H(T), beta=B_T(2), and

    B=h+L-s,    V=R+beta.

Operation O adds a global sink and retains v. Operation Z adds a vertex z that beats only v and loses to every old nonroot vertex, again retaining v. Both preserve s: a Hamilton path starting at v must put the new vertex last, and every old such path extends uniquely. The assumption |U|>=2 prevents the rooted-singleton Z exception.

The exact state updates are

    O: (B,e,V,d) -> (B,0,2V,d),
    Z: (B,e,V,d) -> (B+V-e,d,2V-d,d).           (1)

For O, h and L are unchanged, while R and beta double. For Z, the inherited rooted path decomposition gives h'=h+R-e, e'=d and R'=2R-d. A root-ending cover either puts z last in its ordinary path, contributing L, or puts z immediately before v in its designated path, contributing the arbitrary two-bin covers beta of T. Thus L'=L+beta. Deleting v makes z a sink of T, so beta'=2beta and d'=d. These facts prove(1).

All quantities come from actual paths and covers. Equality of these four coordinates is therefore sufficient to preserve the common B under every identical O/Z suffix, while the separately retained start count s remains unchanged.

## 2. Two strong sources with the same padding state

The two six-vertex tournaments, both rooted at vertex0, have outgoing rows

    U6: (10,20,33,38,13,19),
    U7: (10,12,17,52,35,7).

Their actual triples (h,s,L) are (45,6,72) and (43,7,75). In both cases,

    (B,e,V,d)=(111,9,161,13).                   (2)

Each source is strong. For example the directed spanning cycles

    U6: 0,1,4,3,2,5,0,
    U7: 0,1,3,2,4,5,0

witness this directly.

For every O/Z word w, applying its literal births to these two sources yields equal B values, with preserved endpoint choices6 and7. The separately audited four-core readout therefore gives actual child counts

    3B_w+12 and 3B_w+14.                       (3)

This is a true common-base construction, not a subtraction of counts or a choice of unrelated profile coordinates. The two endpoint choices are only6 and7; no wider endpoint interval or closure under combining choices is asserted.

For the all-Z word of length k>=1, (1) gives

    B_k=148*2^k-33,
    H(C(U6_k))=444*2^k-87,
    H(C(U7_k))=444*2^k-85.

At k=0 the counts are345 and347; the displayed exponential formula is not extended across that boundary. Every Z-only source remains strong because the old strong source reaches z and z reaches the old root. Its four-core child is strong as well. General words may add neutral sink components.

## 3. Exact controls and first missing certificates

The standalone verify_rows.cpp builds every O/Z suffix of length at most4 on both compact sources, recomputes all needed path and cover coordinates directly, and counts both literal four-core children. All 31 words pass. The same checker also tests every labelled tournament of orders 1 through 4 with root fixed at vertex 0, 75 inputs, including the singleton boundary, for the arbitrary amalgamation formulas below. In total 424 literal graph calculations pass in optimized and debug builds. The common-base identity for arbitrary word length follows from(1), not from finite extrapolation.

The same start-preserving padding supplies two concrete profile certificates at the first base omitted from the published conditional compiler, B=999:

- Residue1 uses j=5, s=16. Start from root6 in(2,68,9,19,39,15,61), then append O,Z. The final source has (h,s,L)=(485,16,520), so its base is989=999-10, and the four-core child has H=2999.
- Residue2 uses j=0, s=2. Start from root3 in(2,20,9,3,13), then append O,O,Z,O,Z. The final source has (h,s,L)=(341,2,660), so its base is999, and the child has H=3001.

BASE999_WITNESSES.json contains every original/source/compiled row and the independent exact subset-DP counts. These are new physical certificates for two specified profiles, not an all-order extension of F(8,137). The published compiler remains unchanged, so its missing-certificate controls remain truthful about its frozen data.

## 4. Exact retained-profile interpretation and arbitrary right inputs

The row is not a newly invented counting state. It is exactly the standard retained-root profile (H,e,R,d) of E(U), with the OLD root still marked. Indeed E(U) is the ordered root amalgamation of the rooted edge v->a with U; that edge has d=1 and R=2. The inherited amalgamation laws give H(E(U))=B, e(E(U))=e, R(E(U))=R+beta=V and d(E(U))=d. Thus the two sources yield equal-profile edge extensions, while their marked deletion-endpoint counts start_v(E(U)-a) differ.

For any nontrivial rooted right input W, let its actual coordinates be h_W,s_W,e_W,L_W,R_W,d_W and beta_W=B_(W-root)(2). Ordered root amalgamation gives the following full update of the left row:

    B' = d*h_W +(V-d)*L_W +(e-d)*s_W
           +(B-e-V+d)*d_W,
    e' = d*e_W,
    V' = beta_W*V +d*(R_W-beta_W),
    d' = d*d_W,
    s' = s*d_W.

The first four formulas depend only on the shared row, while the final formula retains the selectable endpoint. They follow either by expanding the inherited h/R/L formulas or by the literal associativity E(U*W)=E(U)*W. The singleton right input is the identity and must use e'=e rather than the nontrivial endpoint formula.

Consequently the compact pair remains equal-row after ANY common rooted right input, and its endpoint choices become 6d_W and 7d_W. This is a genuine physical lift to arbitrary nontransitive inputs. It does not keep the endpoint gap equal to one unless d_W=1, and does not provide independent control of the common base.

## 5. Remaining positive task

The actual row collision survives arbitrary common padding but offers only two endpoint choices. A useful stronger construction would let these choices combine or generate the bounded endpoint menu required by F while keeping a controllable common base. The two-port worker is investigating a shared-module lift of the compact pair. No such lift or eventual-coverage theorem is established here.
