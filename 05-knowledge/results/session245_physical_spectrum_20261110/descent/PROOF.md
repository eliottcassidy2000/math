# A physical synchronized-block lift and palindrome defect compiler

Session245. **PROVED derived construction:** every fixed physical source middle word has an exact quadratic lift under synchronized outer sink blocks. The purely binary case has zero linear term exactly at palindromic middles. **INHERITED:** all-Mersenne coverage and its three-back-arc fan, not new coverage. **FINITE-EXACT:** the displayed finite scans and independently recounted graphs. **OPEN:** arbitrary-target coverage, four-block coverage, and a finite closed root-menu induction.

## Inheritance and why the route changed

The closest mechanisms are Session199 `hidden-structure/PROOF.md` (exact binary block coefficients and the open four-block coverage lemma), `physical-adjoint/PROOF.md` (complete C1 pairing and square-root inverse), and `source-chambers/PROOF.md` (physical membership, not cone membership). Session244 `dyadic/PROOF.md` supplies the hostile: the two signs allowing 2D−1 and 2D+1 at count five are not a continuation-closed menu.

A first attempt to eliminate one source from the bilinear inverse recovers the same physical saturation fiber when its coefficient vanishes. Its dense arithmetic search is not obviously faster than the inherited sparse physical-bank inverse. It is not claimed as a breakthrough here.

A second attempt sought a single Mersenne-root sign-plus seed. The standalone exploratory ray.cpp checks all tournaments admitting the fixed Hamilton order at orders 3 through 8 (2, 8, 64, 1024, 32768, 2097152 orientation words, respectively). Every tournament has some Hamilton order, so this is a complete graph universe up to relabeling, with duplicates. No root with H>=15 Mersenne, E=d, U=H+d+1 was found. This is finite evidence only; no universal root exclusion follows. The positive alternative below changes both outer runs, instead of requiring that one fixed-root extension close on itself.

The resulting Mersenne words were immediately compared to THM-4748, `01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`. They are exactly its inherited Fibonacci fans. We retain that recovery as a positive control and extend the cancellation calculation to every physical middle word. There is no claim to have newly proved Mersenne realizability.

Concept board: physical inverses; restored root menus; synchronized binary scales; the ordered pair correction C; root changes inside the middle; and exact target placement. The last remains unsolved.

## 1. A general two-scale physical lift

Use source state (D,E,U,d), seed t=(1,0,2,1), and the inherited actual graph operations

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D−E+U,d,2U−d,d),
    J(D,E,U,d)=(D−E+U,D−E,D+U,D).

Each letter adds a vertex. O adds a sink; Z adds a vertex beating only the current root; J makes that same new vertex the root. Words act chronologically. Every word from t is an actual start-one rooted tournament. A terminal J gives a strong tournament with exactly one Hamilton cycle.

Fix any word W over O,Z,J and any integer s>=1. Let its linear action on the two formal vectors be

    W(1,0,0,1)=(b_D,b_E,b_U,b_d),
    W(0,0,2,0)=(a_D,a_E,a_U,a_d).

These formal vectors are only a coefficient extraction device. Neither is declared physical. For every integer k>=0, put x=2^k. Then the actual graph word

    O^k W O^(k+s) J

has exactly

    H_W,s(k) = 2^s a_U x² + (a_D+2^s b_U)x + b_D.       (1)

**Proof.** The actual state after O^k is (1,0,2x,1), including k=0. By linearity, its state after W is b+x a. Since k+s>=1, the last O block resets E to zero, leaves D unchanged and multiplies U by 2^(k+s)=2^s x. The final J therefore has count D+U. Expanding gives (1). All graph steps are lawful, with no subtraction or unverified inverse realization.

The graph order is 2k+|W|+s+3. Equation (1) is a symbolic identity for all k, not a fitted recurrence. In particular the linear term vanishes **if and only if**

    a_D+2^s b_U=0.                                    (2)

Thus a finite word satisfying (2) supplies the infinite exact family 2^s a_U 4^k+b_D. The test of (2) is finite exact integer matrix arithmetic, while the word itself supplies the physical certificate. No positivity test on formal coefficient vectors is substituted for physical construction.

## 2. Exact palindrome cancellation for retained-root middles

Now take W over O,Z only, of length l, and let epsilon_i=1 for Z, 0 for O, in chronological order. Define the inherited coefficients

    A=sum epsilon_i 2^i,
    R=sum epsilon_i 2^(l−1−i),
    C=sum_(j<i) epsilon_j epsilon_i 2^(i−j−1)
       +sum_(i>=1) epsilon_(i−1) epsilon_i.

The empty middle has l=A=R=C=0. Substitution into the inherited block formula yields

    D_after_W=1+2Ax−C,
    U_after_W=2^(l+1)x−R.

With s=1, (1) becomes

    H_W,1(k)=2^(2k+l+2)+2(A−R)2^k+1−C.                (3)

The cross term vanishes exactly when A=R. Binary positional uniqueness says this holds exactly when W reads identically in reverse, including leading zeros and the empty word. Hence every palindromic middle gives

    H=2^(2k+l+2)+1−C(W),  k>=0,                       (4)

on 2k+l+4 vertices. C is even: adjacent pair contributions are 2, and every other pair contributes an even power of two. Thus the output is odd, as required. This is a derived specialization of Session199's binary compiler; no priority is claimed for that compiler or the one-apex counting law.

Equation (4) isolates a positive cancellation mechanism beyond the three-arc fan: any symmetric set of inner backward arcs cancels its entire scale-linear contribution, leaving a fixed exact pair-interaction defect. The adjacency correction must be retained: replacing it with an unweighted distance count gives the wrong defect already for ZZ.

Examples of actual strong C1 families, all k>=0:

    W=ZZ:       H=16·4^k−1.
    W=ZOZ:      H=32·4^k−1.
    W=ZZZ:      H=32·4^k−5.
    W=ZOZOZ:    H=128·4^k−11.
    W=ZZOZZ:    H=128·4^k−21.
    W=ZOZZZOZ:  H=512·4^k−65.

The first two are precisely the already proved all-Mersenne fan. For m>=4 use ZZ when m is even and ZOZ when m is odd. The graph has m+2 vertices. No new Mersenne theorem is asserted.

### Half-word norm formulas and a complete target inverse

The mirror construction gives a compact exact defect compiler. For any binary word V let R be its reversed binary integer, C=C(V), and e be its final bit (e=0 for the empty word). Then

    C(V reverse(V))   =2C+R²+e,
    C(V O reverse(V)) =2C+2R²,
    C(V Z reverse(V)) =2C+2R²+2R+2e.                   (5)

Within the two halves the pair contributions total2C, since reversal preserves distances. Across the two halves they factor as R² or2R². In the even case the adjacent central pair adds e; with central O there is no extra contribution; with central Z its two interactions contribute2R plus the two adjacency corrections2e. This proves (5), including empty halves. These are exact weighted bit norms, not a declaration that arbitrary integers are represented by those norms.

For every binary word of length l>=1,

    0<=C<=2^l−2.                                      (6)

Every summand is nonnegative, so the all-Z word maximizes C. Its first value is C_1=0; adjoining a Z increases C_l by2^l, proving the upper bound inductively. A palindrome ray with count h>=9 and m=2k+l+2 consequently satisfies

    h−1 <=2^m <=4(h−3)/3.                             (7)

For l>=1 this follows from h>=2^m−2^l+3 and2^l<=2^(m−2). For the empty middle h=2^m+1 and m>=4 when h>=9, so (7) also holds. The exceptional small value5 is generated by the empty middle at k=0. No other output below9 occurs.

The interval in (7) has ratio strictly less than2. Thus h uniquely determines m if a qualifying power exists. Set C_target=2^m+1−h. Enumerate l=0,...,m−2 of the same parity as m, and all palindromes of length l, testing C(W)=C_target. On acceptance k=(m−l−2)/2 and the original word is a literal graph witness. If all tests fail, h is absent from the **entire synchronized palindrome family**, at every graph order. This is a complete inverse for this family only. The code's `--target h` mode implements it.

There are O(2^(m/2)) half-words to test. Since (7) bounds2^m by O(h), this is O(sqrt(h)) candidates, each with a polynomial-in-log(h) exact arithmetic cost. No input-bit-polynomial bound is asserted. This efficient finite test does not become an all-target existence theorem.

## 3. Middle root changes genuinely extend the cancellation inventory

The general criterion (2) is not restricted to palindromic binary words. The following exact coefficient triples (quadratic,linear,constant) are obtained by the displayed physical matrix rules:

    W=ZZOJOOZ: (340,0,−35),
    W=ZZJOOOZ: (436,0,−37),
    W=ZZJJOOZ: (420,0,−43).

Thus for every k>=0 these give strong C1 counts

    340·4^k−35, 436·4^k−37, 420·4^k−43.

Their middle length is seven and their graph order is 2k+11. No primality or factorization assumption is required. The first and third families may have obvious common factors; the theorem is physical exact realization, not a new prime-supply guarantee.

In contrast to the Mersenne construction these words include root changes inside W. They demonstrate a concrete extension of the symmetric one-apex mechanism. They do not demonstrate universal representability of the constant term or leading coefficient. The auxiliary exact search at length at most12 obtains147 distinct constant terms under (2), but that finite value does not bound the all-word inventory.

## 4. What still prevents an all-target conclusion

The family is not a freely adjustable anchor: W determines the leading coefficient, constant term and graph order together. Even in the palindrome case C is a coupled quadratic bit statistic. All even integers are not automatically values of C, and a specified target need not have a compatible exponent parity. For instance offset9 is absent from the tested palindrome inventory through length12; this is deliberately stated only as a bounded check.

Equation (1) supplies a finite certificate for each entire quadratic ray. It does not supply a finite set of rays covering all odd targets, nor a failure-preserving inverse from an arbitrary missing target. The Mersenne recovery therefore answers a canonical hostile lane but cannot be used to announce Chan–Pak or the four-block coverage lemma.

A precise next positive obligation is to find a physically closed family of middle words whose triples (2^s a_U,a_D+2^s b_U,b_D) solve the target equation for every remaining primitive target, with well-founded parameters. The present computation exposes those three coordinates without discarding their common word certificate.

## 5. Verification and reproduction

Run `python3 verify.py` and `python3 -O verify.py`; compare CHECK.json and CHECK_OPT.json byte-for-byte. The verifier is standalone and uses explicit exceptions. `python3 verify.py --fail` must fail.

The receipt checks8191 binary words through length12 against (3), every palindrome iff test, 1093 general O/Z/J middles through length6 at three s-values and four k-values (13116 substitutions), and sixteen literal graph counts by subset dynamic programming. The latter also independently count Hamilton cycles, obtaining exactly one. All three half-word norm identities are checked on511 half-words; the target-first inverse is compared to an independent forward image for every odd target through5001, with accepted recipes reapplied. The large-parameter substitutions for the inherited Mersenne family are checks, not the proof. The proof of every all-k assertion is Section1's coefficient identity.

To get an explicit graph use, for example, `python3 verify.py --middle ZZOJOOZ --k 3`. It emits row masks, marked root, exact count and construction word. Read each row mask as the outgoing-neighbor set. No tool access to another repository file is needed.

The auxiliary searches ray.cpp, template.py, offsets.py, and general_offsets.py are exploratory records, not required to trust the proof or replay the final verifier. The proof-facing compiler and receipt are self-contained.
