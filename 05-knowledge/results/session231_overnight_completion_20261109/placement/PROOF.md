# Exact interval placement and the distinguished-root graft collapse

Status: PROVED, independently audited with the k=0 singleton boundary clarified; finite literal controls below. These are restricted-language placement/obstruction theorems, not Chan–Pak completeness. No literature-priority claim.

## Inheritance and live concepts

Inherited: Session229 dynamic/PROOF.md (seven-coordinate composition and binary-context scalar formula), inverse/PROOF.md (all-gadget scalar graft row and its independently audited strict descent), and Session231 SHORT_INTERVALS.md (the twelve-bit moving core). Hostile: scalar density or gcd does not imply occupancy. Corrected near miss: an O(log N)-wide repair is not automatically at the needed O(log N) error. Least-used sidecar: root-ending count E equals deletion count d for the star root. Board: unit windows, context placement, full rooted states, graft collapse, logarithmic repair budget.

Write the inherited deletion interval state as (h,q,A,G,T_j,S,E), where h=S=E=45 and T_j=T_0+2j, 0<=j<=k. Its star closure has count D_k-2j. Put R_k=4096^k and J_k=(R_k-1)/4095. Then

 q=11520R_k,
 A=7516R_k+343260J_k+45,
 G=2353+9532800J_k+5760R_k,
 D_k=(245051264R_k+1602510k+1118023)/24843,
 T_0=2A+45-D_k.

The graphs have 12k+13 vertices, not logarithmic width at a freely chosen center.

## 1. Arbitrary ordinal contexts

Let the left and right nonempty context states be

 (p,r,B,J,U,V,W), (s,t,F,K,Z,M,N).

Apply the seven-state composition formula twice. The resulting star count is

 (Bq+pA)(2t-K)+ph(2F+s-Z)
 -sBG-sUh-spT_j-sWS-pEM.                         (1)

In particular its interval spacing is 2ps. Empty contexts are omitted and contribute Hamilton multiplier 1; no false seven-coordinate identity element is used. Unit spacing for this replacement is equivalent to p=s=1. A tournament has exactly one Hamilton path iff it is transitive: transitive tournaments plainly do; a nontrivial strong component has a directed cycle spanning that component, giving at least three Hamilton paths, and SCC counts multiply. For a nontrivial window (k>=1), only transitive ordinal contexts preserve its consecutive-odd spacing. At k=0 the window is a singleton, so this necessity assertion does not apply. This does not rule out different constructions filling holes between the multiplied outputs.

## 2. Explicit binary target selector and finite search

For prefix binary integer P (empty if P=0) and b-bit suffix integer R, let rev_b(R) be digit reversal, e_b(R) its first digit, and

 f(word)=sum over selected positions i<j of (2^(j-i-1)+[j=i+1]).

Empty words have f=rev=e=0. The exact upper endpoint after padding is

 B_k(P,R;b)=(qP+A)(2^(b+1)-rev_b(R))
  +45(2R+1-f(R))-PG-45f(P)-T_0
  -45(P mod 2)-45e_b(R).                           (2)

The realizable interval is exactly B_k(P,R;b)-2j, 0<=j<=k. For odd N, the exact membership test is

                 0 <= B_k(P,R;b)-N <= 2k.          (3)

Parity is automatic for physical contexts; for an implementation it can also be checked explicitly. The context begins at the first selected prefix bit, or at the already selected first vertex of the interval body, so the closure is strong. Leading unselected prefix vertices form initial singleton SCCs and can be deleted without changing the Hamilton count. They do not enlarge the scalar image.

This selector has an explicit LOGARITHMIC finite search bound. In every physical deletion state T<=A, since the covers counted by T are a subset of those counted by A. Hence every star count is at least A_total+h. For a prefix of length a and suffix of length b,

 A_total=2^b(qP+A)+45R.

If a>=1 then P>=2^(a-1), so A_total>=5760*4096^k*2^(a+b). If a=0, the same inequality follows from A>=7516*4096^k. Thus every successful context with m=a+b obeys

              5760*4096^k*2^m <= N-45.            (4)

One may enumerate all k,m satisfying (4), all splits a+b=m, all prefixes starting 1 when a>0, and all suffixes, then apply (3). The algorithm is complete for this entire arbitrary-singleton-context language, returns physical witnesses, and has O(N log N) candidate contexts (a loose bound), unlike an unspecified physical-state oracle. It still need not accept every large odd N.

A sharper separate k bound comes from monotonicity under vertex addition. Given any directed Hamilton path, insert the new vertex immediately before the first vertex it beats, or at the end if it beats none. This gives an injection from old paths to new paths. Therefore padding cannot decrease a member's count, and N>=D_k-2k. In particular

             (245051264/24843)*4096^k <= N.         (5)

The repaired half-count error is at most k, so its necessary budget is at most (log_2 N)/12+O(1); the count error is at most (log_2 N)/6+O(1). A previous O(log N) approximation theorem cannot be imported merely by matching big-O classes: its constant, direction, actual interface, and physical center must satisfy (3).

## 3. The new rooted-state bridge

Root a star-closure tournament at its closing root r. Its rooted four-state in Session229 inverse conventions is

              x_j=(D_k-2j, h, 2A+h, h).            (6)

Indeed deletion of r leaves L=>v, with h paths. A full path ending at r must end with v,r, since v beats only r. Thus E=h=d. For the two-bin count U, a root-starting bin either consists only of r (h choices in the other bin), or first enters the cut in L. The vertex v can finish either bin, giving 2A further choices. Equivalently use the inherited deletion-SCC recursion and the trailing unselected singleton.

Consequently the interval becomes a genuine PURE-D four-state fiber, not merely a family of scalar counts. This permits an exact test of unrestricted graft amplification.

## 4. All gadgets at this same root collapse

For any nonempty tournament gadget B with H Hamilton paths and Q ordered two-bin covers, graft B=>r and old vertices other than r=>B. The inherited scalar graft row yields, from any state (D,h,V+h,h),

 D'=HD+(Q-H)V, E'=d'=Hh, U'=QV+Hh.                 (7)

The latter three identities can be checked without the scalar row: the deletion is an ordinal sum; all of B beats r, so every deletion Hamilton path can end at r; and the active root-starting cover count V=U-d multiplies by Q when the inactive final block B is appended. Our original source has V=2A.

For ANY finite sequence of arbitrary nonempty gadgets at this same distinguished root, write H_* and Q_* for the products of their H and Q. Induction in (7) gives

 D_out=H_*D+(Q_*-H_*)V,
 E_out=d_out=H_*h, U_out=Q_*V+H_*h.                 (8)

The sequence is equivalent at this full rooted-state level to one gadget formed by their ordinal sum. Thus the source interval becomes an arithmetic progression of step 2H_*.

For k>=1, to retain the actual consecutive-odd window one must have H_*=1. When k=0 every multiplier preserves singleton consecutiveness; the algebraic formulas still hold. Every gadget is then transitive, Q_*=2^m for total gadget order m, and

              D_out=D+2A(2^m-1).                  (9)

But (9) is exactly the star count obtained by appending m unselected singleton vertices to the original deletion body: A becomes 2^m A and T stays fixed. Hence ARBITRARY ITERATED SAME-ROOT GRAFTS add no unit-window placements beyond singleton padding already present in Section 2.

This obstruction is sharply restricted. Root changes, other substitution sites, non-ordinal changes, or new auxiliary fibers may escape it. Grafts with H_*>1 can still contribute useful arithmetic progressions; this theorem does not forbid interleaving progressions from different sources. A useful missing lemma is therefore either a physical baseline selector satisfying (3) for every sufficiently large odd N, or a root-change/interface construction that restores unit spacing after leaving (6). Scalar-only inverse descent does neither.

## 5. Verification and scope

verify.py independently counts the actual k=0 star tournament (13 vertices), using literal endpoint subset DP rather than seven-state composition. It verifies its four-state (9909,45,15167,45), then literally grafts singleton, transitive pair, and directed triangle gadgets and checks all four coordinates against (7). The triangle is a hostile: it scales the gap by 3, not 1. The script also checks the all-depth collapse algebra with mixed gadget types. The inherited short-interval independent checker already covers all j through k15 and the all-order proof; this script does not claim to re-enumerate those larger graphs.

Reproduce: python3 -B verify.py > CHECK.json, repeat with -O, and compare. None of these claims closes Chan–Pak completeness or asserts an unrestricted attainable-count cutoff.

`select.py` implements the exact bounded selector (2)--(4), returning k,j,prefix,suffix and order, or no witness in this language. `check_selector.py` compares its scalar formula against independently composed seven-states in 560 contexts (k=0..3, padding<=4, j=0 or k), verifies both resource bounds, accepts the literal graft targets 9909,25031,55275, and rejects 9907. Normal and optimized receipts agree in `SELECTOR_CHECK*.json`. This cross-check is algebraic; the distinct literal path-DP controls remain in `verify.py`.
