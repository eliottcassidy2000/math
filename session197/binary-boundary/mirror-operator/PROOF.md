# The full binary-mirror bin operator, including the deletion boundary

Session197. This is an all-order identity for a named physical context. It is not a total-positivity claim or a classification of arbitrary substitution contexts. The positive-length scalar and its exceptional empty boundary were proved in the reviewed binary-boundary packet. The new statement here retains the entire labelled-bin polynomial and gives a concrete ratio-collision obstruction to inferring primitive content from cancellation.

## 1. Definitions and exact operator

For a tournament U, let B_U(x) be the polynomial whose value at a nonnegative integer x counts partitions of the vertices into x labelled bins, each bin carrying a directed Hamilton path on its assigned vertices. Empty bins are permitted, B_empty=1, and a singleton path is permitted. Thus B_U(1)=H(U), B_U(0)=0 for nonempty U, and B_U(-x)=(-1)^|U| B_U(x).

Let S be a nonempty tournament with marked vertex a. Form O=S and I=S^op on disjoint corresponding vertex sets. Orient every cross edge I->O except the single edge O_a->I_a. Insert any tournament T, including the empty tournament, with I->T->O. Denote this actual tournament by F_T. Define the endpoint polynomial

    E_S(x) = sum over nonempty directed paths P ending at a of B_(S-P)(x).

Then, with f=B_T,

    B_(F_T)(x) = B_S(x)^2 f(x)
                 + x E_S(x)^2 f(x+1)
                 - x E_S(-x)^2 f(x-1).                         (1)

In particular this is an exact linear operator on bin polynomials, not just an equality at x=1. Since transitive tournaments have B_Tr_m=x^m, including m=0, the formula also specifies its unique linear extension to all polynomials; such virtual inputs need not be physical tournaments.

### Proof by odd-cycle collections

Use the inherited odd-cycle collection identity

    B_U(x) = sum over collections C of vertex-disjoint directed odd cycles
             (2x)^|C| x^(|U|-sum cycle lengths).

Cycles are unrooted directed cycles, not vertex sets or based words. In F_T, every cycle meeting more than one of I,T,O uses the unique reverse cross edge O_a->I_a. Hence a collection contains at most one such cycle. Collections without it are internal to the three modules and contribute B_S(x)^2 f(x).

Anchor a crossing cycle at that reverse edge. It consists uniquely of a nonempty I-path starting at I_a, a T-path R which may be empty, and a nonempty O-path P ending at O_a. The I-path corresponds by reversal to a path Q in S ending at a. There cannot be two T-runs, since a second traversal from O back towards I would reuse the sole reverse cross edge. Let ell=|P|+|Q|. For fixed P,Q, sum over R, including its empty choice, and over internal odd-cycle collections on the remaining T vertices.

The bin convolution identity is

    B_T(x+y) = sum over subsets U of V(T) B_U(y) B_(T-U)(x).

At y=1, B_U(1) counts precisely the paths on U; at y=-1 it counts them with sign (-1)^|U|. Consequently the parity filter for the crossing odd cycle, including its factor 2x, is

    x [ f(x+1) - (-1)^ell f(x-1) ].

Every odd-total-length choice gives one actual directed cycle. The only invalid empty-R boundary would have both side paths singleton, producing a formal length-two walk across the reversed edge. Its length is even, so its contribution is zero in this parity filter. For all other empty-R choices the closing I->O edge is present.

The residual O- and I-side collections give E_S(x)^2 for the positive term. For the signed term, residual parity gives

    sum_P (-1)^|P| B_(S-P)(x) = (-1)^|S| E_S(-x).

The two sides square this expression, canceling the global sign. This proves (1), including the empty input, without extrapolating from positive clone lengths. QED.

## 2. Scalar boundary and a positive difference form

Write h=H(S), q=E_S(1), and e=H_start_at_a(S). The separately proved endpoint reciprocity identity gives

    E_S(-1)=(-1)^(|S|+1)e.

Evaluating (1) at1 yields, for every input including empty,

    H(F_T)=h^2 H(T)+q^2 B_T(2)-e^2 B_T(0).                    (2)

For nonempty T the last term vanishes. For empty T the exact value is h^2+q^2-e^2. In particular the formal base-zero term is real and cannot be discarded when deleting the input module.

All coefficients of B_S and E_S are nonnegative integers, by their odd-cycle descriptions. Write E_S(x)^2=U(x)+V(x), where U is its even-degree part and V its odd-degree part. With

    delta f=f(x+1)-f(x-1),
    Delta f=f(x+1)+f(x-1)-2f(x),

operator (1) becomes

    A_S f = P f + Q delta f + R Delta f,
    P=B_S^2+2xV=B_(F_empty), Q=xU, R=xV.                    (3)

Each of P,Q,R has nonnegative integer coefficients. This supplies explicit coefficients for this mirror subfamily of the previously reviewed positive context calculus. It does not claim that arbitrary products have totally nonnegative matrices or Hankel minors.

For example S=C3 has B_S=x^3+2x and E_S=x^2+x+1, giving

    P=x^6+8x^4+8x^2,
    Q=x^5+3x^3+x,
    R=2x^4+2x^2.

The singleton side gives A=x^2 I+x delta and B_(F_empty)=x^2, checking the smallest deletion boundary.

## 3. Explicit family polynomials

For k>=1, let S_k have vertices0,...,k+1, every edge forward except k+1->0, marked at0. Its only cycles use that backward edge. Selecting an odd number of the k internal vertices gives its odd cycles, so

    B_(S_k)(x)=x^(k+2)+x[(x+1)^k-(x-1)^k].

A path ending at0 is either the singleton0, or a path through a chosen subset of1,...,k followed by k+1,0. Its complement is transitive. Therefore

    E_(S_k)(x)=x^(k+1)+(x+1)^k.

Let Q_k=F_empty for this side. Substitution into (1), with f=1, gives

    B_(Q_k)(x)=x^(2k+4)
      +2x^(k+2)[(x+1)^(k+1)-(x-1)^(k+1)]
      +x[(x+1)^(2k+1)+(x-1)^(2k+1)-2x(x^2-1)^k].          (4)

In particular

    H(Q_k)=2(2^k+1)^2-1,
    B_(Q_k)(2)=6*9^k+24*6^k+16*4^k-8*3^k-8*2^k+2.

The reviewed boundary packet proves these actual graphs strong, self-converse and with exactly one Hamilton cycle; those facts are not inferred merely from the displayed polynomial.

## 4. A physical ratio collision

The graph Q_3 has ten vertices and

    B_(Q_3)=x^10+16x^8+64x^6+64x^4+16x^2,
    H(Q_3)=161=7*23, B_(Q_3)(2)=10304=64*161,
    B_(Q_3)(3)=216009.

The six-vertex transitive tournament has H=1, B(2)=64 and B(3)=729. Thus both have the same ratio B(2)/H=64, but their complete bin polynomials are not proportional:216009 is not161*729.

Consequently, in a scalar context whose coefficients are homogeneous linear functions of a leaf's pair (H,B(2)), replacing Tr6 by Q_3 multiplies every coefficient by161 and preserves every ratio-based cancellation condition. It does not preserve primitive integer content. For a literal elementary example, the cyclic wrapper C3[U,1,1] has H=H(U)+B_U(2), so the two inputs give65 and10465=161*65. The wrapper identity follows by its one- or two-run decomposition through U. This guardrail concerns contexts depending only on that pair; it is not an assertion that arbitrary full-profile contexts respond proportionally.

## 5. Reproduction, dependencies and scope

The standalone C++17 verifier uses no imported research code. It computes all subset Hamilton path totals, partitions each vertex set into unlabelled path blocks, then applies the falling-factorial transform to obtain the full labelled-bin polynomial. This is independent of the odd-cycle proof and of fitting finitely many clone values.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify > receipt.txt
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug > receipt-debug.txt

It checks all285 rooted labelled side tournaments through order4 plus16 deterministic order5 sides. For each it substitutes empty, K1, Tr2, Tr3 and C3:1505 physical full-polynomial checks, including endpoint signs and the scalar boundary. It also checks (4) for k=1,...,4 and the ratio collision. The proof is all-order; these finite checks are separate controls.

Named inherited sources read for this packet:

- `session197-odd-cycle-path-polynomials/UNIVERSAL_CONTEXTS.md`, SHA2568615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51: odd-cycle/bin identities and general positive context calculus. The exact closed mirror operator (1) is derived here, not claimed as a new general positivity theorem.
- `session197-odd-cycle-path-polynomials/ENDPOINT_RECIPROCITY.md`, SHA2564761f31465aabbe5e079c73131cabbece27dabca299a5900e5065baeb221a74c: endpoint sign in Section2. Both named path-polynomial sources were published at69a982bad394b2bac3cdf93f308b1458d0d5fc79 under session197/path-polynomials/.
- Reviewed binary-boundary proof, SHA256a6c358fc78204c266d05ad04d458877abe2aa980e25403740511d1af2ab0abf7, published at01177cf7bf2e95de57bc40f9f5b7686e7ee2331c under session197/binary-boundary/PROOF.md: physical scalar boundary and Q_k structure.

These references delimit inherited scope; no literature-priority claim is made.
