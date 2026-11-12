# Physical signed corrections are a finite memory hierarchy

**Session278. PROVED relative to the inherited physical join law.** This packet does not prove Chan–Pak, a new coverage threshold, or a faster counting algorithm. It extracts a new all-length correction structure from the S25 moving-root law, rather than re-auditing that law.

## Inheritance and novelty boundary

Input: S25, `kps_s25_euler_field_bridge_law_20261111/README.md`, Proposition B3 (itself from THM-4876/4866), at e07f05b9fee231e59ba7c474fd52c59bc63d3980. The 3-state physical law is inherited, including its physical construction and root-path interpretation. Its S25 packet was pre-audit at this source pin. Our statements are algebraic consequences of that displayed law; small physical checks independently support its application here. The new object is its **signed boundary-matching expansion, exact positive remainder, and boundary-matching-number termination criterion**. We do not claim to have invented the general Neumann-series or Bonferroni methods.

An older, nonobvious inspiration is `07-reflections/the-lever-zoo-is-one-moment-hierarchy.md` (June 30): many different correction methods can be levels of one hierarchy. `07-reflections/lrc14-sector-cover-is-a-moment-problem-macmini-S7.md` motivates replacing absolute-value tail estimates by signed certificates. We use these as idea sources, not their historical extremal assertions as theorem inputs. Here the rigorous mechanism is a **positive nilpotent Volterra operator**, not an assumed probability/moment representation.

## 1. Physical objects and their positive state

Take nonempty binary words w_1,...,w_m. A block is a one-apex tournament: the spine is transitive in word order, and bit 1 means apex beats spine vertex. Order the blocks by ordinal sum; reverse the edge between each consecutive pair of apexes. Root at the final apex. This always constructs an actual tournament on sum_i(1+|w_i|) vertices.

For w of length n, let P be its ordinary binary value, Q=sum_{j=0}^{n-1}(1-w_{j+1})2^j, t=2^n. Let

    g(0)=1; g(2v+1)=g(v); g(2v)=g(v)+v+(v mod 2), v>=1.

The inherited law uses A=H, W=o+e, d=o-e, where o,e are the weighted sums of rooted paths with respectively odd/even vertex count. Set initial formal state (A_0,W_0,d_0)=(1,0,0), and

    B_i = [[g(P_i), Q_i+1], [P_i+1, 2^n_i]],
    epsilon_i = [Q_i even], eta_i=[P_i even].

Then x_i=(A_i,W_i)^t satisfies

    x_i = B_i x_(i-1) - epsilon_i d_(i-1) e_1,
    d_i = eta_i A_(i-1).

Every B_i is strictly positive. Every actual x_i is nonnegative; after the first block both coordinates are positive. The initial formal vector is only an initialization device, not a fictitious nonempty tournament. Nonempty words are essential: singleton blocks have d_i=A_(i-1)-d_(i-1), so the next theorem does not apply unchanged to them.

## 2. Eligible boundaries and exact matching expansion

Put D_i=epsilon_i eta_(i-1) e_1 e_1^t for i>=2. Elimination gives

    x_1=B_1 x_0,
    x_i=B_i x_(i-1)-D_i x_(i-2), i>=2.                 (1)

Call the boundary between blocks i-1 and i **eligible** if D_i is nonzero. It is eligible exactly when w_(i-1) ends with 0 and w_i begins with 1. Thus the memory error is controlled by a binary interface language, not by total tournament order.

Define positive coefficient vectors c_(i,r) by

    c_(0,0)=x_0,
    c_(i,r)=B_i c_(i-1,r)+D_i c_(i-2,r-1),           (2)

with unavailable indices zero and no D_1 term. Then

    x_m = sum_r (-1)^r c_(m,r).                      (3)

Proof: substitute (2) in (1), including its initialization. Equivalently, tile the block path with monomers and disjoint eligible dimers; replace monomer i by B_i and dimer (i-1,i) by D_i, multiply from right to left in chronological order, and sum. The sign is (-1)^(number of dimers). Matrix multiplication is not commutative, so these are weighted chronological matchings, not an ordinary scalar matching polynomial.

Let nu be the maximum matching size of the graph consisting of eligible edges on the m-block path. **The exact correction degree is nu**: all first-coordinate coefficients c_(m,r) are positive for 0<=r<=nu and vanish beyond nu. Indeed each matching contributes positively: the initial vector has positive first coordinate, monomers are strictly positive, and dimers project onto a positive first coordinate. A matching of each size <=nu exists. In particular nu is the sum of ceil(ell/2) over runs of ell consecutive eligible edges. No eligible edge means the d-blind product is already exact at every length. Arbitrarily large blocks can therefore require zero correction; a chain of tiny blocks can require arbitrarily high correction order.

## 3. Rigorous sign of every truncation error

The subtle assertion is not merely (3): arbitrary alternating sums of positive coefficients need not give valid alternating bounds. Physical positivity supplies the missing certificate.

Let b_j=B_j...B_1 x_0 be the d-blind history. On vector histories y=(y_0,...,y_m), define

    (K y)_j = sum_{i=2}^j B_j...B_(i+1) D_i y_(i-2),
    (K y)_0=(K y)_1=0.

An empty suffix product is the identity. K is entrywise nonnegative and strictly history-lowering by at least two indices. Variation of constants in (1) gives **x=b-Kx**. Iterating,

    x = sum_{r=0}^s (-K)^r b + (-1)^(s+1) K^(s+1)x. (4)

Moreover (K^r b)_m=c_(m,r). Consequently for S_s=sum_{r=0}^s(-1)^r c_(m,r),

    S_(2r+1) <= x_m <= S_(2r)                       (5)

coordinatewise, and the **exact signed error** is

    S_s-x_m = (-1)^s (K^(s+1)x)_m.                 (6)

This is an all-length identity, not a fitted inequality. The coefficient bound follows from 0<=x<=b:

    0 <= (-1)^s(S_s-x_m) <= c_(m,s+1).             (7)

The lower endpoint is attained whenever no further matching exists; if s<nu the first-coordinate error is strictly positive. The upper estimate is the adjacent opposite-parity partial-sum bracket; additional previously computed brackets may be intersected without discarding correctness. One more zero coefficient after degree nu closes both sides exactly. Thus a memory-forgetting overestimate has a known downward first correction, upward second correction, etc., with physically certified remainder at every stage.

For any partial upper/lower integer bracket [L,U], an exact residue modulo M>U-L identifies H uniquely. This is Session277's interval/residue recovery principle with a new physical interval backend, not a new modular reconstruction theorem. Positivity and odd parity provide L>=1 and H odd. No computational speedup over the original 3-state recurrence is claimed.

## 4. A persistent physical example

Repeat the word 10. Each block is the directed 3-cycle, and every boundary is eligible:

    B=[[3,3],[3,4]], D=e_1e_1^t,
    M=[[3,3,-1],[3,4,0],[1,0,0]].

The counts begin H_0=1,H_1=3,H_2=17,H_3=111,H_4=721,H_5=4671. Cayley–Hamilton gives

    H_(m+3)=7H_(m+2)-4H_(m+1)+4H_m,
    sum H_m z^m=(1-4z)/(1-7z+4z^2-4z^3).

The d-blind count follows the 2-state matrix B instead; it is strictly too large from m=2 onward. At 10 blocks (30 vertices), the actual count is 53,277,425 while the blind product is 59,932,665. At 100 blocks (300 vertices), the blind count is about 2.708 times the actual count: forgetting one tiny binary memory at each boundary becomes a macroscopic error. The exact correction degree is floor(m/2), not a uniform constant.

`results.json` records exact integers at m=2,3,4,5,10,20,50,100, signed errors, and progressively corrected certified intervals. For m=20 the sufficient residue width drops from 53 bits at order 0 to 44 at order 3 and 35 at order 5. These are guaranteed sufficient bit counts, not optimal bit-minimal claims. The full recurrence already evaluates H efficiently; the gain is an explicit explanation of the error and a usable signed enclosure, not an algorithmic race.

## 5. Constructive use and the missing coverage step

For target inversion, retain each proposed block word and compute the signed bounds before expensive physical witness construction. A target outside the certified interval is impossible for that word; an integer/residue singleton certificate identifies its exact count, and the word itself supplies the tournament. This is a sound selector **within a prescribed physical word language**. It does not show that every target has such a word. The remaining global task is to make the set of exact anchored values cover all positive odds except 7,21; neither matching-degree termination nor the existence of arbitrarily many correction layers proves that.

A sharper obstruction is available even in this one-word language. Write beta=(7+sqrt(37))/2, the Perron root of B. In the positive coordinates (A,o,e), the exact matrix is

    N=[[3,2,4],[2,2,2],[1,2,2]],

which is strictly positive and has characteristic polynomial p(t)=t^3-7t^2+4t-4. Its Perron root lambda lies in (6,beta): p(6)=-16, p(beta)=beta-4>0, and p is strictly increasing beyond 6. There is exactly one real root (its discriminant is -3376), so this root is the Perron root. Perron–Frobenius with the nonzero nonnegative initial vectors gives blind/exact ~C (beta/lambda)^m with C>0. Thus the memory-blind relative error actually diverges exponentially in this explicit physical family. The difference is not merely a huge absolute error on a relatively accurate approximation. This supplies an all-order hostile control against omitting the signed correction coordinate.

There is a stronger fixed-depth obstruction. Mark each correction dimer by u. The exact bivariate series is

    sum_(m>=0) H_m(u) z^m = (1-4z)/(1-7z+(3+u)z^2-4u z^3),
    H_m(u)=sum_r (-u)^r c_(m,r).

Indeed the vector generating matrix is I-zB+u z^2D; invert its 2x2 matrix and take the first entry. Expanding in u yields the **individual error-layer series**

    C_r(z)=sum_m c_(m,r)z^m
          = z^(2r)(1-4z)^(r+1)/(1-7z+3z^2)^(r+1).  (8)

Let gamma=(7-sqrt(37))/2 and rho=1/beta. Since 1-4rho>0 and 1-gamma*rho>0, the unique smallest-modulus pole in (8) has order r+1. Elementary partial fractions give, for each fixed r,

    c_(m,r) ~ kappa_r m^r beta^m,
    kappa_r = rho^(2r)(1-4rho)^(r+1)
              / ((1-gamma*rho)^(r+1) r!) > 0.

Thus for every fixed correction depth s, S_s/H_m is asymptotic to a nonzero constant of sign (-1)^s times m^s(beta/lambda)^m. **Every fixed truncation level fails even in relative error** on this family: even upper bounds diverge relative to H, and odd lower bounds eventually become negative. Extra layers cannot be judged by the small size of one local correction. All the signed layers jointly move the dominant pole from 1/beta to 1/lambda, while any finitely many layers retain the wrong pole. This is the concrete bridge to Session277's pole correction: an infinite-system growth-rate error can be a resummation effect, not a small tail constant. It also explains precisely why the finite-length terminating expansion is not a uniform bounded-depth algorithm.

The transferable principle is **interface-memory correction**: quotient a positive physical state only when the lost coordinate can be retained as a positive causal remainder. This connects the old signed-overlap hierarchy to tournament composition with an explicit map: source is a rooted block chain; target is a chronological matching expansion; preserved quantity is the exact Hamilton count; lost coordinate is d; the restoring sidecar is boundary eligibility plus earlier A-values. Generalizing to longer-lived interfaces asks for a positive finite-memory resolvent, not an unsupported claim that every signed formal word is realizable.

## Reproduction and scope

Run `python3 verify.py > results.json` in this directory, or call the script by absolute path from any directory. Standard library only. Tests: all 1,554 words of 1–4 blocks drawn from {0,1,00,01,10,11}, every coefficient/sign bracket; independent subset-DP on all 258 words of 1–3 such blocks; exact long examples as above. The computations test the new claims, not a new finite coverage range. Proofs (1)–(7) hold for every finite chain of nonempty words.
