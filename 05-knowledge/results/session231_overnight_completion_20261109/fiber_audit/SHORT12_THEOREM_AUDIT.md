# Independent audit: twelve-position unit loop

**PASS.** The shorter certificate proves k+1 consecutive odd Hamilton-path counts on **12k+13 vertices**, for every k>=0. The family is saturated with gcd45. Its baselines remain prescribed and exponentially separated; there is no all-target coverage conclusion. The older h37 theorem and all its artifacts are retained unchanged.

## Independent checks

`short12_audit.py` is a standalone stdlib-only checker, importing no producer modules and reading the owned SHORT12_WITNESS.json. Normal and optimized runs match. It verifies:

- Literal permutation/split enumeration of the physical core [4,1,18,7,75,31,175,63], cut124, giving (h,q,A,G,T,S,E)=(45,1440,878,577,332,25,14).
- Every forward edge of X=101000111000, Y=111011100110, with integral next G and strict finite-state bounds. The augmented root is (878,577,0), and the walk returns there after12 edges. No inverse edge or negative walk multiplicity is used.
- The twelve edge weights are −108,406,−492,0,278,646,−72,−682,−296,76,242,0, summing to−2.
- Complete singleton-word states, computed directly from binary words rather than from producer composition:

      X: (1,4096,2616,453,288,1,0)
      Y: (1,4096,3814,1655,1540,1,0).

- Actual flanked states:

      s0 C Y s1: (45,23592960,7535881,16564034,5290654,0,45)
      s0 X C s1: (45,23592960,7535881,16564034,5290656,0,45).

  Thus the two full interfaces agree, while left T minus right T is exactly−2.
- Exact placement checks for k=0,...,8,20,100; these support but do not replace the all-k proof below.
- Direct forward Hamilton-path DP on the actual13-vertex k=0 star, yielding9909, and strong connectivity from every vertex. No generic enumeration of the larger stars is claimed.

## Boundary cancellation and all-k proof

For a closed moving-core word, C Y and X C have equal h,q,A,G. If singleton-word endpoint bits are S_X,E_X,S_Y,E_Y, their unflanked T difference is

    δT = A G_Y − A_X G + h(T_Y−T_X) + E S_Y − E_X S.

Choose left boundary u=E_X and right boundary v=S_Y. The flanked difference becomes

    Δ = δT + u(S−hS_X) + v(hE_Y−E)
      = A G_Y − A_X G + h(T_Y−T_X) − huS_X + hvE_Y.

Here u=0,v=1,E_Y=0. The independent global calculation is

    Δ = 878·1655 −2616·577 +45·(1540−288) = −2.

The boundary is essential: the bare C Y versus X C difference is12, not−2. The selected trailing singleton contributes−14 and gives the certified unit identity.

Now define actual ordinal deletion graphs

    L_(k,j) = s1 ⇒ s0 ⇒ X^j ⇒ C ⇒ Y^(k−j) ⇒ s1,
    0 <= j <= k.

They have12k+11 vertices. Each adjacent placement replaces C Y by X C, preserving h,q,A,G. The outer selected singleton vertices force S=E=h=45, so the full interfaces coincide. Every shared left context ends in0, either the explicit s0 or last(X). Every shared right context starts in1, either first(Y) or the explicit last s1. Both contexts have Hamilton count1. The same boundary calculation therefore applies for every j,k without a multiplier:

    T(L_(k,j+1)) = T(L_(k,j)) + 2.

The common full interface satisfies h=45, q=1440·2^(12k+3), S=E=45, hence gcd(h,q)=45.

Star closure adds v,r with L beating v, v beating r, r beating the selected cut, and the unselected L vertices beating r. Its Hamilton-path count is D=2A+h−T. Consequently the k+1 stars on12k+13 vertices have counts

    N_k, N_k−2, ..., N_k−2k.

They are strong: r reaches the first selected singleton, which reaches every later L vertex; every L vertex reaches v; and v reaches r. Positivity also follows constructively by taking a Hamilton path through C, with all ordinal singleton vertices in order, then v,r. Oddness follows from the exact baseline recurrence below. This proves the all-k theorem, not merely the tested range.

## Exact baseline

For W=Y^k, write A_W,G_W,T_W for its singleton-word scalars, taking these three scalars as0 at k=0 without asserting a false nonempty endpoint identity. Direct composition gives

    N_k =11274·4096^k +135A_W−3758G_W−45T_W−1365.

Equivalently,

    N_0=9909,
    N_(k+1)=4096N_k−264150k−184225,

or

    N_k=(245051264·4096^k +1602510k+1118023)/24843.

These formulas follow from A_Y=3814,G_Y=1655,T_Y=1540 and the concatenation recurrences

    A_(Y^(k+1))=4096A_(Y^k)+3814,
    G_(Y^(k+1))=1655·4096^k+G_(Y^k),
    T_(Y^(k+1))=T_(Y^k)+1655A_(Y^k)+1540,

where the endpoint correction is zero because last(Y)=0. The baseline recurrence immediately preserves oddness. All three baseline formulas are independently checked.

Examples:

- k=0: one count9909 on13 vertices.
- k=1: counts40,403,037 and40,403,039 on25 vertices.
- k=2: counts165,490,399,365 through165,490,399,369 by steps of2 on37 vertices.
- k=3: four odd counts677,848,675,102,893 through677,848,675,102,899 on49 vertices.

## Context and target scope

After the full interfaces have been matched, a common ordinal context with Hamilton counts h_left,h_right scales the T gaps by h_left·h_right. Common singleton contexts preserve gap2; arbitrary tournament contexts generally do not. Strongness under extra context also needs its own check. Keeping an outermost selected leading singleton is a sufficient condition, whereas inserting an arbitrary unselected prefix need not preserve the present reachability argument.

The theorem supplies arbitrarily long consecutive-odd runs at the explicit baselines N_k. It does not prescribe their starting points, fill the gaps between different k, cover every odd target, establish an optimal vertex cost, or establish a twelve-step loop as shortest possible. The small-order census is not contradicted: these deletion graphs have order11 already at k=0 and order23 for the first nontrivial pair.

## Artifacts

Run `python3 -B short12_audit.py` and `python3 -O -B short12_audit.py` in this audit directory. SHORT12_AUDIT.json contains the complete12-edge trace, both word states, both flanked states, all tested placement intervals, and the literal base-star adjacency/count. Matching normal/optimized logs are retained. No producer files were changed and no commit was made by this audit worker.
