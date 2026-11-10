# Round3: a pure-count difference channel with bounded anchors

Status: PROVED universal physical compiler. Exact target placement remains OPEN.

## 1. Inherited scalar-preserving rebases recovered first

Session199 `physical-interpolation/PROOF.md`, Section5, already gives equal-G source pairs with difference(2B,0,−B,0), and a universal six-letter equality. Their OJ closures have the same full count and response, with endpoint/deletion coordinates differing by2B. Thus target-preserving rebasing itself is inherited. The unit local interchange families and every-even-gap existence are also inherited there. Their parent scalar is not freely prescribed.

If two physical states differ by(0,q,0,q), applying J then O gives difference(−q,0,0,0). A subsequent Z preserves that difference and makes reachable-root sources strong. This elementary conversion shows why equal-scalar rebases can lead to a different useful channel: a pure total-count difference with the entire rooted sidecar unchanged.

## 2. Universal seven-letter compiler

Execute the actual O/Z/J moves from left to right. For every source state(D,E,U,d), direct multiplication gives

    OJJJJOZ: (17D+40U, 3D+6U, 21D+50U, 3D+6U),
    OJOZJOZ: (17D+42U, 3D+6U, 21D+50U, 3D+6U).

The difference is exactly(2U,0,0,0). This is a full four-coordinate identity; it does not suppress changes in endpoint, deletion or response coordinates.

If the source root reaches all vertices, the same is true after each O or J. A final Z then makes the result strong. The opening O ensures enough vertices for this argument even from a singleton. Both words have length7 and hence the output orders agree. Their final roots have outdegree2: the last J root beats its preceding root, the O after it adds its second outneighbor, and the final Z beats it. In particular this compiler deliberately leaves the saturated outdegree-one face.

Choose the inherited response selector w(N), with D=g(N)≤N,U=N, for N≥2, and the singleton for N=1. Then every positive N has a strong pair

    L_N=17g(N)+40N, R_N=17g(N)+42N,
    R_N−L_N=2N, R_N≤59N,

with common E=d=3g(N)+6N and common U=21g(N)+50N. The order is at most ceil(log2 N)+8 for N≥2 and8 for N=1. For example N=2 gives the prime pair97,101 with common sidecar(E,U,d)=(15,121,15). This finite prime example is not an infinitude statement.

The quantitative comparison is specific: Session243's one-vertex pair has the better5N scalar bound but does not preserve this full sidecar. THM-5000 preserves a more elaborate interface on different moving-core families with large baselines. The present compiler supplies universal prescribed gaps on a small rooted sidecar and a linear absolute bound. No global minimality is asserted.

## 3. Arbitrary fixed-port continuation preserves the chosen gap

If two rooted states differ by(q,0,0,0), then both O and Z leave this difference unchanged, because

    O(q,0,0,0)=(q,0,0,0), Z(q,0,0,0)=(q,0,0,0).

Therefore **every common O/Z word**, of any length, preserves both the gap2N and equality of all three sidecar coordinates. A nonempty word endingZ gives strong outputs. The anchors can grow without bound, for example under repeated Z because U−d>0 for the constructed sources. This is a stable even-gap channel, unlike the prior synchronized face's difference(q,0,q,0), whose synchronized common continuations could only dilate the gap.

More explicitly, for the lower source(H,E,R,d), the inherited binary-tail formula is

    H(m)=H+R m−d Q(m)−E(m mod2),

with the physical bit order and inherited pair-energy Q. The upper source gives exactly H(m)+2N at the same m. The common tail still supplies this particular nonlinear set of offsets, not every even offset. Fixed modulus surjectivity, an infinite gap fibre, and a linear initial anchor bound do not prove that an arbitrary target appears.

## 4. What a successful rebase still has to do

The new channel allows one to vary the anchor while holding the gap and rooted sidecar relation fixed. It does not yet prescribe that anchor. For example the initial lower target equation is

    M=17g(N)+40N,

and after a chosen tail it becomes the displayed exact pair-energy equation. Neither equation has been solved for every requested M. Allowing J changes the difference vector to(q,q,q,q); other asymmetric moves can retain or alter q, but a physical inverse selector is still needed. No completeness follows merely by matching the number of free parameters to the number of equations.

## 5. Checks

`verify_round3.py` checks both full matrices on a basis,32 literal graph pairs including forbidden full counts7,21 used lawfully as response parameters, all four output coordinates and strongness, and arbitrary common O/Z continuations on finite test ranges. The unbounded continuation claim follows from the two displayed vector identities, not from the tests.

`PURE_DIFFERENCE_EXTENDED.jsonl` is an exploratory matrix search through length11, terminal-Z words. The discovered seven-letter pair is verified separately. The search does not establish a global minimality theorem for graph operations or pair constructions.
