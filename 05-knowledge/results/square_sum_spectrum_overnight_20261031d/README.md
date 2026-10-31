# Readable branching, retained cuts and exact digital inverses

**Status: PROVED scoped constructions and inverse/identity theorems; FINITE-EXACT tournament coverage through99,999,999; universal raw spectrum OPEN.** Same-agent proof review with separate permutation, interface, adjacency, and arithmetic consumers. No independent subagent review or literature-priority claim.

## Inheritance and board

The frozen [previous packet](../square_sum_spectrum_overnight_20261031c/README.md) gives600 independent code channels, the seven-parent digital forms, exact coverage through9,999,999, and the sharp **finite** prefix cap56. Its original source-index threshold1062 came from using A/B code prefixes as the two choices. We now replace that leaf supply, make branching recursive, and recover the previously unused zero-offset channels. All older theorems remain valid.

The closest proved mechanisms are the incoming weak-B pair interface, the positive marked-block action, and the exact singleton digital formula. The canonical hostile is assigning different source parity to opposite offsets: this duplicates295 and omits293. The corrected near miss is treating a deleted prefix edge as if it destroyed the entire ordered source. The least-used sidecar is the complete order of the separated tail, together with its known head and integer atom boundaries.

| Lane | Object / operation | Retained predicate | Cheapest test and outcome |
|---|---|---|---|
| Anchor | Fixed seven-parent tournament bank; digit inversion | Actual source and chronological word | Independent capped support maps agree through100 million |
| Niche | Weak pairs with independently chosen absolute-offset children | Common-label parity and recoverability from **each** member | Complete98 paired and98 terminal interfaces;625-leaf physical decoders |
| Niche | Local odd reversal or even-block exchange | Every label's position parity, prefix and endpoint | A complete950-pair switch bank lowers the uniform leaf index to112 |
| Wildcard | Split source head/tail | Ordered history even when the source edge is absent | All625 native channels decode;25 formerly fixed channels become free |
| Wildcard | Finite-width binary reflection | Exact source address and its width | Recovers the earlier Fibonacci fan, preventing a duplicate Mersenne claim |

The META-PATTERNS card used is **Canonicalize mechanics; retain the predicate**. A decoder proves a count of distinct objects; it does not prove that every admissible numerical target enters the representation.

## 1. Two individually distinguishable weak pairs at every index n>=112

A weak B pair consists of square-sum Hamiltonian paths on[1,n] and[1,n+1], both beginning1,3, ending at8 for an even-length member and10 for an odd-length member, with identical position parity for each common label.

For each112<=n<=1061 the finite file `leaf_switches.json` gives one of two transformations for **each member** of the inherited base pair:

* Reverse an internal segment of odd length, retaining the first two and last positions. Reversal moves every label an even number of positions. Two new boundary sums are checked to be squares; every internal edge is retained in reverse order.
* Exchange two consecutive internal blocks, each of even length, without reversing them. Every moved label again changes position by an even number. Three new boundary sums are checked; all internal edges survive.

The unchanged pair is choice0; applying the two listed moves gives choice1. Each output member is different from its original and retains every individual label parity. Hence the changed pair remains weak B, and either member alone distinguishes the choice. `audit_leaf_switches.py` independently converts the certificates into permutations and checks all1900 paths,1,115,300 vertices,1867 reversals and33 exchanges.

For n>=1062 use the already proved incoming A and B operators: their common prefixes have third vertices13 and6 respectively. Thus each member again distinguishes the two choices. This gives a uniform two-choice supply at **every n>=112**.

The number112 is sufficient, not proved optimal among all possible path constructions. At n=111 the supplied first path admits neither of the two searched moves; that is only a boundary of this local move bank. At n=42 a single segment reversal is structurally unavailable because the square-sum graph has no4-cycle. The finite nested-graph check places the first4-cycle at order46, namely1–3–46–35–1. This explains one obstruction without claiming uniqueness of the base paths.

## 2. A readable25-ary weak-pair tree

Set

    T_l=(225*49^l-1)/2,     l>=0.                (T_0=112)

**Theorem.** For every n>=T_l, there is a family of2^(25^l) weak B pairs of index n such that projection to **either member** is injective on the family. Every member starts with the incoming fixed ray prefix P_(C^l), of length4l+2.

**Proof.** At l=0 use Section1. For l>=1 put m=floor((n-24)/49), so m>=T_(l-1), and d=n-49m lies in24..72. Independently choose one child pair for each a=0,...,24. Use that pair in both signed offsets+a,-a of the fixed C paired template; use the a=0 pair in offset0.

For a common source label x, the two signed copies occupy49x-a and49x+a exactly once because the two selected members have the same position parity. Different source labels have disjoint bands. The exceptional last source label has the parity forced by its length, so the inherited top partial-band rule also applies. Every connector uses only the weak-B prefix and endpoints. Separate head, tail and singleton parity profiles agree across the two outputs at every offset. Thus the same complete paired catalogs still produce a weak pair. These assertions are checked symbolically in all98 states by `audit_tree_interfaces.py`; they do not assume equal colors in unrelated channels.

To prove injectivity, fix either output member, n and l. Its template and all atom lengths are known. For each nonnegative offset a, read its block or tail at the declared positions, undo its orientation, and round each lifted label v to floor((v+24)/49). This recovers the selected child path; if its head was split off, prepend the known label1. The selected child-member index is part of the template. By induction, that one member recovers its entire child choice. Every one of the25 channels has such a readable tail in each output. Consequently all25 independently chosen child codes are recovered, giving

    (2^(25^(l-1)))^25=2^(25^l)

distinct pairs, with injectivity in each member. The fixed C head followed by the intact positive offset4 block gives P_(C^l) in both outputs. This completes the induction.

The decoder retains n, l, the template, orientations and cut positions. It does not claim that an arbitrary weak pair has a unique history under every grammar.

### Cycles with the same prescribed prefix

For l>=1, replace only the top paired template by its terminal C template. Each absolute-offset child remains readable, and the cycle's unique oriented traversal beginning1,3 fixes rotation and reversal. The complete98 terminal catalogs have every required tail. Therefore, for every N>=T_l, there are at least

    2^(25^l)

distinct undirected square-sum Hamiltonian cycles on the fixed labelled interval[1,N], **all beginning with the same P_(C^l)**. The difference is stored in the branching interiors, not in incompatible initial prefixes.

Examples: at least2^25 cycles at every N>=5512, and2^625 at every N>=270112. For beta=log(25)/log(49), approximately0.8270874753, choosing the largest admissible l gives a stretched-exponential lower bound2^(c*N^beta) with an explicit positive constant. This is a lower bound, not an asymptotic formula or upper bound for all cycles.

The exact controls include the complete32-choice cube on five marked channels at n=5512, independently decoded from both members; opposite input parity at5513;625 leaf choices decoded at270112; and terminal cycles at5512,5513,270112. Normal and optimized receipts agree. The induction rests on the complete interfaces and decoder, not the number of sampled choices.

## 3. Retaining the native forest costs none of the625 channel choices

The older native outer construction used24 variable nonzero inner channels at each of25 absolute outer offsets, giving600 choices. Its inner zero-offset source was fixed because the source edge1–3 is cut there.

That edge loss does **not** erase the ordered source. Its label1 is a known separated head, and the complete remaining tail is readable. Restoring the head in the decoder recovers the entire selected source member. The head's location and the named2–7 A ports depend on atom lengths, not on the source choice. Therefore the zero-offset input can vary independently as well: there are25*25=**625 channels**.

For every l>=0 and N>=2401*T_l+1200=T_(l+2), choose a depth-l pair independently in each of these625 channels. The two old protected outer layers produce at least

    2^(625*25^l)=2^(25^(l+2))

distinct cycles, while retaining:

1. the native path(1,8,17,19,6,10,15,21,4,12,24);
2. the seed482–3362–6242;
3. all2303 pairwise disjoint ray prefixes

       (2401*p_j+(-1)^(j+1)*(49u+v))_j,
       (p_j)=P_(C^l), u in[-24,24] except0,2, v in[-24,24].

The prefix cuts still explain why only2303 **intact ray paths** are asserted. They do not reduce the625 independent **recoverable histories**. These are different predicates.

For injectivity, the outer native edge1–8 fixes cycle orientation. Both pieces of every nonnegative-offset A path occur, with the same member index and known lengths. Reassemble them, then recover all25 inner child paths, including the split zero channel. Section2's member decoder recovers every leaf choice. `audit_native_readability.py` checks all196 inner and392 outer states. `native_tree_cycle.py` checks complete cycles at270112 and270113 and recovers all625 choices **from each final cycle**, without consulting the original intermediate banks during decoding. A further control changes only one zero-channel bit, holds the other624 fixed, and obtains a different canonically oriented cycle while retaining the same prescribed forest.

In particular the full protected forest admits at least2^625 cycles for every N>=270112 and2^15625 for every N>=13235512. No disjoint source ports are reused and no labels are relabelled to claim the native embedding.

## 4. A square-root-size exact inverse for each digital source

Inherit f(0)=f(1)=1, f(2a)=f(a)+a+(a mod2), f(2a+1)=f(a). An actual singleton source has output

    N=F0+q*y+g*(y mod2)+h*(f(y)-1).

Put t=N-F0+h. Inverting a fixed source asks for

    t=q*y+g*(y mod2)+h*f(y).                       (1)

The exact two branches are

    y=2a:   (t,q,g,h) -> (t,2q+h,h,h),
    y=2a+1: (t,q,g,h) -> (t-q-g,2q,0,h).          (2)

The terminal y=0 has t=h. For q>0, every candidate lies in

    max(0,ceil((t-h-g)/(q+h))) <= y <= floor((t-h)/q),

because1<=f(y)<=y+1. Modulo gcd(2q,2h), the two parity branches also require t-h congruent to0 or q+g. The implementation tests a small final interval directly, otherwise applies(2). An explicit stack retains the assembled address, avoiding a programming-language recursion-depth cutoff.

For the actual q=0 source one has(F0,g)=(h,0). Repeated odd arguments leave f unchanged. Strip those redundant outer odd digits and take the first nontrivial argument even, reducing to the positive-slope state(t,h,h,h). This retains all attainable counts; it chooses one canonical representative of possibly many words. All numerical targets and reconstructed words still retain the source and integer clock.

**Complexity theorem, arithmetic-operation model.** For0<=g<=h and q>0, the algorithm uses O((1+sqrt(h*t)/q)*log(t+2)) arithmetic operations when t>=h, and logarithmic stack depth. At binary depth d the slope is at least2^d*q, the remaining target is at most t, and the candidate interval width is at most

    h*t/q_d^2 + h/q_d + 1.

Once q_d>=sqrt(h*t), this is at most3, so direct enumeration finishes. A full binary tree to that depth has O(1+sqrt(h*t)/q) nodes. Each leaf makes at most eight O(log(t+2)) evaluations of f; gcd tests have the same arithmetic-operation upper bound. For q=0 the preliminary even representative gives O((1+sqrt(t/h))*log(t+2)). Integer bit costs must also be charged; this is **not polynomial in the binary input length**.

`source_inverse.py` agrees with direct finite enumeration on80498 typed source queries. The former linear-interval inverse remains correct; this is an algorithmic improvement in how its complete search is performed. Nine isolated large targets, including999999999999, have independently counted actual tournament witnesses. Isolated successes do not prove interval or universal coverage.

## 5. Finite-width reflection recovers an earlier stronger fan theorem

For k>=1,0<=x<2^k, let rev_k(x) reverse its k binary digits, including leading zeroes. Then

    f(2^k-1-x)=f(x)+rev_k(x)-x+(x mod2)-floor(x/2^(k-1)).       (3)

Expand f as the sum of its even-branch charges. Complement exchanges the zero and one positions; summing all higher-digit charges gives x-(x mod2), and the finite top bit contributes the last correction. This proves(3). The width k is essential: bit reversal without its width loses the identity's coordinate.

For k>=4, define

    a_k=2^(k+1)-3*2^(k/2)-2                 if k is even,
    a_k=2^(k+1)-5*2^((k-1)/2)-2             if k is odd.

Then f(a_k)=2^k-1. This is **not a new Mersenne existence theorem**. The literal singleton-prefix graph has one dominating source vertex. Remove that explicitly named padding vertex: the apex loses precisely to chain positions

    {floor(k/2), k+1-floor(k/2), k+1}.

This is exactly the previously proved Fibonacci fan in [THM-4748, two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity](../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md). That theorem already proves its strong realization and the minimal three-feedback-arc boundary for k>=5. The new calculation identifies the digital address and the padding map, retaining the stronger prior result.

The reflection identity is exhaustively checked on131070 finite-width addresses. Actual fan adjacency counts agree for k=4..128, and the integer identity is checked through k=4096. Reflection turns the original positive-slope source formula into one with a negative x coefficient and an explicit rev_k(x) term; it is an exact representation, not permission to treat a negative-slope formal state as an attainable tournament source.

The companion [power-of-two slope classification](POWER_TWO_SLOPES.md) gives an exact obstruction to replacing that changing fan word by one fixed pure-zero tail. Every actual singleton/D state with q a power of two has positive F-q+g, equal to h or h/5. Therefore no fixed such source followed only by singleton0 steps supplies infinitely many values2^k-c for a fixed positive c. This is a controller obstruction, not a raw-spectrum exclusion; the existing Mersenne fan itself illustrates the difference.

## 6. Exact finite tournament coverage through one hundred million

Retain the previous seven-parent bank: D and D ordinal-sum T_k ordinal-sum D for k=0..5, all16144 actual cuts, plus singleton-only words. Its fixed-block forms are

    F+(p+Ey)x+qy+(H-s)(x mod2)+(H-t)(y mod2)
      +H*(f(x)+f(y)-2).

With x<=56 and arbitrary suffix length, the exact image through100000000 contains every positive odd except7,21,43,93,173,247:49,999,994 counts. `capped_scan.cpp` uses a source-state table and recursively streams the singleton axis; `audit_capped_scan.cpp` uses the quadratic form and independently updates f by binary carries. Their entire6,250,001-byte packed maps agree, with SHA256

    34530aa67d56eb84e89943e9b525f37e2baaecb2b1ccce39ad094f8f4ef6fa7d.

The producer checks1,261,579,639 parameter pairs; the consumer checks1,522,019,810 using different valid bounds. The latter's carry formula is independently checked against the direct digit expression through one million. If z=y+1, k=v_2(z), a=floor(z/2^(k+1)), then

    f(z)-f(z-1) = -a-(a mod2)                         if k=0,
                = z-z/2^k+1-a-(a mod2)               if k>0.

Only the changed trailing bits contribute. Both scans include the singular product axes via the bank's actual all-zero cut; they do not assume such an axis for an arbitrary isolated profile.

The four saved actual witnesses for43,93,173,247 therefore give **unrestricted coverage of every positive odd through99,999,999 except7 and21**. These remain raw Hamiltonian vertex-word counts, never H/Aut. THM-1370-h-spectrum-omits-7-21-all-n.md already excludes7 and21 at every order. This computation does not exclude any other integer globally.

The cap56 is still sharp in this finite interval, because the previous complete target test for154711 excludes x<=55. Nothing here proves a uniform source cap for all integers.

## Remaining targets and reproduction

**OPEN anchor:** prove the seven-parent digital cover for every odd target>=249, or find a hostile beyond the verified interval. The exact inverse supplies a finite rejection certificate for a fixed source bank, not universal acceptance. The bit-reflection identity suggests retaining both address and reversed address when studying nearly all-one blocks.

**OPEN niche:** lower the sufficient112 leaf boundary with different moves, or enlarge the recursive choice family while retaining member-wise injectivity. The local failure at111 does not exclude another weak-pair family. Arbitrary prescribed forests and a globally coherent infinite Hamiltonian ordering are not proved by the finite-cycle multiplicity.

**Closed shortcut:** independently varying opposite signed sources remains invalid. Cutting a source edge may preserve its decodable history, but does not preserve the intact source-path embedding predicate. The625/2303 distinction is the exact boundary.

Run the Python files from this directory; their standard output names match the committed JSON receipts. `leaf_switches.py` also writes its certificate bank. Add `-O` for the optimization-control replay. Main commands:

```sh
python3 leaf_switches.py
python3 audit_leaf_switches.py
python3 audit_tree_interfaces.py
python3 audit_native_readability.py
python3 tree_pairs.py
python3 native_tree_cycle.py
python3 source_inverse.py
python3 digital_reflection.py
python3 local_switch_obstruction.py
python3 bank_query.py 10000001 10000003 10000019 100000001 999999999 2147483647 1099511627775 1099511627777 999999999999
c++ -std=c++17 -O2 capped_scan.cpp -o /tmp/capped_scan
c++ -std=c++17 -O2 audit_capped_scan.cpp -o /tmp/audit_capped_scan
```

Pass each C++ consumer `100000000 56 /tmp/<output>.bin`, followed by the seven previous-packet TSVs D,DD,DT1D,...,DT5D. `VERIFICATION.json` gives their exact commands, mode agreement and consumer scopes. `MANIFEST.json` pins the packet and inherited proof roots. The ordinary square-sum existence theorem remains inherited classical work; no new claim of that theorem is made.
