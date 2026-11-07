# Session220: exact binary envelopes and a constructive-selector obstruction

**COMPLETE: PROVED scoped theorems, an explicit infinite ribbon family and an exact finite spectral census, with independent audits. Chan–Pak completeness, the three-core inverse conjecture, and eventual F2 coverage remain OPEN.** No unrestricted coverage extension or literature-priority claim.

This proof-focused session connects two inherited arithmetic construction carriers. It obtains exact optimization, a limitation on how a coverage proof can encode witnesses, and a new explicit ribbon family. The positive-density and global induction attempts did not close.

## 1. One binary function controls both carriers

Define f(0)=0, f(2u)=f(u), and f(2u+1)=f(u)+u+(u mod 2). For the inherited one-apex function g and e(v)=1−(v mod 2),

    g(v)+f(v)=v+e(v).

Consequently every fixed-bank join expression is exactly

    a g(v)−b e(v)+c(v+1)
      =(a+c)v−a f(v)−(a−b)(v mod 2)+(a−b+c).

This is the same quadratic carrier as a fixed suffix of Session216's three-core construction. The identity preserves the exact count and the binary witness; it does not identify the two physical families or transfer their coverage statements.

## 2. Exact maxima and minima without enumerating words

For K,h>0 and S>=0, the maximum of Kp−hf(p)−S(p mod 2) on 0<=p<2^m occurs at p=0 or p=2^m−2^j. An internal zero followed by a one can be rounded upward while weakly decreasing f and strictly decreasing binary reversal. The candidate list reduces further to a constant number near an explicit square-root balance.

This gives exact, parity-periodic square-root deficits, including all ties. For the full two-word three-core formula at fixed padding lengths, both maximizing words have at most one run of ones. Thus the maximum needs only polynomially many evaluations, instead of exponentially many binary words.

Every quadratic pair coefficient in the same formula is nonpositive. Therefore its minimum on **any fixed-bit cylinder** is an exact minimum-cut problem with polynomially many vertices and edges. The verifier includes capacities larger than 10^100; augmentation uses a proved capacity bound, not a numeric sentinel.

[Proof, sharp formulas, scope and reproduction](extrema/PROOF.md).

These are exact interval hulls, not occupancy certificates. A concrete legal rectangle has support {517,561,635,669,763,797,871,895}: its hull contains519 but the rectangle does not. This does not exclude519 from the full grammar.

## 3. Finite-state witness selectors are necessarily sparse

Fix any finite bank of a,c>0 and b>=0. A finite-state transducer selecting valid ordinary binary pairs

    (v,N),  N=a g(v)−b e(v)+c(v+1),

can produce only polylogarithmically many distinct outputs N<=X. The transducer may be asynchronous, nondeterministic and partial, and need not report which bank row is used.

The proof is structural: every mixed repeated binary period contributes a strictly positive linear drift to g in addition to its exponential term. A pumped finite-state output cannot reproduce that drift. Thus every useful input loop is monochromatic, forcing a uniformly bounded number of binary runs.

This excludes a finite-state construction retaining the actual F2 witness as a route to all-target coverage. It does **not** exclude a proof recognizing only output support, a different non-rational witness encoding, or an algorithm with unbounded integer registers. Fixed radix changes do not evade the theorem. The known streaming compiler with registers evaluates witnesses but still does not select one for every target.

[Selector theorem and exact boundary](selectors/PROOF.md) · [Independent audit](audit/AUDIT.md).

## 4. An explicit ribbon family, and a bounded induction screen

Let B be either span-two ribbon word 02032333 or 02133232, where each digit is the bitmask of backward arcs to the preceding one or two vertices. For every k>=1, the literal tournament on 8k vertices has

    H(B^k)=(63*265^k+200*2^k)/263.

The first values are65,16825,4457825,1181322025. Exact rank-two matrix algebra proves the formula, and an independent literal Hamilton-path DP verifies both eight- and sixteen-vertex seeds. The physical seed has both a265-mode and a2-mode; the presence of a doubling eigenvalue does not isolate a doubling construction. No strong-connectivity claim is needed here.

The complete finite screen examines349,524 span-two blocks of length<=9 and37,448 span-three blocks of length<=5. Only sixteen span-two blocks have eigenvalue2, the cyclic rotations of these two words; all lack eigenvalue1. Span three has none in its box. Thus a nonconstant repeatedly iterated fixed-block rule2N+c, c!=0, is impossible **within these boxes**, even on a constrained prefix-state family. Larger blocks, state-dependent block choices and nonlocal edits are outside the conclusion.

[Explicit construction, census and exact scope](ribbon/RESULT.md) · [Independent algebra/graph audit and census limits](audit/ribbon/AUDIT.md).

## 5. Inheritance and remaining proof obligation

Source snapshot: acbad6c8636cf8cbb675b4c947fa19e33ac0df23. Session214's induced-growth and surgery barriers, Session215's fixed-cap composition, Session216's three-core formula and density theorem, Session217's exact inverse and run obstruction, and Session219's Fermat/composition work remain in their stated scopes. The curvature obstruction to uniform target halving is credited to Session216, not promoted as new.

During publication, the [root-change chamber packet](https://github.com/eliottcassidy2000/math/tree/58a5a8a3d1d031cc15af397cad7971fb1e658112/05-knowledge/results/root_change_carmichael_chamber_20261107) arrived and was preserved. That pre-audit checkpoint reports a two-root-change census through10^9 and a theorem about singleton word fibers. Its O/Z/J word multiplicities are different objects from the literal binary witness relations studied here. Neither abundance of word descriptions nor finite coverage supplies the missing target-selection theorem. We read its synthesis, did not independently rerun or audit its claims, and use none as a proof dependency.

The concept board is exact occupancy, binary curvature, repeated-word drift, nonuniform composition, and ribbon state spaces. The successful connection preserves the common nonlinear coordinate; removing it loses the very drift that blocks finite-state selection. Sharp extrema discard interior occupancy, which must be restored by a genuinely different invariant.

Neither a physical all-odd constructor nor the needed target-dependent inverse has been found. Theorems here improve exact optimization, rule out a concrete witness architecture and provide an explicit thin family; they are not a solution of Chan–Pak.

The [density investigation](exploration/DENSITY.md) retains the prior bound X^(2/3)/(60(ln X)^(2/3)) for X>=4096. Reversal-product, correlation and multiplicative-closure attempts supplied no stronger unconditional bound. The [cross-edge investigation](exploration/REPORT.md) recovers an exact marked transfer, checked on841 small rooted pairs, but it is already contained in inherited root amalgamation and is not presented as new. These stopping records prevent failed or duplicate routes being mistaken for progress.

## Verified checkpoints

- [Checkpoint1](https://github.com/eliottcassidy2000/math/commit/b792886266a56269422998882a9e34d437bc5acc): extrema, minimum-cut theorem and selector obstruction; all23 remote files read back exactly. The concurrent sevensphere/signature-simplex commit was preserved.
- Checkpoint2: adds the audited ribbon construction and finite screen, records the other routes' stopping points, and replays the complete portable packet. This file's immutable commit URL is the final research index.

## Reproduce

    python3 verify_all.py

Requires Python3, NumPy and SymPy. The portable runner executes seven normal-mode processes, checks exact receipts, and verifies SHA256SUMS, including a full replay of both finite ribbon census boxes. Assertions in normal Python runs are load-bearing; optimized runs are only smoke tests. The independent ribbon audit rechecks all80 first-prime candidates and all16 exact survivors, but does not independently enumerate every noncandidate; that limit is explicit. Other audits record the repaired zero-length cylinder boundary and check arbitrary fixed-bit domains and very large capacities. No coverage scan or coverage-record extension is claimed.
