# Session281: renewable banks, complete experiments, and the cost of residue selection

**Chan–Pak remains OPEN.** This session distinguishes actual constructors from arithmetic shadows. This recovered publication has passed final scoped independent review; see audit/REVIEW.md and RESTORATION.md for receipt provenance after the runtime filesystem rollback.

## 1. The positive construction

Six actual eleven-vertex rooted tournaments have Hamiltonian counts

    979, 981, 983, 985, 987, 989

and share deletion count h=15 and receiver coordinate V=1203. They need not share e. Left fusion with the three-vertex one-apex block C1 (word10) satisfies

    H'=H+2V, h'=h, V'=4V-h, e'=h.

Thus the six actual thirteen-vertex tournaments with counts3385,3387,...,3395 share the full reusable receiver(h,V,e)=(15,4797,15). This repairs the earlier six-bank defect by finding a different raw bank and preparing it, not by pretending the earlier defective union was renewable.

After k more C1 preparations,

    V_k=4792*4^k+5,
    A_k=3385+(9584/3)*(4^k-1)+10k.

For even k the paid carry construction supplies twelve consecutive odd counts. At k=0 they start at3,051,967 and have order23. These are infinitely many finite intervals at explicit locations, NOT an eventual tail. The two halves of the twelve-bank have different future receiver coordinates. Renewable six does not mean renewable twelve.

More generally an actual r-bank sharing(h,V), with odd h, can be prepared for paid doubling whenever r lies in the orbit of V under multiplication by2 modulo h. This is a sufficient preparation criterion for this construction, not a classification of all possible preparations. See receiver/PROOF.md.

## 2. Complete observation is relative to an operation language

The rooted-fusion carrier has seven nontrivial coordinates(H,h,Q,U,V,d,e). Seven actual Hamiltonian-count experiments, each adding at most two vertices, recover every coordinate. Consequently equality of those observations is equivalent to equality under ALL two-sided same-root fusion continuations.

For one-sided renewal, four experiments recover(H,h,V,e). If one allows the common preparation C1 first, only three experiments are required:

    r0=H(T), r1=H(C1*T), r2=H(C1*C1*T),
    V=(r1-r0)/2,
    h=(5r1-4r0-r2)/2.

Equal count offsets in all three observations characterize a bank with common(h,V), precisely the information the preparation consumes. This turns a vague demand for 'compatible seeds' into a finite physical test.

This compression is NOT valid for arbitrary rerooting. Two six-vertex rooted tournaments have identical full fusion carriers but different deletion-count menus. Their strong cores can be concatenated in arbitrary binary orders, producing exponentially many nonisomorphic objects with the same fusion response. The forgotten graph information can matter as soon as the allowed operations change. See wildcard/PROOF.md.

## 3. An exact local compiler can still miss almost every integer

Seven physical active blocks, words10,100,110,1000,1010,1100,1110, reach every admissible scalar residue modulo every modulus. At odd moduli their full packet orbit is exactly the primitive triples. The proof uses elementary matrices away from3, an explicit mod3 certificate and Lie-algebra lifting at3, and a rank-one dyadic attracting line with two seed cosets. Sandwich permutations synchronize odd and dyadic conditions.

These are positive physical programs: finite-ring inverses are replaced by positive powers. They are not inverse tournament operations. The cost is their actual integer height.

The same language is power-sparse. With per-letter growth floors

    G=(5,11,9,23,19,23,17)

and u=(276,147,173,82,95,82,104), exact inequalities u_i^5 G_i^4>=1000^5 and sum u_i=959 give at most(1000/41)X^(4/5) programs of count<=X. Hence infinitely many odd integers are omitted despite every finite congruence test passing. Large dyadic windows contain runs of omitted odds of order at least X^(1/5). Its infinite scalar-support generating function is not rational, though word-counting generating functions can be rational. No claim of nonalgebraicity follows.

The structural lesson is operational: residue, interface and absolute height are three distinct resources. Full control of the first two does not supply the third. See selection/PROOF.md.

## 4. A new central switch, with an exact limitation

A finite classification at deletion order6 finds coprime central switches of gap4, including h17 and h29. Common fusion contexts multiply such gaps by odd deletion counts; they cannot divide4 into2. Combining the h29 gap4 switch with the renewable six-bank gives a common-receiver twelve-set with offsets

    0,58,60,116,118,174,176,232,234,290,292,350.

It contains five adjacent odd pairs but is not a twelve-interval. Coefficient gcd2 is insufficient for binary subset-sum coverage. The finite absence of coprime gap2 switches at this order is not a global obstruction.

## 5. Principles and next constructive obligations

- **Context-complete observation:** compress only after naming all allowed continuations.
- **Preparation before renewal:** an imperfect bank may become reusable if a common physical context erases exactly the unstable coordinate.
- **Local controllability with height debt:** finite-ring transitivity is a selector only when its absolute cost is controlled.
- **Interface-aware addition:** differences, gcds and interval widths must retain their common receiver and absolute baseline.
- **Operation-relative indistinguishability:** a complete quotient for fusion need not survive rerooting or arbitrary switching.

The strongest next targets are a renewable wider bank, a constructor that repairs the twelve-bank defect without losing its overlap, or a height-sensitive selector that anchors an interval at a prescribed target. None is supplied by this session. The links to recognizable series, finite automata, local/global arithmetic and information loss are mechanisms, not arguments by analogy.
