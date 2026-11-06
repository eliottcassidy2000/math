# Positive profile defects and exact halving in the start-two grammar

Status: RESERVED / proof candidate under final independent audit.

## Scope and inheritance

This packet proves a universal rooted-cover inequality and derives strictly positive scalar defects for every non-neutral operation in the complete start-one/start-two grammar. It also sharpens the complete start-two membership search order cap by one vertex. It does **not** prove that an inverse branch exists at an arbitrary scalar target, does not prove staggered profile occupancy F(8,137), and gives no new eventual odd tournament-spectrum threshold.

The inherited source is `session203/start-two-gates/PROOF.md`, frozen here as `dependencies/START_TWO_GATES.md`, published at commit `46af9e96b8681662ae41dd316df7f91cee6ef0e9`. Its proof supplies the complete physical grammar and the exact nine-coordinate transfer maps. Its former cap `max(2, floor(log2 B)+3)` is valid but is improved below. No prior claim is retracted.

For a nonempty tournament U with root v put

    p = (h,s,e,L,R,d,C,b,T),       B = h+L-s.

Here h counts Hamilton paths of U, s those starting at v, e those ending at v, and d counts Hamilton paths of U-v, with the empty count equal to one. C and b count ordered two-path covers of U and U-v, respectively; bins may be empty. L counts covers with one designated nonempty path ending at v and one ordinary bin. R counts covers with a designated nonempty path starting at v and one ordinary bin. T similarly has a designated root-start path and two ordinary labeled bins. B is the odd Hamilton-path count of the inherited edge extension E(U).

The seed is the root-first directed edge, with state

    (1,1,0,1,2,1,4,2,3), B=1.

O appends an out-sink while retaining the root. Z adds z with v losing to z and every other old vertex beating z, retaining the root. J is this same physical extension re-marked at z. K1 and K2 are the two physical gates defined in the dependency; they add two or three vertices and double the start count on nontrivial inputs. Complete start-two sources have an O/Z/J prefix, exactly one gate K1 or K2, and an O/Z/J suffix. Singleton exceptions remain as recorded in the dependency; all grammar inputs here are nontrivial.

## 1. Elementary physical cover inequalities

For every nonempty rooted tournament:

    L >= h, R >= h, d >= e, h >= s,
    f0 := R-h+e-d >= 0,
    C >= 2L, C >= 2R, T >= R, T >= d.

For the first two, cut a Hamilton path immediately after or before its root. For f0, the R covers contain both the h images obtained by cutting a Hamilton path before the root and the d covers with singleton designated root path. Their overlap consists exactly of the e Hamilton paths ending at the root. The inequalities d>=e and h>=s follow by deletion or inclusion. In a labeled two-cover the root lies in one of the two bins, giving C>=2L and C>=2R. Emptying the extra ordinary bin gives T>=R; a singleton designated root path gives T>=d.

If U is nontrivial, two disjoint subfamilies of T have respectively a singleton designated root path and a spanning designated root path. Therefore

    w0 := T-b-s >= 0.

If the root has an outneighbor, also

    q0 := b-L-d >= 0,
    a0 := C-2L-b >= 0.

For q0, partition the ordered two-covers of U-v by the endpoint of the first bin. Precisely,

    b-L = sum_{i : v->i} L_i(U-v),

where L_i is the designated-end count at i in U-v. Every summand is at least d by cutting after i, and the index set is nonempty.

For a0, take any labeled two-cover of U-v. Choose the first labeled bin that contains an outneighbor of v, and insert v immediately before its first such outneighbor. The preceding vertex, if present, is an inneighbor of v, so this is legal. The root is nonterminal in its bin. Deleting the root recovers the original labeled cover. This injects the b input covers into the C-2L covers where the root is nonterminal.

Every nontrivial positive-start root has an outneighbor, so all these inequalities apply to every grammar input.

## 2. A universal forest inequality

**Theorem 1 (PROOF CANDIDATE).** Every nonempty rooted tournament satisfies

    tau := 2T+2L-C-2R-b+2d >= 0.

The singleton is included, and no positive-start assumption is required.

Let A=U-v. Work in the commutative squarefree algebra on the vertices of A, and extract the coefficient of the product of all its vertex variables, denoted [A]. Let P enumerate all nonempty directed paths on A; X enumerate those whose first vertex is beaten by v; and Y enumerate those whose last vertex beats v. Distinct directed orders are counted separately. Put f=1+P, S=1+X and E=1+Y. Root cutting gives

    T=[A]f^2 S, L=[A]f E, C=2[A]f E S,
    R=[A]f S, b=[A]f^2, d=[A]f.

Consequently tau is the full coefficient of

    1+2P+P^2-2X-2XY+2P^2 X-2PXY.

Group its terms by unordered spanning directed-path forests of A. For each path i let x_i be one if v beats its first vertex, and y_i be one if its last vertex beats v. A forest with ell nonempty paths has coefficient

    ell=0: 1;
    ell=1: 2(1-x_1);
    ell=2: 2[1-(y_1 x_2+y_2 x_1)];
    ell=3: 2 sum_j x_j(2-sum_{i != j} y_i);
    ell>3: 0.

The only negative coefficient is -2 for two paths both of type (x,y)=(1,1). Choose one of these two paths by a fixed deterministic ordering. It starts at a root outneighbor and ends at a root inneighbor. Split it at its first consecutive transition from an outneighbor to an inneighbor. The resulting three paths have respective types (1,0), (0,1), (1,1), and their forest has coefficient +2.

This map is injective: concatenate the unique (1,0) path with the unique (0,1) path to recover the split path; the unique (1,1) path is the unchanged path. Thus every negative weight cancels a distinct positive weight of equal size. All remaining coefficients are nonnegative. The empty A case contributes one, so the singleton is also covered. This proves the theorem.

## 3. A positive cone on the complete grammar

Set

    g := C-3L-h.

This packet does not assert g>=0 for every arbitrary positive-start rooted tournament. It proves the following sufficient closed class.

**Theorem 2 (PROOF CANDIDATE).** Starting at the rooted edge, every sequence of O,Z,J,K1,K2 has g>=0. In particular it holds on the complete start-one and start-two grammars.

The inherited transfer maps give the exact identities

    g(OU)=2g+3L+h,
    g(JU)=tau+f0+q0,
    g(ZU)=g(JU)+3a0,
    g(K1 U)=tau+2(C-L)+(2T-R-h)+e,
    g(K2 U)=2tau+5(C-2L)+(12T-3h-2R-2d)+2e.

For J and Z the input root has an outneighbor. The gate right sides are nonnegative for every rooted input, using T>=R>=h, T>=d and C>=2L. The seed has g=0. Every grammar operation preserves positive start and nontriviality. Induction proves the claim. The assertion about arbitrary repeated gates is only a sufficient construction class; the completeness assertion is confined to the inherited one-gate start-two grammar.

## 4. Exact scalar defects and halving

Write p0=d-e>=0, in addition to f0,q0,a0,w0 above. Exact subtraction in the transfer maps gives

    B_O=B,
    B_Z-2B=s+f0+2p0+q0,
    B_J-2B=g+(R-e)+s,
    B_K1-2B=2q0+a0+2s+2w0+3R-h-e,
    B_K2-2B=h+2(R-e)+4(R-d)+2d
               +C+2(C-2L)+2T+2(T-b).

The last identity uses R>=d, which follows from the singleton designated-root subfamily. For grammar inputs, the four non-O defects are strictly positive. In particular the Z defect is at least s, the J defect at least 2s (R-e>=h-e>=s because the disjoint start/end Hamilton-path families have no intersection for nontrivial U), and the K1 defect is positive since 3R-h-e>=R>0. The K2 expression is positive since h>0, T>=b+s, and all its other terms are nonnegative.

Thus every existing non-O grammar step has an inverse whose B value is strictly less than half the output B. Every non-O defect is odd, since both source and output B are odd. If a representation has q non-O steps, then

    B >= 2^(q+1)-1,
    q <= floor(log2(B+1))-1.

Neutral O steps are not bounded by this inequality. They are bounded in target search by the order cap below.

For an exact Z ray define delta=s+f0+q0. Z preserves s,f0,q0, resets p0 to zero, and satisfies B_Z=2B+delta+2p0. Therefore for every k>=1,

    B(Z^k U)=2^k(B+delta+p0)-delta.

This is a formula for existing rays, not an all-target occupancy statement.

## 5. Sharper complete order cap

**Theorem 3 (PROOF CANDIDATE).** Every represented start-two pair (B,2) has a physical representative of order at most

    max(2, floor(log2 B)+2).

Remove terminal O operations, which leave B and s unchanged. If nothing remains, the source is the seed (and has start one). Otherwise consider the last non-O operation, with input order n. Any ordered two-cover count of an m-vertex tournament is at least 2^m: every vertex subset and its complement each admit a Hamilton path. Thus C>=2^n and b>=2^(n-1). Also B_out>=L_out, since h_out>=s_out.

- For Z, L_out=L+b>=2^(n-1); output order is n+1.
- For J, L_out=C-L>=C/2>=2^(n-1); output order is n+1.
- For K1, L_out=C-2L+2T-b>=2T>=2b+2s>=2^n; output order is n+2. This uses a0>=0 and w0>=0 and is the one-vertex improvement over the inherited argument.
- For K2, L_out=3C-2L+4T-2b>=2C+2b>=2^(n+1); output order is n+3.

Every case implies output order <=floor(log2 B)+2. Removing terminal O is a valid operation on the physical grammar word. The retained word still contains its unique gate when the start count is two. This proves the result.

## 6. Exact compiler, checks, and limits

`compiler.py` is the inherited complete start-two compiler with exactly one functional change: its order cap is now `max(2, bound.bit_length()+1)`. It enumerates the complete one-gate grammar, merges only equal full nine-coordinate states, and retains the least order. The preceding theorems justify monotone B pruning and the sharper order cap. Resource exhaustion returns INCOMPLETE and permits no absence inference. `dependencies/COMPILER_BEFORE_CAP.py` preserves the prior compiler.

`verify.py` independently enumerates directed-path forests and their cancellation injection, checks direct subset-DP counts and every transfer identity, and compares support against the frozen prior compiler on fixed finite boxes. The computation is a check of the all-order arguments, not their source. Exact commands and receipts are recorded with the final audit.

## Remaining occupancy obligation

The inverse halving theorem says: given an existing grammar representation, one can descend rapidly through its non-neutral steps. It does not say that every prescribed B admits one of those inverse steps. Even complete start-two membership is only one of the endpoint slices needed by F(8,137). The all-target staggered profile property and the unconditional eventual odd-spectrum conclusion remain OPEN.
