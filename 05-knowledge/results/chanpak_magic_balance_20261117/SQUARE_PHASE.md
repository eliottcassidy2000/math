# Magic-square contact grids with arbitrary prescribed square-residue precision

A 3-by-3 magic square with center c has the Lo Shu affine shape

    c+a     c−a−b   c+b
    c−a+b   c       c+a−b
    c−b     c+a+b   c−a.

The eight line sums are 3c, and the nine signed offsets are exactly `xa+yb` for x,y in {-1,0,1}. The arithmetic squares problem then requires all nine c+xa+yb to be distinct integer squares at the **same** c,a,b. See the recent primary [algebraic-geometry preprint](https://arxiv.org/abs/2609.09351) for the six-quadric formulation and [Cain's Gaussian-integer treatment](https://arxiv.org/abs/1908.03236) for a different exact formulation. Neither is a proof that an integer square exists. The current repository's Session57 and Session143 contain the corresponding decorated-frame and cubic-arithmetic work; see the linked provenance in the session index.

**Theorem.** Fix j>=3 and any positive odd M coprime to 389·521. For every return height k>=3 one can build nine strong tournaments whose H counts lie at the above nine distinct magic-grid positions, all congruent to 1 modulo 2^j M. The two balanced steps a,b are independent of k. Yet at all sufficiently large k (explicit threshold below) at most one of the nine counts can be an actual square. Thus simultaneous, arbitrarily high *specified finite* square-residue agreement does not produce an absolute square grid, even though the affine grid itself is exactly physical. No independent claim about the magic square of integer squares is made.

## Physical transfer (inherited, restated)

Take an ordinal stack R of marked blocks and a vertex q beating exactly the tags of R; add TT_k losing to R and beating q. For R use scalars J,P,W,Z,E,S,U: ordinary H, optional two-slot count, weighted starts at a tag, weighted ends at a tag, full ends at a tag, full starts at a tag, and two nonempty slots with first ending tagged and second starting tagged. Put K=J−U. Direct classification around q gives H=W2^k+K: if q is fed by the return block, the nonempty initial return run and optional final return run have 2^k−1 possible splits times W; the other terms are J+W−U. Indeed the term J counts q terminal; if q is followed by a tag and not fed by the return block, the preceding R path may be empty or must end in a non-tag, giving W−U.

For an appended ordinal block with profile (h,p,w,z,e,s,k_0), using K in place of U in this tuple, the transfer is

    J'=hJ, W'=pW+Jw, E'=Je,
    K'=hK−J(h−k_0)−Es−Wz.

To verify the last equality, partition the two marked endpoints by left or right block. The both-left and both-right cases give hU and J(h−k_0), respectively. Left-end/right-start gives E s, and right-end/left-start gives W z; the latter has the optional restrictions internal to each side. The W law follows by placing the marked start in the left or right block. Hence changing only a local k_0 changes final K by its increment times the product of all other h's, without any mixed terms. This is the exact additive geometry that the source-bank assumptions must justify; formally placing independent offsets is not enough.

## Two switches and one residue repair

Each row is an outgoing mask and the marked-set masks are listed in switch order:

| tile | outgoing rows | tags | (h,p,w,z,e,s,k_0) |
|---|---|---|---|
| A | (158,164,56,242,98,1,167,52) | 144,6 | (389,5600,1168,1376,102,88,217→225) |
| B | (10,52,81,70,105,141,163,31) | 40,20 | (521,7376,1662,1662,124,124,331→267) |

These are the two stackable fixtures in the earlier ordinal-digit packet. Both blocks are strong and every tag set is nonempty. The whole tournament is strong: q reaches the first strong block through a chosen tag, that block reaches every subsequent block, and all R reaches the return block and thus q.

Let L=2^(j−3)M. Use the order A^(2L), B^(2R), A^(L−1), taking R an odd multiple of M at least L (specifically L if L is odd, L+M if L is even). The first 2L copies of each type divide into two groups of size L, and any extra B copies stay unswitched; the x parameter selects 0,1 or 2 full A groups switched, and y similarly for B. The final L−1 A copies repair the common phase and do not belong to either axis. Let J be the whole stack H. The balanced steps are

    a = 8L J/389,        b = −64L J/521.

Their magnitudes have ratio 8·389/521>2, so none of the nine offsets collide. Both steps vanish modulo 2^j M.

At k>=3 the baseline is 1 modulo 8 before any switches. Start an empty stack as (J,W,E,K)=(1,0,0,1). Appending A preserves K mod8 because h≡5, h−k_0≡4, s=z≡0 and J,K odd. The A prefix is even in length, so J≡1, W≡0, E even before B begins. Each B has h≡1, p≡0, w=z≡6, s≡4 and h−k_0≡6. Thus W stays 6 after the first B and K shifts +2 initially, then +6 for every subsequent B. Over the 2R B copies the total shift 2+6(2R−1)=12R−4 vanishes modulo8 since R is odd. The trailing A helper bank preserves K modulo8, as does any switch, and W2^k vanishes modulo8.

Each helper switched adds 8J/389, of valuation exactly 3 and odd part coprime to M. Thus 0,...,L−1 helper switches run through all L residues that are 1 mod8 modulo 2^j M, and select the unique desired phase for the entire grid. One may use R=L for odd L or R=L+M for even L. Setting j=3,M=1 gives R=L=1, so no extra copies are needed in the concrete example below. W is a positive constant independent of all switches.

For fixed j,M the nine outputs lie in [c−d,c+d], where d=|a|+|b| does not grow with k, and c is W2^k plus one of finitely many constants (the helper settings). Once c−d>d² no two integer squares in that interval can be distinct, since the next square gap at its left reaches strictly more than 2d. This is the absolute obstruction despite all finite local tests encoded in the chosen modulus.

## Smallest concrete grid

For j=3,M=1, stack A,A,B,B and read with TT_3. Then c=2967778437612897, a=844724392, b=−5045647424. The common W is 382559934716446 and all cells are 1 modulo 8. Larger-k grids stay exact, independently of whether any entry happens to be square.
