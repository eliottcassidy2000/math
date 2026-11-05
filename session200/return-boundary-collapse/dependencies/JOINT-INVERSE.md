# A synchronized square-root inverse for the joint source image

Session 200. **PROVED, independently audited in the stated scope.** The pair
`(D+2U, 2U-d)` has an exact bounded inverse on actual physical sources.
The auxiliary boundary need not be a source, but is transported by the
same word as a physical boundary and is coordinatewise bounded by it.
A physical Mersenne boundary then gives a necessary-and-sufficient pair
criterion for every prescribed C1 target. **OPEN:** universal target
coverage. No assertion about Carmichael's conjecture or rings of integers
is used as a premise.

## 1. Inheritance, objects, and the missing synchronization

Let S be the complete physical orbit of `t=(1,0,2,1)` under chronological
operations

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

A profile records total Hamilton paths, root-ending Hamilton paths,
root-start reserved two-path covers, and root-deletion Hamilton paths.
Every source has exactly one root-start Hamilton path. The frozen physical
adjoint dependency proves that these words are actual rooted tournaments
and that every nontrivial source response has such a word. Only J changes
the root. Use

    F(s)=D-E+U, G(s)=D+2U, H(s)=2U-d,
    b=(3,1,3,1), z=(1,0,2,0),
    Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]],
    K(v,w)=v^T Q w.

The physical adjoint reverses a word and exchanges O and Z while fixing
J. Write this involution as dagger. Then

    K(v,b)=G(v), K(v,z)=H(v),
    K(G_alpha v,w)=K(v,G_(alpha-dagger) w).

The known C1 image is F(S); appending a closing J realizes the count.
The balanced image is `{3} union G(S)`. The known joint descent reads

    F(S)={3} union G(S)
          union {m+2(2^k-1)h : (m,h) in P, k>=1},
    P={(G(s),H(s)):s in S}, 3<=h<m<2h.             (1)

The parent `(m,h)` must be realized by one source. The canonical hostile
443 has parent choices `(9,7),k=5` and `(65,63),k=2`, and none with k=1.
The least-used sidecar is z: it has deletion zero and is not a source.
Concept board: physical pair inverse; joint image; auxiliary boundary;
Mersenne boundary; root-change budget; false fiber occupancy.

The present construction turns the final synchronization obligation in
`CARMICHAEL.md`, Section 7, into an exact, smaller finite inverse. It does
not replace that obligation by separate membership of m and h.

## 2. A bounded decorated orbit

Define the synchronized orbit

    R={(w,q)=(G_beta b,G_beta z): beta a source word}.             (2)

The word beta is shared. The first member is physical; the second is an
auxiliary integer vector. Retain the pair, not just w.

**Lemma 1 (bounded sidecar).** For every `(w,q)` in R, all raw coordinates
satisfy `0<=q_i<=w_i`. In particular they are at most `U_w`.

**Proof.** The extended positive-source cone is given by

    (a,f,e,p)=(D-E-d,U-D+E-d,E,d-E)>=0.

In these coordinates the three generators act as

    O(a,f,e,p)=(a+e,a+2f+e+2p,0,e+p),
    Z(a,f,e,p)=(2a+f+e+p,f,e+p,0),
    J(a,f,e,p)=(f+p,0,a+e+p,e).

All entries are nonnegative. Since z has coordinates `(1,1,0,0)`, every
q has nonnegative raw coordinates. For the difference `Delta=w-q`, use
another cone

    D>=d>=E>=0, U>=d.                            (3)

It is preserved by O, Z, and J: substitute the formulas, using `U>=E`.
The initial difference `b-z=(2,1,1,1)` lies in (3), so all its raw
coordinates remain nonnegative. This proves q<=w. Finally physical w
satisfies `E<=d<=D<=U`, proving the stated height bound. QED.

This is a bound on an auxiliary vector, not a physical realization theorem
for q. For example `z.Z.J=(7,3,7,3)` must not be added to S merely because
it satisfies relaxed source-cone inequalities.

## 3. Exact square-root criterion for simultaneous realizability

Put `B(m)=floor(2 sqrt(m-1))`. Let S[C] be the complete physical orbit
through `U<=C`, and R[C] the synchronized orbit through `U_w<=C`.
Both inventories are finite: U strictly increases at every physical letter
and all source coordinates are bounded by U; Lemma 1 bounds q as well.
Every stored item retains its predecessor recipe.

**Theorem 2 (joint inverse).** For integers m,h, `(m,h)` lies in P if and
only if m>=5 and there exist

    v in S[B(m)], (w,q) in R[B(m)]
    with K(v,w)=m and K(v,q)=h.                  (4)

The same equivalence holds with a root-change budget r: retain the minimum
J cost of each v and of each synchronized pair, and require that the two
costs sum to at most r. It is the pair's minimum cost, not merely w's.

**Proof of sufficiency.** If `v=t.alpha` and
`(w,q)=(b,z).beta`, take the actual source

    s=t.alpha beta-dagger.

The adjoint gives `G(s)=K(v,w)` and `H(s)=K(v,q)`. Its J cost is the sum
of those of alpha and beta. This proves realizability with a literal word.

**Proof of necessity.** Start with the physical word for s in
`m=K(b,s), h=K(z,s)`. Split it as `alpha gamma`. Moving gamma to both
boundaries gives

    v=t.alpha, (w,q)=(b,z).gamma-dagger,
    m=K(v,w), h=K(v,q).                          (5)

The physical adjoint's product bound gives

    U_v U_w<=2(m-1).

As the cut moves right, v gains a letter, its cover count grows strictly
by a factor at most two, and w loses its last physical letter. Choose the
first cut with `U_v>=sqrt(m-1)`. If the initial cut qualifies, `U_v=2`.
Otherwise the crossing height is less than `2 sqrt(m-1)`. The product
bound controls `U_w<=2 sqrt(m-1)`. If no crossing occurs, use the terminal
cut: `U_v<sqrt(m-1)` and w=b has height3, which is at most
`2 sqrt(m-1)` since m>=5. Thus both heights are at most B(m). Moving
both boundaries along the same gamma keeps h exact. J costs add to the
original cost, proving the budget version. QED.

No independent integer choices for D and U were introduced. The compiler
is an exact physical inverse, not a solution of a relaxed cone or an
appeal to Diophantine definability. Its finite search is polynomial in a
crude numerical height, not asserted polynomial in the target's bit length.

## 4. Physical Mersenne boundaries and a target criterion

For k>=1 set `L_k=2(2^k-1)` and

    c_k=b+L_k z=(2^(k+1)+1,1,2^(k+2)-1,1).

**Lemma 3.** The boundary c_k is the actual source `t.O^k Z`, and

    K(c_k,s)=G(s)+L_k H(s)=F(s.O.Z^k).           (6)

This follows directly from the source operations or the adjoint of
`O Z^k`. The linear expression `b+L_k z` is not being declared physical
by cone membership; its displayed word supplies physicality. Likewise
`w+L_k q=G_beta c_k` is physical for every synchronized pair in R.

Every F(s) is odd by the source parity invariants, so even targets are
rejected. Here is an exact necessary-and-sufficient criterion for a
prescribed odd integer n>=3 to lie in F(S):

1. n=3; or
2. there are `v in S[B(n)]` and `(w,q) in R[B(n)]` with `K(v,w)=n`; or
3. there are k>=1, v in S[C_k(n)], and `(w,q)` in R[C_k(n)] such that

       m=K(v,w), h=K(v,q), n=m+L_k h,
       2^k<=(n+1)/6,
       C_k(n)=floor(2 sqrt(floor(n/2^k)-1)).      (7)

A k with `floor(n/2^k)<2` is omitted. In case 2, use the balanced closing
on `t.alpha beta-dagger`. In case 3, use the source word
`alpha beta-dagger O Z^k` and then close by J. Thus accepted arithmetic
pairs produce actual C1 graphs, not merely compatible residues.

**Proof.** Cases 1 and 2 are the balanced image. The projection of R onto
w is exactly the physical orbit from b, so the known balanced inverse
applies. Case 3 follows from (6) and the constructive direction of Theorem 2.
Conversely (1) supplies a joint parent for every remaining representation.
Since `h>m/2`, its target satisfies `n>2^k m`, whence
`m<=floor(n/2^k)` (n is odd). Thus B(m)<=C_k(n), and Theorem 2 supplies
the required pair. Also `m>=5,h>=3` gives `n>=6*2^k-1`, proving the finite
k bound. QED.

The larger k is, the smaller the unweighted inventories become. The full
weighted boundary `w+L_k q` may be much larger; no small-height claim is
made for it. For k=1 the physical inventories have height at most
`floor(sqrt(2n-6))`. This improves the old linear-height parent search while
retaining both readouts. It does not prove that every admissible n passes.

The criterion also preserves a prescribed precursor J budget r. In (1),
removing a terminal `O Z^k` changes no J count; a last J instead produces
a balanced recipe with no larger J count. All-Z words have budget zero.
Apply Theorem 2's budget form in cases 2 and 3. This is an exact restricted
inverse, not a claim that any fixed budget covers all large n.

## 5. Two indispensable hostile controls

**Independent boundary choices give false positives.** Put

    v=t.O=(1,0,4,1), w=b.J=(5,2,6,3), q=z.O=(1,0,4,0).

Then `K(v,w)=15`, `K(v,q)=13`, and `3<=13<15<26`. Both chosen words fit
the square-root cutoff B(15)=7. But `(15,13)` is not in P. Indeed
`D+2U=15`, `2U-d=13` and positive D,d force U=7 and D=d=1. The physical inequality
D-E>=d then forces E=0, contrary to the invariant parity E=U modulo2. Equivalently the complete physical fiber at m=15 has only
h=9 and h=11. The mistake is using J for w and O for q.

**Even equal physical boundaries do not determine q.** From (b,z),

    OJJZ -> ((27,9,27,9),(15,5,15,5)),
    ZZOJ -> ((27,9,27,9),(23,7,23,7)).

Both are valid synchronized states. With v=t they give `(m,h)=(45,25)`
and `(45,39)`. Merging the states only by w discards genuine continuation
data. Keep the full eight-coordinate pair, or an equivalent exact sidecar.

The Carmichael lesson is that separate nonempty fibers do not supply a
common member. The Diophantine lesson is to encode the shared physical
witness, not merely compatible integer equations. Here the common word,
the simultaneous adjoint, and the bounded sidecar supply that witness.

## 6. Verification and limits

`verify.py` checks the two cone invariances on their integer generating
rays, all bilinear identities on a coordinate basis, common-word transfers,
and both hostile controls. It compares the joint inverse against a complete
forward physical source inventory, including minimum J costs, and compares
(7) with the complete C1 target image in a declared finite universe. It
also counts Hamilton paths and cycles in literal reconstructed small graphs.
Normal and optimized Python receipts must agree. Finite agreement verifies
the implementation; the all-order claims rest on the proofs above.

Frozen dependencies are `dependencies/ADJOINT.md` and
`dependencies/CARMICHAEL.md`. Source-pair existence for every odd target
outside the global exceptional set remains open. Neither local congruence
solvability, plentiful signed-response witnesses, nor this finite inverse
by itself settles it.
