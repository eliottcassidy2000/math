# The ternary sum palette has an aligned all-size decoder

October 23, 2026. **PROVED, elementary:** the identities below for every
finite tournament. **FINITE-EXACT:** all labeled tournaments of orders 2
through 6 and all 6,880 order-eight class representatives. This is a local
modular-incidence decoder with alignment, not reconstruction from two
unlabeled outer multisets and not an all-size metagraph classification.

Let A be the 0/1 adjacency matrix, with A_uv=1 for u->v. Let c_v count cyclic
triangles through v and c_uv count cyclic triangles through the pair {u,v}.
All formulas below with residues are in F_3. Give each actual arc u->v the
ternary sum color

    chi_T(u->v) = c_v + c_uv mod 3.

Converse preserves c_v and c_uv, but the head of the reversed arc changes:

    chi_Top(v->u) = c_u + c_uv mod 3.

## 1. Vertex decoding from aligned forward and converse rows

For each same labeled vertex v, retain the forward and converse outgoing
color histograms. Their weighted sums are

    r_v  = sum_{v->w} chi_T(v->w),
    r*_v = sum_{w->v} chi_Top(v->w).

A histogram (N0,N1,N2) has weighted sum N1+2*N2 in F_3. The triangle-incidence
balance gives

    sum_{v->w} c_vw = c_v = sum_{w->v} c_vw.

Consequently

    r = (I+A)c,
    r* = (I+A^t)c.                                    (1)

Since A+A^t=J-I and sum_v c_v=3*c3(T),

    r+r* = (I+J)c = c  in F_3.                        (2)

Thus the TWO ALIGNED vertex histograms recover every c_v mod3, with no
unknown additive constant. The alignment is essential: the outer multiset
of forward histograms and the separate outer multiset of converse
histograms do not say which rows belong to the same vertex.

## 2. Recovering the edge residue needs actual edge association

If the color of each actual arc is also retained with its head identity,
then (2) gives

    c_uv = chi_T(u->v) - c_v  in F_3.                 (3)

This recovers the complete modular pair (c_v mod3,c_uv mod3) attached to
that arc. Vertex histograms alone do NOT supply which colored outgoing
edge reaches which head; (3) requires that additional association. No
unlabeled inversion or improved asymptotic census is claimed.

The difference of the two colors attached to the same reversed edge is

    chi_T(u->v)-chi_Top(v->u) = c_v-c_u mod3.           (4)

It is the exact coboundary of the vertex triangle-count potential. Its sum
around any directed closed walk is zero. This is a genuine gradient on the
TOURNAMENT vertex carrier; it is not an electrical-current reconstruction
on the metagraph's class vertices. Arbitrary colored arrays satisfying
(2)--(4) need not be realizable triangle incidences of a tournament.

## 3. Why the plus sign has an advantage over the minus sign

If instead one uses psi_T(u->v)=c_v-c_uv, the corresponding weighted sums
are

    d=(A-I)c,  d*=(A^t-I)c,
    d+d*=(J-3I)c=0 mod3.                              (5)

Their sum loses c rather than decoding it. On a cyclic triangle all c_v
and c_uv are one, so both minus palettes are identically zero, whereas
the plus weighted rows are all two and (2) recovers one.

This does not by itself prove any unlabeled-signature completeness or
failure. The separate finite quotient search proves that the sum palette
has a reflection-marginal forest at order eight, while the difference
palette does not. The algebra explains a compatible mechanism: forward
and converse plus observations retain a vertex potential that the same
combination of minus observations annihilates.

## 4. General linear palettes and a necessary caution

For chi_(alpha,beta)(u->v)=alpha*c_v+beta*c_uv in F_3, the aligned
weighted-row sum is

    r+r*=(2*beta-alpha)c.                              (6)

Before reducing modulo three, the same expression has the additional
term 3*alpha*c3(T) times the all-ones vector. Thus this paired-row decoder
works exactly when 2*beta-alpha is nonzero in F_3. The plus palette has
coefficient one; the minus palette has coefficient zero.

This local condition is NOT a predictor of unlabeled completeness. The
edge-only palette (alpha,beta)=(0,1) also decodes the aligned c vector,
yet its residue-three anonymous signatures identify only 6,693 order-eight
classes, and their paired converse labels identify only 6,846. The local
decoder and the global reflection forest are independent requirements: one
describes retained modular information under alignment, the other certifies
uniqueness after anonymous signature bucketing.

## Reproduction

    python -B incidence/three_color_check.py --classes wiring/classes_input_n8.tsv
    python -O -B incidence/three_color_check.py --classes wiring/classes_input_n8.tsv

The consumer enumerates cyclic triples directly, independently of the
bit-intersection implementation in the quotient/converse checker. It
checks both operators, decoding, edge residues, gradient differences,
and minus collapse. Normal and optimized runs reproduce
`three_color_receipt.json` exactly. No removable assertions or external
packages are used.
