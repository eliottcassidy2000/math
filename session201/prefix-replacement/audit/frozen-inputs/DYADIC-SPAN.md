# Actual unit-index points span each dyadic residue class

**PROVED, pending independent audit; separate addendum to DIGIT-ANSATZ.md.**
The formal coefficient obstruction also rules out the same even-scaling
identity on an entire product of two independent dyadic residue classes.
Coupled indices and special subsets remain outside the statement.

## 1. A three-dimensional affine span in every tail

Use the recurrence for g from DIGIT-ANSATZ.md. Let B be any fixed O/Z index
word of length k. Its actual index map has the form

    u'=2^k h-r, 0<=r<2^k,
    g(u')=g(h)+a h+b(h mod2)+c, b in {0,-1},                  (1)

for fixed integers a,c. The first formula follows by iterating h->2h or
2h-1. These choices give every r, hence every residue class modulo2^k.
For the second, O leaves g unchanged and Z adds the current index minus
its parity. Only an initial Z can retain the old parity h mod2; after the
first operation the current index parity is fixed by its last letter.
Induction proves (1), including the empty word with r=a=b=c=0.

The linear map from (g(h)+b(h mod2),h,1) to (g(u'),u',1) has determinant
2^k and is invertible. For any arbitrarily large m, take

    h1=2^m, h2=2^(m+1), h3=2^(m+2)-1.

The first two g-values are1. The third is2^(m+1)+1, and its parity is odd.
Thus the first two points lie at first coordinate1, while the third has
first coordinate2^(m+1)+1+b, unequal to1 for m>=1 and b in {0,-1}.
Their second coordinates are distinct, so the three vectors are linearly
independent. Increasing m puts all three above any prescribed threshold.

Consequently, for every dyadic residue class C and every lower bound H,

    span_Q{(g(u),u,1):u in C, u>=H}=Q^3.                    (2)

This is a statement about actual g-values, not independent formal variables.

## 2. Consequence for independent-index recurrence identities

After fixing all original and transformed parities, each central expression
and each separate substitution in DIGIT-ANSATZ.md is bilinear in

    X(u)=(g(u),u,1), Y(v)=(g(v),v,1).

If an identity holds for all sufficiently large independent u,v in a product
of dyadic residue classes, its difference is X(u)^T C Y(v)=0 for a fixed
3-by-3 rational matrix C. Choose three independent vectors from each class
using (2). Their invertibility forces C=0, so the identity would have to be
coefficient-wise. If the original classes do not fix the required parities,
refine them to dyadic subclasses; an identity on the original product also
holds on that refinement.

Therefore the no-even-scaling proposition in DIGIT-ANSATZ.md remains valid
for actual identities holding throughout any such product of independent
dyadic classes, not only formal coefficient identities. It still does not
rule out an identity with u and v coupled, with finitely many fixed indices,
or on some thinner target-dependent subset. No target-occupancy or spectrum
nonexistence assertion follows.
