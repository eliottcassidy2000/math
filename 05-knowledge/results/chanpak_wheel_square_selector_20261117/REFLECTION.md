# Reflected flanks: an exact gap with any prescribed congruence

A much stronger local control is possible than independently realizing
two residues. **For any positive even gap D, modulus M, and odd integer
N, a pair of strong tournaments of the same order has counts H and H+D
with H≡N (mod M).** Both arise from the same transitive spine and one
apex; they have one common protected center vertex. The receiver swaps
the two reflected address roles. This works with an arbitrary prescribed
leading binary prefix of the address variable used below.

This couples the two archive mechanisms in a precise way. Square-sum
reflection transports one port without changing strand length. Here the
two addresses may swap because the intersection term is symmetric;
the singleton port reads their difference precisely. The role of the
prime wheel is then entirely in the middle-word choice, with no effect
on the *exact* difference.

## Literal switch

Use the marked singleton center that beats the apex. Its profile in the
signed flank law is (1,2,0,1,0,-1,1). Thus
`H(u,w)=2uw+u+[u odd]+F(u)+F(w)-1`.
For w=u+D with D even, swapping u and w adds exactly D. These are
physical strong tournaments: there is a Hamiltonian cycle following
the transitive spine with the apex before its first vertex and after
its last, because each positive address starts with a 1. Swapping the
addresses preserves the total number of vertices. No claim is made
that this swap is a single arc edit.

## Gap-preserving congruence steering

Choose t>=2 with D<=2^(t-2), put r=2^(t-2)-D, and write
`u=2^t x+r, w=2^t x+r+D`, where x is any positive binary
address. Both r and r+D are even. Their t-bit suffixes start in 0;
the next bit is 0 for r, 1 for r+D. For a t-bit s, expansion of F
gives, with rev the bit reversal of a t-bit word,

    F(2^t x+s) = F(x) + (2^t-1-rev(s))x
                  + [leading bit of s is 0][x odd] + F(s)-1.

(The coefficient of x is the sum of its contributions from appended
zero bits; only the first new bit sees x's old parity.) So with
epsilon=[x odd], the lower count H has the form

    H = 2 F(x) + J(x) + 2 epsilon,

where J is an explicit integer quadratic polynomial. Its nonconstant
coefficients are multiples of 4: the product term and u coefficient
clearly have this property; the two reverse-complement coefficients
sum to 0 modulo 4 since rev(r) and rev(r+D) are 0 and 2 modulo 4.
The constant J(0) is odd. The formula at x=0 need not be viewed as
a physical word; J(0) here means its polynomial constant.

On every odd modulus the coefficient of F(x) is the unit 2.
For a dyadic modulus at least 4, first choose x's parity so
J(0)+2+2 epsilon matches the desired residue mod 4. Then solve
for the corresponding odd residue of F(x) modulo the dyadic part
(one bit less suffices). Modulo 2 the lower count is already odd.
The protected-fringe address theorem lets x and F(x) attain these
compatible residues, while carrying any chosen binary prefix of x.
It also supplies arbitrarily long joint witnesses in two equal-length
returning address blocks. **The difference D stays exact throughout**,
not approximate and not only modulo M.

## Exact-target audit: a real scale obstruction

For fixed D this does not fill every absolute target. Indeed induction
on x gives 1<=F(x)<=x+1 (for odd x, the stronger F(x)<=(x+1)/2
handles the odd predecessor in the even recurrence). The lower
count at any positive u lies in the disjoint integer strips

    B(u)+[1,2u+D+2],   B(u)=2u^2+(2D+1)u.

Since B(u+1)-B(u)=4u+2D+3, a prescribed target picks at most one
candidate u, by integer square root or binary search and then the
strip test. The exact F values decide equality directly; both
directions are constructive. There are at most O(sqrt(X)) such
lower counts below X for any fixed D, despite perfect congruence
coverage for exact pairs. So even a single tiny core can make the
stronger local test pass everywhere that parity permits. A universal
selector needs a richer choice of absolute scales, not merely better
access to local residues.
