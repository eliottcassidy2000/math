# Balanced C1 witnesses for every 3*2^b-1 with b>=5

Session198. **PROVED:** every integer 3*2^b-1, b>=5, is the Hamilton-path count of an actual balanced C1 source graph on b+4 vertices. This supplies an infinite scalar rebase of a transitive-source residual family. It does not prove universal or cofinite odd coverage.

## 1. Physical operations and the balanced readout

Use the reviewed complete O/Z/J physical source grammar at
https://github.com/eliottcassidy2000/math/blob/6f3f93a3915ef2a7eba9a611f67a1bb75bd392c9/session198/bilinear-descent/PROOF.md .
A source has an actual root with one root-start Hamilton path and response (D,E,U,d). Begin at the two-vertex source tau=(1,0,2,1). All words below are chronological: the leftmost letter is applied first.

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

O adds a universal sink, retaining the root. Z and J add a vertex beating only the old root; Z retains the old root and J marks the newborn. Their applicability is to nontrivial sources, as always here.

The chronological word JZ sends any source q to

    (D+2U,D,D+2U,D).                         (1)

It is the actual M2 gadget extension from the complete source theorem. The resulting graph is strong C1, and its retained root has U=H; call this a balanced source. Thus

    B(q)=D+2U

is always physically realized by a balanced C1 graph obtained by adding two vertices to q.

## 2. Even exponents

Let i>=1 and x=2^i. Apply to tau the word

    O^i J O O Z^i J.                       (2)

After O^i J the profile is (2x+1,1,2x+1,1). Writing h=2x+1, the next OO gives (h,0,4h,1). For y=2^i=x, the Z^i block gives

    D_z=(4x-3)h-x+2=8x^2-3x-1,
    E_z=1,
    U_z=4xh-x+1=8x^2+3x+1,
    d_z=1.

The final J gives the actual source

    q_even(i)=(16x^2-1, 8x^2-3x-2, 16x^2, 8x^2-3x-1).

Consequently

    B(q_even(i))=48x^2-1=3*2^(2i+4)-1.    (3)

Appending JZ realizes this count as a balanced C1 graph. The source word (2) has 2i+4 letters, so the final graph has 2+(2i+4)+2=2i+8 vertices. With b=2i+4 this is b+4.

## 3. Odd exponents

Again let i>=1 and x=2^i. Use instead

    O^i J O Z O Z^i J.                     (4)

After O^i J, the profile is (h,1,h,1), h=2x+1. The next OZO gives (3h,0,8h-2,1). The Z^i block gives

    D_z=(8x-5)h-3x+4=16x^2-5x-1,
    E_z=1,
    U_z=x(8h-2)-x+1=16x^2+5x+1,
    d_z=1.

Its final J is therefore

    q_odd(i)=(32x^2-1, 16x^2-5x-2, 32x^2, 16x^2-5x-1).

Hence

    B(q_odd(i))=96x^2-1=3*2^(2i+5)-1.      (5)

Appending JZ again yields a balanced C1 graph, now of order 2i+9=b+4 for b=2i+5.

## 4. The remaining seed and arithmetic mechanism

For b=5, the chronological source word O Z J Z O applied to tau gives

    q_5=(19,0,38,5),  B(q_5)=95=3*2^5-1.

Appending JZ gives a nine-vertex balanced C1 witness. Equations(3) and(5) cover every even b>=6 and odd b>=7 respectively, completing all b>=5.

The crucial cancellation is visible before the two variable runs are tied together. In (2), replace the second Z^i by Z^j and write y=2^j. The balanced readout is

    48xy-18x+18y-1.

In (4), the corresponding readout is

    96xy-30x+30y-1.

Choosing j=i cancels the intermediate term. This is a physical target-dependent rebase with two linked run lengths, not a universal word-transport identity or an assumption that a formal inverse tuple is realizable. Both q families are primitive, since their D and U are consecutive integers.

The transitive residual family identified independently in the balanced-source normal-form work is

    (2^a-1)(2^b-1)+2, a,b>=2.

The theorem eliminates its entire a=2 row beyond the three small values 11,23,47, and by symmetry also the b=2 column. It makes no claim about whether those small values are balanced counts or about the other rows. No theorem from the still-separate residual normal-form work is needed for the present construction: the equality of this displayed row with 3*2^b-1 is arithmetic.

## 5. Verification and scope

The standalone checker verifies both unrestricted two-run polynomial identities, the closed source profiles, the balanced operation(1), primitivity, and the b+4 vertex count. It also constructs the actual graphs and counts paths and cycles directly for b=5,...,10, independently checking source U and root-start counts. Normal and optimized Python receipts must agree. The displayed source computations prove the infinite families; finite graph checks are controls, not the quantified proof.
