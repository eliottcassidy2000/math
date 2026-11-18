# A physical Bezout family with exact unit-gain difference control

**PROVED construction, not universal absolute placement.** The finite-edit obstruction can be bypassed by growing the edit footprint. This is not hypothetical source supply: two order-five sources suffice.

Use outgoing masks `(20,1,2,7,14)` and `(16,1,3,7,14)`, both rooted at vertex 3. Literal rooted profiles are respectively

    L=(9,1,23,5), R=(9,1,25,5).

For every r=4j+2 (j>=0), put

    q=(2^(r+1)+2)/5,
    u=reverse_r_bits(q),
    WL=2^r+2^(r+2)V, WR=u+2^(r+2)V, V>=0.

q is an even integer less than 2^r, so 0<=u<2^(r-1). The roots' common binary transfer constructs literal tournaments for every V. Direct expansion of the checksum gives

    K(WL)=K(V)+2V,
    K(WR)=K(u)+K(V)+4qV.

Therefore their exact counts are

    HL=9+23*2^r+(23*2^(r+2)-10)V-5K(V),
    HR=9+25u-5K(u)-(u mod2)+(25*2^(r+2)-20q)V-5K(V).

Their difference is exactly

    HL-HR = C_r-2V,
    C_r=23*2^r-25u+5K(u)+(u mod2).

C_r is positive and even. Positivity follows from u<2^(r-1), giving C_r>(23-25/2)2^r. Evenness follows either from the literal tournament counts or the displayed checksum parity. Thus V=C_r/2 produces an exact physical recontact at every level. More generally every even requested signed gap <=C_r is obtained by the explicit nonnegative integer V=(C_r-gap)/2. The family is **congruence-complete for differences**, with no modulus beyond parity.

The mask is z_r=1 and z_j equal to minus the jth bit of u for j<r. Its reverse-weighted sum is 1-2q, and the high-tail coefficient from the finite-edit decoder is

    A=2^(r+2)+5(1-2q)=1.

This is the exact Bezout mechanism: d=5 and the phase r=2 mod4 let the two actual weighted-start reads differing by 2 cancel all but one unit of the remote-bit response. The periodicity here is the order of 2 modulo 5, not a conjectured numerical analogy.

## Smallest family, entirely explicit

At r=2 one has u=1, and

    HL(V)=101+358V-5K(V),
    HR(V)=33+360V-5K(V),
    HL-HR=68-2V.

At V=34 both counts equal 12,233. The shared nonlinear term does not need to be approximated: it cancels exactly from the difference. All actual outputs are odd, and the replay checks them against independent subset DP for small V.

## What this does and does not solve

This provides a physically supplied, periodically growing family of coordinated contacts and a true unit-step controller for half the signed difference. It avoids asking that each large read remain frozen individually. It does **not** make either output equal a prescribed N. Fixing the desired gap fixes V; the remaining absolute baseline must still be selected by r, the low-word freedom, or additional source changes. Consequently a useful next target is to retain this exact coefficient A=1 while supplying a wide independent baseline palette. Merely proving more modular surjectivity would not address that target.

## The growing mask itself has a finite recursion

Writing r=4j+2, the selected low address and its checksum obey

    u_0=1, u_(j+1)=16u_j+6,
    K(u_0)=0, K(u_(j+1))=K(u_j)+6u_j+2.

Equivalently,

    u_j=(7*16^j-2)/5,
    K(u_j)=[14(16^j-1)-10j]/25,
    C_j=(299*16^j+36-10j)/5+[j=0].

Thus the exact collision uses V_j=C_j/2, with no search for its high continuation. The first four common counts are 12,233; 2,827,529; 721,199,909; 184,602,568,309. The checksum recurrence follows by appending the four low bits of 6: their mutual contribution is 2 and their cross contribution to the old address is 6u. Induction gives the closed forms. Distinctness also follows from crude bounds: K(v)<=v, V_j is of order 2^r, and the common count is bounded between 400*(2^r)^2 and 1100*(2^r)^2; increasing j multiplies 2^r by 16, separating these bounds.

Replay: `python3 finite_edit_recontact/unit_gain_replay.py`. It verifies the literal sources and small literal appended graphs, including the first recontact, and exact bigint identities through r=42 plus a large test continuation. These checks accompany the all-order proof, not replace it.
