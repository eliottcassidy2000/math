# Independent audit: endpoint common roots and content growth

Date: 20 October 2026. Status: **PROVED within the stated endpoint scope; FINITE-EXACT independent replay.** No publication was performed. No general two-edit exclusion or Collatz conclusion is asserted.

## Verdict and exact scope

Both proposed content identities and the subsequent short-relation growth bounds are correct for

    p >= 3, p/2 < r < p, gcd(p,r)=1, q=p+r.

Use the zero-intercept lower mechanical word of length q and weight p, labelled by source positions

    s_j = j + ceil(jr/p) - 1, a_j = -jr mod p,
    s_(j+p)=s_j+q, a_j in {0,...,p-1}.

The plus endpoint path moves the one at residue p-1-r right, then the one at residue p-1 right. The minus endpoint path moves residue r left, then residue 0 left. With B the chronological carry and D=2^q-3^p, their exact contents are

    g_plus = gcd(2*13^p-10^p, 3*13^r-2*10^r),
    g_minus = gcd(9^p-2*5^p, 2*9^r-3*5^r).

These are full integer gcd identities, not radical equalities or necessary-only gates. In particular, they retain all prime-power depth.

The inherited mechanism is the moving-anchor formula in `two-edit/symbolic/TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md`, Sections 2–5. The hostile is precisely that packet's endpoint norm barrier: its large norm does not identify actual content. The least-used sidecar retained here is the common clock root together with the moved anchor. The live concept board consists of physical legality, cut transport, whole-modulus units, common-root uniqueness, and short lattice relations. This audit makes no literature-priority claim.

## 1. Legality and the cut convention

At k=1 an initial right move is legal exactly for a<r; an initial left move is legal exactly for a>=p-r. Consecutive labels satisfy a_(j+1)=a_j-r mod p.

For plus, the residue p-1 label has successor residue p-1-r. The latter is less than r because r>p/2, hence it can move right into a zero. Its source gap from the follower is one. The leading move unlocks the follower. The opposite ordering is blocked. When r=p-1, the leader is label 0 at q-1, and its move is q-1 to 0. This crossing is not excluded or mistaken for a new label.

For minus, residue r is immediately before residue 0. Residue r can move left since r>p-r; residue 0 cannot initially move left. Moving r first unlocks 0. The opposite ordering is blocked. In the ordinary cut these moves are q-2 to q-3, then q-1 to q-2.

There are distinct moved labels because 0<r<p. The assumption p>=3 includes the minimal instance p=3,r=2. No undeclared unchanged-anchor assumption is used.

For a left rotation by one bit e, the exact carry relation is

    2 B(rotated word) = 3^e B(word) + e D.

Since gcd(D,6)=1, this proves rotation preserves gcd(|D|,B), with either sign of D. The independent verifier checks this identity at every cut of every tested final word.

## 2. Common-root transport over a whole modulus

Work modulo any m dividing D, with m>1. Then 2 and 3 are units. Choose u,v with pu+rv=1 and define

    z=2^u (3/2)^v.

The compatibility 2^r=(3/2)^p follows from D=0, and directly implies

    z^p=2, z^r=3/2.

This root is unique among units satisfying both equations: the ratio of two roots has both p-th and r-th powers equal to 1, so its pu+rv-th power is 1. This is valid in every commutative unit group, not only over a field.

Moreover

    (z-1)(1+z+...+z^(p-1))=z^p-1=1,

so z-1 is a unit without any separability or prime-field assumption.

Let d_j be the actual labelled displacements and anchor the final word at the final position of label 0. Its ordered offsets are

    j+ceil(jr/p)+d_j-d_0, 0<=j<p.

They lie in [0,q) by preserved cyclic label order. Directly substituting the two clock equations into the chronological carry sum gives

    B_anchor = 2^(-d_0) 3^(p-1)
               [S(z)+sum_j (2^(d_j)-1) z^(a_j)] mod m,
    S(z)=1+z+...+z^(p-1).

Thus the carry vanishes if and only if

    Qraw(z)=1+(z-1)sum_j(2^(d_j)-1)z^(a_j)=0.

All omitted factors, including 2^(-d_0), are units. This proves both directions of carry transport at the whole modulus. The special case m=1 is trivial and needs no modular inverse operation.

## 3. Exact identities and denominator primes

### Plus

For the actual plus vector,

    Qraw = 3-z^(p-1)+z^(p-r)-z^(p-r-1) = Qplus.

The clocks imply z^(p-r)=4/3, and consequently

    z Qplus = (13z-10)/3.

At the actual content g, carry vanishing therefore forces 13z=10. This relation implies 13 is a unit modulo g: any common prime divisor of 13 and g would divide 10. It also implies 5 is a unit, since z and 13 are units. Thus z=10/13, and the two clock equations imply g divides each displayed plus gcd operand.

Conversely let N be that gcd. Its two operands immediately show

    gcd(N,2*3*5*13)=1.

For 2 use the second operand, for 3 and 5 use the first, and for 13 either operand suffices. Define z=10/13 modulo N. The operand congruences are exactly z^p=2 and z^r=3/2. Raising to cross powers gives 2^(p+r)=3^p, so N divides D. The same root has z-1=-3/13 a unit; it gives Qplus=0, hence carry transport gives N divides B. Therefore N divides g as well.

### Minus: notation warning

The integral minus polynomial is

    Hminus = 3-z-z^(r+1)+z^r = 2 Qraw,
    2 Hminus = 9-5z.

If this integral polynomial is named Qminus, the correct identity is 2 Qminus=9-5z. It must not be confused with the unscaled universal Qraw.

At the actual content, Hminus=0 forces 5z=9. Hence 5 is a unit modulo g, z=9/5, and g divides both minus gcd operands. Conversely their gcd N satisfies gcd(N,2*3*5)=1, z=9/5 gives both clocks and z-1=4/5 a unit, and the same transport proves N divides both D and B.

This completes both gcd equalities at every prime-power depth. No passage from a common root modulo a prime to a root modulo a prime power is used; that would have been an unnecessary and potentially invalid extra step.

## 4. Audited short-relation bound for actual content

Set N=ceil(sqrt(p)), so N<p. Put the N+1 residues jr mod p, 0<=j<=N, into the N bins indexed by floor(Nx/p). Two fall in the same bin. Order their indices i>j and put

    b=i-j, a=-(floor(ir/p)-floor(jr/p)), m=ap+br.

Then

    1<=b<=N, -N<=a<=0, 0<|m|<p/N<=N.

The nonzero assertion follows because b<p and gcd(p,r)=1. The two common-root clocks imply

    z^m = 2^(a-b) 3^b mod g.

For the plus root z=10/13, clearing denominators gives a nonzero integer divisible by the full g:

    m>0: 10^m 2^(b-a) - 13^m 3^b;
    m<0: 13^(-m) 2^(b-a) - 10^(-m) 3^b.

The respective pairs of positive terms are bounded by (40^N,39^N) and (52^N,30^N). Prime 13 occurs on exactly one side, proving nonvanishing. Since the absolute difference of two positive numbers is strictly below their maximum,

    g_plus < 52^ceil(sqrt(p)).

For the minus root z=9/5 the integers are

    m>0: 9^m 2^(b-a) - 5^m 3^b;
    m<0: 5^(-m) 2^(b-a) - 9^(-m) 3^b.

The term bounds are (36^N,15^N) and (20^N,27^N). Prime 5 proves nonvanishing, giving

    g_minus < 36^ceil(sqrt(p)).

These are actual-content bounds on these two endpoint families, not norm bounds and not bounds for every two-edit word.

For p>=25,

    ceil(sqrt(p)) <= sqrt(p)+1 <= 6p/25,
    log(52)<4, log(3)>13/12,

so both families satisfy

    g/3^p < exp(-37p/300) < exp(-p/9).

The logarithm lower bound follows from the first two terms of 2 artanh(1/2), namely 1+1/12, with positive omitted terms. This asymptotic separation from 3^p alone is not a comparison with |2^q-3^p| and is not an integer-cycle exclusion.

## 5. Independent exact evidence

`audit_endpoint_common_root.py` imports no repository mathematical implementation and uses no Python assert statements. It builds words with floor differences, applies the actual adjacent bit swaps, and computes B separately by affine branch composition and the positional power sum.

Its complete batch contains all 1,521 coprime shapes with 3<=p<=100 and p/2<r<p, for both endpoint families:

- 3,042 endpoint paths and 6,084 physical swaps;
- 3,042 independently blocked reverse orders;
- 98 original-cut crossing swaps;
- 355,426 exact rotated-carry checks;
- 3,042 nonzero short-relation certificates divisible by full content;
- 2,528 positive-D paths and 514 negative-D paths;
- nontrivial content on 30 plus paths and 32 minus paths.

Selected ordinary-cut examples:

    plus  (p,r)=(3,2): 01011 -> 11010 -> 11001;
           B=31, D=5, g=1, edges 4 then 3.
    minus (p,r)=(3,2): 01011 -> 01101 -> 01110;
           B=38, D=5, g=1, edges 2 then 3.
    plus  (p,r)=(17,12): B=4020218554, D=407730749, g=29.
    minus (p,r)=(11,9): B=3406858, D=871429, g=13.
    plus  (p,r)=(28,17): B=612071803798318,
           D=12307579633871, g=1.
    minus (p,r)=(54,29): B=1031102309005476706302493790,
           D=-48478330446123026292740761, g=127.

Additional exact prime-power controls, outside the finite batch:

    plus  (418,237): g=17^2*41=11849;
           a=-4, b=7, m=-13,
           divisible certificate=598418218300934144.
    minus (203,165): g=13^2=169;
           a=-4, b=5, m=13,
           divisible certificate=1301138673245073.

Both extra controls also check every cut, the canonical common root, unit conditions, carry transport, and the short-relation bound. Their counters are intentionally not mixed into the stated p<=100 batch counters.

Reproduce from this audit directory:

    python3 audit_endpoint_common_root.py > receipt.replay.json
    python3 -O audit_endpoint_common_root.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

Normal and optimized production runs were byte-identical. SHA-256:

    script: fea9f463c4703dbe8adc458de87b966947ffd272a1601463d0dbd11b1f05bfd7
    receipt.json: e54b73f1de176a9ba05256b35fc5c2c1b175b580efccdea9b1db8fbb57eea9de
    receipt.optimized.json: same
    full finite-case trace:
      d5682e924748ec70fc575dc62eb901918cfb176e4afbc49f5ea024ba7fc6a25d

The trace digest includes the ordinary words, complete direct carries, gcd operands, every case's cut-check count, and the exact short-relation certificate, in the deterministic script order. The compact receipt retains selected complete examples and every nontrivial-content shape.
