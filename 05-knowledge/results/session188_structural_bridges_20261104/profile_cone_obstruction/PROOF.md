# Positive same-parity Stirling cones do not prove ambient cover domination

**PROVED OBSTRUCTION:** The proposed nonpositive ordered-cover contrast fails on positive Stirling-profile cones in both parities, even with leading generator coefficient one and an even lower coefficient. **NOT CLAIMED:** a counterexample from an actual tournament profile, an ambient tournament collision, or equality of full named responses.

## 1. The proposed cone argument

Write t^(r)_k=k!S(r,k) for the complete ordered-cover profile of TT_r, with coordinate zero at k=0. The odd-cycle expansion of Grinberg and Stanley, [The Rédei–Berge symmetric function of a directed graph](https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf), Theorem1.39, implies

    q_k(A)=sum_(odd-cycle covers Gamma) 2^(number of nontrivial cycles)
                                      k!S(number of cycles,k).

All cycle lengths are odd, so the generator index has the parity of |A|. The all-singleton cover gives t^(|A|) with coefficient one. Thus actual tournament profiles lie in a positive cone of transitive profiles of one parity. This input is classical, not a new result here.

The hoped-for argument was to polarize the cubic ordered-cover contrast from `covered_ambient/PROOF.md` and show that every same-parity generator triple has the domination sign. This is false. More strongly, the entire positive cone contains profiles with the wrong sign, so no argument using only cone membership can prove that sign for every actual tournament.

## 2. Exact functional and polarized failure

Use K,E0,R and L_q from the companion ordered-cover formula, and define

    D_h(q)=(h/3)L_q(K E0^(h-1)/R^(h+1)).

For actual equal-profile cyclic children, this is q_h(T+)-q_h(U+), where the physical root query is a=V1,c=V2. The desired domination sign is D_h<=0.

For h=2, the three odd generator profiles t^(1),t^(1),t^(5) give

    L_(t1,t1,t5)(K E0/R^3)=30>0.

The mixed coefficient of w1^2 w5 in D_2(w1 t^(1)+w5 t^(5)) is therefore +60. The three possible placements of the TT5 child in an actual cyclic triple, with the other two children TT1, have h=2 contrasts summing to exactly 60; the checker verifies these by literal full-cover counting. This is not an inference from floating-point fitting.

## 3. Positive cone counterexamples in both parities

For the odd cone, set

    q(w)=t^(5)+w t^(1)=(0,1+w,30,150,240,120).

Exact coefficient extraction gives

    D_2(q(w))=60w^2+1746120w-2910041940.                  (1)

At the positive even weight w=1582,

    q=(0,1583,30,150,240,120),
    D_2(q)=2483340>0.

Within this one-parameter family, 1582 is the first nonnegative even weight with positive D_2: (1) is strictly increasing for w>=0, its value at1580 is negative and its value at1582 is positive.

For the even cone, set

    p(w)=t^(6)+w t^(2)=(0,1+w,62+2w,540,1560,1800,720).

The h=4 contrast is

    D_4(p(w))=-192w^3+834144w^2
                         +41383812864w-223898369213472.  (2)

At the positive even weight w=8192,

    p=(0,8193,16446,540,1560,1800,720),
    D_4(p)=65543165754336>0.

The coefficient +834144 of w^2 is independently checked as the sum of the h=4 contrasts for the three actual placements of TT6 among two TT2 children. These are literal unequal-child port-core controls, not equal-profile tournament counterexamples.

Both artificial profiles have zero initial coordinate, positive finite support, one unit of the largest same-parity generator, and a nonnegative even coefficient on the lower generator. Nevertheless neither is a tournament profile: its Hamilton coordinate exceeds n!, the number of all vertex permutations. Additional relations among realizable odd-cycle counts are therefore essential. This obstruction does not refute the still-open domination conjecture restricted to actual tournament profiles.

The scalar h=1 positive-square theorem is unaffected. It says D_1(q)<0 for every nonzero finite nonnegative profile. The wrong signs here first concern higher ordered-cover coordinates.

## 4. Exact certificate and scope

`verify.py` expands each q-coordinate as an affine polynomial in w and multiplies the three profile factors directly. It records all four coefficients of the cubic functional, including zero cubic coefficient in (1), and checks the displayed evaluations exactly. It does not infer a cubic identity from four sampled values. Six actual cyclic substitutions verify the two stated mixed coefficients by independent full subset-cover counts.

Run `python verify.py > receipt.json` and `python -O verify.py > receipt.optimized.json`. The sibling `covered_response` and `response` primitives must be available. All arithmetic is integral and explicit exceptions remain enabled with -O. No full response is truncated to derive a collision: this packet concerns only the sufficiency of a proposed cone argument for ordered-cover domination.
