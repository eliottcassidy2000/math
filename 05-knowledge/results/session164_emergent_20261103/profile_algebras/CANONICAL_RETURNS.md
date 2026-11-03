# Canonical guard returns, multiplier 91, and every 7-power residue

**Status.** Exact classification within the C/D word language started at a directed triangle, explicit return words, and a prime-power modular realization theorem. Multipliers and Hamilton counts are distinguished throughout. This does not classify arbitrary tournament counts or supply a complete ternary observer bank.

The positive one-vertex calculus and its proofs are in [GUARD_ALGEBRA.md](GUARD_ALGEBRA.md). In endpoint coordinates, the two physical moves are

    C(s,u,d)=(u,s+u+d,s),
    D(s,u,d)=(s+d,2u,0).

Here s counts Hamilton paths starting at the current guard, u counts the others, d is the unique-in-neighbor defect, and H=s+u. Both moves preserve strongness and cost one vertex. Words below are processed from left to right.

The directed triangle has state v0=(1,2,0) and H=3. A canonical state is t v0. A word has **multiplier m** if it takes t v0 to mt v0 for every positive t. Its output Hamilton count from the triangle is 3m, not m.

## 1. Exact return language

**Theorem.** A word returns v0 to its positive canonical ray if and only if it is empty, is CD, or has the form

    W D C D,                 W any finite C/D word.            (1)

**Proof.** Every reachable state has s,u>0 and slack

    ell=u−s−d >=0.

The seed has ell=1. Direct calculation gives ell(Cv)=d>=0 and ell(Dv)=u+ell>0, proving the invariant. A nonempty return must end in D, because the final defect after C is s>0. For D(s,u,d)=(s+d,2u,0) to be canonical, its input must satisfy u=s+d, that is, ell=0. Such an input cannot be the seed or an output of D, so its last move is C. Since ell(Cv)=d, the state just before that C has d=0. Only the seed and outputs of D have d=0, because C has strictly positive final defect. This proves necessity of CD or WDCD. Conversely, CD sends v0 to 3v0, and for any state with d=0,

    D(C(s,u,0))=(s+u)v0.                                    (2)

Every WD ends in defect zero, proving sufficiency. Linearity gives the same multiplier for every scalar t. ∎

There is a useful equivalent characterization of the multiplier set:

    {canonical multipliers} = {1} union {H(W v0): W any word}. (3)

Indeed CD has multiplier H(v0)=3. For every nonempty descendant word, replacing its last letter by D leaves its Hamilton count unchanged, because C and D both change H to s+2u+d. Appending CD then returns canonically with that count as multiplier by (2). Conversely WDCD has multiplier H(WD v0).

**Multiplier 11 is impossible in this seed language.** The Hamilton counts at depths 0,1,2 are respectively 3,5,9; at depth 3 they are 15 or 17. Every move strictly increases H by u+d>0. Hence no descendant has H=11, and (3) rules out multiplier 11. The same argument excludes multiplier 13. This is not a claim that tournaments with H=11 or H=13 do not exist, nor that a noncanonical C/D descendant cannot have H=33.

## 2. Two parametric return families

From the seed,

    D C² v0=(5,10,4),
    D C³ v0=(10,19,5).

For every m>=1, repeated D followed by the canonical reset CD therefore gives

    D C² D^m C D: multiplier 9+10·2^m,   length m+5,
    D C³ D^m C D: multiplier 15+19·2^m,  length m+6.           (4)

These follow because D^m(s,u,d)=(s+d,2^m u,0) when m>=1, then (2) applies. They supply, among others,

| Word | Multiplier | Output H from C3 | Output order |
|---|---:|---:|---:|
| DCCDCD | 29 | 87 | 9 |
| DCCDDCD | 49 | 147 | 10 |
| DCCCDDCD | 91=7·13 | 273 | 11 |

There is also the length-six word CCDDCD with multiplier 29. No primality assertion is made for the general expressions in (4).

For the useful multiplier 91, the complete state certificate is

    (1,2,0)
      D→(1,4,0)
      C→(4,5,1)
      C→(5,10,4)
      C→(10,19,5)
      D→(15,38,0)
      D→(15,76,0)
      C→(76,91,15)
      D→(91,182,0)=91v0.

The intermediate nine-vertex source has H=91; the final eleven-vertex canonical source has H=273. This is precisely the distinction between a source count and a reusable return multiplier.

## 3. Every residue modulo every power of seven

Combining CD (multiplier 3, two new vertices) with DCCCDDCD (multiplier 91, eight new vertices) gives canonical strong sources with

    H=3^(a+1)91^b,       order=3+2a+8b,       a,b>=0.          (5)

**Theorem.** For every k>=1 and every residue c modulo 7^k, some source in (5) has H congruent to c modulo 7^k.

**Proof.** The residue 3 has order 6 modulo 7, and

    3^6−1=728=7·104,       v_7(728)=1.

The elementary lifting identity v_7(3^(6·7^j)−1)=j+1 follows by induction using the binomial expansion of (1+7^e a)^7 with 7 not dividing a. Thus the order of 3 modulo 7^k is exactly 6·7^(k−1), the size of the unit group. Consequently its nonnegative powers represent every unit.

If c is nonzero modulo 7^k, write c=7^b u with 0<=b<k and u a unit modulo 7^(k−b). Since 91^b=7^b13^b, choose a>=0 satisfying

    3^(a+1)13^b = u mod 7^(k−b).

Such a exists by unit generation. Then (5) has the required residue. For c=0, choose b=k. All sources remain strong and canonical by the literal return identities. ∎

The earlier multipliers 1+2^j are never divisible by 7, because powers of 2 modulo 7 are 1,2,4. The multiplier 91 therefore supplies a genuinely missing valuation mechanism for that earlier restricted bank. Nevertheless (5) is only a modular universality theorem for a sparse multiplicative family, not exact realization of every integer.

## Reproduction and scope

Run `python3 return_check.py` and `python3 -O return_check.py` beside `check.py`. The focused checker verifies the return-language criterion for all 511 words of length at most 8, the depth-three count boundary, the parametric identities for m=1,…,10, and the three displayed words. It independently builds their tournaments and checks actual endpoint states by subset DP. It verifies unit orders and complete scalar residue coverage through modulus 7^5. The language classification, parametric formulas, and prime-power theorem are proved above; the computations are finite controls.

A remaining concrete question is whether a useful different seed or interface admits a canonical multiplier 11, or an affine observer return that avoids the restrictions of this particular C/D seed language.
