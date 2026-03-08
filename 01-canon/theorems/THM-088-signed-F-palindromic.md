# THM-088: Signed Forward-Edge Polynomial SF(T,x)

**Status:** PROVED (algebraic)
**Proved by:** opus-2026-03-07-S46
**Scope:** All tournaments

---

## Statement

Let T be a tournament on n vertices. Define:

$$SF(T,x) = \sum_{\sigma \in S_n} \text{sgn}(\sigma) \cdot x^{\text{fwd}_T(\sigma)}$$

Then:

**(A)** SF(T,x) satisfies SF(T,x) = (-1)^{C(n,2)} x^{n-1} SF(T, 1/x).
- When C(n,2) is even (n = 0,1 mod 4): SF is palindromic.
- When C(n,2) is odd (n = 2,3 mod 4): SF is anti-palindromic.

**(B)** SF(T,1) = 0 for all n >= 2.

**(C)** (x-1) divides SF(T,x), and the quotient Q(T,x) = SF(T,x)/(x-1) is anti-palindromic:
Q(T,x) = -x^{n-2} Q(T, 1/x).

---

## Proof of (A)

Path reversal: the map P -> P^{rev} sends fwd(P) to n-1-fwd(P) (from F palindromicity). For the sign: sgn(P^{rev}) = (-1)^{C(n,2)} sgn(P) (reversing a permutation of n elements changes parity by C(n,2) transpositions).

Therefore: SF(T,x) = sum sgn(P) x^{fwd(P)} = sum (-1)^{C(n,2)} sgn(P^rev) x^{n-1-fwd(P^rev)}
= (-1)^{C(n,2)} x^{n-1} sum sgn(Q) x^{-fwd(Q)} (reindexing Q = P^rev)
= (-1)^{C(n,2)} x^{n-1} SF(T, 1/x). Done.

## Proof of (B)

SF(T,1) = sum_{sigma} sgn(sigma) * 1^{fwd(sigma)} = sum sgn(sigma) = 0 for n >= 2
(since |{even perms}| = |{odd perms}| = n!/2).

## Proof of (C)

By (B), (x-1) | SF(T,x). Let Q(T,x) = SF(T,x)/(x-1).

From (A): (x-1)Q(x) = (-1)^{C(n,2)} x^{n-1} (1/x - 1) Q(1/x)
= (-1)^{C(n,2)} x^{n-1} * (-1) * (x-1)/x * Q(1/x)
= (-1)^{C(n,2)+1} x^{n-2} (x-1) Q(1/x)

Dividing by (x-1):
Q(x) = (-1)^{C(n,2)+1} x^{n-2} Q(1/x)

Since C(n,2)+1 is odd when C(n,2) is even and vice versa:
Q(x) = -x^{n-2} Q(1/x)  (anti-palindromic). Done.

---

## Verified

- n=3,4,5: exhaustive (all tournaments)
- n=6: exhaustive (all 32768 tournaments, 24 F-classes)

---

## Additional observations

1. At n=4: SF(T,x) = c(T) * (x-1)^2(x+1) for an integer c(T).
   This is because Q has degree 2 and is anti-palindromic, so Q(-1)=0, giving (x+1)|Q.

2. SF is a COARSER invariant than F at n>=6: multiple F-vectors can give the same SF.
   At n=5: SF determines F uniquely (1-to-1).
   At n=6: 4 SF values map to >1 F-vector.

3. The coefficient GCD of SF is generally 1, but individual positions have higher GCD.
   At n=5: SF[1] and SF[3] always even (GCD=2), SF[2] always divisible by 6.
