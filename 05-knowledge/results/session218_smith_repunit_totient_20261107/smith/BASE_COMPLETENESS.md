# Complete existence of a Smith base for composite integers

**Status: PROVED**, with independent finite controls. No historical novelty claim. This extends the earlier Smith/repunit lane without changing that frozen artifact.

## Theorem

A composite positive integer n is a Smith number in at least one integer base b>=2 **if and only if**

`n not in {6,8,9,12}`.

An explicit base is:

- n=4: any b>=5;
- every other composite n outside the four exceptions: `b=n-F(n)+1`, where `F(n)=sum_p v_p(n)p` counts prime factors with multiplicity. For n=2p with p>=5 prime this is b=p-1; for all other eligible n it is larger than n/2.

These cases cover all composites except the four excluded integers. Prime integers are not Smith numbers by definition, even though they satisfy digit-sum equality.

## Proof

Write s_b(n) for digit sum and A_b(n)=sum_p v_p(n)s_b(p). Smith means composite and s_b(n)=A_b(n).

### Factor-weight bounds

For all integers m>=2, `F(m)<=m`, with equality for composite m only when m=4. Indeed, the claim is equality for primes. For a composite split m=uv with u,v>=2, induction gives F(m)=F(u)+F(v)<=u+v<=uv. The last inequality is equivalent to `(u-1)(v-1)>=1`, and equality forces u=v=2.

For composite m>=4, the stronger bound `F(m)<=m/2+2` holds: use the same split and `(u-2)(v-2)>=0`, which gives u+v<=uv/2+2.

Now suppose n>12 is composite and not twice a prime.

- If n is odd, choose a nontrivial factorization n=uv. Both u,v>=3, so `(u-3)(v-3)>=0` gives u+v<=n/3+3. Consequently `F(n)<=n/3+3<(n+2)/2`; the strict inequality is exactly n>12.
- If n is even, write n=2m. Since n is not twice a prime, m is composite. Thus `F(n)=2+F(m)<=2+m/2+2=n/4+4<(n+2)/2`; again strictness is exactly n>12.

### The constructed two-digit representation

Set `b=n-F(n)+1`. The preceding bound gives b>n/2. Also b<n, since F(n)>1. Therefore b>=2, every prime factor p of n satisfies p<=n/2<b, and

`n=b+(F(n)-1)`, with `0<=F(n)-1<b`.

The latter strict inequality is equivalent to `2F(n)<n+2`, already proved. Hence the base-b representation of n is the two-digit word `(1,F(n)-1)`, and

`s_b(n)=1+(F(n)-1)=F(n)=A_b(n)`.

The equality A_b(n)=F(n) holds because every prime factor is a single base-b digit. No primality test for b is required.

### Remaining constructive cases

For n=4 and b>=5, all relevant numbers are single digits, so s_b(4)=4=2+2=A_b(4).

If n=2p with p>=5 prime, take b=p-1>=4. Then p has digits11, 2p has digits22, and 2 is a single digit. Thus s_b(2p)=4=2+2=A_b(2p). Here b=p-1=n-F(n)+1 is the same formula as before. The proof is separate because b>n/2 and the digit representation (1,F(n)-1) fail for twice-prime n; the valid representation is22.

The composite integers <=12 are exactly 4,6,8,9,10,12. The integer10=2*5 is already constructed, leaving precisely the four claimed exceptions.

### Complete exclusion of the exceptions

For any base-b Smith number, reduction modulo b-1 gives

`b-1 | n-F(n)`.

This follows from s_b(x)=x mod(b-1), applied both to n and to each prime factor with multiplicity. For the four candidates n-F(n)>0, so its positive divisors enumerate every possible b, with no upper cutoff needed.

| n | F(n) | n-F(n) | possible b | (s_b(n),A_b(n)) in that base order |
|---|---:|---:|---|---|
| 6 | 5 | 1 | 2 | (2,3) |
| 8 | 6 | 2 | 2,3 | (1,3), (4,6) |
| 9 | 6 | 3 | 2,4 | (2,4), (3,6) |
| 12 | 7 | 5 | 2,6 | (2,4), (2,7) |

None has equality. This proves both directions of the theorem.

## What this changes and does not change

The unqualified question “can these composites be Smith numbers in some bases?” is now completely settled by this elementary classification. Moving the base is a very weak compatibility requirement: it accepts every composite except four small exceptions. It cannot distinguish Camion pseudoprime exponents, classical Carmichael numbers, nontotients, tournament spectra, or an arbitrary composite sequence beyond those exceptions.

Conditional application to any independently verified composite even nontotient n: it has a Smith base unless n is one of6,8,12. (The actual nontotient predicate must still be checked; this theorem does not prove any integer is a nontotient.) Every composite n>12, whether even or odd and regardless of its other arithmetic predicates, is covered.

This theorem varies the base with n. It proves neither infinitude of Smith numbers in a specified fixed base nor any fixed-base assertion about repunit primes. The b=n-F(n)+1 construction also does not require n to contain a repunit prime factor.

## Verification

Run `python3 check_base_completeness.py` in this artifact directory. The script independently builds smallest-prime-factor data, constructs a Smith base for every eligible composite n<=100000, and checks every possible base from the divisibility sieve for composites n<=5000. It also exhaustively checks every base 2<=b<=n+1 for composites n<=300 to verify the sieve through a distinct enumeration path. These finite checks support, but do not replace, the proof.
