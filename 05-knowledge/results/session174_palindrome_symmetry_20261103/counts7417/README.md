# 7417: decimal clocks, mirrored traces, and a residue automaton

Status: elementary propositions with exact certificates; computational census is explicitly finite. No novelty claim. All code uses the Python standard library. Run `python check.py` (also `python -O check.py`); checks remain active under optimization.

## 1. Exact arithmetic and group sizes

Trial division through floor(sqrt(7417))=86 proves p=7417 prime. Its unit group has order 7416=2^3*3^2*103 and is cyclic. The element 10 has order h=2472=2^3*3*103: 10^2472=1, while 10^1236=-1=7416, 10^824=5104, and 10^24=4895 modulo p. Checking h divided by each of its prime factors proves minimality.

Thus H=<10> is the unique index-three subgroup of F_p^*, equivalently its cube subgroup. The cyclic exponent clock has 2472 states. Its dihedral symmetry group (translations of exponents and inversion) has order 4944. The full affine group AGL(1,p) has order p(p-1)=55,004,472. These are different objects: the 1237-element trace image below is a quotient set, not a subgroup of that order.

Banks–Hart–Sakata Corollary 4.4 applies: p>10 and ord_p(10)>=3 sqrt(p), certified without rounding by 2472^2=6,110,784>9p=66,753. It supplies asymptotic equidistribution over residue classes for all decimal palindromes up to x; it does not prove infinitely many prime palindromes or equidistribution within a one-parameter sparse family such as Belphegor numbers. Their general-modulus condition gcd(q,10(10^2-1))=1 fails for each of 240,1001,2222, so these clocks need separate exact handling.

Primary source: https://arxiv.org/html/math/0405056v1 (Corollary 4.4, Corollary 4.5, Theorem 5.1); bibliographic abstract https://arxiv.org/abs/math/0405056 .

## 2. The trace-coordinate theorem

Let w_r=10^r+10^(-r) in F_p. Then:

1. w_r=w_s iff r=s or r=-s modulo2472.
2. Consequently there are (2472+2)/2=1237 distinct values; r=0,...,1236 gives one of each.
3. w_r=0 iff r=618 modulo1236; these are r=618,1854 in one period.
4. w_(r+1236)=-w_r, and the sequence has least period2472.
5. w_(r+1)=(10+10^(-1))w_r-w_(r-1), with w_0=2.

Proof of (1): write x=10^r,y=10^s. Multiplying x+x^(-1)=y+y^(-1) by xy gives (x-y)(xy-1)=0. Since this is a field, x=y or xy=1, giving the exponent alternatives. Inversion on the cyclic group H fixes exactly 1 and -1, proving (2). For (3), x+x^(-1)=0 iff x^2=-1=10^1236; hence 2r=1236 mod2472, or r=618 mod1236. Parts (4) and (5) follow directly; least period follows from w_t=w_0 implying t=0 mod2472 by (1).

This is a Chebyshev/Dickson trace recurrence on a cyclic subgroup, not an asserted new theorem.

## 3. Half-digit automaton and practical counts

For a decimal palindrome of length2m+1 with m>=1 write
N=c*10^m+sum_(j=0)^(m-1) a_j(10^j+10^(2m-j)),
with a_0 in1,...,9, all other digits in0,...,9. Start C_0(0)=1. For each mirrored pair with weight v_j=10^j+10^(2m-j) modq, update
C_(j+1)(s)=sum_(allowed a) C_j(s-a*v_j).
Finally convolve with central digit weight10^m. The separate one-digit case m=0 uses central digits1,...,9, not0,...,9. This counts each palindrome exactly once in O(m*q*10) arithmetic operations and O(q) storage. It is valid even when10 is not a unit moduloq.

For q=p, scaling by10^(-m) gives normalized pair weights w_(m-j), and center weight1. Thus the trace theorem describes every local digit transition, including repeated, opposite, and zero coefficients.

Exact census (digits,total,occupied residues,min count,max count,count at residue0):
- 3,90,90,0,1,0
- 5,900,900,0,1,0
- 7,9000,6467,0,3,1
- 9,90000,7417,8,16,13
- 11,900000,7417,110,133,128
- 13,9000000,7417,1193,1233,1219

These are finite computations, not extrapolated asymptotics. The9-digit row supplies a finite all-residue coverage certificate when independently replayed. Prime counting would additionally reject divisibility by other primes and certify survivors; residue counts alone do not do that.

At distance618 a pair coefficient vanishes. The1237-digit palindromes
A=10^1236+10^618+1,
B=2*10^1236+10^618+2
are distinct but congruent to2737 modulo7417, since B-A=10^1236+1=0 modulo7417. This is the first length with an individually invisible mirrored digit pair. It is NOT the first collision among entire palindromes; those already occur at7 digits in the census. Invisible digits give a factor10 for an internal unrestricted pair, or9 for an outer nonzero pair, in the count after all other coordinates are fixed.

## 4. Three smaller clocks: transient versus unit group

CRT determines the exact groups:
- U(240) ~= C4 x C4 x C2 x C2: size64, exponent4.
- U(1001) ~= C6 x C10 x C12: size720, exponent60.
- U(2222) ~= C10 x C100: size1000, exponent100.

These group exponents must not be confused with decimal power periods:
- mod240: 10^k follows1,10,100,40,160,160,...; preperiod4, eventual period1. The stable value comes from zero mod16 and5, and1 mod3. Ten is not a unit.
- mod1001: period6 from orders6,2,6 mod7,11,13; powers1,10,100,1000,991,901. Here10 is a unit.
- mod2222: preperiod1 then period4; powers1,10,100,1000,1112,10,... . Mod2 the power becomes0, and mod1111=11*101 the order is lcm(2,4)=4. Ten is not a unit modulo2222.

The general automaton in section3 handles these transient cases without illegal inverses. Over101, normalized pair weights vanish exactly at odd distance since ord_101(10)=4. Over11 no mirrored pair weight vanishes in an odd-length palindrome:10^r+10^(-r)=2(-1)^r. Hence vanishing modulo101 does not imply vanishing modulo1111 or2222.

Repository connection: THM-583's half-system reconstruction is structurally analogous: independent half-data plus a seam determines the mirror. THM-582's reversal involution gives a different fixed-point parity statement, not a prime sieve. Prior reflection `07-reflections/the-three-sixties-and-the-erdos-592-dictionary-boxeph-S156.md` already records ord_1001(2)=60. Decimal order6 and binary order60 are intentionally distinguished here.

## Scope

The prime/order/trace proofs are elementary. The automaton is standard exact residue dynamic programming. The presented combination and examples are a research synthesis with no unsupported novelty claim. A listed Belphegor padding index2472 matching ord_7417(10) is an interesting numeric coincidence only; nothing here proves a causal or primality connection.
