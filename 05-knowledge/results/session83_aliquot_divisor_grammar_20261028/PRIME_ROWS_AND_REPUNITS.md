# Prime-target rows, divisor counts, and repunit collisions

Status: PROVED elementary reductions; CITED published finite bound; FINITE-EXACT local check. No new solution of a repunit conjecture is claimed.

## 1. Prime-power sources of a prime target

Write F(q)={n>=1:s(n)=q}, with s the sum of proper positive divisors. A class-I source is n=p^a with p prime, a>=1. For q prime, a>=2 and

\[
q=1+p+\cdots+p^{a-1}.
\]

**Proposition.** Such an exponent a is prime and a divides q-1. Apart from q=3,p=2,a=2, the exponent a is odd.

If a=bc with b,c>1, factor the geometric sum into (1+p+...+p^{b-1})(1+p^b+...+p^{b(c-1)}), both factors >1. Hence a is prime. Also q>p and q>a. If p=1 mod q, the sum would give q|a, impossible. Since p^a=1 mod q, its multiplicative order is exactly a. Lagrange's theorem gives a|q-1. If a=2, q=p+1 can be prime only when p=2.

This supplies an exact finite inverse algorithm: enumerate prime divisors a of q-1 with 2^a-1<=q; for each, solve the strictly increasing equation 1+p+...+p^{a-1}=q by integer bisection and retain prime p only. There is at most one p per exponent a. In particular the number of class-I sources is at most the number of distinct prime divisors of q-1. This bound need not be sharp.

Examples: F(3) contains4, F(7) contains8, F(13) contains27, while F(11) has no prime-power source. This is a structural explanation for the prime-power rows, not a statement that other classes are absent from those rows.

## 2. The first collision connects 3 and 5

\[
s(2^5)=31=s(5^3).
\]

These are two *different* divisor-chain lengths: the proper-divisor sum for32 has five terms, and that for125 has three. Both remain class I under the user's priority partition. Thus the same target already forgets term count even within one class.

Equality of two geometric sums is the classical Goormaghtigh equation. Prime-target collisions cannot involve exponent2 except at3, whose source is unique. They therefore fall within its nontrivial domain (both lengths at least3). The second classical repunit example is

\[
8191=(2^{13}-1)/(2-1)=(90^3-1)/(90-1).
\]

The base90 is composite, so the second expression is **not** s(90^3). A geometric polynomial in a composite base does not list every divisor. This prime-base guard is essential.

Jon Grantham's published computation proves there are no Goormaghtigh primes below 10^700 other than31 and8191. Combining that cited result with the prime-base guard implies: **31 is the only prime target below10^700 with two distinct prime-power aliquot predecessors.** At31 these are32 and125. At8191 only8192 qualifies: the prime-exponent divisors of8190 allowed by 2^a-1<=8191 are2,3,5,7,13. Bisection gives bases8190,90,no integer,no integer,2 respectively; only2 is prime. This exhausts the possible prime-power sources, rather than merely discarding the known base90 representation. This is a consequence of a cited computation, not a computation repeated here and not an unrestricted uniqueness theorem.

Source: Jon Grantham, *No New Goormaghtigh Primes up to 10^700*, Integers24(2024), A98, §§1–2, https://math.colgate.edu/~integers/y98/y98.pdf . The paper also gives an abc-conditional finiteness result for Goormaghtigh primes. That supplies a genuine conditional connection with prior abc work; it does not prove abc, classify all fibers, or imply that an untouchable-prime conjecture is equivalent to abc.

## 3. Distinct questions at a prime row

A class-II predecessor pq, p<q prime, has s(pq)=p+q+1. For a prime target r, it is therefore a distinct-prime Goldbach representation of r-1. A class-I predecessor instead requires the prime-base repunit equation above. Class III supplies the remaining exponent shapes. Missing one mechanism does not imply an empty whole fiber.

For example11 has no class-I predecessor but21=3*7 gives s(21)=11. The source343 in the verified1001 chain is class I, but its target57 is composite; the prime-target proposition does not apply to57. Its square-sum property343+57=400 is instead the square-sigma condition inherited from Session82.

## 4. Reproduction and limits

Run `python3 check_prime_rows.py` and `python3 -O check_prime_rows.py`.

The stdlib-only checker enumerates all prime-power geometric sums whose target is at most100000 and compares every prime target with the independent exponent/order-and-bisection algorithm. All9592 prime rows agree;31 is the sole collision in that finite window. It checks the composite-base8191 hostile example. No enumeration of all targets up to10^700 is performed locally.
