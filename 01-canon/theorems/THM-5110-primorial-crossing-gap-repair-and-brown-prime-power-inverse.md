---
id: THM-5110
title: "Primorial crossing-gap repair and the exact Brown prime-power inverse"
status: "PROVED, independently audited; finite arithmetic controls separate from all-order proofs. Named prime conjectures and Chan–Pak remain OPEN."
session: Session243-prime-gap-bridges-20261110
source: 05-knowledge/results/session243_prime_gap_bridges_20261110/provenance/PROOF.md
audit: 05-knowledge/results/session243_prime_gap_bridges_20261110/arithmetic_audit/REPORT.md
---

# Exact arithmetic bridges, with the missing predicates retained

Let P_n be the product of the first n primes, n>=2. Let F_n and L_n be the
least offsets d>1 for which P_n+d and P_n-d, respectively, are prime.

1. Restore the two skipped offset-one bits: the nearest prime above P_n is
   at distance1 if P_n+1 is prime, and F_n otherwise; the nearest below is
   at distance1 if P_n-1 is prime, and L_n otherwise. Their distance sum is
   the actual consecutive-prime gap crossing P_n.
2. If F_n=L_n=d, the two offset-d primes are consecutive exactly when both
   P_n-1 and P_n+1 are composite. These are the only possible intervening
   primes. In particular30030±17 are not consecutive:30029 intervenes and
   the crossing gap is18. In contrast9699690±23 are consecutive, gap46.
3. For each fixed integer d>1 and fixed positive integer c, cP_n±d are both composite for
   every sufficiently large n. A prime divisor r of d eventually divides
   P_n, and cP_n-d>r eventually makes it a proper divisor of both endpoints.
   Thus this primorial-center template cannot establish Polignac at any
   fixed non-twin symmetric gap2d. The case d=1 is expressly excluded.
4. For n>=4 put C=n!/4. Positive Brown solutions n!+1=m² are in bijection
   with partitions of the full prime-power factors of C into coprime
   products A,B with |A-B|=1. Name A as the even factor. Then m=A+B, and
   exactly2^(pi(n)-1) partitions are candidates. This is a finite exact
   inverse at each n, not a proof that only n=4,5,7 work.

Proofs are elementary. Parts1–2 follow by partitioning interior offsets
into0,1 and2,...,d-1. Part3 is the proper-divisor argument just given. For
part4, set A,B to the even and odd members of (m-1)/2,(m+1)/2. Conversely,
(A+B)²-4AB=(A-B)²=1. Coprimality places each full prime power on one side.
The complete arguments, boundary cases and connection contracts are in the
[prime-gap proof](../../05-knowledge/results/session243_prime_gap_bridges_20261110/provenance/PROOF.md)
and [Brown inverse](../../05-knowledge/results/session243_prime_gap_bridges_20261110/brown/PROOF.md).

The Brown sign conditions retain odd prime-power exponents and the exact
2-adic exponent e-1, where e=v_2(n!), but imply only n! | m²-1. They lose
the quotient needed for equality. The positive examples5,11,71 are roots,
not their squares. Prime-square Brocard is a different conjecture.

The inherited source is Session51's Fortunate wheel and the September27
factorial divisor involution. No literature-firstness is claimed. The new
synthesis isolates the exact two-bit repair and the complete packet inverse.
Primality of a value does not supply a tournament: Euclid number7 is the
canonical hostile. Nor can a tournament witness for49 be divided or rooted
arithmetically to produce the forbidden count7.

Independent controls cover all484,682 canonical partitions for4<=n<=60,
complete sign residue classes in the declared small moduli, primorial
indices2..8, and fixed-offset proper-factor certificates. Full finite scopes,
source pins, normal/optimized replays and failing controls are in the audit.
No prime-gap infinitude or density result, named prime conjecture, or
Chan–Pak completeness follows.
