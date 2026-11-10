---
id: THM-5111
title: "Synchronized strong graft windows, exact twin-gap filters and the prime-Euclid boundary"
status: "PROVED, independently audited. All-order physical families and finite exact prime controls; prime infinitude and Chan–Pak remain OPEN."
session: Session243-prime-gap-bridges-20261110
source: 05-knowledge/results/session243_prime_gap_bridges_20261110/constructions/PROOF.md
audit: 05-knowledge/results/session243_prime_gap_bridges_20261110/construction_audit/AUDIT.md
---

# Four physical counts from a single gadget

For a nonempty tournament B, put H=H(B) and
Q=sum_X H(B[X])H(B[V(B) minus X]), with H(empty)=1. Then

    6Q-H, 6Q+H, 6Q+3H, 6Q+5H

are all Hamiltonian-path counts of explicit strong tournaments. Their orders
are |B|+4, |B|+4, |B|+5, |B|+5. No same-order claim is made for all four.

**Literal proof certificate.** Use row-bit tournaments A4=[8,1,3,6] and
A5=[20,1,2,7,14]. Their marked states (D,E,U,d), with root respectively
3,2,1,2, are

    (5,1,7,1), (5,1,9,3), (9,3,9,3), (9,3,11,5).

The inherited proper-graft formula is Ha+Qb, where
(a,b)=(D+2d-E-U,U-d). These four states have b=6 and a=-1,1,3,5.
Direct path enumeration certifies the finite states; their directed Hamilton
cycles certify strongness. Orient A-r -> B -> r. A strong source makes this
graft strong even if B is not strong, proving the all-order statement.

A transitive m-vertex gadget gives H=1,Q=2^m, hence the four consecutive odd
counts6·2^m-1,6·2^m+1,6·2^m+3,6·2^m+5, for every m>=1. Identical transitive
padding of an arbitrary fixed gadget gives infinitely many distinct strong
count pairs at gaps2H,4H,6H. Universal even-gap existence was already proved
in Sessions199/221; it is not claimed new here.

## Exact arithmetic consequences and scope

- For two fixed nonsink sources and common gadget (H,Q), synchronized padding
  gives difference H(a2-a1)+2^t Q(b2-b1). A fixed gap occurs at infinitely
  many t iff b1=b2 and H(a2-a1) is that gap. A positive twin gap forces H=1
  and a2-a1=2. This is an architecture-specific iff.
- For6·2^m±1, survival of every prime divisor through31 is exactly
  m=5 mod6, m!=1 mod5, m!=2 mod11. The wheel has40 classes modulo330.
  Simultaneous primality implies this for m>=2; m=1 is the equality exception
 11,13. At m=23 both candidates survive and are composite. The wheel is not
  a prime oracle or an infinitude theorem.
- None of the four transitive-gadget window counts is an ordinary Euclid
  number P_n+1, for n>=1. Modulo3 forces the +1 offset for n>=2, but then
  P_n=6·2^m contradicts v2(P_n)=1. This does not exclude other constructors:
  A5 rooted at3 has (a,b)=(-5,18), whose singleton graft gives31.
- For a nonempty fixed-root gadget word, the inherited product identity is
  a·product H_i+b·product Q_i.
  Requiring its symbolic intercept a·product H_i=1 forces a=1 and all
  gadgets transitive. A fixed source bank then has only finitely many
  possible primorial centers. This restriction is not an exclusion of
  incidental Euclid outputs with nonunit intercepts.

For the symmetric general pair6Q·2^t±H, gcd(H,6Q)>1 forces compositeness.
If that gcd is1, both outputs avoid the entire prime support of6QH. This
is Euclid-style escape of prime factors, not primality of the outputs.

[Full proof and exact controls](../../05-knowledge/results/session243_prime_gap_bridges_20261110/constructions/PROOF.md).
The physical graft and ordinal path-cover diagonalization are inherited from
Sessions229/227. Independent audit checks literal grafts, common-root state
laws, source boundaries, local sieve, Brown graph suppliers and all stated
hostiles. No literature-firstness, general target placement or conjecture
closure is claimed.
