---
id: THM-5123
title: "Common gap flags, bounded Jordan memory and universal periodic-apex recurrences"
status: "PROVED universal upper bounds and physically sharp positive-mode examples for one and two apexes; independently audited."
session: Session246, 2026-11-10
source: 05-knowledge/results/session246_unifying_principles_20261110/periodic/PROOF.md
audit: 05-knowledge/results/session246_unifying_principles_20261110/audit/AUDIT.md
---

# Shared triangular structure retains ordered switching memory

Fix any tournament A on k apexes and let remaining vertices form a transitive
tape. A letter is its neighborhood mask against A. Use the inherited positive
permutation/gap automaton, dimension k!2·3^k, from Sessions215/216/221.
For any nonempty block W of length l, let M be its full transfer matrix.

**Universal annihilator.** For every k>=0 and l>=1,

    M^(k+1) product_(j=1)^(k+1)(M−j^l I)^(k−j+2)=0.

For any fixed physical prefix u and suffix v, therefore,

    H(T_A(u W^m v))=sum_(j=1)^(k+1)P_j(m)j^(lm),
    degree P_j<=k−j+1,   m>=k+1.

The coefficients are rational; a zero-mode transient may precede m=k+1.
The corresponding universal recurrence has order at most (k+1)(k+4)/2,
independent of l. These are upper bounds, not a minimal recurrence claim.

**Mechanism.** In every nonfinal gap, use the basis
(e1−e2,e2,e0), for empty/bad/good occupancy states e0,e1,e2. All letter
assignment operators preserve its flag and have diagonal (0,1,0). In the
final gap use (e2,e0), with diagonal (1,0). Tensor rank is the sum of indices;
every off-diagonal transfer strictly lowers it. The diagonal label j counts
nonfinal index-one factors plus a final index-zero factor, independently of
letter. A length-l product consequently has diagonal j^l.

At fixed label j>=1, the possible ranks lie from j−1 to 2k−j+1 in steps two,
giving k−j+2 nonzero grades. Label zero has k+1 grades. The generalized
j^l eigenspace inherits this filtration; M−j^l I lowers its nonzero grades.
This proves the displayed nilpotent exponents and all-context readout law.

Different letters share a flag, not necessarily an eigenbasis. This extends
THM-5119 (`THM-5119-commuting-gap-projections-and-fixed-apex-run-budget.md`):
individual letters are diagonalizable, but mixed products can have Jordan
blocks and polynomial factors.

**Minimal physical mixed example.** For one apex and repeated masks (1,0),
where bit one means the apex beats the tape vertex,

    H(T((1,0)^m))=(4·4^m+6m+5)/9, m>=0.

The graph has 2m+1 vertices. Its nonzero linear term shows the distinction
is observable, not merely an unused ambient Jordan block. Zero apexes and
period length one cannot supply such a term. At two apexes with arc 0→1,
block (3,2,1) realizes maximal degrees 2,1,0 at bases 1,8,27; the leading
coefficients at bases one/eight are 13/49 and 1684/6517. The full exact
formula and recurrence/initial-value certificate are in the proof.

Physical sharpness at every positive mode is proved for k=1,2 only.
General-k sharpness and sharp transient lengths are not asserted. The
universal law makes periodic families symbolically analyzable; it does not
choose a period or context for an arbitrary desired odd count. Chan–Pak
coverage remains OPEN.

[Proof, standalone matrix certificates and literal graph controls](../../05-knowledge/results/session246_unifying_principles_20261110/periodic/PROOF.md).
