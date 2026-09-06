---
id: THM-4455
title: "Three-atom minimizing-sequence rigidity for signed-root stability"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
depends_on:
  - THM-4454 / sharp-global-signed-root-duplication-stability
source: long-frontier-sep06 minimizing-sequence classification
primary_script: 04-computation/long_frontier_sep06_minimizers.py
primary_output: 05-knowledge/results/long_frontier_sep06_minimizers.out
primary_script_sha256: 6d3be5780138c24f58312733270981d3375eae93b9281b26710529bf955dc7cd
primary_output_sha256: 19d6c54228a78db3b7231e30e9d6897bdc9a71ba904dd5159b9c861ab2040b20
report: 05-knowledge/results/long_frontier_sep06_minimizers.md
audit: 05-knowledge/results/long_frontier_sep06_minimizers_audit.md
hash_basis: raw LF repository bytes
---

# THM-4455 — three-atom minimizing-sequence rigidity

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[complete proof and exact boundary analysis](../../05-knowledge/results/long_frontier_sep06_minimizers.md)
and [independent proof/source audit](../../05-knowledge/results/long_frontier_sep06_minimizers_audit.md)
are part of this theorem. The sole mathematical dependency is
[THM-4454 / sharp global signed-root duplication stability](THM-4454-sharp-global-signed-root-duplication-stability.md).

Let r^(n) be arbitrary finite real lists, of possibly varying lengths,
with sum r_i=sum r_i^2=1 and E=(1-p4)/2>0. Let a>=b>=c>=0 be the
three largest positive entries, padded with zeros. Put

    J=(5-8p3+3p4)/(6E), c_*=(13-8sqrt(2))/3,
    d2=2-sqrt(2)(a+b), R=(J-c_*)/d2,
    K3=4sqrt(3)/[3(1+sqrt(2))(1+sqrt(3))],
    Delta3=2-(2/sqrt(3))(a+b+c).

Then R(r^(n))->K3 **if and only if** Delta3(r^(n))->0. Equivalently,
after permutation and zero padding, the lists converge in square norm
to (1/sqrt(3),1/sqrt(3),1/sqrt(3),0,...). In particular the three
leading positive entries converge to 1/sqrt(3), and the square mass
of every remaining coordinate tends to zero.

There are no other minimizing boundary shapes. Near one positive unit
atom, R tends to sqrt(2)-2/3. Near two equal positive unit-norm atoms,
liminf R is at least (64-44sqrt(2))/3. Both constants exceed1/2>K3.
These denominator-sensitive expansions exclude the two extraneous
zeros of the unnormalized objective in THM-4454. Its remaining envelope
zero forces the top two roots to 1/sqrt(3); a strictly increasing signed
secant then forces the third positive root to capture the entire
remaining square mass. Continuity proves the converse.

If N_minus is the number of negative coordinates, every eligible list
satisfies the independent count consequence

    N_minus*Delta3 >=
      max(sqrt(3)-1-(sqrt(3)/2)*Delta3,0)^2.

Thus N_minus tends to infinity along every minimizing sequence. The
signed sum outside the three leading roots tends to 1-sqrt(3).
The positive and negative sums separately need not converge: the proof
constructs normalized minimizing lists with n^4 positive and n^4
negative small entries whose separate magnitudes both diverge.
Square-norm convergence must not be promoted to separate first-moment
convergence.

The exact checker verifies91 always-active gates, including the uniform
boundary identities and seven declared actual mixed-dust controls.
Normal and optimized output agree, as independently replayed. The
universal argument is analytic and does not extrapolate a finite bank.
This settles the minimizing-sequence question left in THM-4454; it
does not settle finite-length sharp constants or supply actual Laurent
coefficient transport. No external priority or proof-assistant claim
is made.
