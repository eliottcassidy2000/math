# Independent audit of the sharp local distance and moment moduli

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT — PASS.** This audits the
complete [local-modulus note](long_frontier_sep06_local_modulus.md),
[exact source](../../04-computation/long_frontier_sep06_local_modulus.py),
and [frozen output](long_frontier_sep06_local_modulus.out). No correction
is needed. The distinction between convergence to the sharp ratio and
equality of a liminf, required by the initial audit, is now stated correctly.

## 1. Types and the exact global identity

The accepted domain is a finite real list with p1=p2=1 and E>0, using the
actual ordered positive coordinates and zero padding. The three matched
target atoms have value z=1/sqrt(3). When fewer than three positive
coordinates occur, additional zero slots may be inserted before matching;
negative coordinates remain in the tail. With e_i=r_i-z on the three
matched slots and e_i=r_i on the tail, p2=1 gives

    sum e_i^2=1+3z^2-2z(a+b+c)=Delta.

The moment is M=sum r_i^2(r_i-z)^2. Its coordinate expansions are

    leading: (z+e)^2 e^2-z^2 e^2=2z e^3+e^4,
    tail:    e^2(e-z)^2-z^2 e^2=-2z e^3+e^4.

Summing proves the exact signed identity

    M-Delta/3
      =2z(sum_leading e_i^3-sum_tail e_i^3)+sum_all e_i^4.

The norm inequalities sum|e_i|^3<=Delta^(3/2) and sum e_i^4<=Delta^2
then give the claimed global bound. Neither a sign assumption on the
tail nor a bound on list length is used. This audit independently checked
the sign difference between leading and tail cubic errors; replacing it
by a common sign would give an incorrect identity.

If Delta<=1/256, the relative error is at most

    1/(8sqrt(3))+1/256 <1/12.

The final strict comparison is equivalent to 96/sqrt(3)<61; both sides
are positive and 96^2<3*61^2. Subtracting and adding 1/12 from 1/3
therefore proves Delta/4<=M<=5Delta/12. The stated weak inequalities
are safe even though the actual positive-Delta comparison is strict.

M=0 forces every coordinate to be zero or z. The second moment then
forces exactly three z entries, contradicting the first moment. Delta=0
forces the same list. Thus both denominators are strictly positive on
the actual domain; the target list itself is not silently admitted as
an eligible equality case.

## 2. Local constants and uniformity

The proved local expansion was checked against the current canonical
statement **THM-4456 / sharp finite-length signed-root stability
asymptotics**, at
`01-canon/theorems/THM-4456-sharp-finite-length-signed-root-stability-asymptotics.md`.
It has status PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.
The present result uses its actual local expansion, not a reserved number
or a finite-size optimizer assumption.

By **THM-4455 / three-atom minimizing-sequence rigidity**, at
`01-canon/theorems/THM-4455-three-atom-minimizing-sequence-rigidity.md`,
R tending to K3 forces the three leading coordinates to z and the
remaining square mass m to zero. Write ell=(a+b+c)/3, h=ell-c, and
v=sum_leading(r_i-ell)^2. Ordering gives h>=0 and v<=6h^2, while
3ell^2=1-m-v. Hence

    Delta=2-2sqrt(1-m-v)=m+O(m^2+h^2).

The inherited response expansion has constants independent of the
length and tail signs:

    R-K3=A m+B h+O(m^(3/2)+h^2),
    A=-2-sqrt(6)/3+2sqrt(2)+(7/3)sqrt(3)>0,
    B=8-5sqrt(2)-4sqrt(3)+3sqrt(6)>0.

For any fixed positive epsilon, subtracting (A-epsilon)Delta leaves
epsilon*m+B*h+O(m^(3/2)+h^2). Choose a sufficiently small common
neighborhood so that the two error terms are absorbed separately into
epsilon*m and B*h. This gives the first sharp liminf bound. The exact
global remainder already gives M/Delta->1/3, yielding the second
liminf bound with constant 3A. This argument covers h^2-dominated
approaches as well as radial ones; no restriction h=O(m) is imposed
to establish the lower bounds.

The local formulation with 0<R-K3<eta follows from the minimizing-
sequence classification by contradiction: if no common neighborhood
worked as the quotient gap tends to zero, a violating sequence would
contradict the proved concentration. This supplies existence of eta,
not an explicit globally optimal threshold.

## 3. Sharpness and the correct equality quantifier

On the actual equal-three/equal-negative family, h=0, m tends to zero,
and Delta/m tends to one. Thus the quotient over Delta tends to A.
The quotient over M tends to 3A. The source independently checks the
first derivatives of the actual normalized family, with epsilon=1/n:

    Delta'(0)=(sqrt(3)-1)^2,
    M'(0)=(sqrt(3)-1)^2/3,
    (R-K3)'(0)=A*(sqrt(3)-1)^2.

The first derivative is positive, so these derivative ratios establish
actual sharp limits rather than a formal zero-over-zero substitution.

Near a minimizing sequence the three leading coordinates sum to more
than one, so p1=1 forces negative tail and m>0. If either ratio remains
bounded along a subsequence, positivity and local coercivity give

    c0(m+h)<=R-K3<=c1(m+h^2).

For the moment ratio, the same upper bound follows using M/Delta->1/3.
Absorbing the small h^2 term proves h=O(m). Division by m is then
legitimate and the error becomes o(1), giving

    (R-K3)/Delta=A+B*h/m+o(1),
    (R-K3)/M=3A+3B*h/m+o(1).

Convergence to the respective sharp constant therefore holds if and
only if h=o(m). The two-way statement is correctly restricted to actual
minimizing sequences. The alternating equal and h=1/n split families
have sharp liminf but a positive h/m limit on the split subsequence.
They are a valid normalized hostile to replacing convergence by mere
liminf equality. The final note preserves this distinction explicitly.

## 4. Exact replay and accepted scope

The entire source was read. Its fifteen explicit runtime gates check
the two coordinate identities, rational constant comparison, moment
constant, radial and variance derivatives of Delta, and the actual
normalized family's base values and first derivatives. The trace digest
records the gate labels; the raw source hash binds the formulas checked.
No disabled assertion, producer import, floating comparison, or sampled
minimizer supplies a universal quantifier.

The following independent normal and optimized replays were compared
byte for byte with the frozen output:

```
python3 -B 04-computation/long_frontier_sep06_local_modulus.py > /tmp/long-frontier-local-modulus-audit.normal.out
python3 -B -O 04-computation/long_frontier_sep06_local_modulus.py > /tmp/long-frontier-local-modulus-audit.optimized.out
cmp /tmp/long-frontier-local-modulus-audit.normal.out /tmp/long-frontier-local-modulus-audit.optimized.out
cmp /tmp/long-frontier-local-modulus-audit.normal.out 05-knowledge/results/long_frontier_sep06_local_modulus.out
```

All three outputs agree, with 15 explicit gates. Freshly computed hashes:

* Source: `4757c9d537a3d73ad3093fd94a695fc32eb915264e2b93dcb3c5a3e1e5033366`.
* Output: `4ada2b28cb6d07de820419b5f65c5b2bc6da863eb1a4c5086345713cfd793206`.
* Trace: `8bc555dee04e5483a9798da2d28b7642965a588e2f4670d092e7de67675f38f8`.

The accepted claims are the global explicit remainder, the stated small-
distance comparison, and the sharp local asymptotic modulus. No sharp
global coefficient, entire-product rate, or actual Laurent-row transport
is inferred. This audit is complete and frozen.
