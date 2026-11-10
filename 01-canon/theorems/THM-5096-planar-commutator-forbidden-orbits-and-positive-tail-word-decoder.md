---
id: THM-5096
title: "A planar Hamiltonian commutator reaches forbidden counts only after inversion; positive tail words have an exact decoder"
status: "PROVED, independently audited. Group-orbit points need not be physical; Chan-Pak completeness OPEN."
session: Session238, 2026-11-10
source: 05-knowledge/results/session238_planar_chanpak_20261110/planar_actions/PROOF.md
audit: 05-knowledge/results/session238_planar_chanpak_20261110/planar_actions_audit/AUDIT.md
---

# THM-5096: exact planar cancellation exposes the missing positive-word constraint

[Full proof](../../05-knowledge/results/session238_planar_chanpak_20261110/planar_actions/PROOF.md) · [Independent audit](../../05-knowledge/results/session238_planar_chanpak_20261110/planar_actions_audit/AUDIT.md)

For x>0, put h=I(Omega(T),x), d=I(Omega(T-u),x)>0 and R=E+O in the inherited
THM-5085 fixed-port state. After a sink, Delta=E-O=0. On the fixed-d plane,
A=S and C=S composed with B (chronological bridge then sink) act by

    A(h,R)=(h,2R), C(h,R)=(h+xR/2,4R-2d).

Positive words are actual constructions and preserve R>=2d. Adjoining
algebraic inverses instead gives, with rightmost-first products,

    D=C A^-2, K=A D A^-1 D^-1,
    J=D K D^-1 K^-1:(h,R)->(h-3xd/8,R).

D,K,J are explicit time-one polynomial Hamiltonian flows for dh wedge dR;
A,C themselves have determinants2,4 and are not symplectic. This transfers
the October19 planar cancellation operation, not an identification of its
source-adic flow group with this affine group.

At d=1, J^-8 maps the physical sink-extended singleton (1,2) to(1+3x,2),
and J^-24 maps the sink-extended triangle(1+x,6) to(1+10x,6).
Further sink padding gives targets(1+3x,16),(1+10x,24), which at x=2,3
satisfy even the inherited necessary cone inequality h<=xR/2.
Thus group orbits contain forbidden binary7,21 and ternary10,31, despite
positive integral terminal states, correct base congruence and these cone
conditions. The permanent exclusions are inherited THM-1370/4940 inputs.

Positive A,C words have a faithful normal form. If A has weight1 and C weight2,
a word of weight n sends R to2^n R-2dm, where m's binary expansion has no
adjacent ones and its highest one is at most position n-2. Every such pair
(n,m) comes from one word, by concatenating binary codes A=0,C=01.
The positive semigroup is therefore free for d!=0. From the seed R=2d,
the final R value alone determines n,m and the complete word, since

    d[(4/3)2^n+2/3,2^(n+1)]

are disjoint possible-value intervals for successive n. This concerns only
the specified sink-ending alphabet, not all S/B maps or all tournaments.
Normal/-O independent physical counts, exact rational identities and deliberate
failure controls supplement the proofs. Group completion, positive growth,
actual symmetry and modular observation are different operations.
