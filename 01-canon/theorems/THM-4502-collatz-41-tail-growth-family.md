---
id: THM-4502
title: "Certified long-growth Collatz family containing 27, logarithmic occurrence, and zero-dimensional dyadic closure"
status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED
source: crossroads-poset-20260926
depends_on: []
related:
  - 05-knowledge/results/crossroads_poset_20260926_bridge.md
  - 01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md
scripts:
  - 04-computation/experiments/crossroads_poset_20260926_bridge.py
outputs:
  - 05-knowledge/results/crossroads_poset_20260926_bridge.out
audit: "Root proof; independent flow and geometry count/closure audits; automata orbit/formula audit. Exceptional n=27 retained."
---

# THM-4502 — A certified family beyond 27

For the shortcut Collatz map T(n)=n/2 or (3n+1)/2, let

    p_k=2*3^(k-1),
    41*2^(r_k)=-1 mod 3^k,  0<=r_k<p_k,
    n_(k,t)=2^k(41*2^(r_k+p_k t)+1)/3^k-1,
    k>=1, t>=0.

The phase is unique and has the unique lift
`r_(k+1)=r_k+p_k d_k`, d_k in {0,1,2}. These phases tend to infinity.
Every n_(k,t) is a positive integer with k initial strictly growing shortcut
steps, followed by r_k+p_k t halvings to 41 and its certified tail to 1.
Its exact total stopping time is `69+k+r_k+p_k t`. In particular
`n_(1,0)=27`, `n_(2,0)=291`, `n_(3,0)=12439`.

The parameter map is injective. The union H satisfies

    #(H intersect [1,X])=Theta(log X).

Its 2-adic closure is exactly

    H union {(2/3)^k-1:k>=1} union {-1}.

At most O(J) cylinders modulo 2^J cover this closure, so its upper box
and Hausdorff dimensions are zero, as is its Haar measure. For fixed k,
letting B_k=2^(p_k), the family obeys the integer affine recursion

    n_(k,t+1)=B_k n_(k,t)+(B_k-1)(1-(2/3)^k).

Full elementary proofs, generalization to a known odd tail coprime to 3,
and exact controls are in [the bridge note, sections 2–3](../../05-knowledge/results/crossroads_poset_20260926_bridge.md).

**Boundary:** n=27 is the only r=0 case; it has two initial odd steps even
though k=1 prescribes only one. For all other parameters the initial odd
run has exactly k steps. Least-phase starts need not grow with k. These
are certified convergent examples, not record-holder or counterexample claims.
Neither the family nor its zero-dimensional closure proves Collatz.
