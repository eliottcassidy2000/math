# Validation and scope

**Same-agent exact validation, not an independent external audit.**
The derivations in PROOFS.md were checked separately against concrete
objects. No theorem ID or audited-canon promotion is requested.

Reproduce from the repository root:

```sh
python3 -B 05-knowledge/results/number_bridges_57_27_16_20261028/check.py > /tmp/bridge-normal.json
python3 -B -O 05-knowledge/results/number_bridges_57_27_16_20261028/check.py > /tmp/bridge-optimized.json
cmp /tmp/bridge-normal.json /tmp/bridge-optimized.json
cmp /tmp/bridge-normal.json 05-knowledge/results/number_bridges_57_27_16_20261028/receipt.json
```

Both modes passed and produced byte-identical output: **80,520 active
checks**, including arithmetic/elimination checks. The script uses explicit
exceptions, so optimization removes no test. Python standard library only.

Declared universes and independent computational paths:

* All512 subsets of the labelled square-sum graph Q9; literal component
  and conditional-branch bijections. Full polynomial multiplication agrees
  with enumeration. The eight-state toggle graph has10 edges, not12.
* Every simple graph through five vertices and every edge in each graph:
  5,325 edge-deletion cases. Direct set differences agree with the
  closed-neighbourhood bijection, including the dominating-edge iff.
* All25 double stars with1≤r,s≤5: exhaustive independent sets agree with
  the formula. The all-size power-of-two boundary is proved separately.
* Odd cycles of sizes3,5,9,15,27 and integer weights−3 through3:
  fraction-free determinants agree with the two-cycle-cover expression.
* 64 weighted port-elimination cases on cycles of sizes3,6,9,12, with
  heterogeneous signed edges and varying integer sources. Full rational
  Gaussian solves agree with reduced solves at named ports and with
  internal reconstruction; determinants agree independently by Bareiss.
* The golden five-cycle operator is checked as a full matrix identity and
  by the explicit five-by-two intertwiner. The height cycle is computed
  from divisibility of the declared integer rather than from a transition
  table asserted by the note.
* All1,331 triples in[−5,5]^3 pass the transform/inverse/norm identities;
  all1,331 candidate output triples pass the integer-image parity test.
  Allowing the omitted input coordinate to vary produces the explicit
  collision(5,7,8,0)/(6,8,9,1).
* The13-by13 circulant determinant is159. Among nontrivial13th roots
  modulo53, the root13 kills1+X²+X^10, proving the characteristic hostile.
* All16 elements of the Eisenstein order modulo4 are checked for the
  golden relation, yielding exactly the roots2+3ω and3+ω. The displayed
  generator change is bijective on16 coefficient pairs. Compatible lifts
  through2^16 verify the relation, invertible basis and exact multiplicative
  order by the prime divisors2 and3 of the predicted order.

Positive controls do not replace scope checks. In particular:

* Sparse cycles beyond order3 are not complete tournament supports.
* Independent-set equivalence does not certify Hamiltonian extension or
  six-colour boundary forcing.
* Determinants do not retain the transported source norm. The sharp local
  amplification21 and the inherited ROOT57/27/16 hostile remain visible.
* A16-element ring of characteristic4 is not the field F16.
* Matching reductions modulo2 do not preserve marked generator clocks at
  higher powers; the compatible ring isomorphism explicitly changes them.
* The longer owner tuples remain undefined as mathematical observables.

The repo documentation checker passed. Relative links are checked against
the Git tree as well as local files because this is a sparse worktree.
No demonstrated false previously published claim was left live; the
six-cycle correction and conductor40 correction were inherited with their
existing lineage, rather than logged again as discoveries.
