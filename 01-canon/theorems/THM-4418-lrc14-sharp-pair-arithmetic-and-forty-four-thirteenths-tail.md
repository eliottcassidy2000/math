---
id: THM-4418
title: "LRC14 sharp pair arithmetic and forty-four-thirteenths tail"
status: >
  PROVED ANALYTICALLY with FINITE-EXACT rational residue obligations and
  independent proof audit. For positive distinct odd a<b<c, all prime to
  three, ordered-pair mass T_ab<=12/77 and component count N_ab<=24b/13
  are sharp. Every such triple with c/b>=44/13 has strict degree-zero
  network certificate below 6/77. Comparable-speed universal height,
  physical entry, synchronization, and LRC(14) remain OPEN.
source: synthesis-sep05 / Bernoulli pair arithmetic and separated-component tail
depends_on:
  - THM-4409-lrc14-third-sheet-component-network-certificate
proof: 05-knowledge/results/lrc14-component-separation-synthesis-sep05.md
script: 04-computation/lrc14_component_separation_synthesis_sep05.py
output: 05-knowledge/results/lrc14_component_separation_synthesis_sep05.out
script_sha256: 3164de8aeab6ae2376fd057e86460fd26bc71574100669299d0e2c372ecc72dc
output_sha256: e52cc962ce15528a276d5bfa84019b10a622ff069fb487c43d54cf86fbb7d8c3
hash_basis: raw LF bytes
---

# THM-4418 -- sharp pair bounds and an infinite local LRC tail

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [full proof and verification](../../05-knowledge/results/lrc14-component-separation-synthesis-sep05.md)
also credits concurrent independent recovery of the separated-contact
identity in [the incoming synthesis](../../05-knowledge/results/synthesis_20260905_lrc_sparse_transport.md).
That independent identity is now also canonical in
[THM-4414](THM-4414-lrc14-six-separated-contact-capacity-collapse.md).
The proof linked here retains its own elementary derivation.

Use D_(v,s)={x in R/Z: ||v(x+s/3)||<1/14}. For distinct positive odd
a<b, neither divisible by three, let T_ab and N_ab be the total mass and
number of connected components of the six distinct-sheet pair intersections.
If p=a/g,q=b/g,g=gcd(a,b), then

```text
T_ab <=12/77,     equality iff (p,q)=(1,11);
N_ab <=24b/13,    equality iff (p,q)=(11,13).
```

For any odd c>b not divisible by three, THM-4409's degree-zero certificate
selecting pair (a,b) satisfies

```text
c/b>=44/13  =>  U_net(a,b;0,0)<6/77.
```

The proof retains the exact sheet phase 1/3. The ordered-pair count is
6g times the number of integers in the open interval of radius (p+q)/14
centered at 1/3. A seven-residue calculation bounds this count sharply.
Periodized interval roofs yield a Bernoulli formula for T_ab; its 49
residue cases plus four small reduced pairs give the sharp mass constant.

Separation makes each contact graph a star forest and its degree-zero
capacity the sum of edgewise minimum lengths. For b<c<6b,
U<=N_ab/(7c)<=24b/(91c); equality at c/b=44/13 is excluded by oddness.
For c>=6b, tooth-center counting gives
U<=T_ab/7+8N_ab/(49c)<=508/7007<6/77.

The independent arithmetic/geometric replay covers all 2,910 primitive
odd ternary-unit triples through height 79, retains the equality comb
(1,5,11), and retains strict quotient loss on (1,19,79).
It performs 216,039 explicit exact checks, also under Python optimization.

The residual local region c/b<44/13 still has unbounded primitive height.
This theorem supplies no global LRC entry map or fourteen-runner closure.
