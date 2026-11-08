# Universal separated binary-compatible defect law

Status: symbolic exact proof; independent audit requested.

A local module transfer is determined by eight counts. Write s for covers with P empty and R starting at an outneighbor, h−s for P empty and R starting at an inneighbor, e for R empty and P ending at an inneighbor, h−e for R empty and P ending at an outneighbor. For both paths nonempty write a,b,c,d for (P ends in,R starts out), (in,in), (out,out), (out,in). Then

H=h, Q=2h+a+b+c+d, alpha=s+a+c.

Compare a physical module M with h times the transfer of a binary word of length k and value m, assuming Q=h2^k and alpha=hm. Their differences satisfy deltaH=deltaQ=deltaalpha=0. Put x=delta s, y=delta e, a=delta a, b=delta b. The remaining differences are delta c=−x−a and delta d=x−b. Thus every such defect is

F(x,y,a,b)=[
0 0 0 0 0 0;
x 0 0 0 0 0;
y 0 y 0 y 0;
a a+b+y a+x a+b+y a a+b+y;
−y 0 −y 0 −y 0;
−a−x −a−b−y −a−x −a−b−y −a −a−b−y].

Let W=span(e1−e5,e3−e5), using0-indexed six-state coordinates. Direct multiplication proves:

F(W)=0 for every x,y,a,b;
S0(image F), S1(image F) are subsets of W;
S0(W), S1(W) are subsets of W.

These are symbolic identities, not sample-module tests. Therefore for ANY two binary-compatible module defects Fi,Fj, and ANY nonempty singleton word with transfer C,

Fi C Fj=0.

Consequently, separated replacements of binary subwords by arbitrary binary-compatible modules satisfy the exact scalar formula

D=Hprod·D_base + sum_i (Hprod/h_i)·delta_i,

where delta_i is the scalar defect of the i-th replacement in its binary context and Hprod=product_i h_i. The same additive formula holds for the entire transfer matrix, before applying root-start and terminal vectors. Sidecar U/d remains that of the original binary word by the matched H,Q,alpha conditions.

Arithmetic boundary: since all delta_i are integers, every resulting D is divisible by

Hprod / lcm(h_1,...,h_t).

This holds for any collection of separated modules, independent of their placements and internal structure. In particular a fixed finite bank of nontrivial H-values cannot remove growing common divisibility while repeating modules. Escaping this obstruction within the separated construction needs unbounded lcm growth (for example new coprime Hamilton counts), or abandoning the separating-singleton hypothesis so that interaction terms survive. The theorem does not exclude other constructions and makes no claim of available modules with arbitrary Hamilton counts.

universal.py constructs the4-parameter transfer family from the literal category update rule and verifies the common invariant space symbolically.
