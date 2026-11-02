# The repaired order eleven graph target

Status: **PROVED** necessary wedge charge and isolated-vertex subcase;
**FINITE-EXACT** seeded probe; **OPEN** general graph inequality and the
order-eleven weighted tournament theorem.

The bare extension fails on P3 disjoint union 4K2, as calculated in
[STABILITY_10.md](STABILITY_10.md). The smallest possible constant wedge
repair for that witness is alpha=7/40. Its candidate endpoint polynomial is

    J11(G)=499m-320m2+80m3-20m4+5m5+80p21+56w.         (1)

Here J11=5120[P11+P21+(7/40)W-(781/1280)V], where
P11=V-M22+M222-M2222+M22222. The wedge and matching counts use exactly
the conventions of the order-ten proof.

## A proved subcase

If an eleven-vertex graph G has an isolated vertex, delete it to get G0.
The counts in (1) are unchanged. Let J10 be the integer polynomial in
the proved order-ten certificate. Direct coefficient comparison gives

    J11(G)=20J10(G0)-m+5m5+56w.                       (2)

If w>=1, then m<=45, so 56w-m>0 and J11>0. If w=0, G is a matching
with r<=5 edges, and (1) has values

    r=0,1,2,3,4,5:   0,499,678,617,376,0.

Thus J11>=0 for EVERY graph with an isolated vertex, with equality only
for the empty graph and five disjoint edges plus an isolated vertex.
Separate affinity also proves the corresponding q-box inequality whenever
one specified vertex has all incident q values zero. This is a proved
intermediate graph inequality; it does not by itself pay the entire
order-eleven cycle ledger.

## Probe and the next certificate algorithm

[graph11_scout.cpp](graph11_scout.cpp) evaluated (1) on 100,000 seeded
graphs over 19 edge densities. None was negative. Its receipt is
[graph11_scout.json](graph11_scout.json). This sample is not an exhaustive
certificate and does not prove the proposed general inequality.

There is an exact pruning rule for a future complete certificate. The
all-size vertex-addition compiler in the order-ten proof gives

    J11(G0+S)=J11(G0)+sum_(i in S) Li
                         +sum_({i,j} in S)(80Qij+56),

where

    Li=499-320m1(G0-i)+80m2(G0-i)-20m3(G0-i)
           +5m4(G0-i)+80ell_i+56d_i.

All quadratic coefficients are positive. A minimum is therefore attained
on a subset of the vertices with Li<0: deleting a selected vertex with
Li>=0 cannot increase the value. One may enumerate only those subsets for
each ten-vertex base, retaining its complete deletion responses. Equation
(2) already certifies the empty subset. This is an exact reduction of the
search, not a claim that the number of negative coordinates is uniformly
small or that the remaining enumeration has been run.

After a graph certificate, the other required step is a complete signed
order-eleven cycle ledger that can pay the 7/40 wedge term. The order-ten
reserve cannot simply be reused at a different size. Strong-target closure,
the 81/27 port control and unmarked A000568 enumeration remain separate tasks.

Reproduce the bounded probe with `g++ -O3 -std=c++17 graph11_scout.cpp -o
graph11_scout` and run the executable. No `geng` or graph atlas is used by
this scout, and no completeness claim is attached to it.
