#!/usr/bin/env python3
"""coloring_tutte_atlas_s682.py — everything as a coloring; tie-induction = the
deletion-contraction (Tutte/Potts) recursion.

Prompt: pursue tie-induction; see everything as graph colorings (some nodes, some
edges, some both); reframe as many problems as possible.

THE UNIFICATION: the deletion-contraction / Tutte-Potts polynomial is the common
frame. 'some NODES' = chromatic polynomial (vertex colorings); 'some EDGES' = flow
polynomial (edge flows); 'BOTH' = the full Tutte plane / Potts model. Tie-induction
(delete an edge => introduce a TIE => oriented graph) IS the deletion step of the
Tutte recursion.

VERIFIED HERE: Hamiltonian-path deletion-contraction H(T) = H(T del e) + H(T/e)
(paths split by whether they use the directed edge e). Deleting e makes it a tie;
contracting merges its endpoints. So computing H by repeatedly DELETING edges
(introducing ties) is the Tutte-shaped recursion bottoming out at tie-graphs.

THE COLORING ROSETTA (repo problems as colorings; node / edge / both):
  problem            colored elements        coloring object            DC/Tutte tie
  -----------------  ----------------------  -------------------------  ------------------
  LRC                NODES (runners)->R/Z    circular chromatic /        tension<->flow =
                                             tension; dual = nowhere-     Tutte chromatic-
                                             zero FLOW on edges (S537o)   flow duality
  Tournament / H     EDGES -> {->,<-}        Ham-path count (Redei-      H=H(del e)+H(/e)
                     (ties = 3rd color)      Berge)                       [verified]
  metagraph G_n/Z2   NODES -> chi=n-1        chromatic number            chromatic poly
  conflict graph Om  NODES -> independent    H=I(Omega,2) hard-core /     Potts/Tutte
                     sets                    independence polynomial      specialization
  unit distance      NODES (plane pts)       Hadwiger-Nelson chi in       chromatic (plane)
                                             {5,6,7}; hexagonal 7-color
  Collatz            NODES (integers)->      parity 2-coloring + map      (weak fit)
                     {even,odd}
  partition fns      (HYP-2183)              Potts Z = q-coloring count   Tutte (x,y) plane

THE GEM: tie-induction (request 1) and everything-as-coloring (request 2) are ONE
object -- the Tutte/Potts polynomial under deletion-contraction. node/edge/both =
chromatic/flow/full Tutte; LRC tension-flow duality (S537o) = Tutte chromatic-flow
duality; H's DC = the Tutte recursion; 'partition functions everywhere' (HYP-2183)
= the Potts partition function. The 'useful tie-family' (HYP-2261) = the
intermediate oriented graphs of the DC recursion.

Session: claude-2026-06-06-S682 (coloring-tutte-atlas)."""
import sys; sys.stdout.reconfigure(line_buffering=True)
import random
def H(n,out):
    if n<=1: return 1
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if c:
                for w in out[v]:
                    if not(mask>>w&1): dp[mask|1<<w][w]+=c
    return sum(dp[size-1])
def tour(n,bits):
    out=[set() for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if bits[idx]: out[i].add(j)
            else: out[j].add(i)
            idx+=1
    return out
rng=random.Random(5); n=6; ok=tot=0
for _ in range(3000):
    bits=[rng.randint(0,1) for _ in range(n*(n-1)//2)]; out=tour(n,bits); a=0
    if not out[a]: continue
    b=min(out[a]); outd=[set(s) for s in out]; outd[a].discard(b); Hdel=H(n,outd)
    keep=[v for v in range(n) if v!=b]; idx={v:i for i,v in enumerate(keep)}; m=len(keep)
    outm=[set() for _ in range(m)]
    for v in keep:
        if v==a:
            for w in out[b]:
                if w!=a: outm[idx[a]].add(idx[w])
        else:
            for w in out[v]:
                if w==a: outm[idx[v]].add(idx[a])
                elif w!=b: outm[idx[v]].add(idx[w])
    if H(n,out)==Hdel+H(m,outm): ok+=1
    tot+=1
print("TIE-INDUCTION = DELETION-CONTRACTION (Tutte recursion):")
print(f"  H(T) == H(T del e) + H(T contract e):  {ok}/{tot} hold  [delete e = introduce a TIE]")
print("\nEVERYTHING AS COLORING -- node / edge / both (see docstring Rosetta):")
print("  LRC = node circular-coloring (tension) <-> edge nowhere-zero flow (Tutte duality, S537o)")
print("  Tournament/H = edge {->,<-} coloring; H deletion-contraction [verified]")
print("  metagraph chi=n-1, unit-distance Hadwiger-Nelson chi in {5,6,7}: node chromatic")
print("  conflict graph Omega / partition functions: Potts/Tutte (q-coloring) specializations")
print("\nGEM: tie-induction + everything-as-coloring = ONE object: the Tutte/Potts polynomial")
print("  under deletion-contraction. node/edge/both = chromatic/flow/full Tutte.")
