"""Exact raw-frontier growth envelopes. Transition from the archived THM-4880 reference implementation."""
from itertools import combinations
import itertools
from fractions import Fraction

def letters(w):
    return [frozenset(c) for r in range(w + 1) for c in itertools.combinations(range(1, w + 1), r)]

NEW = -10  # sentinel for the endpoint at the new vertex

def step(state, letter, w):
    """all successor states (with multiplicity, as a list) of a profile state under a letter."""
    frags = list(state)
    results = []
    def finish(newfrags):
        out = []
        initial = 0
        for (h, t) in newfrags:
            if h == NEW:
                h2 = 0
            elif h == -1:
                h2 = -1
            else:
                h2 = h + 1
                if h2 >= w:
                    h2 = -1
            if t == NEW:
                t2 = 0
            elif t == -1:
                t2 = -1
            else:
                t2 = t + 1
                if t2 >= w:
                    t2 = -1
            if h2 == -1:
                initial += 1
            out.append((h2, t2))
        if initial >= 2:
            return None
        return tuple(sorted(out))
    # option A: new singleton
    s = finish(frags + [(NEW, NEW)])
    if s is not None:
        results.append(s)
    tail_ok = [(t == -1) or ((t + 1) not in letter) for (h, t) in frags]
    head_ok = [(h != -1) and ((h + 1) in letter) for (h, t) in frags]
    m = len(frags)
    # option B: append after fragment i
    for i in range(m):
        if tail_ok[i]:
            nf = frags[:i] + [(frags[i][0], NEW)] + frags[i + 1:]
            s = finish(nf)
            if s is not None:
                results.append(s)
    # option C: prepend before fragment j
    for j in range(m):
        if head_ok[j]:
            nf = frags[:j] + [(NEW, frags[j][1])] + frags[j + 1:]
            s = finish(nf)
            if s is not None:
                results.append(s)
    # option D: join i -> new -> j
    for i in range(m):
        if not tail_ok[i]:
            continue
        for j in range(m):
            if i == j or not head_ok[j]:
                continue
            nf = [frags[k] for k in range(m) if k != i and k != j] + [(frags[i][0], frags[j][1])]
            s = finish(nf)
            if s is not None:
                results.append(s)
    return results

def transfer_system(w):
    """reachable profile states, their index, the matrices M_a (dict (i,j)->mult meaning j -> i),
    the initial index and the final (one-fragment) indices."""
    L = letters(w)
    states = {(): 0}
    order = [()]
    mats = {a: {} for a in L}
    frontier = [()]
    while frontier:
        nf = []
        for st in frontier:
            i = states[st]
            for a in L:
                for s2 in step(st, a, w):
                    if s2 not in states:
                        states[s2] = len(order)
                        order.append(s2)
                        nf.append(s2)
                    j = states[s2]
                    mats[a][(j, i)] = mats[a].get((j, i), 0) + 1
        frontier = nf
    final = [states[st] for st in order if len(st) == 1]
    return order, states, mats, final


def need(x):
 if not x:raise AssertionError
def rowmult(g,M):
 w=[0]*len(g)
 for (i,j),c in M.items():w[j]+=c*g[i]
 return w
def main():
 o,_,M,F=transfer_system(3);n=len(o)
 need((n,len(F))==(119,14))
 f=[int(i in F) for i in range(n)]
 g=f[:];v=f[:]
 layers=[sum(x>0 for x in g)]
 for k in range(3):
  rows=[rowmult(v,A) for A in M.values()]
  v=[sum(row[i] for row in rows) for i in range(n)]
  g=[a+b for a,b in zip(g,v)]
  layers.append(sum(x>0 for x in g))
 need(layers==[14,67,111,119] and (min(g),max(g))==(208,5512))
 R=max(Fraction(y,x) for A in M.values() for x,y in zip(g,rowmult(g,A)))
 need(R==Fraction(1809,214))
 p,q=39,10
 w=[1]*n
 for k in range(200):
  rows=[rowmult(w,A) for A in M.values()]
  v=[max([old]+[(q*row[i]+p-1)//p for row in rows]) for i,old in enumerate(w)]
  if v==w:break
  w=v
 need(k==64 and min(w)==5 and max(w)==78)
 for A in M.values():need(all(q*y<=p*x for x,y in zip(w,rowmult(w,A))))
 need(all(x>=y for x,y in zip(w,f)))
 # Use the literal transitive-prefix sharpness witness, specified by actual arcs.
 graph=(118,100,96,55,38,64,24)
 def ham(G):
  n=len(G);d=[[0]*n for _ in range(1<<n)]
  for u in range(n):d[1<<u][u]=1
  for z,row in enumerate(d):
   for u,v in enumerate(row):
    if v:
     for k in range(n):
      if not z>>k&1 and G[u]>>k&1:d[z|(1<<k)][k]+=v
  return sum(d[-1])
 need(ham(tuple(r&63 for r in graph[:-1]))==1 and ham(graph)==39)
 print("119-state exact positive envelope: upper factor 39/10; three-step H-sum upper factor 1809/214.")
 best=0; scanned=0
 for n in range(1,7):
  for p in itertools.permutations(range(n)):
   at={x:i for i,x in enumerate(p)}
   if any(at[j]<at[i] for i in range(n) for j in range(i+4,n)):continue
   old=[sum(1<<p[j] for j in range(i+1,n)) for i in range(n)]
   old=[old[at[i]] for i in range(n)]
   for mask in range(1<<min(n,3)):
    A=old[:]+[0]
    for i in range(n):
     if n-i<=3 and mask>>(n-1-i)&1:A[n]|=1<<i
     else:A[i]|=1<<n
    best=max(best,ham(A));scanned+=1
 need(best==39)
 print(f"Independent literal subset-DP sharpness witness 1 -> 39 and {scanned} transitive-prefix cases verified.")
if __name__=='__main__':main()
