#!/usr/bin/env python3
"""Cyclic-triple hypergraph cycle index from recovered prime equivariant indices."""
import recover_converse as r
import json

def recover_c3():
 *_,pe,pt=r.recover()
 root={p:v/2 for p,v in r.add(pe,pt).items()}
 c={():r.F(1)};d={}
 for n in range(1,r.N+1):
  candidate=r.add({(1,):r.F(1)},r.comp(root,r.add(c,{():1},-1)))
  d.update({p:v for p,v in candidate.items() if sum(p)==n})
  logarithm={}
  for j in range(1,n+1):logarithm=r.add(logarithm,r.psi(d,j),r.F(1,j))
  # Weighted Euler derivation: n C_n = [ (sum_k k L_k) C ]_n.
  dc={p:sum(p)*v for p,v in logarithm.items()}
  update=r.mul(dc,c)
  c.update({p:v/n for p,v in update.items() if sum(p)==n})
 if r.mul(dc,c)!={p:sum(p)*v for p,v in c.items() if p}:raise RuntimeError('SET exponential equation')
 if d!=r.add({(1,):1},r.comp(root,r.add(c,{():1},-1))):raise RuntimeError('prime root equation')
 if any(v<0 for a in (c,d,root) for v in a.values()):raise RuntimeError('negative')
 if any(sum(v for p,v in a.items() if sum(p)==n).denominator!=1 for a in (c,d,root) for n in range(r.N+1)):raise RuntimeError('nonintegral')
 expected=[1,1,1,2,3,9,32,247,3328,93612,4774417]
 if r.counts(c)[:min(r.N+1,11)]!=expected[:min(r.N+1,11)]:raise RuntimeError('prior counts')
 return c,d,root
if __name__=='__main__':
 c,d,root=recover_c3()
 print(json.dumps({'N':r.N,'C3_isotypes':r.counts(c),'prime_root_or_singleton_C3_isotypes':r.counts(d),'prime_C3_isotypes':r.counts(root),'C3_cycle_index':r.serial(c),'D_cycle_index':r.serial(d),'checks':'PASS nonnegative coefficients, integral totals, prime root equation, weighted Euler exponential, prior supplied-input counts through10'},indent=2))
