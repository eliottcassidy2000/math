#!/usr/bin/env python3
"""Shared-bit edge/order constraints, weighted cycle sums, exact minimization."""
import itertools,json
from order_atlas import Atlas,Orders,edges,need

def compile_joint(o,edge_queries=(),order_queries=()):
 out=o.cube.compile(edge_queries)
 if out is None:return None
 pins=dict(out['pins'])
 for order in order_queries:
  q=o.compile(order)
  if q is None:return None
  for b,z in q['pins'].items():
   if b in pins and pins[b]!=z:return None
   pins[b]=z
 return {'pins':pins,'free_bits':o.cube.bits-len(pins),'one_code':''.join(str(pins.get(b,0)) for b in range(o.cube.bits))}

def partition(cube,pins,weights):
 def product(es):
  ans=1
  for e in es:ans*=weights.get(e,1)
  return ans
 fixed=cube.base-set().union(*(r for r,a in cube.moves.values()))
 out=product(fixed)
 for b,(r,a) in cube.moves.items():
  pair=(product(r),product(a));out*=pair[pins[b]] if b in pins else sum(pair)
 return out

def cheapest(cube,pins,costs):
 code=['0']*cube.bits
 for b,(r,a) in cube.moves.items():
  costs01=[sum(costs.get(e,0) for e in part) for part in (r,a)]
  code[b]=str(pins.get(b,int(costs01[1]<costs01[0])))
 code=''.join(code)
 return code,sum(costs.get(e,0) for e in cube.edge_set(code))

def inversion(r):
 p={v:j for j,v in enumerate(r['new'])}
 return next(list(q) for q in itertools.combinations(r['old'],2) if p[q[0]]>p[q[1]])

def main():
 c=Atlas().cycle(2812,'A');o=Orders(c)
 r=o.regions[0];pair=inversion(r);b=r['bit'];added=min(c.moves[b][1])
 need(o.compile(pair) is not None and c.compile([(added,True)]) is not None,'separate contexts feasible')
 need(compile_joint(o,[(added,True)],[pair]) is None,'shared bit rejects jointly incompatible edge and order')
 other=next(s for s in o.regions if s['bit']!=b);order=inversion(other)
 need(compile_joint(o,[(added,True)],[order]) is not None,'independent source contexts compose')
 # A 64-cycle fibre, independently summed from its actual edge sets.
 free=list(range(6));pins={b:0 for b in range(c.bits) if b not in free};weights={};costs={}
 for b in free:
  for e in c.moves[b][0]|c.moves[b][1]:weights[e]=1+(e[0]+3*e[1])%4;costs[e]=(7*e[0]+e[1])%17-8
 actual=0;minimum=None;records=[]
 for word in itertools.product('01',repeat=len(free)):
  bits=['0']*c.bits
  for b,z in zip(free,word):bits[b]=z
  bits=''.join(bits);es=edges(c.walk(bits),True);w=1
  for e in es:w*=weights.get(e,1)
  actual+=w;cost=sum(costs.get(e,0) for e in es);minimum=cost if minimum is None else min(minimum,cost)
 expected=partition(c,pins,weights);need(actual==expected,'weighted Hamilton polynomial specialization')
 code,cost=cheapest(c,pins,costs);need(cost==minimum,'exact independent-bit minimum')
 # Two conflicting local orders in a tied coarse source become compatible after refinement.
 a=Atlas();coarse=Orders(a.cycle(70312,'A'));fine=Orders(a.cycle(70312,'AA'))
 first=coarse.regions[0];second=next(r for r in coarse.regions if r['bit']==first['bit'] and r is not first)
 q0=inversion(first);q1=inversion(second)[::-1]
 need(compile_joint(coarse,order_queries=[q0,q1]) is None,'tied order obstruction')
 out=compile_joint(fine,order_queries=[q0,q1]);need(out is not None and out['free_bits']==623,'refined joint-order extension')
 for q in [q0,q1]:need(fine.observed(out['one_code'],q)==q,'actual refined order witness')
 print(json.dumps({'status':'passed','mixed_edge_order_hostile':{'order':pair,'required_added_edge':added},'weighted_complete_fibre':64,'weighted_sum':actual,'minimum_signed_edge_cost':minimum,'coarse_fine_order_hostile':{'N':70312,'ordered_pairs':[q0,q1],'coarse_bits':25,'refined_bits':625,'refined_free_bits':623},'scope':'Joint edge/rooted-order relation on the exact family; Hamilton-polynomial factorization and optimization are graph-complete after primitive refinement.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
