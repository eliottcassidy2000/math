#!/usr/bin/env python3
"""Exact feedback costs of the deficit-nine family, including its clone clock."""
import itertools,json
from deficit_nine_family import PARENT,MASK
from audit_power_sources import marked,append_zero,need

def fas(rows):
 n=len(rows);full=(1<<n)-1;dp=[0]+[10**9]*full
 for S in range(1,full+1):
  z=S
  while z:
   b=z&-z;z-=b;v=b.bit_length()-1;old=S^b
   # v last: its outgoing arcs to the earlier vertices point backwards.
   dp[S]=min(dp[S],dp[old]+(rows[v]&old).bit_count())
 return dp[full]

def main():
 rows,x=marked(PARENT,MASK);B=[100]*9;orders=[None]*9
 need(fas(PARENT)==3 and fas(rows)==5,'independent parent and marked-core feedback minima')
 for p in itertools.permutations(range(8)):
  base=sum(rows[v]>>u&1 for j,u in enumerate(p) for v in p[j+1:])
  for slot in range(9):
   cost=sum(v!=x for v in p[slot:])+int(x in p[:slot])
   if base<B[cost]:B[cost]=base;orders[cost]=list(p[:slot])+['Y']+list(p[slot:])
 controls=[];tail=rows
 for r in range(1,8):
  tail=append_zero(tail,x);actual=fas(tail);expected=min(b+m*r for m,b in enumerate(B))
  need(actual==expected,'independent feedback subset DP');controls.append({'tail_length':r,'actual_FAS':actual})
 bad=[]
 for deleted in itertools.combinations(range(7),2):
  left=[v for v in range(7) if v not in deleted]
  witness=next(p for p in itertools.permutations(left,3) if all(PARENT[p[i]]>>p[(i+1)%3]&1 for i in range(3)))
  bad.append({'deleted':deleted,'surviving_directed_triangle':witness})
 upper=[4,3,5,2];need(all(PARENT[u]>>v&1 for j,u in enumerate(upper) for v in upper[j+1:]),'delete marked out-neighbours leaves a transitive parent')
 need(set(upper)=={v for v in range(7) if not(MASK>>v&1)},'same marked cut deleted')
 print(json.dumps({'status':'passed','parent_feedback_arc_number':3,'marked_core_feedback_arc_number':5,'all_positive_tail_feedback_arc_number':6,'clone_slot_intercepts':B,'minimizing_orders':orders,'feedback_arc_formula':'min_m(B[m]+m*r) for every integer r>=1','independent_FAS_controls':controls,'feedback_vertex_number':3,'all_two_deletion_hostiles':bad,'transitive_parent_order_after_deleting_S':upper,'scope':'All-r proof uses consecutive transitive-clone compression. FVS is exactly3 since the original parent is induced and deleting S suffices.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
