#!/usr/bin/env python3
"""Exact partition and linear-observer audit; no dependencies, active under -O."""
from collections import Counter
from fractions import Fraction
import json
CHECKS=0
def need(b,msg='check failed'):
 global CHECKS
 CHECKS+=1
 if not b: raise RuntimeError(msg)
def partition(labels):
 d={};return tuple(d.setdefault(x,len(d)) for x in labels)
def refine(edges,labels):
 p=partition(labels); history=[len(set(p))]
 while True:
  q=partition([(p[i],tuple(sorted(Counter(p[j] for j in edges[i]).items()))) for i in range(len(p))])
  if q==p:return p,history
  p=q;history.append(len(set(p)))
def rank(rows):
 piv={}
 for row in rows:
  v=list(map(Fraction,row))
  for k,b in sorted(piv.items()):
   if v[k]:
    a=v[k];v=[x-a*y for x,y in zip(v,b)]
  if any(v):
   k=next(i for i,x in enumerate(v) if x);a=v[k];piv[k]=[x/a for x in v]
 return len(piv)
def wreath():
 arc=lambda x,y: ((y//3-x//3)%3==1) if x//3!=y//3 else (y-x)%3==1
 tau=lambda v:3*((-(v//3))%3)+(-(v%3))%3
 pairs=[(v,tau(v)) for v in range(9) if v<tau(v)]
 ix={v:i for i,pr in enumerate(pairs) for v in pr}
 states=[(m,v) for m in range(1,16) for v in range(1,9) if m>>ix[v]&1]
 idx={s:i for i,s in enumerate(states)}
 edges=[[idx[(m|1<<ix[w],w)] for w in range(1,9) if not m>>ix[w]&1 and arc(v,w)] for m,v in states]
 base=[(m.bit_count(),int(m==15 and arc(v,0))) for m,v in states]
 return states,edges,base

def main():
 states,edges,base=wreath();need(len(states)==64);need(sum(map(len,edges))==96)
 p,h=refine(edges,base);need(len(set(p))==16,(len(set(p)),h))
 q,k=refine(edges,[(b,v==3) for b,(m,v) in zip(base,states)]);need(len(set(q))==35,k)
 # Exact aggregate propagation, with arbitrary initial masses supported at level1.
 def count(part):
  blocks=[[i for i,x in enumerate(part) if x==a] for a in range(max(part)+1)]
  mass=[sum(states[i][0].bit_count()==1 for i in block) for block in blocks]
  for a in sorted(range(len(blocks)),key=lambda a:base[blocks[a][0]][0]):
   rows=[Counter(part[j] for j in edges[i]) for i in blocks[a]]
   need(all(r==rows[0] for r in rows))
   for b,w in rows[0].items():mass[b]+=mass[a]*w
  return sum(mass[a] for a,block in enumerate(blocks) if base[block[0]][1])
 need(count(p)==33);need(count(q)==33)
 # Linear pullback closure of base colour indicators: DAG means A^4=0.
 cols=[]
 for color in sorted(set(base)):
  f=[int(b==color) for b in base]
  for _ in range(5):
   cols.append(f);f=[sum(f[j] for j in row) for row in edges]
 linear_rank=rank(cols)
 # Terminal-only future-count equivalence can differ from stable partition.
 f=[b[1] for b in base]; seq=[]
 for _ in range(5):seq.append(f);f=[sum(f[j] for j in row) for row in edges]
 count_classes=len(set(zip(*seq)))
 marks={str(v):len(set(refine(edges,[(b,w==v) for b,(m,w) in zip(base,states)])[0])) for v in range(1,9)}
 # Joint-coordinate hostile: target grid states are absorbing (no outgoing).
 # x sees 00+11; y sees 01+10. Single coordinates cannot distinguish x/y.
 e=[[2,5],[3,4],[],[],[],[]]
 first=['s','s',0,0,1,1];second=['s','s',0,1,0,1]
 a,_=refine(e,first);b,_=refine(e,second)
 need(a[0]==a[1] and b[0]==b[1]);meet=partition(list(zip(a,b)))
 joint,jh=refine(e,list(zip(a,b)));need(meet[0]==meet[1] and joint[0]!=joint[1])
 # Repeated backward observation provides explicit depth certificates.
 depth=[]
 for t in range(5):
  labels=[(base[i],states[i][1]==3) for i in range(64)]
  r=partition(labels)
  for _ in range(t):r=partition([(r[i],tuple(sorted(Counter(r[j] for j in edges[i]).items()))) for i in range(64)])
  depth.append(len(set(r)))
 # Quadratic closure and full named-observer menu census.
 need(linear_rank==12)
 g=[sum(base[j][0]==3 for j in row) for row in edges]
 need(rank(cols+[[v*v for v in g]])==13)
 need(len(set(zip(*cols)))==16)
 products=[[x*y for x,y in zip(a,c)] for i,a in enumerate(cols) for c in cols[:i+1]]
 quadratic_rank=rank(cols+products);need(quadratic_rank==16)
 menus={}
 for mask in range(256):
  labels=[(base[i],tuple(v==w for w in range(1,9) if mask>>(w-1)&1)) for i,(m,v) in enumerate(states)]
  menus[mask]=len(set(refine(edges,labels)[0]))
 ranges={str(k):[min(c for m,c in menus.items() if m.bit_count()==k),max(c for m,c in menus.items() if m.bit_count()==k)] for k in range(9)}
 minimal_full=[m for m,c in menus.items() if c==64 and all(menus[m^(1<<i)]<64 for i in range(8) if m>>i&1)]
 need(minimal_full==[126,190,222,237,245,249])
 need(menus[1]+menus[4]<menus[0]+menus[5])
 need(menus[1]+menus[2]>menus[0]+menus[3])
 need(menus[3]==22 and menus[5]==43)
 namedlabels=[(base[i],v==3) for i,(m,v) in enumerate(states)]
 namedcols=[]
 for color in sorted(set(namedlabels)):
  f=[int(b==color) for b in namedlabels]
  for _ in range(5):
   namedcols.append(f);f=[sum(f[j] for j in row) for row in edges]
 need(rank(namedcols)==21);need(len(set(zip(*namedcols)))==34)
 x,y=[states.index(st) for st in [(8,7),(4,8)]]
 need(all(f[x]==f[y] for f in namedcols));need(q[x]!=q[y])
 f=[int(z==(4,0) and v!=3) for z,(m,v) in zip(base,states)]
 for _ in range(2):f=[sum(f[j] for j in row) for row in edges]
 successor_values=[[f[j] for j in edges[i]] for i in [x,y]]
 need(successor_values==[[1,1,1],[0,0,3]])
 need([sum(v*v for v in row) for row in successor_values]==[3,9])
 # Centralizer C3 action, preserving root and anti-involution.
 tau=lambda v:3*((-(v//3))%3)+(-(v%3))%3
 pairs=[(v,tau(v)) for v in range(9) if v<tau(v)]
 pairindex={v:i for i,pr in enumerate(pairs) for v in pr}
 def action(t,state):
  m,v=state
  g=lambda w:3*(w//3)+(w%3+t*(0 if w//3==0 else 1 if w//3==1 else -1))%3
  return sum(1<<pairindex[g(pr[0])] for i,pr in enumerate(pairs) if m>>i&1),g(v)
 orbit=partition([min(action(t,st) for t in range(3)) for st in states])
 need(len(set(orbit))==24)
 need(all(p[i]==p[j] for i in range(64) for j in range(64) if orbit[i]==orbit[j]))
 need(count(orbit)==33)
 result={'states':64,'arcs':96,'centralizer_orbits':24,'quadratic_observer_rank':quadratic_rank,'menu_size_ranges':ranges,'minimal_full_menus':minimal_full,'named3_linear_rank':21,'named3_linear_signature_classes':34,'named3_nonlinear_witness_states':[[8,7],[4,8]],'named3_nonlinear_successor_values':successor_values,'base_partition_history':h,'named3_partition_history':k,'folded_paths':33,'linear_base_observer_rank':linear_rank,'terminal_count_sequence_classes':count_classes,'single_named_vertex_partition_sizes':marks,'named3_depth_sizes':depth,'joint_hostile_separate_classes':[len(set(a)),len(set(b))],'joint_hostile_raw_meet_classes':len(set(meet)),'joint_hostile_stable_classes':len(set(joint)),'checks':CHECKS}
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
