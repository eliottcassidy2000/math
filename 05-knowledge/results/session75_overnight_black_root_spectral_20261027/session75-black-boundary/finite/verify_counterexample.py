#!/usr/bin/env python3
"""Independent exact n10 gate, including every Hamilton word and literal K edges."""
from independent_replay import arcs,connected,core,blue,reverse,opposite,replay,require,cut_parameters
import json,itertools
ROWS=(1008,1001,995,997,494,448,640,768,576,48)
ROUTE=((0,5,6,7,8,9,4,3,2,1),(2,3,7,9,8,6,5,0,4,1))
def every_word(e,n):
 def go(p,unused):
  if not unused:yield tuple(p);return
  for v in sorted(unused):
   if not p or(p[-1],v)in e:yield from go(p+[v],unused-{v})
 yield from go([],set(range(n)))
def main():
 n=10;e=arcs(ROWS);require(connected(e,range(n)),'not strong');delta,lam,ss=cut_parameters(e,n);require((delta,lam,ss)==(2,1,[992]),'wrong exact invariants')
 direct=[]
 for mask in range(1,1023):
  count=sum(bool(mask>>u&1)and not(mask>>v&1)for u,v in e)
  if count==lam:direct.append(mask)
 require(direct==ss,'score formula/direct cuts disagree')
 allp=list(every_word(e,n));allowed=[p for p in allp if 992>>p[0]&1 and not(992>>p[-1]&1)];require(len(allowed)==9,'wrong cone size');require(all(blue(e,p)for p in allowed),'black at lower bound')
 # Stronger direct certificate: exhaust EVERY black neighbor, not just minimum-shore cone.
 blackp=[p for p in allp if not blue(e,p)];require(blackp,'not black-active');require(all(connected(reverse(e,p),range(n))for p in blackp),'black nonstrong neighbor')
 replay(e,n,{'route':ROUTE,'dual':False,'port':1},'B',2)
 moves=[];a=e
 for j,p in enumerate(ROUTE):
  r={p[i]:p[-1-i]for i in range(n)};bad=next((u,v,r[v],r[u])for u,v in sorted(a)if(r[v],r[u])not in a);b=reverse(a,p);s=a if j%2==0 else opposite(a);q=p if j%2==0 else p[::-1];k=opposite(b)if j%2==0 else b
  row=lambda E:tuple(sum(1<<v for u,v in E if u==w)for w in range(n))
  moves.append({'R_path':p,'R_source_rows':row(a),'R_target_rows':row(b),'black_failed_anti_arcs':bad,'K_path':q,'K_source_rows':row(s),'K_target_rows':row(k),'R_target_strong':connected(b,range(n))});a=b
 result={'status':'PASS','rows':ROWS,'scores':[sum(u==v for u,w in e)for v in range(n)],'lambda':lam,'delta':delta,'unique_minimum_shore':ss[0],'all_Hamilton_words':len(allp),'all_black_Hamilton_words':len(blackp),'minimum_cone_words':allowed,'all_minimum_cone_words_blue':True,'every_black_neighbor_strong':True,'exact_black_nonstrong_distance':2,'exact_black_core_distance':2,'terminal_core_strong':core(a,n),'moves':moves}
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
