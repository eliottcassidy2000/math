"""Finite-exact order-five named quartic classification, with full vectors."""
import itertools,json
from collections import defaultdict,Counter
from pathlib import Path
from verify import tournament,quartic_normal_form,hp,extend,polynomial,strong,need

def special(a):
 n=len(a);sources=[u for u in range(n) if a[u].bit_count()==n-1];sinks=[u for u in range(n) if a[u]==0]
 if len(sources)!=1 or len(sinks)!=1:return None
 rest=[u for u in range(n) if u not in sources+sinks]
 if not all(sum(a[u]>>v&1 for v in rest)==1 for u in rest):return None
 return {'source':sources[0],'sink':sinks[0],'triangle':rest}

def main():
 groups=defaultdict(list);basis=None
 for code in range(1024):
  a=tournament(5,code);Q=quartic_normal_form(a);keys=sorted(Q)
  if basis is None:basis=keys
  need(keys==basis,'same complete named basis');groups[tuple(Q[s] for s in basis)].append(code)
 need(len(groups)==502,'quartic class count');hist=Counter(map(len,groups.values()));need(hist=={2:492,4:10},'class-size histogram')
 records=[];exceptional_triangles=set();strong_count=0
 for vector,codes in sorted(groups.items(),key=lambda t:t[1][0]):
  rows=[tournament(5,c) for c in codes];detail=[]
  if len(codes)==2:
   need(codes[0]^codes[1]==1023,'ordinary class is whole-converse pair')
   need(all(special(a) is None for a in rows),'special class unexpectedly size two')
  else:
   forms=[special(a) for a in rows];need(all(x is not None for x in forms),'exception has required ordered shape')
   triangles={tuple(f['triangle']) for f in forms};need(len(triangles)==1,'same named triangle');exceptional_triangles.update(triangles)
   need(len({f['source'] for f in forms})==2,'both exterior orders')
   detail=forms
  for a in rows:
   if strong(a):need(len(codes)==2,'strong quartic reconstruction');strong_count+=1
  records.append({'codes':codes,'rows':rows,'quartic_vector':list(vector),'exceptional_ordered_shapes':detail})
 need(exceptional_triangles==set(itertools.combinations(range(5),3)),'all ten named triangles')
 # An exact quartic-only collision, explicitly separated by the actual full observer.
 a=tournament(5,5);b=tournament(5,131);p=polynomial(a);q=polynomial(b);d={s:p.get(s,0)-q.get(s,0) for s in set(p)|set(q) if p.get(s,0)!=q.get(s,0)}
 need(all(s.bit_count()<=3 for s in d),'same quartic implies cubic difference');need(max(s.bit_count() for s in d)==3,'exception has a genuine cubic')
 fa=hp(extend(a,1,4));fb=hp(extend(b,1,4));need((fa,fb)==(31,29),'literal full-response separator')
 return {'status':'PASS','scope':'Exact finite classification of the full named quartic tensor at old order five only; not full-response equality.','tournaments':1024,'classes':502,'class_sizes':dict(sorted(hist.items())),'strong_tournaments':strong_count,'quartic_basis':[[s&31,s>>5] for s in basis],'basis_convention':'All 100 bidegree-(2,2) coordinates are retained in this fixed physical-name order, including zeros.','records':records,'quartic_only_collision':{'T':a,'U':b,'a_mask':1,'c_mask':4,'F_T':fa,'F_U':fb,'complete_difference':[[s&31,s>>5,v] for s,v in sorted(d.items())],'omitted_coefficients':'Zero in the full Boolean basis.'}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
