"""MITM enumeration and shortest buffered routes for the inherited18/37 circuit."""
import argparse,json
from pathlib import Path
from collections import defaultdict,deque
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,default=Path('session35-publication/clean-repo'));args=ap.parse_args();repo=args.repo
r=sorted(json.loads((repo/'05-knowledge/results/native_packet_decomposition_20261023/native/domains.json').read_text())['sources'],key=lambda t:t['source'])
source=[t['source']for t in r];x=[t['X']for t in r];y=[sum(v<<j for j,v in enumerate(t['code'][1:]))for t in r];w=[sum(a!=b for a,b in zip(t['old_before'],t['old_after']))for t in r]
A={1394,1425,1457,1480,2176,2185,2581,2835,2945};B={1477,1549,1691,1921,1974,2837,2905,3000,3369}
mask=lambda ss:sum(1<<source.index(s)for s in ss)
a,b=mask(B),mask(A);core=a|b;outside=((1<<28)-1)^core
half=[]
for shift in (0,14):
 z=[(0,0,0,0)]
 for m in range(1,1<<14):
  j=(m&-m).bit_length()-1;v=z[m&(m-1)];z.append((v[0]+1,v[1]+x[shift+j],v[2]^y[shift+j],v[3]+w[shift+j]))
 half.append(z)
def info(m):
 p,q=half[0][m&16383],half[1][m>>14]
 return (p[0]+q[0],p[1]+q[1],p[2]^q[2],p[3]+q[3])
k,base,yy,_=info(a)
if info(b)[:3]!=(k,base,yy):raise RuntimeError('endpoint full charge and integer source sum')
lo=defaultdict(list)
for m,(kk,xx,z,ww)in enumerate(half[0]):lo[kk,xx%128,z].append(m)
ss=[]
for m,(kk,xx,z,ww)in enumerate(half[1]):
 for n in lo[k-kk,(base-xx)%128,yy^z]:ss.append(n|(m<<14))
ss.sort()
def route(states,budget):
 left=set(states);left.remove(a);todo=deque([a]);prev={a:None}
 while todo and b not in prev:
  m=todo.popleft();found=[n for n in left if(m^n).bit_count()<=8 and info(m^n)[3]<=budget]
  for n in found:left.remove(n);todo.append(n);prev[n]=m
 if b not in prev:return None
 out=[];m=b
 while m is not None:out.append(m);m=prev[m]
 return out[::-1]
def record(path):
 out=[]
 for m in path:
  out.append({'mask':m,'V_sources':[s for i,s in enumerate(source)if m>>i&1], 'integer_X_sum':info(m)[1], 'U_integer_current_change':2*(info(m)[1]-base), 'active_outside_sources':[s for i,s in enumerate(source)if(m&outside)>>i&1]})
 ev=[]
 for m,n in zip(path,path[1:]):
  dd=m^n
  ev.append({'changed_sources':[s for i,s in enumerate(source)if dd>>i&1],'changed_slots':dd.bit_count(),'ordinary_owner_clock_work':info(dd)[3]})
 return {'states':out,'events':ev,'peak_clock':max(z['ordinary_owner_clock_work']for z in ev),'total_clock_work':sum(z['ordinary_owner_clock_work']for z in ev)}
report={'targets':[1965,2234],'endpoint_margin':k,'endpoint_X_sum':base,'charge':[base%128,yy],'full28_states':len(ss),'outside_sources':[s for i,s in enumerate(source)if outside>>i&1]}
for exact in (False,True):
 domain=[m for m in ss if not exact or info(m)[1]==base]
 answer=None
 tested=[]
 for budget in range(10,23):
  pp=route(domain,budget);tested.append([budget,pp is not None])
  if pp is not None:answer=record(pp);break
 report['exact_X'if exact else'dyadic']={'domain_states':len(domain),'tested_clock_thresholds':tested,'route':answer}
print(json.dumps(report,indent=2));(P/'route18_37.json').write_text(json.dumps(report,indent=2)+'\n')
