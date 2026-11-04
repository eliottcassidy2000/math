"""Explicit ordered-witness and exterior-return checks for the new physical contexts."""
from pathlib import Path
import json,hashlib
from physical_audit import exterior,mul,induced,module,profile
D=Path(__file__).resolve().parent

def need(c,s):
 if not c:raise RuntimeError(s)
def main():
 data=json.loads((D/'all_cover_grouped_controls.json').read_text());out=[]
 for r in data:
  T=r['source_rows'];A,X,Y=r['child_sets'];outside=list(range(len(A),len(T)));R=induced(T,outside);Xr=[outside.index(v) for v in X];Yr=[outside.index(v) for v in Y]
  need(module(R,Xr+Yr),'ordered witness union is a module after distinguished child deletion')
  need(all(T[a]>>x&1 and T[x]>>y&1 and T[y]>>a&1 for a in A for x in X for y in Y),'ordered triangle incidences with complete physical witness sets')
  into=sum(1<<i for i,v in enumerate(outside) if T[v]>>A[0]&1);alpha,beta,tau,*_=exterior(R,into);power=[1]+[0]*((1<<len(R))-1);theta=[]
  for j in range(len(R)+1):theta.append(mul(mul(alpha,beta),power)[-1]);power=mul(power,tau)
  need(theta[0]>0,'positive contracted Hamilton weight')
  expected=sum(1<<v for v in r['exterior_names'] if T[v]>>r['root_names'][0]&1);need(expected==r['base_physical_mask'],'both fixed port bits imitate the source core on each exterior name')
  if r['ambient']:need(theta[3]>0,'genuine third excursion survives in the test exterior')
  out.append(dict(m=r['m'],ambient=r['ambient'],exterior_physical_names=outside,exterior_rows=R,into_mask=into,complete_theta=theta,fixed_ambient_mask=expected,grouped_scalar_difference=r['grouped'][0]['complete_combined_contrast'][0]-r['grouped'][1]['complete_combined_contrast'][0]))
 p=D/'scope_controls.json';txt=json.dumps(out,separators=(',',':'))+'\n'
 if p.exists():need(p.read_text()==txt,'retained scope controls stable')
 else:p.write_text(txt)
 return dict(status='PASS',source_contexts=len(out),genuine_third_excursion_contexts=sum(z['complete_theta'][3]>0 for z in out),source_sha256=hashlib.sha256((D/'scope_audit.py').read_bytes()).hexdigest(),data_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
