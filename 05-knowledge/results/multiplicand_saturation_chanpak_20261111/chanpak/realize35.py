"""Complete literal constructor for every N(35,r,s), r,s positive."""
from pathlib import Path
import argparse,json
from head_lift_probe import inflate,ends,B
from realize35_s2 import tail,seed_metrics,realize as realize_s2
HERE=Path(__file__).resolve().parent
OFFSET_SEEDS={
 3:((26,44,57,48,34,1),(1,2),1,(31,140,22)),
 4:((62,124,120,112,32,64,17),(4,6),3,(27,70,10)),
 5:((30,124,24,112,96,69,5),(2,3),2,(95,280,62)),
}
OLD_SEED=(4,1,18,7,75,31,47)

def realize(r,s):
 if min(r,s)<1:raise ValueError('r,s must both be positive')
 if s==1:return inflate(B,4,r)
 if s==2:return realize_s2(r)
 if s in OFFSET_SEEDS:
  C,O,shift,_=OFFSET_SEEDS[s]
  return tail(C,O,r+shift)
 return tail(inflate(OLD_SEED,5,r+5),(0,3),s-5)

def strong(T):
 reach=[row|(1<<i) for i,row in enumerate(T)]
 for k in range(len(T)):
  for i in range(len(T)):
   if reach[i]>>k&1:reach[i]|=reach[k]
 return all(row==(1<<len(T))-1 for row in reach)

def check():
 seeds=[]
 for s,(C,O,shift,expected) in OFFSET_SEEDS.items():
  actual=seed_metrics(C,O);assert actual==expected
  assert actual[0]-actual[2]==1+2**s
  assert actual[1]*2**shift==35*2**s
  assert strong(C)
  seeds.append({'s':s,'rows':C,'selected':O,'H_U_W':actual,'tail_shift':shift})
 controls=[]
 for r in (1,2,3):
  for s in range(1,8):
   T=realize(r,s);H=sum(ends(T)[-1]);expected=1+2**s+35*2**(r+s)
   assert H==expected
   if s!=2:assert strong(T)
   controls.append({'r':r,'s':s,'H':H,'order':len(T),'strong':strong(T),'rows':T})
 receipt={'scope':'Finite literal controls for the complete all-r,s proof in PROOF.md. Old s1 and s>=6 mechanisms inherited; new s2 product and s3..5 offset seeds.',
          'new_offset_seeds':seeds,'literal_controls':controls}
 (HERE/'realize35_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':'pass','new_seed_stats':[s['H_U_W'] for s in seeds],'literal_instances':len(controls),'maximum_order':max(c['order'] for c in controls)}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--r',type=int);p.add_argument('--s',type=int);a=p.parse_args()
 if a.r is None and a.s is None:check()
 elif a.r is None or a.s is None:p.error('Supply both --r and --s')
 else:
  T=realize(a.r,a.s)
  print(json.dumps({'r':a.r,'s':a.s,'H':1+2**a.s+35*2**(a.r+a.s),'order':len(T),'rows':T},indent=2))
