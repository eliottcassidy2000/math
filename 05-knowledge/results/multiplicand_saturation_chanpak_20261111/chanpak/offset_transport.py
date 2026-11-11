"""Focused physical offset-advance probes, retaining the old selected set."""
from pathlib import Path
import json
from realize35 import OFFSET_SEEDS
from realize35_s2 import seed_metrics
from head_lift_probe import inflate

def born(B,mask):
 n=len(B);R=list(B)+[mask]
 for i in range(n):
  if not mask>>i&1:R[i]|=1<<n
 return R

def main():
 out=[]
 for s,(B,O,shift,stats) in OFFSET_SEEDS.items():
  D=stats[0]-stats[2];U=stats[1];hits=[]
  for mask in range(1<<len(B)):
   R=born(B,mask)
   for selected_new in (0,1):
    Q=tuple(O)+((len(B),) if selected_new else())
    h,u,w=seed_metrics(R,Q)
    if (h-w,u)==(2*D-1,2*U):hits.append({'outmask':mask,'selected_new':selected_new,'H_U_W':[h,u,w],'rows':R,'selected':Q})
  entry={'s':s,'D_U':[D,U],'target':[2*D-1,2*U],'hits':hits};out.append(entry);print(json.dumps(entry),flush=True)
 Path(__file__).with_name('offset_transport_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
