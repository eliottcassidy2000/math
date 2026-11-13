"""Finite arithmetic gate for one constant growing run; unbounded laws are in the note."""
from pathlib import Path
import json
from defect_two_inverse import shell
HERE=Path(__file__).resolve().parent

def Q(x):return 3*x+shell(x)+(x&1)
def reverse(d,r):return int(format(d,'0%db'%r)[::-1],2) if r else 0
def odd(x):return x//(x&-x)

def main():
 base=45157534;controls=[]
 for ell in range(1,13):
  x=base<<ell;h=180630136*2**ell-21377441;q=90315068*2**ell-21377441
  assert Q(x)==h and x+shell(x)+(x&1)==q
  near=320*2**ell+44
  assert Q(near)==1215*2**ell+167
  assert near+shell(near)+(near&1)==575*2**ell+79
  controls.append(dict(ell=ell,zero_prefix_H=h,nearby_H=Q(near)))
 choices=[]
 for bit in (0,1):
  for r in range(9):
   for d in range(1<<r):
    if r and ((d>>(r-1))&1)==bit:continue
    L=4*(1<<r)-reverse(d,r)-bit
    if 1215%odd(L):continue
    z=1215//odd(L);n=z-bit
    B=Q(d)-(n+1-(n&1)-shell(n)) if bit==0 else shell(n)+Q(d)-L
    row=dict(bit=bit,low_length=r,low_value=d,high_value=n,coefficient=L*(n+bit),offset=B)
    for ell in (1,2,5):
     y=(n<<ell) if bit==0 else ((n+1)<<ell)-1
     assert Q((y<<r)+d)==row['coefficient']*2**ell+B
    choices.append(row)
 assert len(choices)==10 and all(z['offset']!=215 for z in choices)
 out=dict(status='PASS',base=dict(x=base,g=shell(base),H=Q(base),q=base+shell(base)),
  scope='All normalized fixed singleton outer words surrounding one growing constant-bit run with coefficient oddpart1215; normalization proves r<=8',
  constant_run_gate=choices,arithmetic_controls=controls)
 (HERE/'singleton_continuation_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','base':out['base'],'offsets_by_bit':{b:[z['offset'] for z in choices if z['bit']==b] for b in (0,1)}}))

if __name__=='__main__':main()
