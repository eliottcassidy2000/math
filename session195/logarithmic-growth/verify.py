#!/usr/bin/env python3
import argparse,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def S(x):
 h,e,R,d=x;return h-e+R,h-e,h+R,h
def K(x):
 h,e,R,d=x;return h-e+R,d,2*R-d,d
def D(x):
 h,e,R,d=x;return h,0,2*R,d
def M(x,p):
 h,e,R,d=x;g,s,L,E,Q,B,delta,dg=p
 return g*h-s*e+L*R,E*h-delta*e,Q*h+B*R,dg*h

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--extension',type=Path,default=HERE.parent/'gadget-extension');args=ap.parse_args()
 bank=json.loads((args.extension/'base_inventory.json').read_text())+json.loads((args.extension/'gadgets.json').read_text())
 profiles=[tuple(g['signature']+[g['H_deleted_root']]) for g in bank]
 for g,p in zip(bank,profiles):
  h,s,L,E,Q,B,delta,d=p
  need(h%2==d%2==1 and L>=h>=s>=0 and Q>=h>=E>=0,'gadget bounds')
  if len(g['rows'])>1:need(B>=2*d>=2,'non-singleton empty-slot bound')
 ops=[K,D]+[(lambda x,p=p:M(x,p)) for p in profiles]
 # Check inequalities on a larger abstract cone, not just physically attained inputs.
 one=two=contractions=0
 for h in range(1,7):
  for d in range(1,h+1):
   for R in range(max(h,2*d),15):
    for e in range(h+1):
     x=h,e,R,d
     for oi,op in enumerate(ops):
      y=op(x);one+=1
      if oi!=1 and e<=d:
       need(2*y[0]>=2*h+R>=3*h,'non-D count contraction');contractions+=1
      need(y[2]>=2*y[3] and y[0]>=h,'invariant or H monotonicity')
      for op2 in ops:
       z=op2(y);two+=1;need(z[2]>=2*R,'two-step doubling')
 seed=(1,1,1,1)
 need(K(seed)==seed,'initial fixed loop')
 for op in ops:
  y=op(seed)
  if y==seed:continue
  need(y[2]>=2 and y[2]>=2*y[3],'first-step bound')
  for op2 in ops:need(op2(y)[2]>=2,'initial pair bound')
 # Exact stationary language, canonical inverse word, and all binary inputs through 2^15.
 f=[1]
 for q in range(1,1<<15):
  n=q//2;f.append(f[n] if q&1 else f[n]+n+(n&1))
 for q,H in enumerate(f):
  x=seed
  if q:
   for bit in bin(q)[2:]:x=D(x) if bit=='1' else K(x)
  need(x==(H,1-(q&1),q+1,1),'canonical stationary formula')
  R=q+1;rev=[]
  while R>1:
   if R%2==0:rev.append('D');R//=2
   else:rev.append('K');R=(R+1)//2
  letters=''.join(reversed(rev))
  need(letters==(''.join('D' if c=='1' else 'K' for c in bin(q)[2:]) if q else ''),'canonical inverse')
 for H in [7,13]:need(not any(f[q]==H for q in range(0,2*H,2)),'stationary omission')
 receipt={'status':'PASS','bank_matrices':len(profiles),'abstract_one_step_checks':one,'abstract_two_step_checks':two,'non_D_contraction_checks':contractions,'stationary_binary_inputs':len(f),'logarithmic_length_bound':'2 floor(log2 H)+1','scope':'All-order proof is algebraic; finite checks verify formulas and explicit bank hypotheses.'}
 (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
