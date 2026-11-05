#!/usr/bin/env python3
import argparse,hashlib,itertools,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=HERE.parent/'session195-eliminate-allorder');p.add_argument('--extension',type=Path,default=HERE.parent/'session195-eliminate-grammar-exceptions');args=p.parse_args();SRC=args.source;EXT=args.extension

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def S(x):
 h,e,r,d=x;return h+r-e,h-e,h+r,h

def K(x):
 h,e,r,d=x;return h+r-e,d,2*r-d,d

def D(x):
 h,e,r,d=x;return h,0,2*r,d

def M(x,y):
 h,e,r,d=x;g,s,L,E,Q,B,z,dg=y
 return g*h-s*e+L*r,E*h-z*e,Q*h+B*r,dg*h

def run(seed,word):
 for op in word:seed=(K if op=='K' else D)(seed)
 return seed

def inverse(seed,target):
 r=target[2];base=seed[2];d=seed[3];rev=[]
 if target[3]!=d or r<base:return None
 while r>base:
  op='K' if r&1 else 'D';prev=(r+d)//2 if op=='K' else r//2
  if prev<base or prev>=r:return None
  rev.append(op);r=prev
 word=''.join(reversed(rev))
 return word if run(seed,word)==target else None

def main():
 hashes={name:hashlib.sha256((SRC/name).read_bytes()).hexdigest() for name in ['PROOF.md','verify.py']}
 bankhashes={name:hashlib.sha256((EXT/name).read_bytes()).hexdigest() for name in ['base_inventory.json','gadgets.json','inventory6.json']}
 base=json.loads((EXT/'base_inventory.json').read_text());new=json.loads((EXT/'gadgets.json').read_text());six=json.loads((EXT/'inventory6.json').read_text())
 records=base+new
 keys=[tuple(g['signature'])+(g['H_deleted_root'],) for g in records]
 for y,g in zip(keys,records):
  h,s,L,E,Q,B,z,dg=y
  need(h>=s>=0 and L>=h and Q>=h and 1<=dg<=h and 0<=E<=dg,'physical gadget bounds')
  need(dg&1 and (z==1 or B>=2*dg>=2),'deletion oddness/empty covers')
 seed=(1,1,1,1);ops=[('K',K,1),('D',D,1)]+[(f'M{i}',lambda x,y=y:M(x,y),len(g['rows'])) for i,(y,g) in enumerate(zip(keys,records))]
 nodes=pairchecks=overshoots=inversechecks=contractions=0
 def walk(x,depth,maxdepth,Rhist,lastnond,vertices,corevertices,alphabet):
  nonlocal nodes,pairchecks,overshoots,inversechecks,contractions
  if depth==maxdepth:return
  for name,op,cost in alphabet:
   if depth==0 and name=='K':continue
   y=op(x);j=depth+1;h,e,r,d=y;nodes+=1
   need(0<=e<=d<=h<=r and d>=1,'state cone')
   need(r>=2*d and r>=1<<((j+1)//2),'invariant/exponential bound')
   need(h>=x[0],'H monotonicity')
   if j>=2:need(r>=2*Rhist[-2],'two operation doubling');pairchecks+=1
   ell=lastnond if name=='D' else j
   v=vertices+cost;cv=corevertices if name=='D' else v
   need(ell<=2*(h.bit_length()-1)+1,'logarithmic reduced length')
   need(cv<=1+6*(2*(h.bit_length()-1)+1),'logarithmic order')
   if name!='D':
    need(h>=x[2],'preterminal R control');overshoots+=r>h
    if depth>0:need(2*h>=2*x[0]+x[2]>=3*x[0],'strict non-D contraction');contractions+=1
   if name.startswith('M'):
    key=keys[int(name[1:])];g,s,L,E,Q,B,z,dg=key
    if z:
     need((d,d-e,r-d)==(x[0],x[1],x[2]),'singleton inverse')
     need(h==e+r-d,'singleton inverse identity')
    else:
     need(d//dg==x[0] and d%dg==0,'gadget inverse h')
     need((r-Q*x[0])//B==x[2] and (r-Q*x[0])%B==0,'gadget inverse R')
     if s:need((g*x[0]+L*x[2]-h)//s==x[1] and (g*x[0]+L*x[2]-h)%s==0,'gadget inverse endpoint')
    inversechecks+=1
   walk(y,j,maxdepth,Rhist+[r],ell,v,cv,alphabet)
 walk(seed,0,3,[1],0,1,1,ops)
 # Long words include all slow singleton/stationary cases, both new pieces,
 # and a non-singleton edge. Repeated D tails deliberately stress normalization.
 chosen=[ops[0],ops[1],ops[2],ops[3],ops[-2],ops[-1]]
 walk(seed,0,7,[1],0,1,1,chosen)
 need(K(seed)==seed and S(seed)==D(seed)==(1,0,2,1),'singleton controls')
 for m in range(1,30):
  x=seed
  for _ in range(m):x=K(x)
  need(x==seed,'initial K loop')
  y=S(seed)
  for _ in range(m):y=D(y)
  need(y[0]==1 and y[2]>y[0],'terminal overshoot control')
 seeds=sorted({(g['signature'][0],g['signature'][3],g['signature'][4],g['H_deleted_root']) for g in base+six})
 stationary=zero_growth=wrong=0
 for x in seeds:
  need(x[3]&1,'seed oddness')
  if x[2]==x[3]:need(x==(x[3],)*4 and K(x)==x,'all physical zero growth');zero_growth+=1
  for k in range(11):
   for wordtuple in itertools.product('KD',repeat=k):
    word=''.join(wordtuple);target=run(x,word);canonical=word.lstrip('K') if x[2]==x[3] else word
    got=inverse(x,target);need(got==canonical,('stationary inverse',x,word,target,got))
    need(len(got)<=max(0,(target[2]-x[3]).bit_length()+1),'stationary logarithmic length')
    bad=(target[0]+2,target[1],target[2],target[3]);need(inverse(x,bad) is None,'full coordinate inverse check');wrong+=1;stationary+=1
 # Independent pair expansion for the singleton binary response.
 binarychecks=0
 for q in range(1<<15):
  w=bin(q)[2:] if q else ''
  pair=1+sum(w[i:i+2]=='10' for i in range(len(w)-1))
  pair+=sum(1<<(j-i-1) for i in range(len(w)) for j in range(i+1,len(w)) if w[i]=='1' and w[j]=='0')
  x=run(seed,''.join('D' if c=='1' else 'K' for c in w))
  need(x==(pair,1-(q&1),q+1,1),'binary pair expansion')
  if q and q%2==0:need(q<=2*x[0]-2,'stationary target bound')
  binarychecks+=1
 need(all(hashlib.sha256((SRC/n).read_bytes()).hexdigest()==v for n,v in hashes.items()),'source changed during audit')
 need(all(hashlib.sha256((EXT/n).read_bytes()).hexdigest()==v for n,v in bankhashes.items()),'inventory changed during audit')
 out={'status':'PASS','source_sha256':hashes,'bank_sha256':bankhashes,'non_D_contraction_checks':contractions,'bank_matrices':len(keys),'reduced_word_steps_checked':nodes,'pair_checks':pairchecks,'terminal_R_exceeds_H_controls':overshoots,'root_moving_inverse_checks':inversechecks,'physical_stationary_seeds':len(seeds),'physical_R_equals_d_seeds':zero_growth,'stationary_word_inversions':stationary,'wrong_coordinate_targets_rejected':wrong,'independent_binary_pair_checks':binarychecks,'scope':'All-order claims proved algebraically; no unbounded coverage assertion.'}
 (HERE/'receipt.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
