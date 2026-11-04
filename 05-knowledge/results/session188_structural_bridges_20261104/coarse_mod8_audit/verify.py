"""Independent run-word certificate for all64 coarse modulo-eight queries."""
from collections import defaultdict
from itertools import product
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'covered_quotient'))
from verify import rows,profile,cycle,need


def words(a,c,op=False):
 Q=[4,1,2] if op else [2,4,1];R=Q+[1<<4,0]
 for u in range(3):
  if a>>u&1:R[u]|=1<<3
  else:R[3]|=1<<u
  if c>>u&1:R[u]|=1<<4
  else:R[4]|=1<<u
 D={(tuple(int(i==u) for i in range(3)),0,u):1 for u in range(3)}
 D[(0,0,0),1,3]=1;D[(0,0,0),2,4]=1;coeff=defaultdict(int)
 for length in range(1,12):
  E=defaultdict(int)
  for (ks,p,u),z in D.items():
   if p==3 and min(ks):coeff[ks]+=z
   vs=R[u]
   while vs:
    bit=vs&-vs;vs-=bit;v=bit.bit_length()-1;cs=list(ks);q=p
    if v<3:
     if cs[v]==3:continue
     cs[v]+=1
    else:
     b=1<<(v-3)
     if q&b:continue
     q|=b
    E[tuple(cs),q,v]+=z
  D=E
 return coeff


def affine_coefficients(W):
 out=[0]*8
 affine={1:(1,2),2:(-2,4),3:(4,0)}
 for ks,z in W.items():
  for mask in range(8):
   t=z
   for i,k in enumerate(ks):t*=affine[k][(mask>>i)&1]
   out[mask]+=t
 return out


def main():
 root=Path(__file__).resolve().parent;producer=root.parent/'coarse_mod8_boundary'
 retained=json.loads((producer/'receipt.json').read_text());lookup={(x['block_a_mask'],x['block_c_mask']):x['multilinear_coefficients_by_subset_mask'] for x in retained['coefficient_records']}
 records=[]
 for a,c in product(range(8),repeat=2):
  P=affine_coefficients(words(a,c));Q=affine_coefficients(words(a,c,True));D=[x-y for x,y in zip(P,Q)]
  need(all(x%8==0 for x in D),'symbolic universal blindness')
  need(D==lookup[a,c],'independent run coefficients differ from certificate')
  records.append({'a_blocks':a,'c_blocks':c,'exact_coefficients':D})
 source_profiles=[]
 for n,code in [(4,c) for c in range(64)]+[(6,c) for c in (0,37,1234,8191,16384,32767)]:
  A=rows(n,code);q=profile(A);H=q[1]
  need(H%2==1 and q[2]%8==(2*H-4)%8 and q[3]%8==4 and all(x%8==0 for x in q[4:]),'even profile reduction')
  source_profiles.append({'order':n,'code':code,'profile':q})
 TT4=[14,12,8,0];S4=[8,1,3,6];TT6=rows(6,32767);A6=rows(6,12345)
 literal=[]
 with tempfile.TemporaryDirectory(prefix='coarse8_audit_') as work:
  exe=Path(work)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(root/'literal.cpp'),'-o',str(exe)],check=True)
  for blocks in ([TT4,TT4,TT4],[TT4,S4,A6],[S4,S4,TT6]):
   T,offs=cycle(blocks);U,_=cycle(blocks,True)
   for a,c in ((1,2),(1,3),(3,5)):
    aa=sum(((1<<len(A))-1)<<off for i,(A,off) in enumerate(zip(blocks,offs)) if a>>i&1)
    cc=sum(((1<<len(A))-1)<<off for i,(A,off) in enumerate(zip(blocks,offs)) if c>>i&1);values=[]
    for V in (T,U):
     data=f'{len(V)} {aa} {cc} '+' '.join(map(str,V))+'\n';values.append(int(subprocess.check_output([str(exe)],input=data.encode())))
    need((values[0]-values[1])%8==0,'literal unequal child blindness')
    literal.append({'children':blocks,'profiles':list(map(profile,blocks)),'a_mask':aa,'c_mask':cc,'counts':values})
  blocks=[[2,0]]*3;T,_=cycle(blocks);U,_=cycle(blocks,True);values=[]
  for V in (T,U):values.append(int(subprocess.check_output([str(exe)],input=('6 3 12 '+' '.join(map(str,V))+'\n').encode())))
  need(values==[509,545] and (values[0]-values[1])%8==4,'order-two exclusion boundary')
 print(json.dumps({'status':'PASS','independent_symbolic_records':records,'symbolic_coefficients':512,'actual_profile_controls':source_profiles,'literal_controls':literal,'order_two_boundary':values,'source_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (producer/'PROOF.md',producer/'verify.py',producer/'receipt.json',root.parent/'mod8_transfer'/'PROOF.md',root.parent/'covered_quotient'/'verify.py',root/'verify.py',root/'literal.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
