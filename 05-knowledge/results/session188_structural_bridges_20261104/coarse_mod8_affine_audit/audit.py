from pathlib import Path
from itertools import combinations,permutations,product
from math import comb
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'coarse_mod8_affine'
PINS={'PROOF.md':'dcbbfa8895f41c4649369bd11faf861672ced4fb010ab047fa9addfcdf969924','verify.py':'65991a25e6fdc8fd492bdfd1333f4fcde3b08cc91c0bcc453bae092ee29fd1e5','literal.cpp':'c7f27b14e27b0ab5f9ba88b7635755bb4bb2eb657f53040f07393df2c2a14b57','normal.json':'07c210dd54abde00805c70cbb06d9ef5947ad0b4d0203c1a8758efce4cb01592','optimized.json':'07c210dd54abde00805c70cbb06d9ef5947ad0b4d0203c1a8758efce4cb01592'}
def need(c,m):
 if not c:raise RuntimeError(m)
def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def op(T):return [((1<<len(T))-1)^(1<<i)^v for i,v in enumerate(T)]
def blow(Q,blocks):
 offsets=[];n=0
 for b in blocks:offsets.append(n);n+=len(b)
 T=[]
 for i,B in enumerate(blocks):
  cross=sum(((1<<len(C))-1)<<offsets[j] for j,C in enumerate(blocks) if Q[i]>>j&1)
  T.extend((v<<offsets[i])|cross for v in B)
 return T,offsets
def extend(T,a,c):
 n=len(T);U=T[:]+[1<<(n+1),0]
 for i in range(n):
  if a>>i&1:U[i]|=1<<n
  else:U[n]|=1<<i
  if c>>i&1:U[i]|=1<<(n+1)
  else:U[n+1]|=1<<i
 return U
def mask(b,B,offs):return sum(((1<<len(C))-1)<<offs[i] for i,C in enumerate(B) if b>>i&1)
def profile(T):
 n=len(T);q=[0]*(n+1)
 for w in permutations(range(n)):
  bad=sum(not(T[u]>>v&1) for u,v in zip(w,w[1:]));good=n-1-bad
  for h in range(1,n+1):
   if bad<=h-1<=n-1:q[h]+=comb(good,h-1-bad)
 return q
def evaluate(C,types):
 states={'S':(0,0),'P':(2,0),'E':(6,4),'O':(6,6),'C':(2,2)}
 v=C['K0']
 for i,t in enumerate(types):v+=states[t][0]*C['K2'][i]+states[t][1]*C['K3'][i]
 for i,j,k in C['K22']:v+=states[types[i]][0]*states[types[j]][0]*k
 return v%8
def enumerate_words(T,k):
 out=[];left=list(k)+[1,1]
 def go(w):
  if not any(left):out.append(tuple(w));return
  for i,v in enumerate(left):
   if v and(not w or T[w[-1]]>>i&1):left[i]-=1;go(w+[i]);left[i]+=1
 go([]);return out
def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'frozen '+f)
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','status')
 for f,h in d['dependency_pins'].items():need(hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()==h,'dependency '+f)
 frozen=json.loads((ROOT.parent/'five_profile_types/normal.json').read_text());jobs={};occurrences=0;lookup={}
 def add(T,v):
  nonlocal occurrences
  key=tuple(T);need(key not in jobs or jobs[key]==v,'duplicate expected count');jobs[key]=v;occurrences+=1
 def queue(Q,B,qs,values):
  for side,P in enumerate((Q,op(Q))):
   T,offs=blow(P,B)
   for (a,c),vv in zip(qs,values):add(extend(T,mask(a,B,offs),mask(c,B,offs)),vv[side])
 for r in d['small_model_controls']:
  Q=r['quotient_rows'];m=len(Q);need(len(r['models'])==1+2*m+comb(m,2),'small-model count')
  for z in r['models']:
   B=[tt(n) for n in z['sizes']];need(sum(z['sizes'])+2<=m+4,'small-model physical bound');queue(Q,B,r['block_queries'],z['physical_values'])
  for qi,(a,c) in enumerate(r['block_queries']):
   want=[]
   for side in (0,1):
    vv=[z['physical_values'][qi][side] for z in r['models']];f=vv[0];K2=[];K3=[];K22=[]
    for i in range(m):
     x,y=vv[1+2*i:3+2*i];need((x-f)%2==0 and (y-3*x+2*f)%12==0,'exact extraction divisibility');K2.append((x-f)//2);K3.append((y-3*x+2*f)//6)
    for j,(i,k) in enumerate(combinations(range(m),2)):
     num=vv[1+2*m+j]-vv[1+2*i]-vv[1+2*k]+f;need(num%4==0,'mixed exact extraction');K22.append([i,k,num//4])
    C=dict(K0=f,K2=K2,K3=K3,K22=K22);need(min([f]+K2+K3+[z[2] for z in K22])>=0,'nonnegative coefficients');want.append(C)
   need(r['extracted_coefficients'][qi]==dict(block_query=[a,c],coefficients=want),'complete coefficient equality');lookup[(tuple(Q),a,c)]=want
 canonical=[]
 for r,old in zip(d['canonical_full_banks'],frozen['canonical_records']):
  need(r['types']==old['types'],'canonical order');B=old['child_rows'];profiles=[profile(b) for b in B];unit=1
  for q in profiles:unit=unit*q[1]%8
  need(unit==r['odd_H_product_mod8'],'canonical scalar restoration');vals=[]
  for row,(a,c),vv in zip(r['records'],product(range(8),repeat=2),old['values']):
   pred=[unit*evaluate(C,r['types'])%8 for C in lookup[((2,4,1),a,c)]];need(row==[a,c,*vv,*pred] and pred==[z%8 for z in vv],'all complete canonical vectors');vals.append(vv)
  queue([2,4,1],B,list(product(range(8),repeat=2)),vals);canonical.append([r['types'],[[(x*pow(unit,-1,8))%8 for x in vv] for vv in vals]])
 for r in d['actual_profile_variants']:
  B=r['child_rows'];qs=[profile(b) for b in B];need(qs==r['profiles'],'actual complete profiles by independent word/cut enumeration');unit=1;types=[]
  for q in qs:
   unit=unit*q[1]%8;p=[(q[k] if k<len(q) else 0)*pow(q[1],-1,8)%8 for k in (2,3)];types.append({(0,0):'S',(2,0):'P',(6,4):'E',(6,6):'O',(2,2):'C'}[tuple(p)])
  need(types==r['types'] and unit==r['odd_H_product_mod8'],'actual normalized five-type domain');T,offs=blow(r['quotient_rows'],B)
  need(r['physical_queries']==[[mask(a,B,offs),mask(c,B,offs)] for a,c in r['block_queries']],'actual physical masks')
  for (a,c),vv,pred in zip(r['block_queries'],r['physical_values'],r['predictions_mod8']):need(pred==[unit*evaluate(C,types)%8 for C in lookup[(tuple(r['quotient_rows']),a,c)]]==[z%8 for z in vv],'actual full response restoration')
  queue(r['quotient_rows'],B,r['block_queries'],r['physical_values'])
 need(occurrences==d['physical_H_counts']+16000,'all source and canonical literal occurrences')
 with tempfile.TemporaryDirectory(prefix='coarse-mod8-independent-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True)
  wire=str(len(jobs))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in jobs);got=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()))
 need(got==list(jobs.values()),'every distinct actual physical integer count')
 involutions=[];mapped=0;names={r['name']:r['quotient_rows'] for r in d['small_model_controls']}
 for name,a,c,i,j,kind,n,sha in d['involution_instances']:
  Q=names[name];k=[1]*len(Q);k[i]=3
  if j is not None:k[j]=2 if kind=='i3j2' else 3
  ww=enumerate_words(extend(Q,a,c),k);need(len(ww)==n and n%2==0,'complete word census');W=set(ww)
  for w in ww:
   p=[z for z,v in enumerate(w) if v==i];P=w[:p[0]];A=w[p[0]+1:p[1]];B=w[p[1]+1:p[2]];S=w[p[2]+1:];v=P+(i,)+B+(i,)+A+(i,)+S;need(A!=B and v in W,'fixed-free exchange including shared repeated label')
  need(hashlib.sha256(json.dumps(ww,separators=(',',':')).encode()).hexdigest()==sha,'additional census digest');involutions.append([name,a,c,i,j,kind,n]);mapped+=n
 need(mapped==d['total_mapped_return_words']==8080 and len(involutions)==7876,'involution accounting')
 norm={tuple(t):v for t,v in canonical};base=[(a-b)%8 for a,b in norm[('P',)*3]];derivatives=[]
 for i in range(3):
  for typ in ('E','C'):
   t=['P']*3;t[i]=typ;v=[(a-b-base[j])%8 for j,(a,b) in enumerate(norm[tuple(t)])];need(set(v)<={0,4},'binary derivatives');derivatives.append([x//4 for x in v])
 need(derivatives==d['cyclic_nontrivial_contrast_derivatives'],'full derivatives');basis={}
 for v in derivatives:
  x=sum(z<<i for i,z in enumerate(v))
  while x:
   p=x.bit_length()-1
   if p not in basis:basis[p]=x;break
   x^=basis[p]
 full=[[(a-b)%8 for a,b in norm[t]] for t in product('PEOC',repeat=3)];need(len(basis)==6 and len(set(map(tuple,full)))==64,'complete64-vector rank and distinctness');need([t for t,v in zip(product('PEOC',repeat=3),full) if not any(v)]==[('E','E','E')],'unique finite cyclic blind state')
 return dict(status='PASS',source_pins=PINS,dependency_pins=d['dependency_pins'],literal_occurrences=occurrences,distinct_physical_graphs=len(jobs),maximum_physical_order=max(map(len,jobs)),complete_physical_counts=got,canonical_normalized_banks=canonical,involution_controls=involutions,total_mapped_words=mapped,full_cyclic_nontrivial_contrast_vectors=full,scope='Independent physical counts and exact algebra audit of all-order coarse mod8 response reduction. Not an oracle query bound or a nonannihilation theorem.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
