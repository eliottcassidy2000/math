import sys
if '--fail' in sys.argv: raise RuntimeError('deliberate producer sentinel')
import contextlib,io,json,math
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()): import verify as v
P=Path(__file__).resolve().parent
cs={h:json.loads((P/f'nonunit{h}.json').read_text()) for h in (13,29)}
for h,c in cs.items():
 v.check(v.literal(c['adj'],c['mask'])==tuple(c['core']),'literal core')
 c['xs']=v.word(c['X']);c['ys']=v.word(c['Y']);u=int(c['X'][-1]);w=int(c['Y'][0])
 a=v.chain([v.single(u),c['core'],c['ys'],v.single(w)])
 b=v.chain([v.single(u),c['xs'],c['core'],v.single(w)])
 expect=[0]*7;expect[4]=-c['cycle_gain']
 v.check(tuple(y-x for x,y in zip(a,b))==tuple(expect),'flanked identity')

def family(h,k,j):
 if h==45:return v.family(45,k,j)
 c=cs[h];u=int(c['X'][-1]);w=int(c['Y'][0]);blocks=[]
 if not u:blocks.append(v.single(1))
 blocks.append(v.single(u))
 if j:blocks.append(v.repeat(c['xs'],j))
 blocks.append(c['core'])
 if k-j:blocks.append(v.repeat(c['ys'],k-j))
 blocks.append(v.single(w))
 return v.chain(blocks)

for h,c in cs.items():
 for k in [0,1,2,5,20,100]:
  base=family(h,k,0)
  for j in sorted({0,k//2,k}):
   cur=family(h,k,j);expect=[0]*7;expect[4]=-c['cycle_gain']*j
   v.check(tuple(y-x for x,y in zip(base,cur))==tuple(expect),'family drift')

def pair(h,K,L,i,j):return v.H(v.mul(family(h,K,i),family(45,L,j)))
records=[]
for h,K,Ls,g,overhead in [(13,12,[224,250,300,1000],5,167),(29,28,[134,150,200,1000],3,808)]:
 a=45*g;b=h;C=(a-1)*(b-1)
 for L in Ls:
  M=a*K+b*L;base=pair(h,K,L,0,L)
  for t in range(C,M-C+1):
   i,j=v.selector(a,b,K,L,t)
   if t in {C,C+1,M//2,M-C-1,M-C}:
    v.check(pair(h,K,L,i,L-j)==base+2*t,'weighted physical selector')
  records.append(dict(h=h,K=K,L=L,vertices=12*L+overhead,consecutive_odd_values=M-2*C+1,lower_bit_length=(base+2*C).bit_length()))
# Inner h377 fiber: a145,b39; K38,L149 gives378 consecutive values.
IK,IL=38,149;IC=144*38
innerbase=v.mul(family(13,IK,0),family(29,IL,0))
def inner(r):
 i,j=v.selector(145,39,IK,IL,IC+r)
 return v.mul(family(13,IK,i),family(29,IL,j))
ib=inner(0)
for r in range(378):
 cur=inner(r);expect=[0]*7;expect[4]=-2*r
 v.check(tuple(y-x for x,y in zip(ib,cur))==tuple(expect),'inner unit fiber')
triple=[]
for L in [44,50,100,1000]:
 C=44*376;M=45*376+377*L
 baseline=v.H(v.mul(inner(0),family(45,L,L)))
 for t in range(C,M-C+1):
  i,j=v.selector(45,377,376,L,t)
  if t in {C,C+1,M//2,M-C-1,M-C}:
   actual=v.H(v.mul(inner(i),family(45,L,L-j)))
   v.check(actual==baseline+2*t,'triple physical selector')
 triple.append(dict(L=L,vertices=12*L+4662,consecutive_odd_values=M-2*C+1,lower_bit_length=(baseline+2*C).bit_length()))
crt=[]
reps={((6525*i+1755*j)%377):(6525*i+1755*j,i,j) for i in range(13) for j in range(29)}
v.check(len(reps)==377 and max(x[0] for x in reps.values())==127440,'CRT residues')
for L in [337,400,1000]:
 base=v.H(v.chain([family(13,12,0),family(29,28,0),family(45,L,L)]))
 lo=127064;hi=377*L+376
 for z in range(lo,hi+1):
  r,i,j=reps[z%377];t=(z-r)//377
  v.check(0<=t<=L and r+377*t==z,'CRT bounded selector')
  if z in {lo,lo+1,(lo+hi)//2,hi-1,hi}:
   actual=v.H(v.chain([family(13,12,i),family(29,28,j),family(45,L,L-t)]))
   v.check(actual==base+2*z,'CRT physical count')
 crt.append(dict(L=L,vertices=12*L+962,consecutive_odd_values=hi-lo+1,lower_bit_length=(base+2*lo).bit_length(),lower_sha256=v.hashlib.sha256((base+2*lo).to_bytes(((base+2*lo).bit_length()+7)//8,'big')).hexdigest()))
print(json.dumps(dict(crt=crt,status='PASS',new_loops=[dict(h=h,length=len(c['X']),gain=c['cycle_gain']) for h,c in cs.items()],pairs=records,triple=triple),indent=2))
