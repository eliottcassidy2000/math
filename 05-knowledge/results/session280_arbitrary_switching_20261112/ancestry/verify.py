"""Exact hidden-itinerary decoder; explicit checks survive python -O."""
from pathlib import Path
import importlib.util,itertools,random,json
P=Path(__file__).resolve().parents[2]/'session279_interval_tightening_20261112/physical/verify.py'
spec=importlib.util.spec_from_file_location('prior279',P); old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
def check(p,m):
 if not p: raise RuntimeError(m)
def thin(x):
 a,w,d=x
 return a>=1 and a%2==1 and 0<=d<=a<=w<2*a and (w-d)%2==0

def decode_step(x):
 """Unique active predecessor IN THE THIN STRIP, else None."""
 if not thin(x):return None
 aout,wout,dout=x
 if dout<=0 or dout%2!=1:return None
 # Select unique t=2^n with 3*t*dout < 2*wout < 6*t*dout.
 t=4
 while 6*t*dout<=2*wout:t*=2
 if not 3*t*dout<2*wout<6*t*dout:return None
 n=t.bit_length()-1
 a=(wout*pow(dout,-1,t))%t
 if not t//2<a<t or a%2!=1:return None
 word=format(a-1,'0%db'%n)
 w=(wout-a*dout)//t
 h,b,_=old.block(word)[0]
 parent=(dout,w,h*dout+b*w-aout)
 if not thin(parent):return None
 check(old.mv(old.block(word),parent)==x,'inverse recomposition')
 return word,parent

def decode_to_seed(x,seed):
 if not thin(seed) or not thin(x):return None
 rev=[]
 while x[0]>seed[0]:
  z=decode_step(x)
  if z is None:return None
  word,x=z;rev.append(word)
 if x!=seed:return None
 return tuple(reversed(rev))

def main():
 rng=random.Random(28011);seed=(3,3,1);checks=0
 words=['1'+''.join(z)+'0' for n in range(2,11) for z in itertools.product('01',repeat=n-2)]
 # Structural inequality uses h >= (a+1)/2 and b >= t/2+1.
 for word in words:
  mat=old.block(word);h,b,_=mat[0];a,t,_=mat[1]
  check(2*h>=a+1 and 2*b>=t+2,'coefficient lower bound')
  for j in range(7):
   A=2*rng.randrange(1,1000)+1;W=rng.randrange(A,2*A)
   possible=list(range(W%2,A+1,2));d=rng.choice(possible);x=(A,W,d)
   y=old.mv(mat,x);check(thin(y),'strip invariance')
   check(decode_step(y)==(word,x),'blind step')
   checks+=1
 # Literal physical sources/small constructions independently counted.
 literal=0
 for ws in [('10',),('100',),('110',),('10','100'),('100','110'),('10','10','10')]:
  x=old.packet(ws);check(x==old.literal_packet(ws),'literal packet');literal+=1
 # Empty prefix and one initial active triangle are separate.
 for ws in itertools.product(words[:7],repeat=3):
  x=seed
  for word in ws:x=old.mv(old.block(word),x)
  check(decode_to_seed(x,seed)==ws,'three-block labels')
 # Very long varying/non-Pisot schedules, no label/depth input to decoder.
 maxbits=0
 for trial in range(120):
  ws=tuple(rng.choice(words) for _ in range(rng.randrange(1,151)))
  x=seed
  for word in ws:x=old.mv(old.block(word),x)
  check(decode_to_seed(x,seed)==ws,'blind full itinerary');maxbits=max(maxbits,x[0].bit_length())
 # Explicit cone/parity formal seed at forbidden A7 stays formal.
 formal=(7,7,1);check(thin(formal),'formal example')
 # Thin restriction cannot simply be omitted: two distinct physical predecessors.
 # Find deterministic stored witness for image collision across lengths from broad cone.
 broad=[]
 short=[''.join(z) for n in range(1,4) for z in itertools.product('01',repeat=n)]
 for u,v in itertools.product(short,repeat=2):
  x=old.packet((u,v));broad.append((x,(u,v)))
 seen={};hostile=None
 for x,source in broad:
  for word in words[:15]:
   y=old.mv(old.block(word),x)
   if y in seen:
    px,pw,ps=seen[y]
    if pw!=word and (px,pw)!=(x,word):hostile={'output':y,'first':{'source':ps,'packet':px,'word':pw},'second':{'source':source,'packet':x,'word':word}};break
   seen[y]=(x,word,source)
  if hostile:break
 check(decode_to_seed(seed,seed)==(),'empty itinerary')
 check(decode_step((7,12,1)) is None,'dyadic boundary is rejected')
 result={'seed':seed,'active_words_tested':len(words),'thin_step_checks':checks,'literal_cases':literal,'all_7_letter_three_block_schedules':7**3,'long_unknown_schedules':120,'max_tested_count_bits':maxbits,'broad_cone_collision_found':hostile,'scope':'All-order thin-strip hidden-itinerary theorem; finite controls only. No scalar coverage.'}
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
