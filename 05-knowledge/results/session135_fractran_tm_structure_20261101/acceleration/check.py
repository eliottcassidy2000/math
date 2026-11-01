from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib,json,random
from accelerate import Program,certify_word,eligible_interval,run,PRIMEGAME

def need(x,m):
 if not x:raise RuntimeError(m)
def literal(rows,N,limit):
 trace=[]
 for _ in range(limit):
  for i,f in enumerate(rows):
   if N%f.denominator==0:N=N//f.denominator*f.numerator;trace.append(i);break
  else:return N,trace,True
 return N,trace,False

def main():
 checks=0;intervals=0;rng=random.Random(135)
 for _ in range(3000):
  b=[rng.randrange(-5,10) for _ in range(3)];D=[rng.randrange(-3,4) for _ in range(3)];needv=[rng.randrange(5) for _ in range(3)]
  got=eligible_interval(b,D,needv)
  for t in range(31):
   expect=all(x+t*d>=v for x,d,v in zip(b,D,needv));actual=got is not None and t>=got[0] and(got[1] is None or t<=got[1])
   need(actual==expect,'guard interval');intervals+=1
 rows_pool=[(1,2),(3,2),(5,3),(2,5),(1,10),(2,3),(3,5),(1,1),(2,1)]
 for _ in range(500):
  rows=rng.sample(rows_pool,rng.randrange(1,6));p=Program(rows);N=rng.randrange(1,250);x,inert=p.encode(N)
  need(p.decode(x,inert)==N,'valuation encoding')
  end,tr,halt=literal(p.rows,N,40)
  for length in range(1,min(len(tr),8)+1):
   w=tr[:length];c=certify_word(p,x,w);R=c['repetitions']
   repetitions=min(12, R if R is not None else 12)
   en,tt,hh=literal(p.rows,N,len(w)*repetitions)
   need(tt==w*repetitions,'accelerated prefix changed priority')
   need(p.decode([v+repetitions*d for v,d in zip(x,c['direction'])],inert)==en,'accelerated endpoint')
   if R is not None and R<=12:
    en,tt,hh=literal(p.rows,N,len(w)*(R+1))
    need(tt!=w*(R+1),'maximal repeat certificate');checks+=1
   checks+=1
  if halt:
   accelerated=run(p,x,8,200)
   need(accelerated['status']=='halted' and accelerated['primitive_steps']==len(tr),'halt/cost mismatch')
   need(p.decode(accelerated['state'],inert)==end,'halt endpoint mismatch');checks+=1
 p=Program([(1,10),(5,2)]);x,r=p.encode(8);c=certify_word(p,x,[1]);out=run(p,x)
 need(c['repetitions']==1 and c['limiter']['kind']=='earlier_instruction','priority activation hostile')
 need(p.decode(out['state'],r)==5 and 5**3!=5,'blind drain hostile')
 # A genuine guarded two-step countdown loop, then terminal control.
 p2=Program([(15,14),(7,5),(11,7)]);x,r=p2.encode(7*2**30)
 c2=certify_word(p2,x,[0,1]);out2=run(p2,x)
 need(c2['repetitions']==30 and c2['primitive_steps']==60,'two-phase drain count')
 need(out2['status']=='halted' and out2['primitive_steps']==61,'two-phase terminal')
 need(p2.decode(out2['state'],r)==11*3**30,'two-phase terminal integer')
 big=[0]*len(p2.primes);big[p2.primes.index(2)]=10**18;big[p2.primes.index(7)]=1
 huge=certify_word(p2,big,[0,1]);need(huge['repetitions']==10**18,'huge symbolic loop')
 prime=Program(PRIMEGAME);periods=[]
 for N,L in [(11,2),(4301,6)]:
  x,r=prime.encode(N);o=run(prime,x,12)
  need(o['status']=='certified_infinite_word','known periodic point')
  c=o['infinite_word'];need(len(c['word'])==L and not any(c['direction']),'known period certificate')
  periods.append({'input':N,'period':L,'word':c['word']})
 drift=Program([(2,1)]);od=run(drift,[0]);need(od['status']=='certified_infinite_word' and od['infinite_word']['direction']==[1],'drifting nontermination')
 # A fixed residue does not repair cap under decrement/halving.
 cap=5;mod=6;a=12;b=18
 need((min(a,cap),a%mod)==(min(b,cap),b%mod),'cap-residue equal')
 need(min(a-8,cap)!=min(b-8,cap),'countdown distinguishes cap-residue states')
 return {'status':'PASS','guard_interval_checks':intervals,'random_word_checks':checks,
  'priority_hostile':{'program':[[1,10],[5,2]],'input':8,'actual_halt':5,'blind_drain_halt':125},
  'huge_repeat_count':huge['repetitions'],'huge_primitive_steps':huge['primitive_steps'],
  'known_primegame_periods':periods,'drifting_infinite_word':od['infinite_word'],
  'scope':'Exact repetition and cost certificates; bounded tests are not a general halting decision procedure.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
