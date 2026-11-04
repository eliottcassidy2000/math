"""Independent arithmetic audit of the reset compiler; no external packages."""
import itertools,json,importlib.util,hashlib
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parent.parent/'collatz'
spec=importlib.util.spec_from_file_location('author_reset',ROOT/'reset_phase.py');author=importlib.util.module_from_spec(spec);spec.loader.exec_module(author)

def need(ok,s):
 if not ok:raise RuntimeError(s)
def carry_direct(w):
 return sum(3**(len(w)-1-i)*2**sum(w[:i]) for i in range(len(w)))
def inverse_route(w,b):
 need(b%2==0,'even reset');x=(2**b-1)//3;back=[x]
 for a in reversed(w):
  y=2**a*x-1;need(y%3==0,'integral inverse intermediate');x=y//3;need(x>0 and x%2==1,'positive odd inverse');back.append(x)
 return x,list(reversed(back))+[1]
def replay(n):
 z=3*n+1;a=0
 while z%2==0:z//=2;a+=1
 return z,a

def main():
 count=0;nonfirst=0;first=0;phase_hist={}
 for L in range(1,5):
  P=2*3**L;M=3**(L+1)
  for w in itertools.product(range(1,7),repeat=L):
   A=sum(w);B=carry_direct(w);phases=[b for b in range(0,P,2) if (pow(2,A+b,M)-pow(2,A,M)-3*B)%M==0]
   need(len(phases)==1,'one full exponent phase');beta=phases[0];need(beta%6 in (2,4),'phase class');need(author.reset_phase(w)==beta,'author phase agrees independent congruence')
   b=beta if beta>=4 else beta+P
   for exponent in (beta,b+P):
    n,states=inverse_route(w,exponent);need(n==(2**(A+exponent)-2**A-3*B)//M,'closed source formula')
    for i,a in enumerate(w+(exponent,)):
     nxt,actual=replay(states[i]);need((nxt,actual)==(states[i+1],a),'literal valuation replay')
    isfirst=all(x>1 for x in states[:-1]);need(isfirst==(exponent>=4),'first-hit iff reset at least four')
    first+=isfirst;nonfirst+=not isfirst
   count+=1
 # Independently validate every phase-address compilation through depth eight.
 controllers=[]
 for L in range(1,9):
  P=2*3**L;seen=set()
  for beta in range(2,P,2):
   if beta%3==0:continue
   w=author.prefix_for_phase(L,beta);need(len(w)==L and set(w)<=set(range(1,5)),'alphabet four')
   A=sum(w);B=carry_direct(w);M=3**(L+1)
   need((pow(2,A+beta,M)-pow(2,A,M)-3*B)%M==0,'independent compiled phase congruence');seen.add(beta)
  need(len(seen)==2*3**(L-1),'full spectrum count');controllers.append([L,len(seen)])
 # Direct modulo-9 dead-end argument at every forbidden phase in these ranges.
 forbidden=0
 for L in range(2,9):
  for beta in range(14,2*3**L,18):
   m=(pow(2,beta,27)-1)//3;need(m%9==7,'forbidden terminal residue');need((4*m-1)//3%3==0,'alphabet three forces a dead end');forbidden+=1
 # Large encoded letters exercise only the genuinely modular phase routine.
 large_words=[(10**30+1,), (10**40+1,10**35+2,7)]
 large=[]
 for w in large_words:
  beta=author.reset_phase(w);L=len(w);M=3**(L+1);A=sum(w)
  B=sum(pow(3,L-1-i,M)*pow(2,sum(w[:i]),M) for i in range(L))%M
  need((pow(2,A+beta,M)-pow(2,A,M)-3*B)%M==0,'huge input modular phase');cert=author.completion(w);need('carry' not in cert and 'source' not in cert and cert['carry_residue']==B,'compact huge-letter certificate');need(author.periodic_phase(w,4)==author.reset_phase(w*4),'compact huge-letter periodic phase');large.append({'letters':[str(x) for x in w],'beta':beta})
 # Center error and mismatch identities, with exact integer cross multiplication.
 periodic=0
 for w in [(1,),(2,),(1,2),(1,1,3),(3,1,2,4),(1,4,2)]:
  A=sum(w);B=carry_direct(w);L=len(w);D=2**A-3**L;old=None
  for r in range(1,11):
   wr=w*r;Ar=sum(wr);Br=carry_direct(wr);M=3**(r*L+1)
   need(D*Br==B*(2**Ar-3**(r*L)),'repeated carry identity')
   target=((D+3*B)*pow(D,-1,M))%M;beta=author.reset_phase(wr)
   need(pow(2,beta,M)==target,'center reset address')
   if old is not None:need(beta%old[1]==old[0],'phase compatibility')
   old=beta,2*3**(r*L);periodic+=1
 rational_tests=0
 for L in range(1,6):
  for w in itertools.product(range(1,7),repeat=L):
   A=sum(w);B=carry_direct(w);D=2**A-3**L;C=Fraction(D+3*B,D)
   k=(abs(C.numerator)&-abs(C.numerator)).bit_length()-1
   need(k==w[0] and C.denominator%2==1,'exact dyadic valuation of center address')
   rational_possible=C.denominator==1 and abs(C.numerator)==2**k
   fixed=(all(a==1 for a in w) or all(a==2 for a in w))
   need(rational_possible==fixed,'rational-address classification arithmetic')
   if fixed:
    beta=author.reset_phase(w)
    need(beta==(3**L+1 if w[0]==1 else 2),'correct half-integer versus integer branch')
   rational_tests+=1
 inspected={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'PROOF.md',ROOT/'reset_phase.py')}
 return {'status':'PASS','independent_words':count,'literal_routes':first+nonfirst,'first_hit_routes':first,'non_first_hit_routes':nonfirst,'compiled_phase_counts':controllers,'compiled_phases':sum(v for _,v in controllers),'forbidden_phase_checks':forbidden,'periodic_center_checks':periodic,'rational_address_word_checks':rational_tests,'large_encoded_letter_modular_checks':large,'inspected_sha256':inspected,'complexity_review':'PASS after repair: reset_phase, prefix_for_phase, completion(expand=False), and periodic_phase use modular carries. The symbolic certificate retains word and cost sidecars. Runtime depends on input bit length as well as L; exact B and source are materialized only with expand=True.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
