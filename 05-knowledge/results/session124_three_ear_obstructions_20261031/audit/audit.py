#!/usr/bin/env python3
"""Independent exact three-ear audit; no producer imports or numerical solvers."""
import itertools,json,math,copy,hashlib
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'code';checks=0
ORDERS=list(itertools.permutations(range(3)))
STATES=list(itertools.product(range(6),range(3),range(3),range(3),range(2)))
INDEX={s:j for j,s in enumerate(STATES)};D=len(STATES)
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ValueError(msg)
def transitions(letter):
 columns=[]
 for order,*flags in STATES:
  ears=ORDERS[order];targets=[]
  for run in range(4):
   if run and flags[run]==0 and not(letter>>ears[run-1]&1):continue
   nxt=flags.copy()
   nxt[run]=1 if run==3 or not(letter>>ears[run]&1) else 2
   targets.append(INDEX[(order,*nxt)])
  columns.append(targets)
 return columns
MATS=[transitions(x) for x in range(8)]
def actrow(row,letter):return [sum(row[t] for t in column) for column in MATS[letter]]
def rowword(row,word):
 for letter in reversed(word):row=actrow(row,letter)
 return row
def actvec(v,letter):
 out=[0]*D
 for j,col in enumerate(MATS[letter]):
  for t in col:out[t]+=v[j]
 return out
def output(core):
 arcs={(0,1),(1,2),(2,0)} if core else {(i,j) for i in range(3) for j in range(i+1,3)}
 out=[]
 for o,p,q,r,s in STATES:
  ears=ORDERS[o];valid=p in [0,1]
  for flag,left,right in [(q,ears[0],ears[1]),(r,ears[1],ears[2])]:valid=valid and (flag==1 or flag==0 and (left,right) in arcs)
  out.append(int(valid))
 return out
INITIAL=[int(flags==[0,0,0,0]) for o,*flags in STATES]
ELL=[[int(o==k and q==r==s==0) for o,p,q,r,s in STATES] for k in range(6)]
def occupancy(o,flags):
 bits=[0]+[e+1 for e in ORDERS[o]]
 return sum(1<<b for b,f in zip(bits,flags) if f)
OCC=[[[int(o==k and occupancy(o,[p,q,r,s])==mask) for o,p,q,r,s in STATES] for mask in range(16)] for k in range(6)]
REL=[[x-y for x,y in zip(OCC[k][mask],OCC[0][mask])] for k in range(1,6) for mask in range(16)]
def primitive(row):
 g=math.gcd(*row)
 return [x//g for x in row] if g else row

def verify_certificate(cert,occupancy_relations):
 core=cert['core'];word=cert['word'];rank=cert['rank']
 need(core in [0,1] and len(word)==2 and all(type(x)==int and x in range(8) for x in word),'certificate state')
 need(type(rank)==int and 1<=rank<=16,'row count')
 c=output(core);once=rowword(c,word);twice=rowword(once,word)
 R=primitive([x-17*y+16*z for x,y,z in zip(twice,once,c)]);rows=[]
 for _ in range(rank):rows.append(R);R=primitive(rowword(R,word))
 rows+=ELL
 if occupancy_relations:rows+=REL
 rational=[Fraction(x) for x in cert['multipliers']];need(len(rational)==len(rows),'multiplier arity')
 denom=math.lcm(*(x.denominator for x in rational));y=[int(x*denom) for x in rational]
 need(sum(y[rank:rank+6])<0,'strict negative normalization constant')
 for j in range(D):need(sum(t*r[j] for t,r in zip(y,rows))>=0,'nonnegative separating coordinate')

def brute_h(core,word):
 n=3+len(word);arcs=[0]*n
 for i in range(3):
  for j in range(i+1,3):
   left,right=(j,i) if core and (i,j)==(0,2) else (i,j);arcs[left]|=1<<right
 for k,x in enumerate(word,3):
  for i in range(3):
   left,right=(i,k) if x>>i&1 else (k,i);arcs[left]|=1<<right
  for j in range(3,k):arcs[j]|=1<<k
 dp=[[0]*n for _ in range(1<<n)]
 for j in range(n):dp[1<<j][j]=1
 for mask in range(1,1<<n):
  for j,count in enumerate(dp[mask]):
   if not count:continue
   available=arcs[j]&~mask
   while available:
    bit=available&-available;k=bit.bit_length()-1;dp[mask|bit][k]+=count;available-=bit
 return sum(dp[-1])

def check_quotient():
 for letter in range(8):
  quotient=[[0]*16 for _ in range(16)]
  for mask in range(16):
   quotient[mask|1][mask]+=1
   for ear in range(3):
    bit=1<<(ear+1)
    if mask&bit or letter>>ear&1:quotient[mask|bit][mask]+=1
  for o in range(6):
   need(actrow(ELL[o],letter)==ELL[o],'prefix normalization invariant')
   for mask in range(16):
    actual=actrow(OCC[o][mask],letter)
    for j,(source_o,*flags) in enumerate(STATES):need(actual[j]==(quotient[mask][occupancy(source_o,flags)] if source_o==o else 0),'occupancy quotient coordinate')
 for e in ELL:need(sum(x*y for x,y in zip(e,INITIAL))==1,'initial normalization')
 for row in REL:need(sum(x*y for x,y in zip(row,INITIAL))==0,'initial occupancy agreement')

def complement_row(core,word):
 r=output(core)
 for root in [0,1,1,1,4,4,4,9,9,16]:r=[x-root*y for x,y in zip(rowword(r,word),r)]
 return r

def main():
 need(D==324,'dimension');check_quotient();comparisons=0
 for L in range(5):
  for word in itertools.product(range(8),repeat=L):
   v=INITIAL
   for x in word:v=actvec(v,x)
   for core in [0,1]:need(sum(x*y for x,y in zip(output(core),v))==brute_h(core,word),'transfer versus independent subset DP');comparisons+=1
 files=['farkas_certificates.json','farkas_occupancy_certificates.json'];certs=[];hashes={}
 for use_occupancy,name in enumerate(files):
  raw=(SOURCE/name).read_bytes();hashes[name]=hashlib.sha256(raw).hexdigest();bank=json.loads(raw)
  expected={(core,x,y) for core in [0,1] for x,y in itertools.product(range(8),repeat=2) if not(x&y) and ((x|y).bit_count()==2 if use_occupancy else (x|y).bit_count()%2==1)}
  need(len(bank)==len(expected) and {(c['core'],*c['word']) for c in bank}==expected,'exact exclusion universe')
  for c in bank:verify_certificate(c,use_occupancy)
  certs.append(bank)
 expected_row=[3840*sum(row[j] for row in ELL) for j in range(D)]
 for core in [0,1]:
  for x in range(8):
   actual=complement_row(core,[x,x^7])
   for u,v in zip(actual,expected_row):need(u==v,'complement polynomial operator identity')
 need(Fraction(23040,6*(1-4)**3*(1-9)**2*(1-16))==Fraction(4,27),'forced cubic coefficient')
 need(complement_row(0,[0,0])!=expected_row,'hostile complement template')
 corrupted=expected_row.copy();corrupted[0]+=1;need(complement_row(0,[0,7])!=corrupted,'hostile complement RHS')
 mutations=[]
 for use_occupancy,bank in enumerate(certs):
  for label,action in [('zero',lambda z:z.__setitem__('multipliers',['0']*len(z['multipliers']))),('negate',lambda z:z.__setitem__('multipliers',[str(-Fraction(x)) for x in z['multipliers']])),('zero_word',lambda z:z.__setitem__('word',[0,0])),('truncate',lambda z:z['multipliers'].pop())]:
   z=copy.deepcopy(bank[0]);action(z)
   try:verify_certificate(z,use_occupancy)
   except ValueError as e:mutations.append([use_occupancy,label,str(e)])
   else:raise ValueError('hostile certificate accepted')
 need(sum(not(x&y) for x in range(8) for y in range(8))==27,'growth eligibility count')
 need(len({tuple(c['word']) for bank in certs for c in bank})==26,'all nonzero eligible pairs excluded')
 report={'status':'passed','checks':checks,'states':D,'subset_DP_comparisons':comparisons,'prefix_invariance_row_identities':48,'occupancy_coordinate_equalities':248832,'cross_order_relations':80,'Farkas_certificates':52,'Farkas_coordinate_inequalities':16848,'excluded_nonzero_templates_per_core':26,'growth_ineligible_templates':37,'complement_identities':16,'complement_coordinate_equalities':5184,'forced_cubic_coefficient':'4/27','rejected_certificate_mutations':mutations,'rejected_identity_mutations':2,'certificate_sha256':hashes}
 print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
