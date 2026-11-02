"""Independent framework audit with explicit producer-API conformance checks.

Expected profiles/direct Hamilton counts come from our previously independent
backward endpoint / anchored-partition kernel. The complete <=6 census comes
from enumerate_profiles.cpp, a different forward endpoint / ordered-first-block
algorithm. Producer certify is imported only as the API being tested, never to
calculate an expected value.
"""
from pathlib import Path
from itertools import product,combinations
import importlib.util,json,hashlib,sys
BASE=Path(__file__).resolve().parents[1];D=Path(__file__).resolve().parent
SOURCE=BASE/'zero_framework_later';KERNEL=BASE/'independent/audit_profile_dichotomy.py'
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
G=load('independent_graph_kernel',KERNEL)
sys.path.insert(0,str(SOURCE));C=load('producer_certificate_under_test',SOURCE/'certify.py')
def need(b,message):
 if not b:raise RuntimeError(message)
def pad(v,n):return list(v)+[0]*max(0,n-len(v))
def word_count(r,s,t):
 return sum(tuple(sum((start+k)%3==j for k in range(r+s+t)) for j in range(3))==(r,s,t) for start in range(3))
def literal_kernel(d):
 return [[sum(d[r-1]*word_count(r,s,t) for r in range(1,len(d)+1)) for t in range(1,len(d)+2)] for s in range(1,len(d)+2)]
def closed_kernel(d):
 z=[0]+list(d)+[0,0];n=len(d)+1
 return [[z[i]+3*z[i+1]+z[i+2] if i==j else z[min(i,j)+1]+z[max(i,j)+1] if abs(i-j)==1 else 0 for j in range(n)] for i in range(n)]
def scalar(K,a,b):
 a=pad(a,len(K));b=pad(b,len(K));return sum(x*a[i]*b[j] for i,row in enumerate(K) for j,x in enumerate(row))
def reference_negative(K,v):
 if not v or min(v)<=0 or any(v[i]**2<v[i-1]*v[i+1] for i in range(1,len(v)-1)):return False
 if scalar(K,v,v)>=0:return False
 for c in range(len(v)):
  if all(x>=0 for row in K[:c] for x in row) and all(x<=0 for row in K[c+1:] for x in row):return True
 return False
def lift(v):
 return [(r+1)*((v[r] if r<len(v) else 0)+(v[r-1] if r else 0)) for r in range(len(v)+1)]
def multiply(v,coeff):
 powers=[list(v)]
 for _ in coeff[1:]:powers.append(lift(powers[-1]))
 return [sum(c*(powers[e][r] if r<len(powers[e]) else 0) for e,c in enumerate(coeff)) for r in range(len(powers[-1]))]
def actual_calls():
 point=(0,);tt=(6,4,0);tri=(2,4,1);v6=G.ordinal(2)
 checks=[]
 for old,new,base,pstat,nstat in [(tt,tri,point,'SUFFICIENT_CERTIFICATE','NOT_CERTIFIED'),(tri,tt,point,'NOT_CERTIFIED','SUFFICIENT_CERTIFICATE'),(tt,tt,v6,'WEAK_NONNEGATIVITY_ONLY','NOT_CERTIFIED'),(tri,tt,v6,'NOT_CERTIFIED','SUFFICIENT_CERTIFICATE')]:
  z=C.certify_actual(old,new,base);need(z['actual_profile_difference']==[b-a for a,b in zip(G.profile(old),G.profile(new))],'actual wrapper difference');need(z['base_profile']==list(G.profile(base)),'actual wrapper base');need(z['positive']['status']==pstat and z['negative']['status']==nstat,'actual wrapper status');checks.append(z)
 invalid=[([],point,None),((1,),point,None),((2,1),(2,0),None),((4,0),(2,0),None),(point,(2,0),None),(point,point,[]),((0.0,),point,None),(point,point,(0.0,))]
 for old,new,base in invalid:
  try:C.certify_actual(old,new,base)
  except ValueError:continue
  raise RuntimeError('invalid literal tournament accepted')
 
 for fn,args in [(C.cyclic_kernel,([1.0],)),(C.cyclic_kernel,([],)),(C.negative_certificate,([-1],[1.0]))]:
  try:fn(*args)
  except ValueError:continue
  raise RuntimeError('noninteger/empty formal input accepted')
 return checks,len(invalid)
def main():
 enum=json.loads((D/'profiles_O2.json').read_text());need((D/'profiles_O2.json').read_bytes()==(D/'profiles_UBSan.json').read_bytes(),'independent census modes differ')
 pairchecks=0;entrychecks=0
 for row in enum['orders']:
  n=row['n'];fs=row['profiles'];need(len(fs)==row['profile_count'],'profile count');need(len(list(combinations(fs,2)))==row['unordered_profile_pairs'],'pair count')
  for f,g in combinations(fs,2):
   d=[b-a for a,b in zip(f,g)];K=literal_kernel(d);need(K==closed_kernel(d),'small literal kernel identity');need(K==[r[1:] for r in C.cyclic_kernel(d)[1:]],'small producer kernel API');values=[x for r in K for x in r];need(min(values)>=0 or max(values)<=0,'small mixed kernel');pairchecks+=1;entrychecks+=(n+1)**2
 need(pairchecks==308,'308 pairs')
 data=json.loads((SOURCE/'actual_row_hostile.json').read_text());B=tuple(data['parent_masks']);old=G.add(B,data['old_cut']);new=G.add(B,data['new_cut'])
 need(old==tuple(data['old_module_masks']) and new==tuple(data['new_module_masks']),'same-card construction');need(G.valid(B) and G.strong(old) and G.strong(new),'literal tournament/strongness')
 f=G.profile(old);g=G.profile(new);need(list(f)==data['old_profile'] and list(g)==data['new_profile'],'hostile full profiles');need(f==G.literal_cover(old) and g==G.literal_cover(new),'hostile literal permutation controls');d=[b-a for a,b in zip(f,g)];need(d==[-12,0,240,480,240,0,0]==data['difference'],'actual difference')
 K=literal_kernel(d);need(K==closed_kernel(d),'hostile kernel identity');need(K[0][0]==-36 and K[1][1]==228,'hostile row crossing');contexts=[]
 sidechoices=[(0,),(2,0),(2,4,1)]
 for A,Bside in product(sidechoices,repeat=2):
  before=G.subset_h(G.c3_context(old,A,Bside))[-1];after=G.subset_h(G.c3_context(new,A,Bside))[-1];expect=scalar(K,G.profile(A),G.profile(Bside));need(after-before==expect,'direct actual cyclic context');contexts.append({'side_orders':[len(A),len(Bside)],'old_H':before,'new_H':after,'gap':expect})
 need(contexts[0]['gap']==-36 and contexts[-1]['gap']==67932 and contexts[4]['gap']==data['two_vertex_side_gap']==828,'hostile actual endpoints');need(multiply([1],[2,0,1])==list(G.profile((2,4,1))),'actual C3 multiplier');need(C.certify_actual(old,new,(0,))['negative']['status']=='NOT_CERTIFIED','hostile rejected')
 controls,invalid_count=actual_calls();formal=accepted=chamber=0;coeffs=[(1,),(0,1),(2,0,1),(1,0,0,1)]
 # These are formal algebra controls, not a claim of tournament realizability.
 for n in range(1,5):
  for delta in product((-1,0,1),repeat=n):
   K=closed_kernel(delta)
   for N in range(1,5):
    for v in product((1,2,3),repeat=N):
     truth=reference_negative(K,v);out=C.negative_certificate(delta,v);need((out['status']=='SUFFICIENT_CERTIFICATE')==truth,'formal N API decision');formal+=1
     if truth:
      accepted+=1;ws=[multiply(v,c) for c in coeffs]
      for a,b in product(ws,repeat=2):need(scalar(K,a,b)<0,'formal certified multiplier chamber');chamber+=1
 result={'status':'PASS','independence':'Expected values from independent kernels; producer imported only as API under test','small_census_tournaments':enum['total_tournaments'],'small_distinct_profiles':enum['total_profiles'],'small_profile_pairs':pairchecks,'small_literal_kernel_entries':entrychecks,'actual_hostile':{'module_order':7,'old_masks':old,'new_masks':new,'difference':d,'row2_diagonal':228,'contexts':contexts,'same_card':True,'both_strong':True},'actual_api_controls':controls,'invalid_actual_api_inputs_rejected':invalid_count,'invalid_formal_api_inputs_rejected':3,'formal_api_decisions':formal,'formal_certificates_accepted':accepted,'formal_multiplier_pair_controls':chamber,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'independent_graph_kernel_sha256':hashlib.sha256(KERNEL.read_bytes()).hexdigest(),'independent_census_sha256':hashlib.sha256((D/'profiles_O2.json').read_bytes()).hexdigest(),'producer_certify_sha256':hashlib.sha256((SOURCE/'certify.py').read_bytes()).hexdigest()}
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
