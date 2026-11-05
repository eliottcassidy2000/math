import importlib.util,json,itertools
from pathlib import Path
spec=importlib.util.spec_from_file_location('inherited',Path(__file__).parent/'dependencies/source_cut_verifier.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
checks={}
def need(ok,label):
 if not ok:raise RuntimeError(label)
 checks[label]=checks.get(label,0)+1
F=(1,-1,1,0);G=(1,0,2,0);tau=(1,0,2,1)
def letter(a,k):return ''.join('Z' if a>>i&1 else 'O' for i in range(k))
def values(s,row,k):return [v.dot(row,v.word(s,letter(a,k))) for a in range(1<<k)]
def unit_components(vals,n):
 seen=bytearray(len(vals));largest=[];maxdegree=0
 for a in range(len(vals)):
  if seen[a]:continue
  comp=[a];seen[a]=1
  for x in comp:
   degree=0
   for i in range(n):
    y=x^(1<<i)
    if abs(vals[x]-vals[y])==2:
     degree+=1
     if not seen[y]:seen[y]=1;comp.append(y)
   maxdegree=max(maxdegree,degree)
  if len({vals[x] for x in comp})>len({vals[x] for x in largest}):largest=comp
 return maxdegree,largest

def unit_word(U):
 n=U;rev=[]
 while n!=2:
  if n%2==0:n//=2;rev.append('O')
  else:n=(n+1)//2;rev.append('Z')
 return ''.join(reversed(rev))
def count_graph(w):
 T,c=v.graphword(w);H,_,_=v.inventories(T)
 return T,c,H[-1]

def main():
 rows=[F,G,v.readout(v.P['N'],tau),v.readout(v.P['J1'],v.word(tau,'J'))]
 matching=[]
 for s in [tau,v.word(tau,'J'),v.word(tau,'OZJ')]:
  for row in rows:
   for k in range(5,11):
    vals=values(s,row,k);deg,comp=unit_components(vals,k)
    need(deg<=1,'single_block_matching');matching.append([list(s),list(row),k,deg])
 exchange=[]
 for k in range(2,14):
  vals=values(tau,F,k);ds=[]
  for a in range(1<<k):
   for i in range(k-1):
    if not(a>>i&1) and a>>(i+1)&1:ds.append(vals[a]-vals[a^(3<<i)])
  expect=2 if k==2 else 4 if k==3 else 2+2**((k-2)//2)+2**((k-1)//2)
  need(min(ds)==expect,'sharp_exchange');exchange.append([k,min(ds)])
 vals=values(tau,F,4);deg,comp=unit_components(vals,4)
 # This particular witness is the sharp threshold control, not a claim about all components.
 ws=['ZOOO','ZZOO','ZZZO'];cts=[v.dot(F,v.word(tau,w)) for w in ws]
 need(cts==[27,25,27],'matching_threshold_hostile')
 # Linear basis verifies all-digit toggle and universal scalar identities independently of physicality.
 basis=[tuple(int(i==j) for i in range(4)) for j in range(4)]
 for k in range(1,8):
  B=1<<k
  for a in range(B):
   w=letter(a,k);p=a&1
   for s in basis:
    D,E,U,d=s;lo=v.word(s,'O'+w);hi=v.word(s,'Z'+w)
    diff=tuple(x-y for x,y in zip(hi,lo))
    need(diff==(U-E-(a+p)*d,0,-B*d,0),'toggle_identity')
 for s in basis:
  lo=v.dot(G,v.word(s,'OJJOOO'));hi=v.dot(G,v.word(s,'OJJOZO'))
  need(lo==hi==33*s[0]+68*s[2],'universal_rebase_identity')
 controls=[]
 for delta in range(-12,13,2):
  B=2
  while B+1+delta<3:B*=2
  m=B.bit_length()-1;U=B+1+delta;pref=unit_word(U)
  s=v.word(tau,pref);need(s==(v.g(U),1,U,1),'unit_prefix')
  lo=v.dot(F,v.word(s,'O'+'O'*m));hi=v.dot(F,v.word(s,'Z'+'O'*m))
  need(lo==v.g(U)+2*B*U and hi-lo==delta,'all_even_offsets')
  controls.append({'delta':delta,'B':B,'U':U,'counts':[lo,hi],'prefix':pref})
 for delta,B in [(-2,4),(0,2),(2,2),(2,4)]:
  m=B.bit_length()-1;U=B+1+delta;pref=unit_word(U)
  w0=pref+'O'+'O'*m+'J';w1=pref+'Z'+'O'*m+'J'
  T0,c0,H0=count_graph(w0);T1,c1,H1=count_graph(w1)
  need(len(T0)==len(T1) and sum((a^b).bit_count() for a,b in zip(T0,T1))==2,'literal_one_edge')
  need(v.cycles(T0)==v.cycles(T1)==1 and H1-H0==delta,'literal_c1_offsets')
 for k in range(1,7):
  B=1<<k
  for a in range(B):
   U=a+(a&1)+2*B+1;s=v.word(tau,unit_word(U));w=letter(a,k)
   lo=v.word(s,'O'+w);hi=v.word(s,'Z'+w)
   need(tuple(x-y for x,y in zip(hi,lo))==(2*B,0,-B,0),'parametric_rebase')
   need(v.dot(G,hi)==v.dot(G,lo),'balanced_rebase')
 two=[]
 for k,l in [(i,i) for i in range(1,10)]+[(3,10),(10,3),(5,10),(10,5),(6,10),(10,6),(4,12),(12,4),(8,10),(10,8)]:
  vals=[]
  for a in range(1<<k):
   s=v.step(v.word(tau,letter(a,k)),'J')
   vals.extend(values(s,G,l))
  deg,comp=unit_components(vals,k+l);cts=sorted({vals[x] for x in comp})
  if k==l==5:need(cts==[6679,6681,6683],'two_block_positive_chain')
  if k==l and 6<=k<=9:need(len(cts)==2,'two_block_matching_observation')
  record={'lengths':[k,l],'max_degree':deg,'max_distinct_values_in_component':len(cts),'example_counts':cts,'example_words':[letter(a>>l,k)+'J'+letter(a&((1<<l)-1),l) for a in comp]}
  two.append(record)
 for k in range(3,10):
  for mask in range(1<<k):
   sig=[-1 if mask>>i&1 else 1 for i in range(k)]
   for i,j,h in itertools.combinations(range(k),3):
    need((sig[i]*sig[j])*(sig[i]*sig[h])*(sig[j]*sig[h])==1,'signed_switching_triangle')
 out={'status':'PASS','checks':checks,'single_block_tests':len(matching),'sharp_exchange_minima':exchange,'threshold_control':{'words':ws,'counts':[27,25,27]},'even_offset_controls':controls,'two_block_probes':two,'scope':'Theorem1 all-order proof; two-block probes finite only; no prescribed-parent or all-target result.'}
 data=json.dumps(out,sort_keys=True,indent=2)+'\n';Path(__file__).with_name('RESULT.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
