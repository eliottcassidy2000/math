from pathlib import Path
import collections,itertools,json,math,random
from model import D,P,profile,need

def decode(code,n):
 rows=[0]*n;b=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>b&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
   b+=1
 return rows

def encode(rows):
 code=0;b=0
 for i in range(len(rows)):
  for j in range(i+1,len(rows)):
   code|=((rows[i]>>j)&1)<<b;b+=1
 return code

def canonical(rows):
 n=len(rows);best=1<<(n*(n-1)//2)
 for p in itertools.permutations(range(n)):
  code=0;b=0
  for i in range(n):
   for j in range(i+1,n):code|=((rows[p[i]]>>p[j])&1)<<b;b+=1
  if code<best:best=code
 return best

def symbolic(rows):
 n=len(rows);coeff=[0]*(1<<(2*n))
 for order in itertools.permutations(range(n+2)):
  terms={0:1};valid=True
  for i,j in zip(order,order[1:]):
   if i<n and j<n:
    if not(rows[i]>>j&1):valid=False;break
    continue
   if i>=n and j>=n:
    if i>j:valid=False;break
    continue
   v=i if i<n else j;port=j if j>=n else i;b=1<<(v+n*(port-n))
   if i<n:terms={mask|b:c for mask,c in terms.items()}
   else:
    nxt=terms.copy()
    for mask,c in terms.items():nxt[mask|b]=nxt.get(mask|b,0)-c
    terms=nxt
  if valid:
   for mask,c in terms.items():coeff[mask]+=c
 values=coeff[:]
 for bit in range(2*n):
  for mask in range(len(values)):
   if mask>>bit&1:values[mask]+=values[mask^(1<<bit)]
 return coeff,values

def all_subset_paths(rows):
 n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 incoming=[sum(1<<i for i in range(n) if rows[i]>>j&1) for j in range(n)]
 for s in range(1,N):
  todo=s
  while todo:
   bit=todo&-todo;todo-=bit;i=bit.bit_length()-1;prev=s^bit;value=prev==0;m=prev&incoming[i]
   while m:
    b=m&-m;m-=b;value+=dp[prev][b.bit_length()-1]
   dp[s][i]=value
  h[s]=sum(dp[s])
 return h

def cycle_count(rows):
 n=len(rows);count=0
 for tail in itertools.permutations(range(1,n)):
  order=(0,)+tail
  count+=all(rows[order[i]]>>order[(i+1)%n]&1 for i in range(n))
 return count

def main():
 records=json.loads((P/'five_child_records.json').read_text());need(len(records)==12,'twelve five classes')
 qs={};canonical_codes={};multiplicity=collections.Counter()
 for code in range(1024):
  rows=decode(code,5);qs[code]=profile(rows);key=canonical(rows);canonical_codes[code]=key;multiplicity[key]+=1
 need(set(multiplicity)=={r['code'] for r in records},'exhaustive all1024 five-name sources')
 targets=0;nontransitive_targets=0;changing=set();normalized=collections.Counter();module_pairs=0
 for r in records:
  code=r['code'];rows=decode(code,5);q=qs[code]
  need(r['rows']==rows and r['labeled_count']==multiplicity[code] and r['profile']==q,'source five classification')
  need(r['nontransitive']==(q[0]>1),'nontransitive source status')
  modules=[]
  for mask in range(1,32):
   if mask.bit_count()>1 and all(len({rows[v]>>w&1 for w in range(5) if mask>>w&1})==1 for v in range(5) if not(mask>>v&1)):modules.append(mask)
  need(modules==r['modules'],'all original five modules');need(len(r['targets'])==1<<len(modules),'all original families')
  for family,t in enumerate(r['targets']):
   selected=[m for j,m in enumerate(modules) if family>>j&1];actual=rows[:]
   for i in range(5):
    for j in range(i+1,5):
     if sum((m>>i&1)*(m>>j&1) for m in selected)%2:actual[i]^=1<<j;actual[j]^=1<<i
   target=encode(actual);need(t['family']==selected and t['code']==target and t['canonical_code']==canonical_codes[target] and t['profile']==qs[target],'every source-relative original family target')
   need(all(a<=b for a,b in zip(q,qs[target])),'complete profile monotonicity')
   targets+=multiplicity[code];module_pairs+=1
   if r['nontransitive']:
    nontransitive_targets+=multiplicity[code]
    if q!=qs[target]:changing.add(code)
  for modulus in [8,16]:need(r['normalized_mod'+str(modulus)]==[v*pow(q[0],-1,modulus)%modulus for v in q],'five normalized profiles')
  if r['nontransitive']:normalized[tuple(r['normalized_mod8'])]+=multiplicity[code]
 need(targets==130288 and nontransitive_targets==7408 and changing=={2,5,8,11},'family counts and profile-changing classes')
 need(sum(multiplicity[c] for c in multiplicity if qs[c][0]>1)==904,'904 labeled nontransitive sources')
 source=decode(8,5);target=decode(152,5);need(source==[16,1,3,7,14] and target==[16,5,9,3,14],'overlapping-pair counterexample orientation')
 need(qs[8]==[9,66,186,240,120] and qs[152]==[15,78,198,240,120],'profile invariance hostile control')
 named=json.loads((P/'named_controls.json').read_text());need(len(named)==12,'complete local named domain');cubics=quartics=0
 for r,control in zip(records,named):
  for name,rows in [('source',r['rows']),('converse',[31^v^(1<<i) for i,v in enumerate(r['rows'])])]:
   c,v=symbolic(rows);need(control[name]['rows']==rows and control[name]['coefficients']==c and control[name]['values']==v,'all named values and coefficients')
   need(all(not a for m,a in enumerate(c) if m.bit_count()>4),'complete local quartic degree boundary')
   cubics+=sum(a!=0 and m.bit_count()==3 for m,a in enumerate(c));quartics+=sum(a!=0 and m.bit_count()==4 for m,a in enumerate(c))
 need((cubics,quartics)==(2184,2160),'named cubic/quartic counts')
 # Independent exact q2 identity checks, including even orders and both cyclic
 # exterior orientations. These are tests, not a substitute for the proof.
 identities=[];ambient_controls=[]
 rng=random.Random(290030)
 sources=[(n,code) for n in range(3,6) for code in range(1<<(n*(n-1)//2))]
 sources += [(n,rng.randrange(1<<(n*(n-1)//2))) for n in range(6,10) for _ in range(40)]
 for n,code in sources:
  rows=decode(code,n);h=all_subset_paths(rows);N=(1<<n)-1;q2=sum(h[a]*h[N^a] for a in range(1,N));cyc=cycle_count(rows) if n%2 else 0
  need((q2-2*h[-1]-4+4*cyc)%8==0,'q2 Hamilton-cycle parity identity')
  identities.append([n,code,h[-1],q2,cyc])
 # All small three-block quotients. A cyclic quotient gives3H+q2; a
 # transitive one givesH, with actual exterior vertices individually retained.
 tests=[(n,code) for n in [3,4] for code in range(1<<(n*(n-1)//2))]+[(5,r['code']) for r in records]
 for n,code in tests:
  rows=decode(code,n);h=all_subset_paths(rows);N=(1<<n)-1;q2=sum(h[a]*h[N^a] for a in range(1,N))
  for quotient_code in range(8):
   quotient=decode(quotient_code,3);aug=rows[:]+[0,0]
   for i in range(n):
    for j in range(1,3):
     if quotient[0]>>j&1:aug[i]|=1<<(n+j-1)
     else:aug[n+j-1]|=1<<i
   if quotient[1]>>2&1:aug[n]|=1<<(n+1)
   else:aug[n+1]|=1<<n
   iscycle=all(r.bit_count()==1 for r in quotient);actual=all_subset_paths(aug)[-1];expected=3*h[-1]+q2 if iscycle else h[-1]
   need(actual==expected,'exact two-exterior integer formula')
   ambient_controls.append([n,code,quotient_code,actual])
 (D/'local_identity_controls.json').write_text(json.dumps(dict(mod8_identities=identities,exact_exterior=ambient_controls),separators=(',',':'))+'\n')
 result=dict(status='PASS',five_labeled_sources=1024,five_nontransitive_sources=904,five_classes=12,representative_family_outputs=module_pairs,labeled_family_outputs=targets,nontransitive_labeled_family_outputs=nontransitive_targets,profile_changing_classes=sorted(changing),named_coefficients=24576,nonzero_cubics=cubics,nonzero_quartics=quartics,mod8_identity_controls=len(identities),exact_exterior_controls=len(ambient_controls))
 (D/'local_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
