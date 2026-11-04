import json,pathlib,itertools,hashlib
P=pathlib.Path(__file__).parent.parent

def need(x):
 if not x:raise RuntimeError('codec audit failed')
def polynomial_reduce(a):
 a=[v%7 for v in a]
 while len(a)>3:
  degree=len(a)-1;c=a.pop();a[degree-1]=(a[degree-1]-c)%7;a[degree-3]=(a[degree-3]-c)%7
 return tuple((a+[0]*3)[:3])
def multiply(a,b):
 p=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):p[i+j]+=x*y
 return polynomial_reduce(p)
def exponent(a,n):
 # Deliberately linear powering, independent of binary powering implementation.
 v=(1,0,0)
 for _ in range(n):v=multiply(v,a)
 return v
def chart(g):
 h,x,y=[v%7 for v in g];z=(h,x,y);a=(x+4*(h-y))%7;b=y;c=(h-y)%7
 return {'observer_mod7':list(z),'F343_vector_codec':[list(z),list(exponent(z,7)),list(exponent(z,49))],'F343_codec_degree':1 if exponent(z,7)==z else 3,'F49_charge':[a,b,c],'F49_frobenius':[(a+4*b)%7,-b%7,c]}
need(all((t**3+t*t+1)%7 for t in range(7)))
fixed343=fixed49=0
for z in itertools.product(range(7),repeat=3):
 need(exponent(z,343)==z);q=chart(z);a,b,c=q['F49_charge'];aa,bb,cc=q['F49_frobenius']
 need(((aa+4*bb)%7,-bb%7,cc)==(a,b,c));need(((b+c)%7,(a-4*c)%7,b)==z)
 fixed343+=exponent(z,7)==z;fixed49+=(a,b,c)==(aa,bb,cc)
need((fixed343,fixed49)==(7,49))
old=json.loads((P/'ports/input_audit.json').read_text())
need(len(old['states'])==4)
for state in old['states']:need(chart(state['golden'])==state['charts'])
receipts=[]
for stem,jointfile,key in [('candidate','candidate_joint.json','name'),('continuation','continuation_joint.json','physical_name')]:
 out=json.loads((P/f'phases/{stem}_phases.json').read_text());joint=json.loads((P/'ports'/jointfile).read_text())
 for s in out['states']:
  name=s[key];fields=joint['fields'][str(name)];e1,_,_,_,s2=fields;need((s2-e1)%2==0)
  g=[joint['H'],(s2-e1)//2,(s2+e1)//2];need(s['golden']==g);need(s['charts']==chart(g))
  receipts.append({'endpoint':stem,'physical_name':name,'golden':g,'degree':s['charts']['F343_codec_degree']})
 if 'input_sha256' in out:
  for file,digest in out['input_sha256'].items():need(hashlib.sha256((P/file).read_bytes()).hexdigest()==digest)
result={'status':'PASS','historical_controls':4,'all_field_elements_checked':343,'cubic_irreducible_mod7':True,'F343_fixed_vectors':fixed343,'F49_charge_fixed_vectors':fixed49,'current_ports':receipts,'scope':'Exact endpoint source binding and independent algebra; joint endpoint DP is audited and replayed separately, no clock or automorphism inference.'}
(P/'older/codec_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
