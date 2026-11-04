"""Direct physical forward-weighted subset recount of every grouped polynomial."""
from pathlib import Path
from functools import lru_cache
import json,hashlib,subprocess,tempfile,math
D=Path(__file__).resolve().parent;P=D.parent/'skew_polynomial'
if not P.exists():P=D.parent/'session189-skew-polynomial'
def need(c,s):
 if not c:raise RuntimeError(s)
def op(R):return [((1<<len(R))-1)^(1<<i)^r for i,r in enumerate(R)]
def inflate(Q,cap):
 offset=[];n=0
 for c in cap:offset.append(n);n+=c
 R=[0]*n;groups=[]
 for i,c in enumerate(cap):
  groups.extend([i]*c)
  for u in range(offset[i],offset[i]+c):R[u]=sum(1<<v for v in range(u+1,offset[i]+c))|sum(1<<v for j,k in enumerate(cap) if Q[i]>>j&1 for v in range(offset[j],offset[j]+k))
 return R,groups

def starts(R,repeated=-1):
 n=len(R);full=(1<<n)-1;other=full^(1<<repeated) if repeated>=0 else full
 @lru_cache(None)
 def suffix(S,left,last):
  if not S and not left:return 1
  total=0
  if left and R[last]>>repeated&1:total+=suffix(S,left-1,repeated)
  choices=S&R[last]
  while choices:
   b=choices&-choices;choices-=b;total+=suffix(S^b,left,b.bit_length()-1)
  return total
 if repeated<0:return [suffix(full^(1<<i),0,i) for i in range(n)]
 return [suffix(other,1,i) if i==repeated else suffix(other^(1<<i),2,i) for i in range(n)]
def main():
 selected=json.loads((P/'selected.json').read_text());weighted=json.loads((P/'selected.out').read_text());runs=json.loads((P/'runword.out').read_text());small=json.loads((P/'small_controls.json').read_text());sw=json.loads((P/'small_weighted.out').read_text());cases=[];meta=[]
 expected=[[192,289,99,71,15,153,306,94,189],[320,193,83,71,43,15,432,61,190],[448,209,195,71,45,15,432,312,62],[384,177,147,135,105,77,15,368,126]]
 for j,q in enumerate(selected):
  Q=q['rows'];need(Q==expected[j],'fixed four physical quotient rows');need(all(not(Q[i]>>i&1) and all(((Q[i]>>k)&1)+((Q[k]>>i)&1)==1 for k in range(9) if i!=k) for i in range(9)),'actual tournament')
  modules=[]
  for M in range(1,511):
   if M.bit_count()<2:continue
   a=(M&-M).bit_length()-1
   if all(not(M>>b&1) or ((Q[a]^Q[b])&(511^M))==0 for b in range(9)):modules.append(M)
  need(not modules and any(r.bit_count()!=4 for r in Q),'prime and nonregular')
  S=starts(Q);E=starts(op(Q));need(S==q['S'] and E==q['E'] and all((a+b)%2==0 for a,b in zip(S,E)),'complete exact endpoint metadata and p zero');sk=[[0]*9 for _ in Q]
  for k in range(9):
   a,b=starts(Q,k),starts(op(Q),k)
   for i in range(9):sk[i][k]=a[i]-b[i]
  sigma=[(S[i]-E[i]+2*sum(sk[i]))%4 for i in range(9)];need(sigma==q['sigma']==[0]*9,'independent full repeated-run sigma')
  witness=[(i,j) for i in range(1,9) for j in range(1,9) if i!=j and Q[0]>>i&1 and Q[i]>>j&1 and Q[j]&1 and all(((Q[u]>>i)&1)==((Q[u]>>j)&1) for u in range(1,9) if u not in [i,j])];need(not witness,'no displayed ordered-witness pair at distinguished root')
  V,groups=inflate(Q,[1]+[2]*8);need(V==q['contracted_rows'] and groups==q['groups'],'all seventeen physical names and groups');cases.append((V,0,9,groups));meta.append(dict(node=q['node'],quotient_rows=Q,contracted_rows=V,groups=groups,S=S,E=E,start_run_skews=sk,sigma=sigma,selected_group=q['certificate_group']))
 for v in small['virtuals']:cases.append((v['rows'],0,max(v['groups'])+1,v['groups']))
 need(math.factorial(16)*math.comb(16,8)<1<<63,'exact endpoint coefficient bound');need(math.factorial(17)*16*(1<<15)<1<<127,'exact signed convolution bound')
 with tempfile.TemporaryDirectory() as tmp:
  binary=Path(tmp)/'weighted';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'weighted_forward.cpp'),'-o',str(binary)],check=True);inp=str(len(cases))+'\n'+''.join(f'{len(R)} {z} {ng}\n'+' '.join(map(str,R))+'\n'+' '.join(map(str,g))+'\n' for R,z,ng,g in cases);raw=subprocess.check_output([str(binary)],input=inp.encode()).decode()
 lines=[list(map(int,l.split())) for l in raw.splitlines()];need(len(lines)==7,'all large and small polynomial contexts');polys=[]
 for i,(line,(R,z,ng,g)) in enumerate(zip(lines,cases)):
  n=len(R);need(len(line)==ng*(n-1),'complete polynomial degree range');kk=[line[(n-1)*j:(n-1)*(j+1)] for j in range(ng)];expected_skews=weighted[i]['skews'] if i<4 else sw[i-4]['skews'];need(kk==expected_skews,'every independently recounted physical grouped coefficient including zeros')
  if i<4:
   chosen=kk[selected[i]['certificate_group']];need(chosen[0]==[175478400,158246148,248375296,50196404][i] and all(c>=0 for c in chosen),'strict sign-definite selected certificate');need(runs[i]['skew']==chosen+[0],'producer separate compressed derivative agrees with fresh physical recount');meta[i]['all_grouped_skew_polynomials']=kk
  polys.append(dict(physical_rows=R,distinguished=z,groups=g,all_grouped_skew_polynomials=kk))
 outputs={'four_quotient_recount.json':meta,'all_polynomial_recount.json':polys}
 for name,obj in outputs.items():
  p=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if p.exists():need(p.read_text()==txt,'independent exact artifacts stable')
  else:p.write_text(txt)
 return dict(status='PASS',fixed_prime_quotients=4,exact_endpoint_coordinates=72,exact_repeated_run_coordinates=324,large_grouped_polynomials=36,large_exact_coefficients=576,small_contexts=3,total_exact_polynomial_coordinates=sum(len(x) for x in lines),physical_contracted_order=17,source_pins={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['polynomial_audit.py','weighted_forward.cpp']},output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs},producer_proof_sha256=hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
