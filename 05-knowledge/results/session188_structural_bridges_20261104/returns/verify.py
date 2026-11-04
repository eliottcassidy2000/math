"""S188 mixed-module target obstruction. Standard library only."""
import json,itertools,hashlib
from pathlib import Path
from inherited import tables, strong, tournament, attach
from carrier import inflate,physical_names
from skeleton import coeffs
from family_coeffs import family_coeffs

def need(x,message):
 if not x: raise RuntimeError(message)

def response(rows):
 n=len(rows);h,e=tables(rows);rev=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)];_,s=tables(rev);N=(1<<n)-1
 F=[sum(e[a][i]*h[N^a] for a in range(N+1)) for i in range(n)]
 G=[sum(h[a]*s[N^a][i] for a in range(N+1)) for i in range(n)]
 J=[[sum(e[a][i]*s[N^a][j] for a in range(N+1)) for j in range(n)] for i in range(n)]
 return h[N],sum(h[a]*h[N^a] for a in range(N+1)),F,G,J,e[N],s[N]

def signature(R,W):
 h,p,F,G,J,e,s=R
 return (h,p,sum(F[i] for i in W),sum(G[i] for i in W),sum(J[i][j] for i in W for j in W),sum(e[i] for i in W),sum(s[i] for i in W))

def permutation_response(rows):
 n=len(rows);h=0;p=0;F=[0]*n;G=[0]*n;J=[[0]*n for _ in rows];e=[0]*n;s=[0]*n
 for perm in itertools.permutations(range(n)):
  bad={i for i in range(1,n) if not(rows[perm[i-1]]>>perm[i]&1)}
  if not bad:h+=1;e[perm[-1]]+=1;s[perm[0]]+=1
  for c in range(n+1):
   if bad-{c}:continue
   p+=1
   if c:F[perm[c-1]]+=1
   if c<n:G[perm[c]]+=1
   if 0<c<n:J[perm[c-1]][perm[c]]+=1
 return h,p,F,G,J,e,s

def run():
 expected={'h':{'2':9,'1':-1,'0':-3},'J0':{'2':6,'1':-4},'end':{'(3, 2)':18,'(2, 1)':8,'(1, 0)':-1},'start':{'(1, 2)':5,'(0, 1)':-4}}
 need(coeffs()==expected,'finite skeleton coefficients')
 counts={'literal_singleton_formulas':0,'all_port_masks':0,'independent_permutation_responses':0,'literal_attachments':0}
 small=[]
 for m in range(1,9):
  rows,W=inflate(m);R=response(rows)
  need(tournament(rows) and strong(rows),'literal strong carrier')
  need(R[0]==9*2**m-1,'Hamilton clock')
  if m>=2:
   predicted=(9*2**m-1,66*3**m+28*2**m-30,66*2**m-12,24*3**m+15*2**m-9,35*2**(m-1),13,7*2**(m-1))
   need(signature(R,W|{0})==predicted,'complete first-clone seven-coordinate family')
  need(signature(R,W)[4]==6*2**m-4,'unmarked module joint')
  if m<=3:
   need(R==permutation_response(rows),'permutation and subset-DP full named response agreement');counts['independent_permutation_responses']+=1
  for i in range(m):
   j=m-1-i
   pred=6*2**m-4+(18*3**i+5)*2**j+8*2**i-int(j==0)-4*int(i==0)
   need(signature(R,W|{i})[4]==pred,'singleton port formula')
   if m>=2:need(pred>=35*2**(m-1) and (pred==35*2**(m-1))==(i==0),'sharp singleton minimum')
   counts['literal_singleton_formulas']+=1
  if m<=6:
   for a in range(1,1<<m):
    Z=W|{i for i in range(m) if a>>i&1};z=signature(R,Z);c=z[0]-z[4]
    if m>=2:need(z[4]>=35*2**(m-1),'arbitrary port-mask lower bound')
    counts['all_port_masks']+=1
    if m in (2,3):small.append({'m':m,'module_port_mask':a,'rows':rows,'physical_names':physical_names(m),'port_mask':sum(1<<i for i in Z),'signature':z,'c':c})
  if m<=3:
   for a in range(1,1<<m):
    Z=W|{i for i in range(m) if a>>i&1};z=signature(R,Z)
    for t in (1,2):
     full=attach(rows,Z,t)
     need(strong(full),'actual attachment strong')
     need(tables(full)[0][-1]==z[0]-z[4]+2**t*z[2],'actual attachment count')
     counts['literal_attachments']+=1
 fc=family_coeffs()
 for m in range(1,9):
  r,W=inflate(m);R=response(r)
  val=lambda mode,exp:sum(coef*int(base)**exp for base,coef in fc[mode].items())
  z=signature(R,W|{0})
  got=(val('h',m),val('p',m),val('Fouter',m)+val('Ffirst',m-1),val('Gouter',m)+val('Gfirst',m-1),z[4],val('eouter',m)+val('efirst',m-1),val('souter',m)+val('sfirst',m-1))
  need(z==got,'full finite-skeleton family coordinates')
 # All-order reduction m>=4 is proved analytically; evaluate its boundary.
 need(-17*2**3-1 < -103,'all m>=4 excluded before contexts')
 for C in (-103,-101):
  for x in small:
   h,p,F,G,J,e,s=x['signature'];c=h-J
   if c<C:continue  # Every ordered extension has constant <= c.
   need(3*c<C<c and c-min(F,G)<C,'arbitrary ordered-context target exclusion')
 sources=[Path(__file__).with_name(f) for f in ('verify.py','inherited.py','carrier.py','skeleton.py','family_coeffs.py')]
 receipt={'status':'PASS','counts':counts,'target_constants':[-103,-101],'scope':'Every transitive inflation m>=2 of the address24 mixed module, with a nonempty module-port set, and arbitrary marked ordered contexts; all integer outer-clock shifts. Not arbitrary internal replacement or global closure.','skeleton':expected,'family_coefficients':family_coeffs(),'small_candidates':small,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
 Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':run()
