import json,math,argparse,hashlib
from pathlib import Path
from audit_common import need,literal,product,singleton,fold,word,close,count
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=Path(__file__).resolve().parent.parent/'graft');p.add_argument('--sentinel',action='store_true');args=p.parse_args()
if args.sentinel:need(False,'intentional supplement failure')
cs={h:json.loads((args.source/(f'nonunit{h}.json' if h!=45 else 'core45.json')).read_text()) for h in (13,29,45)}
identities=[]
for h,c in cs.items():
 n=len(c['adj']);need(all(not(c['adj'][i]>>i&1) for i in range(n)),'loops')
 need(all(((c['adj'][i]>>j&1)+(c['adj'][j]>>i&1))==1 for i in range(n) for j in range(i)),'tournament')
 need(literal(c['adj'],c['mask'])==tuple(c['core']),'independent permutations')
 c['x']=word(c['X']);c['y']=word(c['Y']);u=int(c['X'][-1]);w=int(c['Y'][0]);c['u']=u;c['w']=w
 a=fold([singleton(u),c['core'],c['y'],singleton(w)]);b=fold([singleton(u),c['x'],c['core'],singleton(w)])
 d=[y-x for x,y in zip(a,b)];need(d==[0,0,0,0,-c['cycle_gain'],0,0],'exact flank identity')
 need(count(close(c['adj'],c['mask']))==2*c['core'][2]+h-c['core'][4],'literal core closure')
 identities.append({'h':h,'length':len(c['X']),'flanks':[u,w],'T_difference':d[4]})
def family(h,k,j):
 c=cs[h];prefix=[singleton(c['u'])]
 if not c['u']:prefix.insert(0,singleton(1))
 return fold(prefix+[c['x']]*j+[c['core']]+[c['y']]*(k-j)+[singleton(c['w'])])
def H(s):return 2*s[2]+s[0]-s[4]
def fp(s):return s[:4]+s[5:]
for h,c in cs.items():
 for k in (0,1,2,10,100):
  base=family(h,k,0)
  for j in range(k+1):
   x=family(h,k,j);need(fp(x)==fp(base) and x[4]-base[4]==-c['cycle_gain']*j,'all family positions')
# Independent bounded selector by residue scan and shifts.
def coin(a,b,K,L,t):
 j=next(j for j in range(a) if (t-b*j)%a==0)
 while t-b*j>a*K:j+=a
 i=(t-b*j)//a;need(0<=i<=K and 0<=j<=L and a*i+b*j==t,'coin selector');return i,j
pairs=[]
for h,K,a,Ls,overhead in [(13,12,225,(224,250,300),167),(29,28,135,(134,150,200),808)]:
 left=[family(h,K,i) for i in range(K+1)]
 for L in Ls:
  right=[family(45,L,j) for j in range(L+1)];base=H(product(left[0],right[L]));C=(a-1)*(h-1);M=a*K+h*L
  for i in range(K+1):
   for u in range(L+1):
    v=product(left[i],right[L-u]);need(H(v)==base+2*(a*i+h*u),'pair rectangle');need(fp(v)==fp(product(left[0],right[0])),'pair common interface')
  for t in range(C,M-C+1):
   i,u=coin(a,h,K,L,t);need(H(product(left[i],right[L-u]))==base+2*t,'pair all selected targets')
  pairs.append({'h':h,'L':L,'width':M-2*C+1,'vertices':12*L+overhead})
# Direct CRT triple, independent of any nested two-coin construction.
R={}
for i in range(13):
 for j in range(29):
  t=6525*i+1755*j;r=t%377;need(r not in R,'CRT unique residue');R[r]=(t,i,j)
need(len(R)==377 and max(t for t,i,j in R.values())==127440,'CRT size and max')
# Boundary omissions are exact, not merely insufficient selector bounds.
low=127064
need(all((low-1-t)%377 or low-1<t for t,i,j in R.values()),'lower predecessor missing')
left=[family(13,12,i) for i in range(13)];mid=[family(29,28,j) for j in range(29)];inner=[[product(x,y) for y in mid] for x in left]
triple=[]
for L in (337,338,400):
 right=[family(45,L,j) for j in range(L+1)];ref=product(inner[0][0],right[L]);base=H(ref);upper=377*L+376
 for t in range(low,upper+1):
  z,i,j=R[t%377];k=(t-z)//377;need(0<=k<=L and z+377*k==t,'CRT selector bounds')
  x=product(inner[i][j],right[L-k]);need(fp(x)==fp(ref) and H(x)==base+2*t,'CRT every physical target')
 need(all((upper+1-t)%377 or upper+1>t+377*L for t,i,j in R.values()),'upper successor missing')
 need((base+2*low)%2==1,'triple oddness')
 triple.append({'L':L,'lower_offset':low,'upper_offset':upper,'width':upper-low+1,'vertices':12*L+962,'lower_bit_length':(base+2*low).bit_length(),'lower_sha256':hashlib.sha256((base+2*low).to_bytes(((base+2*low).bit_length()+7)//8,'big')).hexdigest()})
# Exact final producer text: nested construction and L1000 CRT endpoint.
inner_states=[]
for r in range(378):
 i,j=coin(145,39,38,149,5472+r)
 inner_states.append(product(family(13,38,i),family(29,149,j)))
for r,x in enumerate(inner_states):
 need(fp(x)==fp(inner_states[0]) and x[4]-inner_states[0][4]==-2*r,'nested entire inner unit fiber')
nested=[]
for L in (44,100,1000):
 base=H(product(inner_states[0],family(45,L,L)));C=44*376;M=45*376+377*L
 for t in (C,C+1,M//2,M-C-1,M-C):
  i,k=coin(45,377,376,L,t)
  x=product(inner_states[i],family(45,L,L-k))
  need(H(x)==base+2*t,'nested selected state endpoints')
 nested.append({'L':L,'width':M-2*C+1,'vertices':12*L+4662})
L=1000;base=H(product(inner[0][0],family(45,L,L)));upper=377*L+376
for z in (low,low+1,(low+upper)//2,upper-1,upper):
 r,i,j=R[z%377];k=(z-r)//377
 need(0<=k<=L,'L1000 selector');x=product(inner[i][j],family(45,L,L-k))
 need(H(x)==base+2*z,'L1000 physical target')
endpoint=base+2*low;digest=hashlib.sha256(endpoint.to_bytes((endpoint.bit_length()+7)//8,'big')).hexdigest()
need(digest=='fe78f2c31a532384fc58a3329368bcdf8bcb3cc22d382b7f8ba3bb9f1fc22ac1','L1000 exact lower endpoint hash')
large_crt={'L':L,'width':upper-low+1,'vertices':12*L+962,'lower_bit_length':endpoint.bit_length(),'lower_sha256':digest}
print(json.dumps({'status':'PASS','finite_identities':identities,'pairs':pairs,'CRT_residues':len(R),'CRT_max_representative':127440,'triple':triple,'nested':nested,'large_crt':large_crt},indent=2))
