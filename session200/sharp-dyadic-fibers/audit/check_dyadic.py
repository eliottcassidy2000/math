import json,collections

def need(q,m):
 if not q:raise RuntimeError(m)
def letter(s,a):
 D,E,U,d=s
 if a=='O':return D,0,2*U,d
 if a=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D
T=(1,0,2,1)
def state(w):
 s=T
 for a in w:s=letter(s,a)
 return s
def K(v,w):
 D,E,U,d=v;A,B,V,b=w
 return D*b+d*A+E*B+U*V-E*b-d*B-U*b-d*V+d*b
prefixes=[T,state('Z'),state('JO'),state('ZJZ')]
records=[]
for r in range(1,5):
 critical=2**r-2;mod=2**r
 for prefix in prefixes:
  for extra in (0,1):
   k=critical+extra;nodes=[prefix]
   for _ in range(k):nodes=[letter(s,a) for s in nodes for a in 'OZ']
   for gate in ('F','G'):
    counts=collections.Counter(((D-E+U) if gate=='F' else (D+2*U))%mod for D,E,U,d in nodes)
    need(set(counts)==set(range(1,mod,2)),('support',r,k,prefix,gate))
    need(all(c%2==(1 if extra==0 else 0) for c in counts.values()),('parity',r,k,prefix,gate))
    records.append({'r':r,'k':k,'gate':gate,'prefix':prefix,'residue_counts':dict(sorted(counts.items()))})
# Exact finite differences recover the advertised quadratic coefficient support.
coefficient_tests=0
for prefix in prefixes:
 for k in range(1,9):
  D,E,U,d=prefix
  for gate in ('F','G'):
   def H(S):
    s=prefix
    for i in range(k):s=letter(s,'Z' if i in S else 'O')
    a,b,c,e=s
    return a-b+c if gate=='F' else a+2*c
   h0=H(set())
   for i in range(k):
    need((H({i})-h0)%2==0,'linear parity')
    for j in range(i+1,k):
     coeff=H({i,j})-H({i})-H({j})+h0
     need(coeff==-d*(2**(j-i-1)+int(j==i+1)),'mixed coefficient')
     need((coeff//2)%2==int(j-i in (1,2)),'odd quadratic support')
     coefficient_tests+=1
for r in range(1,14):
 t=2**(r-1)-1;f0,f1=0,1
 for _ in range(t+1):f0,f1=f1,(f0+f1)%2
 need(f0==1,('critical matching parity',r))
for a in range(1,9):
 h=2**a+1;v=state('Z'*a);w=v
 need(v==(h,1,h,1),'balanced source')
 need(K(v,w)==h*h,'sharp pair')
 q=state('Z'*a+'O'*a+'J');need(q[0]==h*h,'literal sharp word')
 for n in (3,7,21):
  m=h*h-n
  if m>n:need(m%2==0 and K(v,w)==n+m,'threshold')
print(json.dumps({'status':'PASS','critical_and_next_depth_fiber_checks':len(records),'mixed_coefficient_tests':coefficient_tests,'critical_matching_parity_through_r':13,'sharp_source_a_range':[1,8],'critical_edge_r1_checked':True,'fiber_records':records},sort_keys=True,indent=2))
