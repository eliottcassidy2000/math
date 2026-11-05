from verify import *

def universe(n,force_path):
 edges=list(it.combinations(range(n),2));active=[e for e in edges if not(force_path and e[1]==e[0]+1)]
 orbit={};groups=[]
 for i,j in active:
  if (i,j) in orbit:continue
  reflected=(n-1-j,n-1-i);group={(i,j),reflected}
  for x in group:orbit[x]=len(groups)
  groups.append(group)
 expected=(n-1)**2//4 if force_path else n*n//4
 gate(len(groups)==expected,'complete reflected pair universe')
 for mask in range(1<<len(groups)):
  R=[0]*n
  for i,j in edges:
   win=True if force_path and j==i+1 else bool(mask>>orbit[i,j]&1)
   if win:R[i]|=1<<j
   else:R[j]|=1<<i
  yield R

def literal(R):return sum(all(R[a]>>b&1 for a,b in zip(w,w[1:])) for w in it.permutations(range(len(R))))
def is_strong(R):
 n=len(R)
 for v in range(n):
  seen=1<<v
  for _ in range(n):
   for u in range(n):
    if seen>>u&1:seen|=R[u]
  if seen!=(1<<n)-1:return False
 return True
expected={1,3,5,9,11,13,15,17,25,27,29,31,33,37,41}
witness={};rows=[];total_normal=total_full=0
for n in range(1,8):
 normal=set();full=set();nc=fc=0
 for R in universe(n,True):
  h=h_all(R)[-1];gate(h==literal(R),'independent literal normalized H');normal.add(h);nc+=1
  tau=list(reversed(range(n)))
  gate(all(((R[i]>>j)&1)==((R[tau[j]]>>tau[i])&1) for i in range(n) for j in range(n)),'actual anti-involution')
  if h<=41:witness.setdefault(str(h),{'rows':R,'anti_involution':tau})
 for R in universe(n,False):full.add(h_all(R)[-1]);fc+=1
 gate(normal==full,'two complete universes support')
 rows.append({'n':n,'normalized_presentations':nc,'full_fixed_involution_presentations':fc,'support_through41':sorted(x for x in normal if x<=41)})
 total_normal+=nc;total_full+=fc
gate({int(h) for h in witness}==expected,'all small self-converse support')
gate(5**7>41**3,'strong core cutoff')
closure={d*d*h for d in (1,3,5) for h in expected if d*d*h<=41}
gate(closure==expected,'all-order SCC square closure')
escapes={19:[2,20,57,51,33,3],23:[38,8,58,49,35,2],35:[54,48,122,115,96,64,3],39:[250,252,249,48,64,208,136,24]}
for target,R in escapes.items():
 validate(R);gate(h_all(R)[-1]==literal(R)==target,'unrestricted escape H')
 if target!=39:gate(is_strong(R),'strong escape')
gate(all(target not in expected for target in escapes),'extra omissions')
(P/'converse-spectrum-witness.json').write_text(json.dumps({'self_converse':witness,'unrestricted_escapes':escapes},sort_keys=True,indent=2)+'\n')
print(json.dumps({'status':'PASS','orders':rows,'normalized_graphs_and_literal_checks':total_normal,'full_fixed_involution_graphs':total_full,'all_order_support_through41':sorted(expected),'unrestricted_escapes':list(escapes),'integer_cutoff':{'5^7':5**7,'41^3':41**3}},sort_keys=True,indent=2))
