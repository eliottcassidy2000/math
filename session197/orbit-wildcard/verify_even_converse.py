from verify import *

def universe(n):
 edges=list(it.combinations(range(n),2));orbit={};groups=[]
 for i,j in edges:
  if j==i+1 or (i,j) in orbit:continue
  group={(i,j),(n-1-j,n-1-i)}
  for e in group:orbit[e]=len(groups)
  groups.append(group)
 gate(len(groups)==(n-1)**2//4,'normalized orbit count')
 for bits in range(1<<len(groups)):
  R=[0]*n
  for i,j in edges:
   win=j==i+1 or bool(bits>>orbit[i,j]&1)
   if win:R[i]|=1<<j
   else:R[j]|=1<<i
  yield R

def literal(R):return sum(all(R[a]>>b&1 for a,b in zip(w,w[1:])) for w in it.permutations(range(len(R))))
def cyc_key(w):
 i=w.index(min(w));return tuple(w[i:]+w[:i])
expected={1,3,5,9,11,13,15,17,25,27,29,31,33,37,41,45,49,51,57,61,65,67,69,71,73}
witness={};counts=[];total=even=cycle_checks=0
for n in range(1,9):
 vals=set();count=0;tau=list(reversed(range(n)))
 for R in universe(n):
  h=h_all(R)[-1];gate(h==h_backward(R),'independent normalized HP');count+=1;total+=1
  if n%2==0:gate(h%4==1,'even self-converse congruence');even+=1
  if n<=6:
   cs=cycles(R);fixed=sum(cyc_key([tau[x] for x in reversed(w)])==w for _,_,w in cs)
   gate(h%4==(1+2*fixed)%4,'invariant-cycle congruence');cycle_checks+=1
  if h<=73:vals.add(h);witness.setdefault(str(h),{'rows':R,'anti_involution':tau})
 counts.append({'n':n,'presentations':count,'support_through73':sorted(vals)})
gate({int(v) for v in witness}==expected,'complete bounded support')
gate(5**8>73**3,'central strong order cutoff')
gate({d*d*h for d in (1,3,5,7) for h in expected if d*d*h<=73}==expected,'square closure')
words={43:(3,1,2,0),47:(1,1,1,2,2),53:(1,0,2,1,0),55:(1,1,2,3,0),59:(1,0,2,0,0),63:(1,0,0,1,2,2)}
rows={43:[14,20,56,50,33,3],47:[30,96,122,114,98,65,1],53:[38,16,122,115,97,66,3],55:[46,48,122,114,97,64,3],59:[6,16,122,115,97,67,3],63:[38,192,250,243,227,194,129,1]}
for target,w in words.items():
 R=[0]*(len(w)+2);R[0]=2
 for i,c in enumerate(w):
  v=i+2
  for a in (0,1):
   if c>>a&1:R[a]|=1<<v
   else:R[v]|=1<<a
  for j in range(i):R[j+2]|=1<<v
 gate(R==rows[target],'physical escape word')
 gate(h_all(R)[-1]==literal(R)==target,'independent unrestricted escape')
 gate(target not in expected,'extra omission')
gate(h_all([2,4,1])[-1]%4==3,'odd-order hostile')
Q=[2,8,3,5];tau=[2,3,0,1]
gate(h_all(Q)[-1]==5 and all(((Q[i]>>j)&1)==((Q[tau[j]]>>tau[i])&1) for i in range(4) for j in range(4)),'mod8 hostile')
(P/'even-converse-witness.json').write_text(json.dumps({'self_converse':witness,'new_unrestricted_escapes':{str(h):{'rows':rows[h],'word':words[h]} for h in words}},sort_keys=True,indent=2)+'\n')
print(json.dumps({'status':'PASS','orders':counts,'total_graphs_two_independent_DP':total,'even_order_mod4_checks':even,'invariant_cycle_checks':cycle_checks,'all_order_support_through73':sorted(expected),'new_unrestricted_escapes':list(words),'cutoff':{'5^8':5**8,'73^3':73**3}},sort_keys=True,indent=2))
