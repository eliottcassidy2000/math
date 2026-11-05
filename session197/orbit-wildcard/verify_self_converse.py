from verify import *

def anti(r,tau):
 return all(((r[i]>>j)&1)==((r[tau[j]]>>tau[i])&1) for i in range(len(r)) for j in range(len(r)) if i!=j)
def fixed_counts(r,tau):
 n=len(r);fixed=[i for i in range(n) if tau[i]==i];pairs=[(i,tau[i]) for i in range(n) if i<tau[i]]
 gate(len(fixed)<=1,'anti involution fixed count');k=len(pairs);literal=0;example=None
 for bits in range(1<<k):
  reps=[pair[bits>>i&1] for i,pair in enumerate(pairs)]
  for w in it.permutations(reps):
   full=list(w)+fixed+[tau[i] for i in reversed(w)]
   if all(r[a]>>b&1 for a,b in zip(full,full[1:])):
    literal+=1
    if example is None:example=full
 # Independent pair-mask path DP, with endpoint boundary readout.
 dp=[[0]*n for _ in range(1<<k)]
 for i,pair in enumerate(pairs):
  for v in pair:dp[1<<i][v]=1
 for mask in range(1,1<<k):
  for v in range(n):
   if not dp[mask][v]:continue
   for j,pair in enumerate(pairs):
    if mask>>j&1:continue
    for u in pair:
     if r[v]>>u&1:dp[mask|1<<j][u]+=dp[mask][v]
 count=sum(dp[-1][v] for v in range(n) if r[v]>>(fixed[0] if fixed else tau[v])&1)
 gate(count==literal,'two fixed-word counts')
 gate(example==[tau[v] for v in reversed(example)],'actual fixed witness')
 return count,example
Q=[6,24,26,1,9];C=cyclic(7,[1,2,3]);T=substitute(Q,0,C)
tau11=[0,6,5,4,3,2,1,10,9,8,7]
a=runs(Q,0);gate(a==[0,11,17,4],'complete rooted run inventory')
gate(h_all(T)[-1]==h_backward(T)==56483,'independent 11-vertex H')
gate(anti(T,tau11) and all(tau11[tau11[i]]==i for i in range(11)),'11 anti-involution')
g=[(i+1)%7 for i in range(7)]+list(range(7,11))
gate(all(((T[i]>>j)&1)==((T[g[i]]>>g[j])&1) for i in range(11) for j in range(11)),'rotation')
for v in range(11):
 seen=1<<v
 for _ in range(11):
  for u in range(11):
   if seen>>u&1:seen|=T[u]
 gate(seen==(1<<11)-1,'strong 11-vertex graph')
cert=prime_cert(8069);check_cert(cert);gate(all(8069%d for d in range(2,math.isqrt(8069)+1)),'independent 8069 primality')
F11,w11=fixed_counts(T,tau11)
old=json.loads((P/'witness.json').read_text());R=old['rows'];tau14=[7,13,12,11,10,9,8,0,6,5,4,3,2,1]
gate(anti(R,tau14) and all(tau14[tau14[i]]==i for i in range(14)),'14 anti-involution')
F14,w14=fixed_counts(R,tau14);gate(F14==6111,'14 fixed count')
g=old['rotation'];ginv=[g.index(i) for i in range(14)]
gate([tau14[g[tau14[i]]] for i in range(14)]==ginv,'dihedral conjugation')
gate((old['H']+7*F14)//14==1084789 and (old['H']+7*F14)%14==0,'Burnside count')
witness={'eleven':{'rows':T,'H':56483,'prime_quotient':8069,'prime_certificate':cert,'anti_involution':tau11,'fixed_word':w11,'fixed_count':F11},'fourteen':{'anti_involution':tau14,'fixed_word':w14,'fixed_count':F14,'dihedral_word_orbits':1084789}}
(P/'self-converse-witness.json').write_text(json.dumps(witness,sort_keys=True,indent=2)+'\n')
print(json.dumps({'status':'PASS','eleven_H':56483,'eleven_prime_quotient':8069,'eleven_fixed_words':F11,'fourteen_fixed_words':F14,'fourteen_dihedral_word_orbits':1084789,'halfword_universe':128*math.factorial(7)},sort_keys=True,indent=2))
