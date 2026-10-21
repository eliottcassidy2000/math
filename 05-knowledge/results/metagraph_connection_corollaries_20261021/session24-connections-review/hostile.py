import itertools, json
checks=0
def need(value, label):
 global checks
 checks+=1
 if not value: raise RuntimeError(label)
def cyclic(n,steps):
 return [[i!=j and (j-i)%n in steps for j in range(n)] for i in range(n)]
def inflate(a):
 n=len(a)
 return [[((j-i)%3==1) if i//3==j//3 else a[i//3][j//3] for j in range(3*n)] for i in range(3*n)]
def module(a,S):
 return all(len({a[x][y] for y in S})==1 for x in range(len(a)) if x not in S)
def cycle(a,S):
 return all(sum(a[x][y] for y in S)==1 for x in S)
def flip_path(a):
 n=len(a)
 need(all(a[i][i+1] for i in range(n-1)), 'valid path')
 return [[False if i==j else a[i][j] if abs(i-j)==1 else not a[i][j] for j in range(n)] for i in range(n)]
# Exact canonical cyclic-triangle module controls, including composite quotients.
c3=cyclic(3,{1}); r5=cyclic(5,{1,2}); p7=cyclic(7,{1,2,4})
controls={}
for label,a in [('C3',c3),('R5',r5),('P7',p7),('C3[C3]',inflate(c3)),('T3',[[i<j for j in range(3)] for i in range(3)]),('K1',[[False]]),('T2',[[False,True],[False,False]]),('C3_to_1',[[c3[i][j] if i<3 and j<3 else i<3 and j==3 for j in range(4)] for i in range(4)])]:
 t=inflate(a)
 found={S for S in itertools.combinations(range(len(t)),3) if cycle(t,S) and module(t,S)}
 expected={tuple(range(3*i,3*i+3)) for i in range(len(a))}
 need(found==expected, 'canonical bottom triangles '+label)
 controls[label]=len(found)
# Every pair of consecutive nonempty blocks: uniformity fails iff rs>1.
# Transitive local gadgets ensure the identity vertex order is a Hamilton path.
ports=[]
for r in range(1,5):
 for s in range(1,5):
  n=r+s; a=[[i<j for j in range(n)] for i in range(n)]; b=flip_path(a)
  uniform=len({b[i][j] for i in range(r) for j in range(r,n)})==1
  need(uniform==(r*s==1),'port threshold')
  need(module(b,set(range(r)))==(r==1),'left module loss')
  need(module(b,set(range(r,n)))==(s==1),'right module loss')
  ports.append([r,s,uniform])
# Independently compute all automorphisms, Hamilton paths and converse orbits
# on C3[C3], establishing 39 ordinary orbits, 11 fixed, 14 pairs.
t=inflate(c3); n=len(t)
autos=[]; paths=[]
for p in itertools.permutations(range(n)):
 if all(t[p[i]][p[i+1]] for i in range(n-1)): paths.append(p)
 if all(t[i][j]==t[p[i]][p[j]] for i in range(n) for j in range(i+1,n)): autos.append(p)
need(len(autos)==81,'full wreath order')
need(all(all(len({g[x]//3 for x in range(3*i,3*i+3)})==1 for i in range(3)) for g in autos),'actual group preserves blocks')
def canon(p): return min(tuple(g[x] for x in p) for g in autos)
reps={canon(p) for p in paths}
# Product of quotient reflection and local reflection is an anti-involution.
tau=tuple(3*(-(i//3)%3)+(-(i%3)%3) for i in range(n))
need(all(t[i][j]!=t[tau[i]][tau[j]] for i in range(n) for j in range(i+1,n)), 'anti-involution')
alpha={p:canon(tuple(tau[x] for x in reversed(p))) for p in reps}
need(all(alpha[alpha[p]]==p for p in reps),'converse involution')
fixed=sum(alpha[p]==p for p in reps)
need((len(paths),len(reps),fixed)==(3159,39,11),'independent orbit counts')
# Only one block-contiguous orbit, fixed by converse, lifts C3's unique orbit.
contiguous={p for p in reps if all(len({p[3*i+j]//3 for j in range(3)})==1 for i in range(3))}
need(len(contiguous)==1 and all(alpha[p]==p for p in contiguous),'canonical lift fixed control')
print(json.dumps({'status':'PASS','checks':checks,'canonical_triangle_modules':controls,'port_controls':ports,'C3_C3':{'H':len(paths),'aut':len(autos),'h':len(reps),'b':fixed,'paired_degree':(len(reps)-fixed)//2}},sort_keys=True,indent=2))
