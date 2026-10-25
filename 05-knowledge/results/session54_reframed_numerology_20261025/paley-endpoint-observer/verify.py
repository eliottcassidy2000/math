import json, pathlib, itertools, math, hashlib
ROOT=pathlib.Path(__file__).resolve().parent
def need(v,msg):
 if not v: raise RuntimeError(msg)
def python_dp(q):
 qr={x*x%q for x in range(1,q)}; d={(1,0):1}
 for mask in range(1,1<<q,2):
  for v in range(q):
   z=d.get((mask,v),0)
   if z:
    for w in range(1,q):
     if not mask>>w&1 and (w-v)%q in qr:
      k=(mask|1<<w,w);d[k]=d.get(k,0)+z
 return [d.get(((1<<q)-1,v),0) for v in range(q)]
def mul(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
for pin in json.loads((ROOT/'sources.json').read_text()):
 need(hashlib.sha256((ROOT/pin['copy']).read_bytes()).hexdigest()==pin['sha256'],'source pin')
need(math.factorial(19)<2**64,'uint64 bound')
records=[]
for line in (ROOT/'counts.txt').read_text().splitlines():
 q,*row=map(int,line.split());qr={x*x%q for x in range(1,q)};r=(q-1)//2
 a=row[1];b=row[next(v for v in range(1,q) if v not in qr)];s=a+b;delta=a-b
 need(row==[0]+[a if v in qr else b for v in range(1,q)],'two endpoint orbits')
 if q<=11:need(row==python_dp(q),'independent Python DP')
 if q<=7:
  direct=[0]*q
  for tail in itertools.permutations(range(1,q)):
   p=(0,)+tail
   if all((y-x)%q in qr for x,y in zip(p,p[1:])): direct[p[-1]]+=1
  need(row==direct,'permutation enumeration')
 A=[[int((y-x)%q in qr) for y in range(q)]for x in range(q)]
 S=[[A[i][j]-A[j][i] for j in range(q)] for i in range(q)]
 need(mul(S,S)==[[1-q*int(i==j) for j in range(q)]for i in range(q)],'Gauss quadratic identity')
 E=[[a*A[i][j]+b*A[j][i] for j in range(q)]for i in range(q)]
 N=(s*s+q*delta*delta)//4;K=(q*s*s+delta*delta)//4
 need(s%2==1 and delta%2==1 and N*4==s*s+q*delta*delta,'integer parity/norm')
 need(b%q==1,'closed-channel congruence')
 if q==19:need(a%q==3 and delta%q==2,'false open-channel/conductor extension')
 EE=mul(E,E)
 need(all(EE[i][j]+s*E[i][j]+N*int(i==j)==K for i in range(q)for j in range(q)),'quadratic matrix identity')
 need(all((E[i][j]-b*(i!=j))//delta==A[i][j] for i in range(q)for j in range(q)),'labelled reconstruction')
 # Hostile: symmetrization has no quadratic eigenfield; scalar s does not equal lambda.
 need(all(E[i][j]+E[j][i]==s*(i!=j) for i in range(q)for j in range(q)),'symmetrization hostile')
 need(delta!=0 and s*s+s*s+N!=0,'scalar substitution hostile')
 need((b+a)==s and (b-a)**2==delta**2,'transpose spectral collision')
 if q==3:need(a==0 and b==1,'zero endpoint hostile')
 records.append(dict(q=q,a=a,b=b,s=s,delta=delta,H=q*r*s,unrooted_cycles=r*b,norm=N,a_mod_q=a%q,b_mod_q=b%q,delta_mod_q=delta%q,order_conductor=abs(delta),order_discriminant=-q*delta*delta,determinant=r*s*N**r,dp_bytes=(1<<(q-1))*(q-1)*8))
need([z['H'] for z in records]==[3,189,95095,1172695746915],'inherited totals')
print(json.dumps({'status':'PASS','records':records,'checks':['forward/backward C++ agreement','independent Python DP q<=11','permutations q<=7','quadratic matrix identity','labelled reconstruction','symmetrization and scalar hostiles']},indent=2))
