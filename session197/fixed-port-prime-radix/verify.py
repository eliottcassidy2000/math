import itertools,json,math
from graphs import require,run_dp,literal_runs,hp,birth,strong
checked=0
for n in range(1,6):
 edges=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(edges)):
  rows=[0]*n
  for i,(u,v)in enumerate(edges):
   if mask>>i&1:rows[u]|=1<<v
   else:rows[v]|=1<<u
  for r in range(n):
   a=run_dp(rows,r);b=literal_runs(rows,r)
   require(a[1:len(b)]==b[1:] and not any(a[len(b):]),'literal/subset run agreement')
   require(all(a[j]%math.factorial(j-1)==0 for j in range(1,len(a))),'factorial divisibility');checked+=1
examples=[]
for rows,r,ports,A,beta in [([2,12,145,149,99,143,175,19],7,[0,1,4],60,-1),([30,64,10,66,174,335,21,367,95],0,[1,2,3,4],72,1)]:
 a=run_dp(rows,r);require(a[1]==3*A+beta and a[2]==3*A and a[3]==A and not any(a[4:]),'clock run vector')
 seq=[];current=rows[:];ps=ports[:]
 for n in range(5):
  h=hp(current);require(h==A*3**(n+1)+beta,'clock count');require(strong(current),'strong clock');require(ps==ports,'fixed ports');seq.append(h)
  current,ps=birth(current,ps,'R',(1<<len(ps))-1,len(ps))
 examples.append({'rows':rows,'root':r,'ports':ports,'runs':a[:4],'counts':seq,'offset':-2*beta})
controls=0
for p in [3,5,7,11,13]:
 for r in range(1,5):
  for i,d in itertools.product(range(1,15),repeat=2):
   require((i*d**r==p)==((i==p and d==1)or(i==1 and r==1 and d==p)),'prime spectral split');controls+=1
print(json.dumps({'status':'PASS','rooted_literal_controls':checked,'prime_factorization_controls':controls,'physical_examples':examples,'necessary_offsets':[{'prime':p,'residue':p-1,'modulus':2*(p-1),'minimum_ports':p-1,'amplitude_divisor':math.factorial(p-1)}for p in [3,5,7,11]]},sort_keys=True,indent=2))
