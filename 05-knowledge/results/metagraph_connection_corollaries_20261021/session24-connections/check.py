import itertools,json
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
def cyc3():return [[i!=j and (j-i)%3==1 for j in range(3)] for i in range(3)]
def sub(a,b):
 n=len(a);m=len(b)
 return [[(b[i%m][j%m] if i//m==j//m else a[i//m][j//m]) for j in range(n*m)] for i in range(n*m)]
def flip(a):
 n=len(a);need(all(a[i][i+1] for i in range(n-1)),'fixed path valid')
 return [[a[i][j] if abs(i-j)==1 else (i!=j and not a[i][j]) for j in range(n)] for i in range(n)]
def hp(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   for j in range(n):
    if a[i][j] and not(s>>j&1):d[s|1<<j][j]+=d[s][i]
 return sum(d[-1])
def permhp(a):
 return sum(all(a[p[i]][p[i+1]] for i in range(len(a)-1)) for p in itertools.permutations(range(len(a))))
a=cyc3();t=sub(a,a)
cases={'lift':t,'flip_lift':flip(t),'lift_flip':sub(flip(a),a),'flip_both_factors':sub(flip(a),flip(a))}
expected={'lift':3159,'flip_lift':2983,'lift_flip':27,'flip_both_factors':1}
for name,g in cases.items():
 need(hp(g)==expected[name],name+' DP')
 need(permhp(g)==expected[name],name+' permutation')
b=cases['flip_lift']
for block in range(3):
 M=set(range(3*block,3*block+3))
 need(any(len({b[v][u] for u in M})>1 for v in range(9) if v not in M),'original block loses modularity')
print(json.dumps({'status':'PASS','checks':checks,'Hamilton_counts':expected,'scope':'Exact9-vertex noncommutation witness and port-module defect; no full wiring classification.'},indent=2,sort_keys=True))
