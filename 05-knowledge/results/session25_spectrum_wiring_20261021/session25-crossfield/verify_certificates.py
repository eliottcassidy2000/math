import json,itertools

def need(x,msg):
    if not x: raise RuntimeError(msg)

def step(t):
    a,b,c=t
    if b==0:
        if c==0:return None
        if c==1:
            return (5,2,a//2+1) if a%2==0 else (2,0,a//2+3)
        return (2+3*(a%2),a//2+3,c-2)
    if c==0:return (a+3,0,b-1)
    return (a+6,b-1,c-1)

def closed(a,c):
    k,e=divmod(a,2)
    if c==0:return None,1
    if c==1:
        if e:return (2,0,k+3),1
        return ((14,0,0) if k==0 else (17,0,k-1)),3
    if c<=k+4:return (6*c-7+3*e,0,k-c+4),c
    return (6*k+20+3*e,0,c-k-5),k+4

checks=0
for a in range(101):
 for c in range(101):
    t=(a,0,c); count=0
    while True:
        t=step(t);count+=1
        if t is None or t[1]==0:break
        need(count<1000,'unexpected nonreturn')
    need((t,count)==closed(a,c),('return mismatch',a,c,t,count))
    checks+=1
trace=[(5,2,0)]
for i in range(20):
 t=step(trace[-1])
 if t is None:break
 trace.append(t)
need(step((5,2,0))==(8,0,1),'initialization')
t=(8,0,1); boundaries=[]
for _ in range(4):
 t,l=closed(t[0],t[2]);boundaries.append({'state':t,'macro_rules':l})

clocks=[]
for j in range(21):
 m=2+12*j;n=(2**m-1)//3
 need(n%105==1 and n%2==1,'guard')
 x=3*n+1;v=0
 while x%2==0:x//=2;v+=1
 need(x==1 and v==m,'clock')
 clocks.append({'j':j,'valuation':v,'raw_steps':v+1})

sharp=[]
for k in range(2,7):
 n=k+1; edges=[(i,j) for i in range(n) for j in range(i) if (i,j)!=(1,0)]
 equal=different=0
 for colours in itertools.product(range(k),repeat=n):
    if any(colours[i]==colours[j] for i,j in edges):continue
    if colours[0]==colours[1]:equal+=1
    else:different+=1
 need(equal>0 and different==0,'sharp skeleton')
 # Two bags of size k, with every edge covered and connected occurrences.
 bags=[set(range(n))-{0},set(range(n))-{1}]
 need(all(any(i in b and j in b for b in bags) for i,j in edges),'bags')
 sharp.append({'k':k,'equal_colourings':equal,'separating_colourings':different,'width_upper':k-1})
# Independent full expanded-graph versus hard-skeleton support tests.
assemblies=0
for k in (3,4):
 for labels in itertools.product(range(3),repeat=3):
    pairs=[(0,1),(0,2),(1,2)];edges=[];hard=[];n=3
    for lab,(u,v) in zip(labels,pairs):
        if lab==1:edges.append((u,v));hard.append((u,v))
        elif lab==2:edges.extend([(u,n),(n,v)]);n+=1
    observed=set()
    for col in itertools.product(range(k),repeat=n):
        if all(col[u]!=col[v] for u,v in edges):observed.add(col[:3])
    expected={col for col in itertools.product(range(k),repeat=3)
              if all(col[u]!=col[v] for u,v in hard)}
    need(observed==expected,('hard skeleton support',k,labels))
    assemblies+=1
print(json.dumps({'status':'PASS','A9_first_returns_checked':checks,'hard_skeleton_assemblies_checked':assemblies,'boundary_prefix':boundaries,'clock_controls':clocks,'sharp_abstract_skeletons':sharp},indent=2))
