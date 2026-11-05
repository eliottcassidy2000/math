"""Independent core projection/inventory compiler audit, without author imports."""
from heapq import heappop, heappush
from itertools import product, permutations
from collections import defaultdict
import json

T=(1,0,2,1)

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def step(s,c):
    D,E,U,d=s
    if c=='O':return (D,0,2*U,d)
    if c=='Z':return (D-E+U,d,2*U-d,d)
    if c=='J':return (D-E+U,D-E,D+U,D)
    raise ValueError(c)

def state(word):
    s=T
    for c in word:s=step(s,c)
    return s

def bank(cap):
    if cap<2:return {}
    out={T:(0,'')};queue=[(2,T)]
    while queue:
        _,s=heappop(queue)
        cost,word=out[s]
        for c in 'OZJ':
            child=step(s,c)
            need(child[2]>s[2],'strict cover increase')
            if child[2]>cap:continue
            candidate=(cost+int(c=='J'),word+c)
            if child not in out:
                out[child]=candidate;heappush(queue,(child[2],child))
            elif candidate<out[child]:out[child]=candidate
    return out

def corebank(R):
    if R<2:return {}
    out={(1,2,1):(0,'')}
    for p,(cost,word) in bank(R-2).items():
        for c in 'ZJ':
            q=step(p,c);D,E,u,d=q
            if D+d<=R:
                triple=(D,u,d);candidate=(cost+int(c=='J'),word+c)
                if triple not in out or candidate<out[triple]:out[triple]=candidate
    return out

def normalize(s):
    D,E,U,d=s
    need((E==0)==(U>=2*D),'endpoint-height dichotomy')
    if E>0:need(D<=U<=2*D-3,'strict core strip')
    if D==1:
        need(E==0 and d==1 and U>=2 and U&(U-1)==0,'transitive seed boundary')
        a=U.bit_length()-2
    else:
        a=0;u=U
        while u>=2*D:
            need(u%2==0,'exact dyadic divisibility')
            u//=2;a+=1
        need(D<=u<2*D,'canonical half-open range')
    return a,(D,U//(1<<a),d)

def generated(cores,N):
    joints={};targets={3:0}
    for (D,u,d),(cost,word) in cores.items():
        a=0
        while D+2*((1<<a)*u)<=N:
            h=(1<<(a+1))*u-d;m=h+D+d
            joints[m,h]=min(cost,joints.get((m,h),10**9))
            k=0
            while ((1<<(k+1))-1)*h+D+d<=N:
                n=((1<<(k+1))-1)*h+D+d
                recipe=word+'O'*(a+1)+'Z'*k
                s=state(recipe)
                need(s[0]-s[1]+s[2]==n,'literal source recipe')
                need(recipe.count('J')==cost,'root cost retention')
                targets[n]=min(cost,targets.get(n,10**9));k+=1
            a+=1
    return joints,targets

N=1001
sources=bank(N-1)
for s,(cost,word) in sources.items():
    a,triple=normalize(s)
    trailing=len(word)-len(word.rstrip('O'))
    need(a==trailing,'canonical suffix count on retained recipe')
    head=state(word[:-a] if a else word)
    need((head[0],head[2],head[3])==triple,'physical normalized core')
    need(head[1]>0 or head==T,'positive endpoint or exceptional seed')
# Check every word through length8, including non-minimal alternative recipes.
word_checks=0
for length in range(9):
    for letters in product('OZJ',repeat=length):
        word=''.join(letters);s=state(word);a,triple=normalize(s)
        need(a==len(word)-len(word.rstrip('O')),'all-word suffix uniqueness')
        word_checks+=1

counts={}
for R in [2,3,4,5,6,8,16,32,64,128]:
    compiled=corebank(R)
    direct={(1,2,1):(0,'')}
    for s,candidate in bank(max(2,2*R-5)).items():
        D,E,U,d=s
        if E>0 and D+d<=R:
            triple=(D,U,d)
            if triple not in direct or candidate<direct[triple]:direct[triple]=candidate
    need(set(compiled)==set(direct),('complete core triples',R))
    need({s:c[0] for s,c in compiled.items()}=={s:c[0] for s,c in direct.items()},('core costs',R))
    need(len(compiled)<=(R**3-R)//6,'cubic upper bound')
    counts[R]=len(compiled)
need(set(corebank(2))=={(1,2,1)},'r2 slice')
need(set(corebank(4))-set(corebank(3))=={(3,3,1)},'r4 slice')
need(set(corebank(6))-set(corebank(5))=={(5,5,1),(5,7,1)},'r6 slice')

# r<m/2<=N/2 is a complete deficit cap for the joint/target universe.
cores=corebank((N-1)//2)
joints,targets=generated(cores,N)
refj={};reft={};balanced={3}
for s,(cost,word) in sources.items():
    D,E,U,d=s;m=D+2*U;h=2*U-d;n=D-E+U
    if m<=N:
        refj[m,h]=min(cost,refj.get((m,h),10**9));balanced.add(m)
    if n<=N:reft[n]=min(cost,reft.get(n,10**9))
need(joints==refj,'complete joint inverse including minimum J costs')
need(targets==reft,'complete C1 normal form including minimum J costs')

canonical_outer_checks=0
for s,(cost,word) in sources.items():
    D,E,U,d=s;n=D-E+U
    if n>N or n in balanced:continue
    value=U-d;k=(value&-value).bit_length()-1
    need(k>=1,'nonbalanced scalar outer exponent')
    h=value//(1<<k);m=n-2*((1<<k)-1)*h
    parent=(m-h-d,(h+d)//2,d)
    need(word.endswith('O'+'Z'*k),'canonical outer physical suffix')
    q=state(word[:-(k+1)])
    need(parent==(q[0],q[2],q[3]),'canonical parent reconstruction')
    normalize(q);canonical_outer_checks+=1

# Scope hostiles: joint pairs need not determine core triples, and scalar
# nonbalancedness is stronger than profile D != U.
x=state('ZJ');need(x==(5,2,6,3) and x[0]-x[1]+x[2]==9,'scalar-boundary hostile')
pairs=[state(w) for w in ['ZZZJJ','ZZJJJ']]
need([(s[0]+2*s[2],2*s[2]-s[3]) for s in pairs]==[(97,53)]*2,'same joint pair')
need([normalize(s)[1] for s in pairs]==[(27,35,17),(29,34,15)],'different core triples')
need((15,13) not in joints and (11,7) not in joints,'nonphysical hostiles')
# Sharp predecessor/output cutoff family, checked beyond the finite bank cap.
for j in range(20):
    p=state('O'*j);q=step(p,'Z');r=q[0]+q[3]
    need(p[2]==r-2 and q[2]==2*r-5,'all-order sharpness family instances')

low_deficit_counts={}
for R in [2,4,6,16,32]:
    _,values=generated(corebank(R),100001)
    bound=1+((R**3-R)//6)*(100001).bit_length()**2
    need(len(values)<=bound,'finite density bound control')
    low_deficit_counts[R]=len(values)
result={
 'status':'PASS','target_cap':N,'forward_cover_cap':N-1,'forward_profiles':len(sources),
 'all_word_checks_through_length8':word_checks,'complete_cumulative_core_counts':counts,
 'compiled_core_deficit_cap':(N-1)//2,'compiled_cores':len(cores),
 'joint_values':len(joints),'C1_targets':len(targets),
 'joint_minimum_J_costs_match':True,'target_minimum_J_costs_match':True,
 'canonical_nonbalanced_scalar_profiles_checked':canonical_outer_checks,
 'low_deficit_target_counts_through100001':low_deficit_counts,
 'same_joint_different_core':{'joint':[97,53],'words':['ZZZJJ','ZZJJJ'],'cores':[(27,35,17),(29,34,15)]},
 'scope':'Canonical projection of a supplied source; finite complete joint/target and core counts; density obstruction only for low-deficit C1 subfamilies, not the whole spectrum.'
}
# Separate exhaustive-permutation graph check for the four published core rays.
def graph(word):
    edges={(0,1)};n=2;root=0
    for c in word:
        for old in range(n):
            edges.add((n,old) if c!='O' and old==root else (old,n))
        if c=='J':root=n
        n+=1
    return n,edges

def path(edges,order):return all((x,y) in edges for x,y in zip(order,order[1:]))
literal=[]
for coreword in ['', 'Z', 'ZZ', 'OZ']:
    q=state(coreword);D,E,u,d=q
    for a,k in product(range(2),repeat=2):
        target=((1<<(k+1))-1)*((1<<(a+1))*u-d)+D+d
        recipe=coreword+'O'*(a+1)+'Z'*k+'J'
        size,edges=graph(recipe)
        hp=sum(path(edges,p) for p in permutations(range(size)))
        hc=sum(path(edges,(0,)+p) and (p[-1],0) in edges for p in permutations(range(1,size)))
        need((hp,hc)==(target,1),'literal core-ray graph')
        literal.append({'word':recipe,'vertices':size,'H':hp,'C':hc})
result['independent_literal_graphs']=literal
print(json.dumps(result,sort_keys=True,indent=2))
