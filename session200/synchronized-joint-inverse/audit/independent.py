"""Independent synchronized-inverse audit; imports no author verifier."""
from collections import deque
from math import isqrt
import json

T=(1,0,2,1)
B=(3,1,3,1)
Z0=(1,0,2,0)


def need(test, message):
    if not test:
        raise RuntimeError(message)


def step(s, c):
    d0,e,u,d=s
    if c=='O': return (d0,0,2*u,d)
    if c=='Z': return (d0-e+u,d,2*u-d,d)
    if c=='J': return (d0-e+u,d0-e,d0+u,d0)
    raise ValueError(c)


def word(s, text):
    for c in text: s=step(s,c)
    return s


def dagger(text):
    return ''.join({'O':'Z','Z':'O','J':'J'}[c] for c in text[::-1])


def bilinear(a,b):
    matrix=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
    return sum(a[i]*matrix[i][j]*b[j] for i in range(4) for j in range(4))


def positive(s):
    d0,e,u,d=s
    return d0-e-d,u-d0+e-d,e,d-e


def from_positive(y):
    a,f,e,p=y
    return a+2*e+p,e,a+f+2*e+2*p,e+p


def physical_inventory(cap):
    recipes={T:''}
    todo=deque([T])
    while todo:
        s=todo.popleft()
        for c in 'OZJ':
            child=step(s,c)
            if child[2]<=cap and child not in recipes:
                recipes[child]=recipes[s]+c
                todo.append(child)
    return recipes


def paired_inventory(cap):
    recipes={(B,Z0):''}
    todo=deque([(B,Z0)])
    while todo:
        w,q=todo.popleft()
        for c in 'OZJ':
            child=(step(w,c),step(q,c))
            if child[0][2]<=cap and child not in recipes:
                recipes[child]=recipes[w,q]+c
                todo.append(child)
    return recipes


# Cone preservation is certified on all four generating rays.
for i in range(4):
    basis=tuple(int(i==j) for j in range(4))
    for letter in 'OZJ':
        need(min(positive(step(from_positive(basis),letter)))>=0,('cone ray',i,letter))
need(positive(Z0)==(1,1,0,0),'auxiliary cone seed')
need(Z0[3]==0,'auxiliary is not a physical source')
# Product inequalities are certified coefficientwise, including d=0.
for i in range(4):
    x=from_positive(tuple(int(i==j) for j in range(4)))
    for j in range(4):
        y=from_positive(tuple(int(j==k) for k in range(4)))
        value=bilinear(x,y)
        need(x[2]*y[2]-value>=0,'upper product cone form')
        need(2*value-x[2]*y[2]-2*x[3]*y[3]>=0,'lower product cone form')

MAX_M=401
CAP=isqrt(4*(MAX_M-1))
physical=physical_inventory(CAP)
paired=paired_inventory(CAP)
reference=physical_inventory((MAX_M-1)//2)
direct={(s[0]+2*s[2],2*s[2]-s[3]) for s in reference if s[0]+2*s[2]<=MAX_M}
accepted=set()
checked=0
for v,alpha in physical.items():
    for (w,q),gamma in paired.items():
        m=bilinear(v,w)
        if m>MAX_M:continue
        h=bilinear(v,q)
        cutoff=isqrt(4*(m-1))
        if max(v[2],w[2])>cutoff:continue
        physical_word=alpha+dagger(gamma)
        s=word(T,physical_word)
        need((s[0]+2*s[2],2*s[2]-s[3])==(m,h),('reconstruction',physical_word))
        need(h<m<2*h and h>=3,'joint cone')
        need(m<=v[2]*w[2]<=2*(m-1),'physical product bound')
        need(h<=v[2]*q[2]<=2*h,'auxiliary product bound')
        need(min(positive(q))>=0,'auxiliary cone')
        accepted.add((m,h));checked+=1
need(accepted==direct,('finite exact inverse',sorted(direct-accepted),sorted(accepted-direct)))

m15=sorted(h for m,h in direct if m==15)
need(m15==[9,11] and (15,13) not in direct,'joint-membership hostile')
need(any(h==13 for m,h in direct),'h=13 separately physical')
# Same physical state, different sidecar.
pairs=[(word(B,s),word(Z0,s)) for s in ('OJOJ','ZZOJ')]
need(pairs[0][0]==pairs[1][0]==(27,9,27,9),'same w')
need(pairs[0][1]==(15,5,15,5) and pairs[1][1]==(23,7,23,7),'different q')
values=[(bilinear(T,w),bilinear(T,q)) for w,q in pairs]
need(values==[(45,25),(45,39)],'collision readout consequences')

result={
 'status':'PASS',
 'scope':'Finite inverse equality for m<=401, exact cone/product algebra, and hostile controls; no all-target claim',
 'maximum_m':MAX_M,'square_root_cap':CAP,
 'physical_inventory_size':len(physical),'paired_inventory_size':len(paired),
 'direct_inventory_cover_cap':(MAX_M-1)//2,'direct_inventory_size':len(reference),
 'joint_pairs':len(direct),'accepted_presentations':checked,
 'm15_h_fiber':m15,
 'same_w_different_q':{'words':['OJOJ','ZZOJ'],'w':pairs[0][0],'q':[p[1] for p in pairs],'readouts_with_t':values},
 'auxiliary_bound':'h <= U_v U_q <= 2h; d_q may be zero',
}
# Independently minimize root changes by topological height order.
def minimum_costs(inventory, decorated=False):
    seed=(B,Z0) if decorated else T
    costs={seed:0}
    height=(lambda x:x[0][2]) if decorated else (lambda x:x[2])
    for state0 in sorted(inventory,key=height):
        need(state0 in costs,'reachable inventory state')
        for c in 'OZJ':
            child=(step(state0[0],c),step(state0[1],c)) if decorated else step(state0,c)
            if child in inventory:
                candidate=costs[state0]+int(c=='J')
                costs[child]=min(candidate,costs.get(child,10**9))
    return costs

full=physical_inventory(MAX_M-1)
cf=minimum_costs(full)
cv=minimum_costs(physical)
cw=minimum_costs(paired,True)
direct_joint_cost={}
direct_target_cost={}
for source,cost in cf.items():
    m,h=source[0]+2*source[2],2*source[2]-source[3]
    n=source[0]-source[1]+source[2]
    if m<=MAX_M:direct_joint_cost[m,h]=min(cost,direct_joint_cost.get((m,h),10**9))
    if n<=MAX_M:direct_target_cost[n]=min(cost,direct_target_cost.get(n,10**9))
actual_joint_cost={}
actual_target_cost={3:0}
for v in physical:
    for w,q in paired:
        m,h=bilinear(v,w),bilinear(v,q)
        height=max(v[2],w[2]);cost=cv[v]+cw[w,q]
        if m<=MAX_M and height<=isqrt(4*(m-1)):
            actual_joint_cost[m,h]=min(cost,actual_joint_cost.get((m,h),10**9))
            actual_target_cost[m]=min(cost,actual_target_cost.get(m,10**9))
        k=1
        while m+2*((1<<k)-1)*h<=MAX_M:
            n=m+2*((1<<k)-1)*h
            cap=isqrt(4*(n//(1<<k)-1))
            if height<=cap and 6*(1<<k)<=n+1:
                actual_target_cost[n]=min(cost,actual_target_cost.get(n,10**9))
            k+=1
need(actual_joint_cost==direct_joint_cost,'joint minimum root-cost equality')
need(actual_target_cost==direct_target_cost,'full Mersenne criterion and minimum root costs')
# Difference-cone extreme rays: D>=d>=E>=0 and U>=d.
def diffcone(s):
    d0,e,u,d=s
    return d0>=d>=e>=0 and u>=d
for ray in [(1,0,0,0),(0,0,1,0),(1,1,1,1),(1,0,1,1)]:
    for c in 'OZJ':need(diffcone(step(ray,c)),'difference cone invariance')
for w,q in paired:need(all(0<=a<=b for a,b in zip(q,w)),'raw sidecar bound')
for k in range(1,13):
    L=2*((1<<k)-1)
    ck=tuple(B[i]+L*Z0[i] for i in range(4))
    need(ck==word(T,'O'*k+'Z'),'literal Mersenne boundary')
    for i in range(4):
        e=tuple(int(i==j) for j in range(4))
        end=word(e,'O'+'Z'*k)
        need(bilinear(ck,e)==end[0]-end[1]+end[2],'boundary scalar identity')
result.update({
 'full_forward_cover_cap':MAX_M-1,
 'full_forward_profiles':len(full),
 'C1_targets':len(direct_target_cost),
 'joint_minimum_J_costs_match':True,
 'C1_minimum_J_costs_match':True,
 'physical_Mersenne_basis_checks':48,
 'difference_cone_generating_ray_checks':12,
 'raw_sidecar_bound_checks':len(paired),
 'scope':'Exact cone and boundary algebra; finite joint/target inverse and minimum root-cost equality through 401; no all-target claim'
})
print(json.dumps(result,sort_keys=True,indent=2))
