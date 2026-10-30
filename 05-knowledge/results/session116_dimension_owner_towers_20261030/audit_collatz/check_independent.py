"""Independent Session 116 audit: no imports from production research code."""
from fractions import Fraction as Q
from collections import Counter
import json

counts=Counter()
def check(ok,label):
    counts[label]+=1
    if not ok: raise RuntimeError(label)

def words(b):
    yield ()
    for a in range(1,b+1):
        for tail in words(b-a): yield (a,)+tail

def affine(w):
    # Direct rational recurrence, independent of integer coefficient routine.
    m,b=Q(1),Q(0)
    for a in w: m,b=3*m/2**a,(3*b+1)/2**a
    return m,b

def multiply(g,h):
    a,b=g; c,d=h
    return a*c,a*d+b

def vp(x):
    if not x: return float('inf')
    n,d=x.numerator,x.denominator
    k=0
    while n%5==0: n//=5; k+=1
    while d%5==0: d//=5; k-=1
    return k

def residue(x,p): return x.numerator*pow(x.denominator,-1,p)%p

def subgroup(g,p):
    a,b=map(lambda x:residue(x,p),g)
    result=set(); h=(1,0)
    while h not in result:
        result.add(h); c,d=h; h=(a*c%p,(a*d+b)%p)
    check(h==(1,0),'orbit closure')
    return result

def parameters(g):
    m,b=g
    t=next(k for k in range(1,5) if residue(m**k,5)==1)
    r=vp(m**t-1) if t!=1 else min(vp(m-1),vp(b))
    return t,r

def depth(g,h):
    m,b=g; n,c=h; t,r=parameters(g)
    gamma=vp((m-1)*c-b*(n-1))
    if t!=1:
        if residue(n**t,5)!=1: return 0
        beta=vp(n**t-1)
        return min(gamma,beta if beta<r else float('inf'))
    alpha=vp(n-1) if vp(m-1)<=vp(b) else vp(c)
    return max(0,min(gamma-r,alpha if alpha<r else float('inf')))

def weight(g,d):
    t,r=parameters(g)
    if d==float('inf'): return Q(4*r+1,4*t)
    if d<=r: return Q(d,t)
    return Q(r,t)+Q(1,4*t)*(1-Q(1,5**(d-r)))

def ds(w,u):
    return [depth(affine(w),affine(w[:i]+u)) for i in range(1,len(w)+1)]
def rate(w,u): return sum((weight(affine(w),d) for d in ds(w,u)),Q(0))
def root(w):
    return next(w[:j] for j in range(1,len(w)+1) if len(w)%j==0 and w[:j]*(len(w)//j)==w)
def in_language(u,z,winners):
    return any(u[:len(s)]==s and len(u[len(s):])%len(z)==0 and u[len(s):]==z*(len(u[len(s):])//len(z)) for s in winners)

controllers=list(words(6))[1:]; entries=list(words(5))
for w in controllers:
    g=affine(w); d=len(w); suffix=[w[j:] for j in range(1,d+1)]
    matrix=[ds(w,s) for s in suffix]
    scores=[rate(w,s) for s in suffix]; best=max(scores)
    for j in range(d):
        check(matrix[j][j]==float('inf'),'infinite diagonal')
        check(rate(w,suffix[j])==rate(suffix[j]+w[:j+1],()),'rotation rate')
        for i in range(d):
            check(matrix[j][i]==matrix[i][j],'matrix symmetry')
            for k in range(d): check(matrix[j][k]>=min(matrix[j][i],matrix[i][k]),'ultrametric')
    z=root(w); zscores=[rate(z,z[j:]) for j in range(1,len(z)+1)]
    winners=[z[j:] for j in range(1,len(z)+1) if zscores[j-1]==max(zscores)]
    for e in range(1,5):
        p=5**e; orbit=subgroup(g,p); t,r=parameters(g)
        check(len(orbit)==t*5**max(0,e-r),'literal order')
        completion_counts=[sum(x>=e for x in row) for row in matrix]
        for u in entries:
            depths=ds(w,u)
            hits=[]
            for i in range(1,d+1):
                # V_i M^n U=1 iff M^n=(U V_i)^-1; enumerate this literal orbit.
                target=multiply(affine(u),affine(w[:i]))
                target_inv=(1/target[0],-target[1]/target[0])
                literal=tuple(residue(x,p) for x in target_inv) in orbit
                check(literal==(depths[i-1]>=e),'literal membership')
                hits.append(literal)
            j=max(range(d),key=lambda j:depths[j])
            check(sum(hits)<=completion_counts[j],'one completion layer domination')
            if e<=depths[j]: check(hits==[x>=e for x in matrix[j]],'same lower-layer hits')
    # Separate direct chronological orbit: start with entry, append each letter,
    # count identity prefixes over a full controller period (no inverse shortcut).
    for e in range(1,4):
        p=5**e; period=len(subgroup(g,p))
        for u in sorted(set(entries[:8]+suffix)):
            a,b=map(lambda x:residue(x,p),affine(u)); literal_count=0
            for n in range(period):
                for letter in w:
                    inverse_two=pow(2**letter,-1,p)
                    a,b=3*a*inverse_two%p,(3*b+1)*inverse_two%p
                    literal_count+=((a,b)==(1,0))
            check(literal_count==sum(x>=e for x in ds(w,u)),'direct chronological period count')
    for u in entries:
        depths=ds(w,u); j=max(range(d),key=lambda j:depths[j]); value=rate(w,u)
        check(value<=scores[j],'exact completion domination')
        if max(depths)!=float('inf'): check(value<scores[j],'finite-depth strictness')
        check((value==best)==in_language(u,z,winners),'all maximizer language')
        check(value==(len(w)//len(z))*rate(z,u),'primitive root scaling')
    if z==w:
        for i in range(d):
            for j in range(d): check((matrix[i][j]==float('inf'))==(i==j),'primitive separation')

for w,empty,best in [((3,1,2),Q(9,16),Q(13,16)),((1,2,1,1,2),Q(1,4),Q(9,20))]:
    check(rate(w,())==empty,'counterexample empty')
    check(max(rate(w,w[j:]) for j in range(1,len(w)+1))==best,'counterexample maximum')
print(json.dumps({'status':'PASS','controllers':len(controllers),'entries_per_controller':len(entries),'moduli':[5,25,125,625],'checks':dict(sorted(counts.items())),'total_checks':sum(counts.values()),'production_imports':False},sort_keys=True,indent=2))
