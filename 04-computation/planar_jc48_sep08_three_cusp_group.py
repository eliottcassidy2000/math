#!/usr/bin/env python3
"""Exact six-word group controls; no geometric path certification."""
import hashlib,json
WORDS={
 'cusp_plus':[-1,-2,1,1,1,2,1],
 'cusp_minus':[2,1,3,3,3,-1,-2],
 'cusp_two':[2,1,1,1,-2],
 'node0':[2,-1,2,2,1,-2],
 'node1':[2,1,2,3,3,-2,-1,-2],
 'node2':[2,3,2,2,-3,-2]}
PREFIX={
 'cusp_plus':([-1,-2],[1,1,1]),
 'cusp_minus':([2,1],[3,3,3]),
 'cusp_two':([2],[1,1,1]),
 'node0':([2,-1],[2,2]),
 'node1':([2,1,2],[3,3]),
 'node2':([2,3],[2,2])}
GATES=0
def check(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
def inverse_word(w):return tuple(-k for k in reversed(w))
def free_mul(a,b):
    out=list(a)
    for k in b:
        if out and out[-1]==-k:out.pop()
        else:out.append(k)
    return tuple(out)
def action(row,word,mul,inv):
    row=list(row)
    for letter in word:
        i=abs(letter)-1;a,b=row[i:i+2]
        row[i:i+2]=[mul(mul(a,b),inv(a)),a] if letter>0 else[b,mul(mul(inv(b),a),b)]
    return tuple(row)
def fprod(*rows):
    out=()
    for row in rows:out=free_mul(out,row)
    return out

base=((1,),(2,),(3,),(4,));aa,bb,cc,dd=base
e=fprod(bb,cc,inverse_word(bb))
pairs={
 'cusp_plus':(bb,cc),
 'cusp_minus':(bb,dd),
 'cusp_two':(aa,e),
 'node0':(fprod(inverse_word(e),aa,e),bb),
 'node1':(aa,dd),
 'node2':(e,fprod(bb,dd,inverse_word(bb)))}
for name,(prefix,core) in PREFIX.items():
    check(WORDS[name]==prefix+core+list(inverse_word(prefix)),'literal chronological inverse '+name)
    row=action(base,prefix,free_mul,inverse_word);i=core[0]-1
    check(row[i:i+2]==pairs[name],'independent relevant free prefix pair '+name)
    check(action(row,inverse_word(prefix),free_mul,inverse_word)==base,'free prefix inverse '+name)

def mm(a,b):
    return(a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
           a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3])
def mi(a):
    check(a[0]*a[3]-a[1]*a[2]==1,'exact determinant one')
    return(a[3],-a[1],-a[2],a[0])
ident=(1,0,0,1);p=(1,1,0,1);q=(1,0,-1,1)
check(mm(mm(p,q),p)==mm(mm(q,p),q),'literal matrix braid relation')
check(mm(p,q)!=mm(q,p),'matrix image nonabelian')
ee=mm(mm(p,q),mi(p))
check(ee==mm(mm(mi(q),p),q),'first exact conjugacy identity')
check(mm(mm(mi(ee),q),ee)==p,'first node collapses to p,p')
row=(q,p,q,q)
for name,w in WORDS.items():check(action(row,w,mm,mi)==row,'whole matrix word fixed '+name)
power=ident
for n in range(1,21):
    power=mm(power,p)
    check(power==(1,n,0,1),'closed unipotent power control')

def pm(a,b):return tuple(a[b[i]] for i in range(len(a)))
def pi(a):return tuple(a.index(i) for i in range(len(a)))
def trans(n,i,j):
    p=list(range(n));p[i-1],p[j-1]=p[j-1],p[i-1];return tuple(p)
records=[]
for n,row,size in [(3,(trans(3,2,3),trans(3,1,2),trans(3,2,3),trans(3,2,3)),6),
                   (4,(trans(4,3,4),trans(4,2,3),trans(4,1,2),trans(4,3,4)),24)]:
    for name,w in WORDS.items():check(action(row,w,pm,pi)==row,'whole permutation word fixed '+name)
    group={tuple(range(n))};todo=list(group)
    while todo:
        g=todo.pop()
        for h in row:
            gh=pm(g,h)
            if gh not in group:group.add(gh);todo.append(gh)
    check(len(group)==size,'named full symmetric image')
    check({g[0] for g in group}==set(range(n)),'named image transitive')
    check(pm(row[0],row[1])!=pm(row[1],row[0]),'named image nonabelian')
    if n==4:check(row[0]!=row[2],'a=c is not forced either')
    if n==4:
        for name in ['cusp_plus','cusp_minus','cusp_two']:
            prefix,core=PREFIX[name];rr=action(row,prefix,pm,pi);i=core[0]-1
            gg={tuple(range(n))};todo=list(gg)
            while todo:
                g=todo.pop()
                for h in rr[i:i+2]:
                    gh=pm(g,h)
                    if gh not in gg:gg.add(gh);todo.append(gh)
            check(len(gg)==6<len(group),'every declared cusp subgroup remains proper')
            check(len([j for j in range(n) if all(g[j]==j for g in gg)])==1,
                  'each cusp subgroup fixes one label')
    records.append([n,len(group),row])

print('Six-word group obstruction PASS; underlying actual paths remain HEURISTIC')
print('Exact consequence: surjection onto B3; no forced cyclicity or equality of all four generators')
print('Controls: full chronological words; infinite SL2Z image; transitive S3 and S4 images')
print('Always-active gates:',GATES)
print('Semantic SHA256:',hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
