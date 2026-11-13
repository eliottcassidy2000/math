# (score sequence, arc-usage multiset) as a joint observer on merged classes, orders 5-9, and the residual collisions.
# inputs in the working directory: cls_n.txt (classgen), arcms_n.txt (kps_s29 scripts/arcms.c output)
import collections
def decode(code,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1: out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
for n in range(5,10):
    cc={}; aut={}; scf={}
    for line in open('cls_%d.txt'%n):
        p=line.split('|')[0].split(); cc[int(p[0])]=int(p[6]); aut[int(p[0])]=int(p[5]); scf[int(p[0])]=int(p[4])
    M={}
    for line in open('arcms_%d.txt'%n):
        left,right=line.split('|'); a=left.split(); M[int(a[0])]=(int(a[2]),tuple(map(int,right.split())))
    reps=sorted({min(c,v) for c,v in cc.items()})
    def scores(r):   # converse-invariant: the unordered pair {scores(T), scores(T^op)}
        out=decode(r,n); s=tuple(sorted(bin(x).count('1') for x in out)); so=tuple(sorted(n-1-x for x in s))
        return tuple(sorted([s,so]))
    G=collections.defaultdict(list)
    for r in reps: G[(scores(r),M[r])].append(r)
    coll=[g for g in G.values() if len(g)>1]
    print('n=%d merged %d, distinct (scores, arc multiset) %d, collision groups %d'%(n,len(reps),len(G),len(coll)))
    for g in coll: print('   ',[(r,'SC' if scf[r] else 'NSC','aut',aut[r]) for r in g])
