from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,math
ROOT=Path(__file__).resolve().parent
counts={}
def check(ok,category):
    counts[category]=counts.get(category,0)+1
    if not ok: raise RuntimeError(category)
I=(F(1),F(0),F(0))
def join(a,b):
    l,d,r=a; m,e,s=b
    return (m*l,m*d+e,max(r,(d+s)/l))
def fold(w):
    s=I
    for l,d in w: s=join(s,(l,d,max(F(0),d/l)))
    return s
def independent(w):
    l,d,r=F(1),F(0),F(0)
    for a,b in w:
        l,d=a*l,a*d+b
        r=max(r,d/l)
    return l,d,r
def states(w,x):
    out=[x]
    for l,d in w:
        x=l*x-d;out.append(x)
    return out
def ceil(q):return -((-q.numerator)//q.denominator)
for row in json.loads((ROOT/'SOURCES.json').read_text())['sources']:
    check(hashlib.sha256((ROOT/row['snapshot']).read_bytes()).hexdigest()==row['sha256'],'source_hash')
alphabet=list(product([F(1,2),F(1),F(2)],[F(-1),F(0),F(1)]))
for k in range(5):
    for w in product(alphabet,repeat=k):
        s=fold(w); l,d,r=s
        check(s==independent(w),'summary')
        check(join(I,s)==s and join(s,I)==s,'identity')
        for cut in range(k+1):check(join(fold(w[:cut]),fold(w[cut:]))==s,'all_cuts')
        for x in [r-F(1,2),r,r+F(1,2),F(ceil(r)-1),F(ceil(r))]:
            ss=states(w,x)
            check((min(ss)>=0)==(x>=r),'direct_feasibility')
            check(ss[-1]==l*x-d,'endpoint')
        check(min(states(w,F(ceil(r))))>=0 and min(states(w,F(ceil(r)-1)))<0,'integer_boundary')
for a,b,c in product(alphabet,repeat=3):
    x,y,z=map(lambda v:fold([v]),(a,b,c))
    check(join(join(x,y),z)==join(x,join(y,z)),'associativity')
for k in range(6):
    for ds in product(range(-2,3),repeat=k):
        p=0; m=0
        for d in ds:p+=d;m=min(m,p)
        check(fold([(F(1),F(-d)) for d in ds])==(F(1),F(-p),F(-m)),'additive')
# Independent exact valuation-cylinder checks, including rejected coarse lifts.
def v2(x):return (x&-x).bit_length()-1
def literal(n,k):
    vals=[];states_=[n]
    for _ in range(k):
        a=v2(3*n+1);n=(3*n+1)//2**a;vals.append(a);states_.append(n)
    return tuple(vals),states_
for k in range(1,5):
    for w in product(range(1,4),repeat=k):
        K=0;C=0
        for a in w:C=3*C+2**K;K+=a
        P=3**k;M=2**(K+1);r=((2**K-C)*pow(P,-1,M))%M
        for t in range(4):
            n=r+M*t
            actual,ss=literal(n,k)
            check(actual==w,'exact_cylinder')
            check(ss[-1]*2**K==P*n+C,'native_affine')
            # Count actual ordinary C and shortcut T steps separately.
            x=n;ordinary=0;shortcut=0
            for a in w:
                x=3*x+1;ordinary+=1
                for _ in range(a):x//=2;ordinary+=1
                shortcut+=a
            check(x==ss[-1] and ordinary==k+K and shortcut==K,'actual_clocks')
            scalar=[(F(3,2**a),F(1)-F(4,2**a)) for a in w]
            check(states(scalar,F(n-1))==[F(v-1) for v in ss],'native_displacement')
        # The other lift of the coarse cylinder changes the last exact valuation.
        bad=(r+2**K)%M
        actual,_=literal(bad,k)
        check(actual!=w,'missing_final_bit')
# Cone closure, with separate c and signs.
for c in [F(0),F(1),F(3,2)]:
    cone=[(l,c*(l-1)*t) for l in [F(1),F(3,2),F(2),F(5)] for t in [F(0),F(1,2),F(1)]]
    for (l,d),(m,e) in product(cone,repeat=2):
        check(0<=m*d+e<=c*(m*l-1),'cone_closure')
check(fold([(F(1),F(1)),(F(1),F(-1))])==(F(1),F(0),F(1)),'hostile_endpoint')
s=fold([(F(1),F(1))]);check(join(s,s)!=s,'hostile_idempotence')
s2=join(s,s);s3=join(s2,s)
check(s3[2]==3 and join(s2,s2)[2]==4,'hostile_overlap')
check(fold([(F(1),F(2)),(F(1),F(-2))])==(F(1),F(0),F(2)),'hostile_internal_cone')
check(min(states([(F(1,2),F(0))],F(1)))>=0 and states([(F(1,2),F(0))],F(1))[-1].denominator==2,'hostile_integrality')
check(v2(3*1+1)==2 and v2(3*3+1)==1,'hostile_guard')
check(F(3,2)*0+1!=F(3,2)*0,'hostile_omitted_P')
check(sum([-1,1])==0 and min([-1,1])<0,'hostile_vector')
# All-demand service by simple exhaustive owner allocations.
def serves(masks,demand):
    for choices in product(*[(None,)+tuple(m) for m in masks]):
        paid=[choices.count(i) for i in range(3)]
        if all(paid[i]>=demand[i] for i in range(3)):return True
    return False
A=[{0,1},{2}];B=[{0,2},{1}]
check(not serves(A,[1,1,0]) and serves(B,[1,1,0]),'hostile_shared_port')
check(all(sum(i in m for m in A)==sum(i in m for m in B) for i in range(3)),'hostile_equal_singletons')
# Concrete disjoint data structures; all query intervals compared to direct fold.
word=[alphabet[(i*i+3*i+2)%len(alphabet)] for i in range(19)]
leaves=[fold([a]) for a in word];n=len(word)
size=1
while size<n:size*=2
tree=[I]*(2*size)
for i,v in enumerate(leaves):tree[size+i]=v
for i in range(size-1,0,-1):tree[i]=join(tree[2*i],tree[2*i+1])
def segment(l,r):
    left=right=I;l+=size;r+=size
    while l<r:
        if l&1:left=join(left,tree[l]);l+=1
        if r&1:r-=1;right=join(tree[r],right)
        l//=2;r//=2
    return join(left,right)
disjoint=[]
for b in range((n-1).bit_length()):
    arr=[I]*n;half=2**b
    for start in range(0,n,2*half):
        mid=min(start+half,n);end=min(start+2*half,n)
        acc=I
        for i in range(mid-1,start-1,-1):acc=join(leaves[i],acc);arr[i]=acc
        acc=I
        for i in range(mid,end):acc=join(acc,leaves[i]);arr[i]=acc
    disjoint.append(arr)
def sparse(l,r):
    if l==r:return I
    if l+1==r:return leaves[l]
    b=(l^(r-1)).bit_length()-1
    return join(disjoint[b][l],disjoint[b][r-1])
lift=[leaves]
for b in range(1,n.bit_length()):
    h=2**(b-1)
    lift.append([join(lift[b-1][i],lift[b-1][i+h]) for i in range(n-2*h+1)])
def lifted(l,r):
    acc=I
    while l<r:
        b=(r-l).bit_length()-1;acc=join(acc,lift[b][l]);l+=2**b
    return acc
for l in range(n+1):
    for r in range(l,n+1):
        expected=fold(word[l:r])
        check(segment(l,r)==expected,'segment_tree_ranges')
        check(sparse(l,r)==expected,'disjoint_sparse_ranges')
        check(lifted(l,r)==expected,'binary_lift_ranges')
print(json.dumps({'status':'PASS','checks':sum(counts.values()),'counts':counts,'scope':'Finite exact rational controls; not a Collatz census or coverage theorem.'},sort_keys=True,indent=2))
