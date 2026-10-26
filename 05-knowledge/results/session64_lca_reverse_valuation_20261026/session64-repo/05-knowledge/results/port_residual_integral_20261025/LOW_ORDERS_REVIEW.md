# Independent review of the order-thirteen endpoint-four certificate

2026-10-25. **ACCEPTED / FINITE-EXACT:** the complete order13 consequence
in [LOW_ORDERS.md](LOW_ORDERS.md). Every directed closed blue marked path
with endpoint outdegrees(8,4) at order13 has three successive directed
black K moves to its exact labelled K target. The proof combines all464
base cells with the16 residual cells covered by the cited structural
repairs. No all-order threshold is lowered by this isolated result.

**ACCEPTED / FINITE-EXACT:** exactly the two order12 and four order14
failed-ledger cells in that document. Each failure is existential
blackness failure for the particular inherited word. It is not a
black-connectivity obstruction.

## Source-universe completeness

The independent test starts with one unrestricted bit for each forward
arc on the literal vertices0,...,n-1. It imposes only the original
path arcs, the closing arc, the full positional anti-reversal equations,
the exact incoming port set I, and the two branch bits T(1,n-2) and
T(2,n-2). It imposes no module, strongness, mixedness, isomorphism,
disjointness or additional arc condition.

A closed endpoint-four source has exactly three incoming interior
ports. Vertex1 is excluded by the path, so I is an arbitrary3-subset
of {2,...,n-2}. Each actual tournament has exactly one such I and one
value for each branch bit. These data therefore partition the source
universe into

    4*binom(n-3,3)

cells. There are 336,480,660 cells at n12,13,14 respectively.
Every one is nonempty; the independent test verifies this including
the deferred cells. The only deferred case is a closed middle with
both y=n-3 and c=n-2 in I. Its third port is h in{2,...,n-4}.
There are n-5 choices of h and two chord bits, hence2(n-5) deferred
cells:14,16,18. The tested base counts are consequently322,464,642.

As a check against accidentally fixed interior arcs, every cell has
affine dimension18,23,28 at n12,13,14 respectively. This also follows
directly. Positional reflection has floor(n/2) fixed unordered pairs,
so it leaves

    (binom(n,2)+floor(n/2))/2

independent arc bits. The path fixes ceil((n-1)/2) reflection orbits;
the closing edge fixes one more. The endpoint row contributes n-3
additional constraints after its path edge is removed, and the two
branch bits contribute two. The remaining dimensions are exactly the
ones above. Thus the order13 statement covers480*2^23=4,026,531,840
literal sources in marked-path coordinates, without listing them.
The count is for labelled sources with the fixed marked path, not
isomorphism classes.

## Independent method and exact consequence

The independent implementation below imports neither low_orders.py
nor either of its inherited audit engines. It uses a signed union-find
over raw forward arcs, rather than Gaussian elimination on reduced
forms. This works because every source and blue-reversal constraint
equates two orientation bits, possibly with a sign or a constant.

For each word at its actual stage, its every path arc is forced true
by the source equations. Its blue condition is then added to a fresh
copy of that source system with all preceding K flips retained.
An inconsistent augmented system proves blackness for every source
completion. A consistent system records only the existence of a blue
completion. Finally, the implementation checks the literal unordered
edge XOR and the complete set of flipped pairs against K_P. Thus
legality, anti-reversal and exact target are separate controls.

Every one of the464 order13 base cells passes all three blackness
tests. The16 residual cells are precisely

    I={h,10,11}, 2<=h<=9, closed middle, either T(2,11).

The reviewed source domains in [PORT_FRONTIER.md section3](PORT_FRONTIER.md)
cover all of these: h2 by
[ports_owners_f1/endpoint.md section2](../ports_owners_f1_20261025/endpoint.md)
from n>=11; h3 by the complementary 1->3 and3->1 repairs; h4 by the
complementary same chord split; and5<=h<=9 by the
[unprotected-port main theorem](../unprotected_ports_cancellation_20261025/ports/README.md)
from n>=12. Their two chord directions are exhaustive; no further
x,c or x,y restriction survives the union of cases. The h3/h4
proofs and their reviews in this packet were read, as were the
unprotected main statement and the triangle-endcap repair.
Consequently the complete order13 claim is justified.

At n12 the failures are exactly I={2,7,9} and{6,7,9}, open middle,
x->c, and the first word of open_penultimate. At n14 they are exactly
I={3,5,10} or{3,7,10}, closed middle and x->c, and
I={5,9,10} or{7,9,10}, closed middle and c->x.
All four involve the third word. No other base cell or step admits
a blue completion. Every word is uniformly directed and every
three-word ledger reaches the specified target even in these cells.

The six explicit hostile rows were read as data and replayed by the
independent literal implementation below. The n12 flags are
(true,false,false); the n14 flags are(false,false,true). These are
actual counterexamples to extending those fixed ledgers, not claims
that a source has no alternate ledger.

The separate second-word budget test was also reconstructed on raw
arcs. The relaxed systems are inconsistent at n12 and n14, for both
branches. At n13 both have exact minimum four incoming interior ports.
A signed component has an independent root choice, so its contribution
is exactly min(number of positive entries, number of negative entries).
Adding the forced entries proves and attains the calculated minimum.
Since the endpoint-four source requires three ports, Q2 cannot be the
missing blackness condition.

## Reproduction

The producer passed in normal and optimized Python with byte-identical
outputs. The following standalone code passed in both modes with
byte-identical outputs. Run it from the repository root. It reads only
the producer JSON hostile rows; all symbolic algebra, source families
and literal replay are implemented here.

```python
from itertools import combinations, product
from collections import Counter
from pathlib import Path
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

class Signs:
    def __init__(self,n):
        self.parent=list(range(n)); self.sign=[0]*n; self.rank=[0]*n
    def copy(self):
        x=Signs(0);x.parent=self.parent[:];x.sign=self.sign[:];x.rank=self.rank[:]
        return x
    def find(self,a):
        if self.parent[a]!=a:
            p,s=self.find(self.parent[a]);self.sign[a]^=s;self.parent[a]=p
        return self.parent[a],self.sign[a]
    def equal(self,a,b,d):
        a,x=self.find(a);b,y=self.find(b);d^=x^y
        if a==b:return d==0
        if self.rank[a]<self.rank[b]:a,b=b,a
        self.parent[b]=a;self.sign[b]=d
        if self.rank[a]==self.rank[b]:self.rank[a]+=1
        return True
    def known(self,a,value):
        r,s=self.find(a);z,t=self.find(0)
        return r==z and (s^t)==value

def edges(q):return {tuple(sorted(e)) for e in zip(q,q[1:])}

def words(n,I,o,c):
    z=n-1;pen=n-2;y=n-3;mid=list(range(2,pen));tail=list(range(3,pen))
    if o:
        if y not in I:
            return 'open_release',[[1]+mid+[pen,z,0],[z]+mid+[pen,1,0],[0,z]+mid+[pen,1]]
        if pen not in I:
            return 'open_penultimate',[[0,1,pen,z]+mid,[1,pen,0,z]+mid,[1]+mid+[pen,0,z]]
        h=next(iter(I-{y,pen}));cut=y-(h==y-1);L=list(range(2,cut));R=list(range(cut,pen))
        return 'open_endrun',[[0,1,pen,z]+L+R,[1,pen,z]+L+[0]+R,[1]+L+[0]+R+[pen,z]]
    if pen in I:
        need(y not in I,'residual must be separated before word choice')
        cut=min({z-v for v in I}-{1});L=list(range(1,cut));R=list(range(cut,pen))
        return 'closed_last_movable',[[0]+L+[z]+R+[pen],L+[z]+R+[0,pen],L+R+[0,pen,z]]
    if not c:
        return 'closed_zero',[[pen]+mid+[1,z,0],[z,pen]+mid+[1,0],[0,z,1]+mid+[pen]]
    if 2 not in I:
        return 'closed_one',[[pen,1,z,0]+mid,[z,pen,1,0]+mid,[0,z,1]+mid+[pen]]
    if y not in I:cut=min({z-v for v in I if z-v>2});name='closed_second_movable'
    elif 3 not in I:cut=y;name='closed_second_fixed'
    else:cut=n-4;name='closed_second_tail4'
    L=list(range(2,cut));R=list(range(cut,pen))
    return name,[[0,1]+L+[z]+R+[pen],[0]+L+[z]+R+[1,pen],[0]+L+R+[1,pen,z]]

totals={};failures=[];cell_dimensions={}
for n in (12,13,14):
    pairs=list(combinations(range(n),2));index={e:i+1 for i,e in enumerate(pairs)}
    all_edges=set(pairs);P=list(range(n));ep=edges(P)
    def literal(a,b,flipped):
        e=tuple(sorted((a,b)));return index[e],int(a>b)^int(e in flipped)
    def force(ds,a,b,value):
        v,s=literal(a,b,set());return ds.equal(v,0,value^s)
    def anti(ds,q,flipped):
        reflection=dict(zip(q,q[::-1]))
        for a,b in pairs:
            v,s=literal(a,b,flipped)
            w,t=literal(reflection[b],reflection[a],flipped)
            if not ds.equal(v,w,s^t):return False
        return True
    count=Counter();dims=Counter()
    for port_tuple in combinations(range(2,n-1),3):
        I=set(port_tuple)
        for o,c in product((False,True),repeat=2):
            ds=Signs(len(pairs)+1)
            need(anti(ds,P,set()),'raw source anti')
            need(all(force(ds,i,i+1,1) for i in range(n-1)),'raw path')
            need(force(ds,n-1,0,1),'raw closed')
            need(all(force(ds,0,v,int(v not in I)) for v in range(1,n-1)),'raw ports')
            need(force(ds,1,n-2,o) and force(ds,2,n-2,c),'raw branch')
            dim=len({ds.find(i)[0] for i in range(len(pairs)+1)})-1
            dims[dim]+=1;count['all']+=1
            if not o and {n-3,n-2}<=I:
                count['residual']+=1
                need(I-{n-3,n-2} <= set(range(2,n-3)),'residual h range')
                continue
            name,qs=words(n,I,o,c);count['base']+=1
            flips=set();path_xor=set();flags=[]
            for k,q in enumerate(qs,1):
                need(sorted(q)==P,('permutation',n,I,name,k))
                for a,b in zip(q,q[1:]):
                    v,s=literal(a,b,flips)
                    need(ds.known(v,1^s),('not uniform directed',n,I,name,k,a,b))
                is_blue_possible=anti(ds.copy(),q,flips)
                flags.append(is_blue_possible)
                if is_blue_possible:failures.append([n,list(port_tuple),o,c,k,name])
                path_xor ^=edges(q);flips ^=all_edges^edges(q)
            need(path_xor==ep and flips==all_edges^ep,'literal XOR target')
            count[name]+=1
    totals[n]=dict(sorted(count.items()));cell_dimensions[n]=dict(sorted(dims.items()))

boundaries=[]
for n in (12,13,14):
    pairs=list(combinations(range(n),2));index={e:i+1 for i,e in enumerate(pairs)}
    all_edges=set(pairs);P=list(range(n));z=n-1;c=n-2;y=n-3
    def literal(a,b,flipped):
        e=tuple(sorted((a,b)));return index[e],int(a>b)^int(e in flipped)
    def force(ds,a,b,value):
        v,s=literal(a,b,set());return ds.equal(v,0,value^s)
    def anti_full(ds,q,flipped):
        r=dict(zip(q,q[::-1]))
        for a,b in pairs:
            v,s=literal(a,b,flipped);w,t=literal(r[b],r[a],flipped)
            if not ds.equal(v,w,s^t):return False
        return True
    for branch in ('open','closed'):
        ds=Signs(len(pairs)+1)
        need(anti_full(ds,P,set()) and all(force(ds,i,i+1,1) for i in range(n-1))
             and force(ds,z,0,1),'relaxed source')
        if branch=='open':
            guards=[(1,c,1),(0,y,0),(0,c,1)]
            q1=[0,1,c,z]+list(range(2,c));q2=[1,c,0,z]+list(range(2,c))
        else:
            guards=[(1,c,0),(2,c,1),(0,2,1),(0,c,1)]
            q1=[c,1,z,0]+list(range(2,c));q2=[z,c,1,0]+list(range(2,c))
        need(all(force(ds,*guard) for guard in guards),'relaxed guards')
        if not anti_full(ds,q2,all_edges^edges(q1)):
            boundaries.append([n,branch,'inconsistent'])
            continue
        zero,zero_sign=ds.find(0);minimum=0;components={}
        for v in range(2,n-1):
            root,sign=ds.find(index[(0,v)])
            if root==zero:minimum+=(sign^zero_sign)==0
            else:components.setdefault(root,[]).append(sign)
        minimum+=sum(min(s.count(0),s.count(1)) for s in components.values())
        need(minimum>=4,'second-word port budget')
        boundaries.append([n,branch,minimum])
need(boundaries==[[12,'open','inconsistent'],[12,'closed','inconsistent'],
                 [13,'open',4],[13,'closed',4],
                 [14,'open','inconsistent'],[14,'closed','inconsistent']],
     'boundary minima')

expected=[
[12,[2,7,9],True,True,1,'open_penultimate'],
[12,[6,7,9],True,True,1,'open_penultimate'],
[14,[3,5,10],False,True,3,'closed_one'],
[14,[3,7,10],False,True,3,'closed_one'],
[14,[5,9,10],False,False,3,'closed_zero'],
[14,[7,9,10],False,False,3,'closed_zero']]
need(failures==expected,('failed ledger cells',failures))
need([totals[n]['base'] for n in (12,13,14)]==[322,464,642],'base counts')
need([totals[n]['residual'] for n in (12,13,14)]==[14,16,18],'residual counts')

# Literal hostile rows are data, not imported executable code.
receipt=json.loads(Path('05-knowledge/results/port_residual_integral_20261025/low_orders.normal.json').read_text())
hostile_replays=[]
for h in receipt['failed_ledger_hostiles']:
    R=h['source_rows'];n=len(R);P=list(range(n));U=set(combinations(range(n),2))
    def arc(R,a,b):return (R[a]>>b)&1
    def K(R,q):
        f=U^edges(q)
        return [sum((arc(R,a,b)^int(tuple(sorted((a,b))) in f))<<b
                    for b in range(n) if b!=a) for a in range(n)]
    def blue(R,q):
        r=dict(zip(q,q[::-1]))
        return all(arc(R,a,b)==arc(R,r[b],r[a]) for a,b in U)
    need(all(not arc(R,a,a) for a in range(n)) and
         all(arc(R,a,b)^arc(R,b,a)==1 for a,b in U),'literal tournament')
    need(blue(R,P) and arc(R,n-1,0) and
         all(arc(R,i,i+1) for i in range(n-1)),'literal blue source')
    need({v for v in range(1,n-1) if arc(R,v,0)}==set(h['ports']),'literal ports')
    need(sum(arc(R,0,v) for v in range(n))==n-5 and
         sum(arc(R,n-1,v) for v in range(n))==4,'literal scores')
    need(arc(R,1,n-2)==h['middle_open'] and arc(R,2,n-2)==h['x_to_c'],'literal guard')
    target=K(R,P);flags=[]
    for q in h['words']:
        need(sorted(q)==P and all(arc(R,a,b) for a,b in zip(q,q[1:])),'literal legality')
        flags.append(blue(R,q));R=K(R,q)
    need(R==target and flags==h['literal_blue_flags'],'literal endpoint/blackness')
    hostile_replays.append([n,h['ports'],flags])
print(json.dumps({'status':'PASS','all_cell_counts':totals,'source_dimensions':cell_dimensions,
                  'failed_cells':failures,'literal_hostiles':hostile_replays,'boundary_minima':boundaries},sort_keys=True))

```

The independent output is PASS. Its source counts are(336,480,660),
base counts(322,464,642), and residual counts(14,16,18). The per-cell
dimensions are(18,23,28). Its failed-cell set and six literal flags
are exactly those listed above; its boundary minima are inconsistent,
inconsistent,4,4,inconsistent,inconsistent in increasing order and
open/closed branch order.

The first failed implication in the nearby orders is absence of a
particular positional anti-reversal. The preserved survivor is legality
and the exact labelled target for every cell. The unused restoration
sidecar is the offending word's actual feasible symmetry; a new repair
could use that symmetry, but the present audit supplies no such repair
and makes no global-connectivity inference.

