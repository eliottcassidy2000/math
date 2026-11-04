"""Exact arbitrary-third-block identities, congruences, and cube controls."""
import itertools,json,random
from collections import defaultdict
from functools import lru_cache
from math import prod
ONE={(0,0,0):1}

def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def opposite(a):
    return tuple(((1 << len(a))-1) ^ row ^ (1 << i) for i,row in enumerate(a))

def tournament(n, bits):
    out = [0]*n
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        u,v = (j,i) if bits >> bit & 1 else (i,j)
        out[u] |= 1 << v
    return tuple(out)

def hamilton_subsets(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1,1 << n):
        for v in range(n):
            count = dp[mask][v]
            if not count:
                continue
            unused = a[v] & ~mask
            while unused:
                bit = unused & -unused
                unused -= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return [1]+[sum(row) for row in dp[1:]]

@lru_cache(None)
def hamilton(a):
    return hamilton_subsets(a)[-1]

@lru_cache(None)
def covers(a):
    n = len(a)
    H = hamilton_subsets(a)
    dp = [[0]*(n+1) for _ in range(1 << n)]
    dp[0][0] = 1
    for mask in range(1,1 << n):
        sub = mask
        while sub:
            rest = mask ^ sub
            for k in range(1,n+1):
                dp[mask][k] += H[sub]*dp[rest][k-1]
            sub = (sub-1)&mask
    return tuple(dp[-1][1:])

def substitute(quotient, blocks):
    offsets = [0]
    for block in blocks:
        offsets.append(offsets[-1]+len(block))
    out = [0]*offsets[-1]
    for i, block in enumerate(blocks):
        for u, row in enumerate(block):
            out[offsets[i]+u] |= row << offsets[i]
            for j in range(len(blocks)):
                if quotient[i] >> j & 1:
                    out[offsets[i]+u] |= ((1 << len(blocks[j]))-1) << offsets[j]
    return tuple(out)

def extend(a, x, c):
    n = len(a)
    out = list(a)+[1 << (n+1),0]
    for v in range(n):
        if x >> v & 1:
            out[v] |= 1 << n
        else:
            out[n] |= 1 << v
        if c >> v & 1:
            out[v] |= 1 << (n+1)
        else:
            out[n+1] |= 1 << v
    return tuple(out)

def add(*polys):
    out = defaultdict(int)
    for poly in polys:
        for powers,value in poly.items():
            out[powers] += value
    return {p:v for p,v in out.items() if v}

def scale(poly, value):
    return {p:v*value for p,v in poly.items() if v*value}

def multiply(*polys):
    out = ONE
    for poly in polys:
        new = defaultdict(int)
        for p,v in out.items():
            for q,w in poly.items():
                new[tuple(a+b for a,b in zip(p,q))] += v*w
        out = {p:v for p,v in new.items() if v}
    return out

def bilinear(left, matrix, right):
    return add(*(scale(matrix[i][j],left[i]*right[j])
                 for i in range(3) for j in range(3)))

def port_numerator(G, d, alpha, gamma):
    e = (1,1,1)
    oa = tuple(1-x for x in alpha)
    oc = tuple(1-x for x in gamma)
    return add(
        multiply(add(d,bilinear(e,G,alpha)),
                 add(d,bilinear(oa,G,gamma)),
                 add(d,bilinear(oc,G,e))),
        multiply(add(d,bilinear(e,G,gamma)),
                 bilinear(oc,G,alpha),
                 add(d,bilinear(oa,G,e))))

def q_at(q, k):
    return q[k] if k < len(q) else 0

def uniform_bits(blocks):
    bits,offset = [],0
    for block in blocks:
        bits.append(((1 << len(block))-1) << offset)
        offset += len(block)
    return bits

def literal_three_queries(blocks):
    T = substitute((2,4,1),blocks)
    U = substitute((4,1,2),blocks)
    masks = uniform_bits(blocks)
    records = []
    for i,j in ((0,1),(1,2),(2,0)):
        a,c = masks[i],masks[j]
        ft,fu = hamilton(extend(T,a,c)),hamilton(extend(U,a,c))
        records.append({"a":a,"c":c,"T":ft,"U":fu,"gap":ft-fu})
    return T,U,records


def coefficient_forms(m):
 a,b,c=({(1,0,0):1},{(0,1,0):1},{(0,0,1):1});p=multiply(a,b,c);d=add(ONE,scale(p,-1))
 G=[[a,multiply(a,b),p],[p,b,multiply(b,c)],[multiply(a,c),p,c]];Gt=list(map(list,zip(*G)))
 D=add(port_numerator(G,d,(1,0,0),(0,1,0)),scale(port_numerator(Gt,d,(1,0,0),(0,1,0)),-1))
 need(all(max(p)-min(p)<=3 for p in D),'numerator run imbalance exceeds three')
 q=(0,)+covers(tournament(m,0));M=[[0]*(m+4) for _ in range(3)]
 for (i,j,l),v in D.items():
  for k in range(m+4):
   I,J,L=i+k,j+k,l+k;w=v*(k+1)*(k+2)//2
   for row,cindex,u,z in ((0,L,I,J),(1,J,I,L),(2,I,J,L)):
    coefficient=w*q_at(q,u)*q_at(q,z)
    if coefficient:
     need(cindex<=m+3,'nonzero run tail outside bound');M[row][cindex]+=coefficient
 return [v[1:] for v in M]


def formula(forms,C):
 q=(0,)+covers(C)
 return [sum(v*q_at(q,k) for k,v in enumerate(row,1)) for row in forms]


def cube(m,C):
 blocks=(tournament(m,0),tournament(m,0),C);active=[i for i,b in enumerate(blocks) if len(b)>1];r=len(active)
 masks=uniform_bits(blocks);offsets=[0,m,2*m];Q=((2,4,1),(4,1,2));states={}
 for theta in range(2):
  for bits in range(1<<r):
   changed={i for j,i in enumerate(active) if bits>>j&1}
   states[theta,bits]=substitute(Q[theta],tuple(opposite(b) if i in changed else b for i,b in enumerate(blocks)))
 root_pairs=((0,1),) if m==2 else ((1,2),(2,0));queries=[(masks[i],masks[j]) for i,j in root_pairs];dec=[]
 for theta in range(2):
  for bit,i in enumerate(active):
   T,U=states[theta,0],states[theta,1<<bit];rep=offsets[i]
   ext=sum(1<<v for v in range(len(T)) if not masks[i]>>v&1 and T[v]>>rep&1);choice=None
   for a,c in itertools.product(range(1<<len(blocks[i])),repeat=2):
    x,y=ext|(a<<offsets[i]),ext|(c<<offsets[i]);ft,fu=hamilton(extend(T,x,y)),hamilton(extend(U,x,y))
    if (ft-fu)%4==2:choice=(x,y);break
   need(choice is not None,'child isolation missing');idx=len(queries);queries.append(choice)
   dec.append({'theta':theta,'child':i,'bit':bit,'query_index':idx,'values':[ft,fu]})
 records=[];roots=[set(),set()];seen=set()
 for(theta,bits),T in states.items():
  values=tuple(hamilton(extend(T,a,c)) for a,c in queries);roots[theta].add(values[:len(root_pairs)])
  need(values not in seen,'fixed cube battery collision');seen.add(values)
  for d in dec:
   if d['theta']==theta:need(values[d['query_index']]==d['values'][(bits>>d['bit'])&1],'other child affects decoder')
  records.append({'theta':theta,'bits':bits,'values':values})
 need(all(len(s)==1 for s in roots),'root depends on child interiors')
 root_sums=[sum(next(iter(s))) for s in roots];need(root_sums[0]!=root_sums[1],'root not separated')
 return {'m':m,'C':C,'states':len(states),'adaptive_bound':len(root_pairs)+r,'fixed_count':len(queries),
         'root_sums':root_sums,'queries':queries,'decoders':dec,'all_fixed_signatures':records}


def main():
 forms={2:coefficient_forms(2),3:coefficient_forms(3)}
 expected={2:[[-8,-14,4,0,0],[56,-46,-98,-24,0],[-58,11,94,24,0]],
           3:[[-44,-270,-120,108,0,0],[1192,1494,-2090,-2484,-432,0],[-950,-1777,1538,2484,432,0]]}
 need(forms==expected,'coefficient forms mismatch')
 rng=random.Random(188117);profiles=[];controls=0
 for n in range(1,10):
  count=(1<<(n*(n-1)//2)) if n<=5 else {6:100,7:50,8:20,9:10}[n]
  codes=range(count) if n<=5 else [rng.randrange(1<<(n*(n-1)//2)) for _ in range(count)]
  for code in codes:
   C=tournament(n,code);q=(0,)+covers(C);gap2=formula(forms[2],C)[0];g3=formula(forms[3],C)
   need(q[1]%2==1,'Redei parity control')
   if n>=2:
    need(q[2]%4==2 and q_at(q,3)%6==0,'cover divisibility control')
    need(q[2]>=2*q[1],'two-cover lower bound control')
    need(gap2%8==4,'uniform mod-eight separator')
    need(g3[1]+g3[2]<0,'TT3 two-query separator')
   else:need(gap2==-8 and g3[1]+g3[2]==242,'singleton signs')
   controls+=1
 # Literal queries on a declared range, including nontransitive arbitrary C.
 tests=[tournament(n,code) for n in range(1,4) for code in range(1<<(n*(n-1)//2))]
 tests += [tournament(4,x) for x in (0,1,4,9,17,29)]
 tests += [tournament(n,rng.randrange(1<<(n*(n-1)//2))) for n in (5,5,6,7)]
 literal=[]
 for m in (2,3):
  for C in tests:
   T,U,queries=literal_three_queries((tournament(m,0),tournament(m,0),C))
   need([r['gap'] for r in queries]==formula(forms[m],C),'literal arbitrary-third identity')
   literal.append({'m':m,'C':C,'scalar':[hamilton(T),hamilton(U)],'queries':queries})
 hostileC=(2,4,1);_,_,hostile=literal_three_queries((tournament(3,0),tournament(3,0),hostileC))
 need([r['gap'] for r in hostile]==[-2472,0,-4284],'one-probe cancellation')
 C5=tuple(sum(1<<((i+j)%5) for j in(1,2)) for i in range(5))
 cubes=[cube(m,C) for m,C in [(2,hostileC),(2,C5),(3,hostileC),(2,(0,)),(3,(0,))]]
 return {'status':'PASS','forms':forms,'profile_controls':controls,'literal_cases':literal,'single_probe_cancellation':hostile,'cubes':cubes,
         'scope':'Two small transitive blocks and arbitrary third tournament; no universal three-block or ambient claim.'}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
