"""Independent directed-word and finite-order controls; no target imports."""
from pathlib import Path
import hashlib,json,itertools
from literal import profile,labelled,join,converse
P=Path(__file__).resolve().parent
def need(x,m):
    if not x: raise RuntimeError(m)
def state(p):return (p['h'],p['R']-p['d'],p['d']-p['e'],p['d'])
def mat(p):
    h,e,s,R,L,d,b=(p[k] for k in ('h','e','s','R','L','d','b'))
    return ((d,L-d,d-s,h-d),(0,b,0,R-d),(0,0,0,d-e),(0,0,0,d))
def apply(M,v):return tuple(sum(a*b for a,b in zip(row,v)) for row in M)
def mul(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)) for i in range(4))
def F(n):return 1+sum((n>>(i+1))+((n>>(i+1))&1) for i in range(n.bit_length()) if not(n>>i&1))
def birth(T,z):
    n=len(T);out=[x|(1<<n) for x in T]+[0]
    if z:out[0]^=1<<n;out[n]=1
    return out
def prefix(q):
    out=[0]
    if q:
        for bit in bin(q)[2:]:out=birth(out,bit=='0')
    return out
def graph(C,c,q,t):
    A,r=join(prefix(q),0,C,c)
    while t:
        A=birth(A,t&1);t//=2
    return A
def scalar(p,q,t):
    h,e,s,R,L,d,b=(p[k] for k in ('h','e','s','R','L','d','b'))
    return h+d*(F(q)+F(t)-2)+(L-d)*q+(d-s)*(q&1)+b*q*t+(R-d)*t+(d-e)*(t&1)
def main():
    bank={};small=[];pairs=associative=binary=matrix_pairs=normalizations=0
    for n in range(1,5):
        for C in labelled(n):
            for c in range(n):
                p=profile(C,c)
                if n<=3:small.append((C,c,p))
                if n>=2:
                    key=tuple(p[k] for k in ('h','e','s','R','L','d','b'))
                    bank.setdefault(key,(C,c,p))
                    need(all(x>=0 for row in mat(p) for x in row),'positive matrix')
                    need((p['L']>p['d'])==(sum(bool(C[v]>>c&1) for v in range(n))>0),'incoming slope iff')
                    need((p['R']>p['d'])==(C[c]!=0),'outgoing slope iff')
                    for A,a,pa in small if n==4 else []:
                        T,r=join(A,a,C,c)
                        need(state(profile(T,r))==apply(mat(p),state(pa)),'positive retained state against literal words')
                        pairs+=1
    # Small factors independently check representation composition as matrices.
    for A,a,pa in small:
        for B,b,pb in small:
            if len(B)==1:continue
            T,r=join(A,a,B,b)
            need(state(profile(T,r))==apply(mat(pb),state(pa)),'small positive retained state')
            pairs+=1
            if len(A)==1:continue
            need(mat(profile(T,r))==mul(mat(pb),mat(pa)),'exact matrix anti-homomorphism')
            matrix_pairs+=1
    tiny=[A for n in range(1,4) for A in labelled(n)]
    for A in tiny:
        for B in tiny:
            for C in tiny:
                X,_=join(A,0,B,0);X,_=join(X,0,C,0)
                Y,_=join(B,0,C,0);Y,_=join(A,0,Y,0)
                need(X==Y,'literal associativity');associative+=1
    inverse=[];M=200;added=M.bit_length()+1
    for key,(C,c,p) in bank.items():
        for q in range(8):
            for t in range(8):
                T=graph(C,c,q,t); got=profile(T,0)
                need(got['h']==scalar(p,q,t),'physical binary history')
                Opp=graph(converse(C),c,t,q)
                # Map prefix vertices to the reversed final suffix and vice versa.
                left=q.bit_length();right=t.bit_length();n=len(C)
                perm={0:0}
                for i in range(1,left+1):perm[i]=n+right+left-i
                for j in range(n-1):perm[left+1+j]=right+1+j
                for j in range(right):perm[n+left+j]=right-j
                need(all(((T[x]>>y)&1)==((Opp[perm[y]]>>perm[x])&1) for x in range(len(T)) for y in range(len(T)) if x!=y),'actual binary converse arcs')
                binary+=1
        # Exhaust all canonical bit-length pairs within the proved order bound,
        # not the source's rectangular index scan or nested slope bounds.
        values=set()
        for q in range(1<<added):
            for t in range(1<<(added-q.bit_length())):
                val=scalar(p,q,t)
                if val<=M:values.add(val)
        inverse.append(dict(signature=list(key),outputs=len(values),support_sha256=hashlib.sha256(json.dumps(sorted(values),separators=(',',':')).encode()).hexdigest()))
    source=json.loads((P/'source'/'rooted_amalgamation_normal.json').read_text())
    need({tuple(x['signature']):(x['outputs'],x['support_sha256']) for x in inverse}=={tuple(x['signature']):(x['outputs'],x['support_sha256']) for x in source['inverse_sets']},'all complete bounded images')
    for n in range(4096):
        q=n
        while q&1:q//=2
        need(F(q)==F(n) and (q==0 or (q%2==0 and q<=2*F(n)-2)),'zero-axis physical index normalization')
        normalizations+=1
    C=[24,1,3,6,14];c=2;p=profile(C,c);q=24084;t=0;T=graph(C,c,q,t)
    need(T==source['repair']['rows'],'independently constructed repair graph')
    need(len(T)==20 and scalar(p,q,t)==462661 and F(q)==20279,'repair arithmetic')
    need(p['e']==3,'repair cycle-product prediction')
    (P/'repair-graph.txt').write_text(' '.join(map(str,[len(T),462661,3]+T))+'\n')
    need(mat(p)[0][1]*mat(p)[1][2]-mat(p)[0][2]*mat(p)[1][1]==-96,'negative minor')
    print(json.dumps(dict(status='PASS',literal_retained_pairs=pairs,matrix_compositions=matrix_pairs,associative_triples=associative,binary_graphs=binary,inverse_sets=len(inverse),zero_axis_normalizations=normalizations,repair_graph_sha256=hashlib.sha256((P/'repair-graph.txt').read_bytes()).hexdigest()),sort_keys=True,indent=2))
if __name__=='__main__':main()
