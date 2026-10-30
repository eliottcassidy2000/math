"""Literal sharpness checks for AP1 and GP. Standard library; assertions unused."""
from itertools import product
import json

def req(ok,why):
    if not ok:raise RuntimeError(why)

def build(M,kind):
    head=4 if kind=='AP1' else 3
    n=head+2*M; gates=[]
    def add(op,a,b):gates.append((op,a,b));return n+len(gates)-1
    g=add('and',0,1)
    if kind=='AP1':h=add('and',2,3);k=add('xor',g,h)
    else:k=add('xor',g,2)
    products=[add('and',head+2*j,head+2*j+1) for j in range(M)]
    ip=products[0]
    for p in products[1:]:ip=add('xor',ip,p)
    output=add('xor',k,ip)
    return n,gates,output

def evaluate(n,gates,mask):
    v=[(mask>>j)&1 for j in range(n)]
    for op,a,b in gates:v.append((v[a]&v[b]) if op=='and' else v[a]^v[b])
    return v[-1]

def simplify(n,gates,kind):
    # GP G is constant by the restricted three-point domain, rather than an
    # affine assignment of either x or y.
    aliases={i:('input',i) for i in range(n)}
    if kind=='AP1':aliases[2]=('constant',1);aliases[1]=('constant',0)
    survivors=[];removed=[]
    for j,(op,a,b) in enumerate(gates):
        x,y=aliases[a],aliases[b]
        if kind=='GP' and j==0:out=('constant',0)
        elif x[0]==y[0]=='constant':out=('constant',(x[1]&y[1]) if op=='and' else x[1]^y[1])
        elif x[0]=='constant' or y[0]=='constant':
            fixed,other=(x,y) if x[0]=='constant' else(y,x)
            if op=='and':out=other if fixed[1] else('constant',0)
            else:
                req(fixed[1]==0,'examples need no unary complement');out=other
        else:
            req(x!=y,'unexpected repeated-pin survivor');out=('gate',n+j);survivors.append(j)
        aliases[n+j]=out
        if out!=('gate',n+j):removed.append(j)
    return removed,survivors

def rank(rows):
    basis={}
    for row in rows:
        while row:
            p=row.bit_length()-1
            if p in basis:row^=basis[p]
            else:basis[p]=row;break
    return len(basis)

def check(M,kind):
    n,gates,out=build(M,kind);head=4 if kind=='AP1' else 3
    def domain(mask):return kind=='GP' or bool((mask&1) or (mask&4))
    vals=[];restricted=0
    for x in range(1<<n):
        bits=[(x>>i)&1 for i in range(n)]
        tail=sum(bits[head+2*j]*bits[head+2*j+1] for j in range(M))%2
        expected=(bits[0]*bits[1])^(bits[2]*bits[3] if kind=='AP1' else bits[2])^tail
        req(evaluate(n,gates,x)==expected,('literal truth table',M,kind,x));vals.append(expected)
        selected=(bits[2]==1 and bits[1]==0) if kind=='AP1' else not(bits[0] and bits[1])
        if selected:
            req(expected==(bits[3] if kind=='AP1' else bits[2])^tail,'restricted identity');restricted+=1
    for i in range(n):
        req(any(domain(x) and domain(x^(1<<i)) and vals[x]!=vals[x^(1<<i)] for x in range(1<<n)),('initial essentiality',kind,M,i))
    def selected(x):
        return (((x>>2)&1)==1 and ((x>>1)&1)==0) if kind=='AP1' else not((x&1) and (x&2))
    residual_indices=[3 if kind=='AP1' else 2,*range(head,n)]
    for i in residual_indices:
        req(any(selected(x) and selected(x^(1<<i)) and vals[x]!=vals[x^(1<<i)] for x in range(1<<n)),('residual essentiality',kind,M,i))
    req(len(gates)==n-1,'initial minimum size')
    removed,survivors=simplify(n,gates,kind)
    expected_removed=3 if kind=='AP1' else 2
    req(len(removed)==expected_removed,'sharp deletion count');req(len(survivors)==2*M,'residual minimum size')
    # Polar coefficients from literal f values, independent of the declared pairs.
    rows=[sum((vals[(1<<i)^(1<<j)]^vals[1<<i]^vals[1<<j]^vals[0])<<j for j in range(n)) for i in range(n)]
    R=rank(rows);req(R==(n if kind=='AP1' else n-1),'polar rank')
    if kind=='GP':req(rows[2]==0 and vals[1<<2]==1,'nonzero quadratic value on radical')
    d=M+2 if kind=='AP1' else M+1
    a=n-1 if kind=='AP1' else n
    req(a>=d+2 and a-1>d,'loop and nonconstant postguard')
    return {'case':kind,'M':M,'inputs':n,'d':d,'safe_core':a,'gate_count':len(gates),'deleted_gate_indices':removed,'surviving_gates':len(survivors),'truth_table_evaluations':1<<n,'restricted_evaluations':restricted,'polar_rank':R}

if __name__=='__main__':print(json.dumps({'status':'PASS','checks':[check(M,k) for k,M in product(['AP1','GP'],range(1,5))]},indent=2,sort_keys=True))
