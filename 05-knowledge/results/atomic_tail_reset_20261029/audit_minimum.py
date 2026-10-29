"""Independent arbitrary-integer affine-cylinder minimum; no producer imports.

An input candidate is authenticated by literal native execution before it
is used as an upper bound. The entire smaller-source universe is searched.
"""
import argparse,json,time

BLOCKS=((1,1,1),(4,6,179),(7,11,16177),(7,11,12289),(7,11,9697))
WORDS=((1,),(3,1,1,1),(3,3,1,1,1,1,1),(3,1,3,1,1,1,1),(3,1,1,3,1,1,1))

def step(n):
    z=3*n+1;a=0
    while z%2==0:z//=2;a+=1
    return z,a

def block(n):
    q=1;word=[]
    while q:
        n,a=step(n);word.append(a);q+=-1+3*(a%3==0)
    return n,tuple(word)

def witness(row):
    n=row['least_source']
    for letter in row['word']:
        h,w=block(n)
        if w!=WORDS[letter] or h<=n:raise ValueError('invalid candidate upper bound')
        n=h
    if n!=int(row['endpoint']):raise ValueError('candidate endpoint mismatch')

def audit(row):
    witness(row);k=row['k'];bound=row['least_source']+1
    # At depth d, b branches give J=d+3b and D=d+5b.
    data={}
    for depth in range(1,k+1):
        for b in range(2*depth+1):
            J=depth+3*b;D=depth+5*b;mod=1<<(D+1)
            data[depth,b]=(D,mod-1,pow(3**J,-1,mod),3**J)
    stack=[(0,0,0,1)] # depth, branch count, affine carry, canonical source
    nodes=0;leaves=0;found=None
    while stack:
        depth,b,C,r=stack.pop();nodes+=1
        if r>=bound:continue
        D=depth+5*b
        if (1<<(D+1))>=bound:
            # This whole cylinder contains only r below the bound.
            # Check the remaining native trajectory, without chart transport.
            n=((3**(depth+3*b))*r+C)>>D
            okay=True
            for _ in range(k-depth):
                h,w=block(n)
                if h<=n:okay=False;break
                n=h
            leaves+=1
            if okay:bound=r;found=r
            continue
        children=[]
        for i,(j,d,c) in enumerate(BLOCKS):
            bb=b+(j-1)//3;CC=3**j*C+(1<<D)*c
            DD,mask,inv,_=data[depth+1,bb]
            rr=(((1<<DD)-CC)*inv)&mask
            if rr>=bound:continue
            if depth+1==k:bound=rr;found=rr
            else:children.append((depth+1,bb,CC,rr))
        # Descending push order means least lower-bound source is explored first.
        stack.extend(sorted(children,key=lambda x:x[3],reverse=True))
    if found!=row['least_source']:raise ValueError(('minimum mismatch',k,found,row['least_source']))
    return dict(k=k,least_source=found,affine_nodes=nodes,native_tail_checks=leaves,status='PASS')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('receipt');p.add_argument('--horizons',default='40,50,55,60');args=p.parse_args()
    source=json.load(open(args.receipt));bank={r['k']:r for r in source['rows']}
    rows=[audit(bank[int(k)]) for k in args.horizons.split(',')]
    print(json.dumps(dict(scope='Independent affine-cylinder audit; each stated horizon searches all positive odd sources',rows=rows),indent=2,sort_keys=True))
