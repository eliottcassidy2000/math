"""Independent determinant witnesses, prime exits, and named embeddings."""
from pathlib import Path
from itertools import permutations
from math import gcd,isqrt
import json
from pair_growth import LEFT,RIGHT,inventory,strong,require
ROOT=Path(__file__).resolve().parent

def extend2(rows,attachment):
    a,b,e=attachment
    return tuple(row|((a>>i&1)<<8)|((b>>i&1)<<9) for i,row in enumerate(rows))+((255^a)|(e<<9),(255^b)|((1-e)<<8))

def ordinal(blocks):
    n=sum(map(len,blocks));rows=[];offset=0
    for block in blocks:
        offset2=offset+len(block);later=((1<<n)-1)^((1<<offset2)-1)
        rows.extend((r<<offset)|later for r in block);offset=offset2
    return tuple(rows)

def birth(rows,r):
    n=len(rows)
    return tuple(x|((1<<n) if i==r else 0) for i,x in enumerate(rows))+(((1<<n)-1)^(1<<r),)

def nr(rows,r):
    x=len(rows);return birth(birth(rows,r),x),x

def source(M,N,prepared):
    rows=list(ordinal([(0,),(0,),M,(0,),N,(0,),(0,)]));n=len(rows)
    names=['r','s']+['M'+str(i) for i in range(len(M))]+['t']+['N'+str(i) for i in range(len(N))]+['u','v']
    beats={'r','s','u'} if prepared else {'r','t'}
    beats.update('N'+str(i) for i in range(len(N)))
    apex=0
    for i,name in enumerate(names):
        if name in beats:apex|=1<<i
        else:rows[i]|=1<<n
    return tuple(rows+[apex]),names+['a']

def hp_permutations(rows):
    return sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(rows)-1)) for p in permutations(range(len(rows))))

def automorphisms(rows,fixed=None):
    return sum((fixed is None or p[fixed]==fixed) and all(((rows[i]>>j)&1)==((rows[p[i]]>>p[j])&1) for i in range(len(rows)) for j in range(i+1,len(rows))) for p in permutations(range(len(rows))))

def delete(rows,r):
    keep=[i for i in range(len(rows)) if i!=r]
    return tuple(sum(((rows[i]>>j)&1)<<k for k,j in enumerate(keep)) for i in keep)

def from_code(n,code):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:rows[i]|=1<<j
            else:rows[j]|=1<<i
            k+=1
    return tuple(rows)

def check_aut_boundary():
    tested=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            rows=from_code(n,code)
            for r in range(n):
                out,x=nr(rows,r);indeg=sum(row>>r&1 for row in rows)
                expected=automorphisms(rows,r) if indeg else 3*automorphisms(delete(rows,r))
                require(automorphisms(out)==expected,'NR automorphism boundary')
                tested+=1
    return tested

def main():
    data=json.loads((ROOT/'pair_growth2.json').read_text());records=[]
    base,names0=source(LEFT,RIGHT,True);old,x0=nr(base,0);names0+=['x','y']
    for i,hit in enumerate(data['hits']):
        M=extend2(LEFT,hit['left_attachment']);N=extend2(RIGHT,hit['right_attachment'])
        require(inventory(M)==tuple(hit['M']),'left independent profile')
        require(inventory(N)==tuple(hit['N']),'right independent profile')
        require(strong(M) and strong(N),'strong modules')
        h,p=hit['M'];j,q=hit['N'];e=j*(6*p+h)-6*h*q
        A=12*p*q+6*h*q+2*p*j-h*j;Z=12*p*q+2*h*q+6*p*j-h*j;Y=3*Z+2*h*j
        require(e==hit['e'] and Y==3*A+2*e,'exact digit')
        require(gcd(Y,Z)==1,'primitive balanced port')
        before,names=source(M,N,False);prepared,_=source(M,N,True);out,root=nr(prepared,0);names+=['x','y']
        changed=[(a,b) for a in range(len(before)) for b in range(a+1,len(before)) if ((before[a]>>b)&1)!=((prepared[a]>>b)&1)]
        require(len(changed)==3,'three preparation flips')
        embedding=[names.index(name) for name in names0]
        require(all(((old[a]>>b)&1)==((out[embedding[a]]>>embedding[b])&1) for a in range(len(old)) for b in range(len(old))),'old output induced')
        require(embedding[x0]==root and strong(out),'retained root and strong output')
        require([v for v in range(len(out)) if sum(row>>v&1 for row in out)==1]==[root,root+1],'canonical newest pair')
        divisor=next((d for d in range(3,isqrt(Y)+1,2) if Y%d==0),None)
        record=dict(M_rows=M,N_rows=N,source=A,prepared=Z,target=Y,e=e,prime=divisor is None,first_divisor=divisor,balanced_state=[Y,Z,Z,Y],output_rows=out,root=root,old_output_embedding=embedding,changed_arcs=changed)
        if i==2:
            counts=[hp_permutations(M),hp_permutations(N)]
            require(counts==[h,j],'independent permutation counts')
            record['permutation_H']=counts
        records.append(record)
    require(sum(z['prime'] for z in records)==4,'prime target count')
    boundary=check_aut_boundary()
    out=dict(status='PASS',scope='Independent witnesses and embeddings, not full large-output path enumeration',records=records,NR_automorphism_controls=boundary)
    (ROOT/'pair_check.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS six physical compensated pairs, four prime exits, induced marked embeddings;',boundary,'NR group controls')

if __name__=='__main__':main()
