#!/usr/bin/env python3
"""Independent graph/path-cover consumer of a retained local-migration matrix bank.
No producer import. Rebuilds all edges and counts direct two-path partitions.
"""
import itertools,json,pathlib

def need(c,msg):
    if not c:raise RuntimeError(msg)
def tournament(n,code):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            a,b=(i,j) if code>>k&1 else (j,i);rows[a]|=1<<b;k+=1
    return rows

def dp(rows):
    n=len(rows);size=1<<n;d=[[0]*n for _ in range(size)];h=[1]+[0]*(size-1)
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,size):
        h[mask]=sum(d[mask])
        for i,v in enumerate(d[mask]):
            if not v:continue
            todo=rows[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1;d[mask|bit][j]+=v
    return h,d

def add(rows,cut):
    n=len(rows);return [rows[i]|(0 if cut>>i&1 else 1<<n) for i in range(n)]+[cut]

def state(rows,cut):
    n=len(rows);full=(1<<n)-1;h,d=dp(rows)
    rev=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
    _,first=dp(rev)
    q=sum(h[full^S]*sum(first[S][i] for i in range(n) if cut>>i&1) for S in range(1,full+1))
    g=sum(d[full][i] for i in range(n) if cut>>i&1)
    F=dp(add(rows,cut))[0][-1]
    return h[-1],F,q,g

def build(A,cut,local,z):
    n=len(A);m=len(local);full=add(A,cut)
    full+= [0]*(m-1)
    # Overwrite every local edge, leaving the A--old-apex cut untouched.
    for i in range(m):full[n+i]|=local[i]<<n
    for a in range(n):
        for j in range(1,m):full[a]|=1<<(n+j)
    root=n+z;keep=[i for i in range(n+m) if i!=root]
    B=[sum(1<<j for j,v in enumerate(keep) if full[u]>>v&1) for u in keep]
    S=sum(1<<j for j,v in enumerate(keep) if full[root]>>v&1)
    return B,S

def literal(rows,cut):
    n=len(rows);h=q=g=0
    for w in itertools.permutations(range(n)):
        e=[bool(rows[a]>>b&1) for a,b in zip(w,w[1:])]
        if all(e):h+=1;g+=bool(cut>>w[-1]&1)
        for split in range(n):
            if cut>>w[split]&1 and all(e[:max(0,split-1)]) and all(e[split:]):q+=1
    expanded=add(rows,cut)
    F=sum(all(expanded[a]>>b&1 for a,b in zip(w,w[1:])) for w in itertools.permutations(range(n+1)))
    return h,F,q,g

def main():
    root=pathlib.Path(__file__).resolve().parent
    bank=json.loads((root/'migration_matrix_bank.json').read_text());records=bank['records']
    need(len(records)==210,'local graph universe changed');cases=literal_cases=0;distinct={}
    for n in range(1,4):
        for code in range(1<<(n*(n-1)//2)):
            A=tournament(n,code)
            for cut in range(1<<n):
                old=state(A,cut)
                for r in records:
                    B,S=build(A,cut,r['local'],r['new_root']);new=state(B,S)
                    want=tuple(sum(a*b for a,b in zip(row,old)) for row in r['matrix'])
                    need(new==want,f'compiler mismatch {n,code,cut,r,new,want}')
                    if n<=2 and len(r['local'])<=3:
                        need(new==literal(B,S),'literal cover control');literal_cases+=1
                    cases+=1
    for r in records:distinct.setdefault(len(r['local']),set()).add(tuple(map(tuple,r['matrix'])))
    print(json.dumps({'status':'PASS','old_marked_parents':74,'local_rooted_tournaments':len(records),'actual_extension_checks':cases,'literal_permutation_and_cover_controls':literal_cases,'local_orders_distinct_matrices':{n:len(v) for n,v in distinct.items()},'matrix_bank_producer_sha256':bank['producer_sha256'],'scope':'All old labelled parents of order1..3 with all cuts; all local tournaments of orders2..4 with every distinct new root. The all-order compiler claim additionally requires its path-decomposition proof.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
