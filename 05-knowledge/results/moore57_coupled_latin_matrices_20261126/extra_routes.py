from pathlib import Path
from fixed_rank import rank_mod
import json,sys

def load(path):
    a=iter(map(int,Path(path).read_text().split()));m,e=next(a),next(a)
    rows=[set() for _ in range(m*m)]
    for _ in range(e):
        u,v=next(a),next(a);rows[u].add(v);rows[v].add(u)
    return m,rows

def main(grid,dom):
    m,D=load(grid)
    T=[list(map(int,l.split())) for l in Path(dom).read_text().splitlines()]
    x,n=T[0];T=T[1:]; assert n==len(T)
    rows=[]; keys=[]
    for k in range(1,5):
        for j in sorted({t[k] for t in T}):
            rows.append([int(t[k]==j) for t in T])
            keys.append((k,j))
    assert len(keys)==4*(m-1-len(D[x]))
    added=[]
    for y in range(m*m):
        if y//m in (0,28) or len(D[y])!=m-1 or y==x:continue
        want=int(y//m!=x//m and y%m!=x%m)-int(y in D[x])-len(D[y]&D[x])
        opts=[int(t[0] in D[y]) for t in T]
        if want:
            assert want==1
            added.append((y,opts))
        else:assert not any(opts)
    checks={}
    for p in [2,1009,1013]:
        r0=rank_mod(rows,p)
        new=[*rows,*[o for y,o in added]]
        q=rank_mod(new,p); b=rank_mod([u+[1] for u in new],p)
        checks[p]=[r0,q,b]
    print(json.dumps({"vertex":x,"candidates":n,"base_rows":len(rows),"extra_exact_routes":[[y//m,y%m,sum(o)] for y,o in added],"ranks_base_new_augmented":checks},indent=2))
if __name__=="__main__":main(*sys.argv[1:])
