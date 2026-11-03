"""All two ears that may distinguish every fixed boundary vertex of a word family."""
from ear_cuts import *
from math import factorial


def run_coefficients(T,z):
    """Coefficient of j! S(k,j) when z is inflated to a transitive k-module."""
    fixed=[i for i in range(len(T)) if i!=z];e=len(fixed);width=e+2
    dp=[[[0]*width for _ in range(e+1)] for _ in range(1<<e)]
    dp[0][e][1]=1
    for i in range(e):dp[1<<i][i][0]=1
    for mask in range(1<<e):
        for i,v in enumerate(fixed):
            if mask>>i&1 and T[v]>>z&1:
                for j,x in enumerate(dp[mask][i][:-1]):dp[mask][e][j+1]+=x
        for i,v in enumerate(fixed):
            if mask>>i&1:continue
            target=dp[mask|1<<i][i]
            for last,u in enumerate(fixed+[z]):
                if T[u]>>v&1:
                    for j,x in enumerate(dp[mask][last]):target[j]+=x
    return tuple(sum(dp[-1][i][j] for i in range(e+1)) for j in range(width))


def evaluate(c,k):
    s=[1]+[0]*(len(c)-1)
    for _ in range(k):s=[0]+[j*(s[j]+s[j-1]) for j in range(1,len(s))]
    return sum(a*b for a,b in zip(c,s))


def inflate(T,z,k):
    addresses=[(i,j) for i in range(len(T)) for j in range(k if i==z else 1)]
    U=tuple(sum(1<<b for b,(v,t) in enumerate(addresses) if
                (u==v and s<t) or (u!=v and T[u]>>v&1)) for u,s in addresses)
    return U,addresses


def identity(old,new,offset):
    return all(v-3*(old[j] if j<len(old) else 0)==(offset if j==1 else 0)
               for j,v in enumerate(new) if j>=1)


def main():
    atlas=json.loads((HERE.parent/'ready_return_pentagon_20261103/languages.json').read_text())
    out=[]
    for row in atlas:
        if row['target'] not in (4,6):continue
        for f in row['families']:
            tail=tuple(map(tuple,f['tail']));T,p,_=source(((1,1),)+tail)
            h=2*f['A']+f['B'];require(tables(T)[0][-1]==h,'literal boundary seed')
            old=run_coefficients(T,2);hits=[]
            for m,q,value in all_two(T,{3*h+2,3*h+4}):
                U=birth(birth(T,m),q);c=run_coefficients(U,2)
                if identity(old,c,value-3*h):
                    hits.append(dict(masks=(m,q),offset=value-3*h,strong=strong(U),coefficients=c))
            record=dict(C=row['target'],tail=tail,A=f['A'],B=f['B'],source=T,source_coefficients=old,hits=hits)
            out.append(record)
            print(row['target'],tail,'hits',len(hits),'offsets',dict(Counter((x['offset'],x['strong']) for x in hits)),flush=True)
            (HERE/'boundary.json').write_text(json.dumps(out,indent=2)+'\n')
    print('complete',len(out),'hit families',sum(bool(r['hits']) for r in out),flush=True)


if __name__=='__main__':main()
