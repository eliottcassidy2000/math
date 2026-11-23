"""Compare all physical neighbors of each near witness with both prescribed counts."""
from argparse import ArgumentParser
from random import Random
from subprocess import run

from near_witnesses import (A, D, EDGE, TARGET, WINNING_SLICES, flipped,
                            literal_receiver, nontags, ordinal, star_rows, strong)


def scan(binary, rows, sizes):
    n=len(rows)
    input_text = f"{n} {len(sizes)}\n"+' '.join(map(str,sizes))+'\n'
    input_text += ''.join(''.join(str((r>>i)&1) for i in range(n))+'\n' for r in rows)
    out=run([binary],input=input_text,text=True,capture_output=True,check=True)
    result={}
    for row in out.stdout.splitlines():
        u,v,d=map(int,row.split())
        assert 0<=u<v<n and (u,v) not in result
        result[u,v]=d
    assert len(result)==n*(n-1)//2
    return result


def brute(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    for u in range(n): dp[1<<u][u]=1
    for mask in range(1,1<<n):
        for u,d in enumerate(dp[mask]):
            if d:
                for v in range(n):
                    if (rows[u]>>v)&1 and (mask>>v)&1==0:
                        dp[mask|(1<<v)][v]+=d
    return sum(dp[-1])


def check(binary):
    rng=Random(682743)
    visits=0
    for sizes in ((1,),(2,),(1,1),(3,), (1,3,1),(4,2),(2,1,3),(2,4,1),(6,1)):
        n=sum(sizes)+1
        for _ in range(3):
            rows=[0]*n
            boundary={}
            for i,c in enumerate(sizes):
                for u in range(sum(sizes[:i]),sum(sizes[:i])+c): boundary[u]=i
            for u in range(n):
                for v in range(u+1,n):
                    a,b=(u,v) if (v!=n-1 and boundary[u]!=boundary[v]) or rng.getrandbits(1) else (v,u)
                    rows[a]|=1<<b
            h=brute(rows)
            deltas=scan(binary,rows,sizes)
            for (u,v),d in deltas.items():
                assert brute(flipped(rows,u,v))==h+d,(rows,u,v,d)
            visits+=len(deltas)
    print("small hosts",27,"literal arc flips",visits,"passed")


def main(binary):
    recv,tags=literal_receiver(nontags(4))
    marks=12|(35<<len(A))|(tags<<(len(A)+len(D)))
    for row in WINNING_SLICES:
        x,y,order,delta=row[-1]
        full=[]
        arc=tuple(x.bit_length()+len(A)+v for v in EDGE)
        print("base lower offset",delta,"original arc",arc,flush=True)
        for i,block in enumerate((flipped(D,*EDGE),D)):
            core=ordinal(ordinal(A,block),recv)
            host=star_rows(core,marks,x,y)
            sizes=(1,)*x.bit_length()+(len(A),len(D),len(recv))+(1,)*y.bit_length()
            assert len(host)==order
            deltas=scan(binary,host,sizes)
            assert deltas[arc]==(2 if i==0 else -2)
            vals={edge:delta+2*i+dh for edge,dh in deltas.items()}
            full.append(vals)
            twins={edge for edge in vals if (host[edge[0]]^host[edge[1]]) &
                   ~((1<<edge[0])|(1<<edge[1]))==0}
            neutral={edge for edge,d in deltas.items() if d==0}
            assert neutral==twins
            in_band={edge for edge,z in vals.items() if abs(z)<=10**9}
            assert in_band==twins|{arc},(arc,sorted(in_band-(twins|{arc})))
            closest_external=min(abs(z) for z in vals.values() if abs(z)>10**9)
            print("neutral arcs",len(neutral),"twins",len(twins),
                  "non-twin neutrals",sorted(neutral-twins),
                  "nearest magnitude outside 1e9",closest_external,flush=True)
            nearest=min(min(abs(z),abs(z-2)) for z in vals.values())
            best=[(edge,z) for edge,z in vals.items() if min(abs(z),abs(z-2))==nearest]
            nearby=sorted(set(z for z in vals.values() if abs(z)<=10**9))
            print("orientation",i,"arcs",len(vals),"distinct heights",len(set(vals.values())),
                  "target hits",sum(z in (0,2) for z in vals.values()),"closest distance",nearest,
                  "attained offsets",sorted({z for _,z in best}),"ties",len(best),"first",best[0],
                  "offsets within 1e9",nearby,flush=True)
        rails={e:(a,full[1][e]) for e,a in full[0].items() if e!=arc and abs(a-full[1][e])==2}
        closest=min(abs((a+b-2)//2) for a,b in rails.values())
        best=[(e,v) for e,v in rails.items() if abs((sum(v)-2)//2)==closest]
        e,endpoints=best[0]
        assert all(strong(flipped(star_rows(ordinal(ordinal(A,block),recv),marks,x,y),*e))
                   for block in (flipped(D,*EDGE),D))
        print("parallel arcs preserving unit flip",len(rails),"closest midpoint distance",closest,
              "ties",len(best),"first",e,endpoints,flush=True)


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("binary")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    if args.check: check(args.binary)
    else: main(args.binary)
