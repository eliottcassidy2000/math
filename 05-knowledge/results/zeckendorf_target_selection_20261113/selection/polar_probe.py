from pathlib import Path
import json


def require(ok,msg):
    if not ok: raise ValueError(msg)


def hp(rows):
    n=len(rows);full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for i in range(n): dp[1<<i][i]=1
    for S in range(1,full+1):
        for v,z in enumerate(dp[S]):
            if not z: continue
            rem=rows[v]&~S
            while rem:
                b=rem&-rem;rem-=b;dp[S|b][b.bit_length()-1]+=z
    return sum(dp[-1])


def birth(rows,r,t,side,move):
    n=len(rows);full=(1<<n)-1
    if side=='R':
        out=tuple(row|((1<<n) if i==r else 0) for i,row in enumerate(rows))+(full^(1<<r),)
        return out,n if move else r,t
    out=tuple(row|((1<<n) if i!=t else 0) for i,row in enumerate(rows))+(1<<t,)
    return out,r,n if move else t


def canon(rows,r,t):
    # Canonical by chronology is intentionally not an isomorphism quotient.
    return rows,r,t


def main():
    layer={((2,4,1),0,t):'' for t in (1,2)}
    supports={};records=[]
    for n in range(3,10):
        for (rows,r,t),word in layer.items():
            h=hp(rows);supports.setdefault(h,(rows,r,t,word))
        record=dict(order=n,named_states=len(layer),attained=len(supports),max=max(supports),missing_below_301=[h for h in range(1,302,2) if h not in supports])
        records.append(record);print({k:record[k]for k in ('order','named_states','attained','max')},flush=True)
        if n==9:break
        nxt={}
        for (rows,r,t),word in layer.items():
            for side in 'RL':
                for move in (False,True):
                    out,u,v=birth(rows,r,t,side,move)
                    require(sum(row>>u&1 for row in out)==1,'entry guard')
                    require(out[v].bit_count()==1,'exit guard')
                    nxt.setdefault((out,u,v),word+side+('1' if move else '0'))
        layer=nxt
    Path(__file__).with_name('polar_probe.json').write_text(json.dumps(dict(scope='All specified two-polar-port birth histories from C3, orders3..9; no other repairs',records=records,witnesses={h:supports[h] for h in (13,17,23,25,27,31,35,39,41,47,67) if h in supports}),indent=2)+'\n')
if __name__=='__main__':main()
