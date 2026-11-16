import csv, sys
from collections import Counter
from itertools import combinations
from check_27_63 import restriction, profile, key, scc_sizes

def entries(path):
    table={}
    all_hosts=[]
    for row in csv.DictReader(open(path)):
        h=int(row['H'])
        if h not in (27,63):continue
        a=tuple((int(row['code'])>>(8*i))&255 for i in range(8))
        if h==63:
            assert row['strong']=="1"
            all_hosts.append(a)
        elif h==27:
            for S,k in enumerate(profile(a)):
                if k==27:
                    b=restriction(a,[i for i in range(8) if S>>i&1])
                    table[key(b)]=b
    for a in sorted(all_hosts):
        h=profile(a)
        print('host',a,'induced',sorted(Counter((S.bit_count(),k) for S,k in enumerate(h) if k>=21).items()))
    return table,all_hosts

def main(path):
    sources,hosts=entries(path)
    core27=(62,124,120,112,32,64,17)
    k0=key(core27)
    assert k0 in sources
    best={}
    for a in sorted(hosts):
        for S in range(1,256):
            J=[i for i in range(8) if S>>i&1]
            if len(J)<6:continue
            b=restriction(a,J)
            for R in range(3):
                for E in combinations(list(combinations(range(len(b)),2)),R):
                    c=list(b)
                    for i,j in E:
                        c[i]^=1<<j;c[j]^=1<<i
                    if profile(c)[-1]!=27:continue
                    kk=key(c)
                    assert kk in sources
                    if kk not in best or best[kk][0]>R:best[kk]=R,(tuple(c),a,J,E)
    for k,c in sorted(sources.items()):
        print('source', len(k),scc_sizes(k),c,'original' if k==k0 else '',best.get(k))
    print('summary',Counter((len(k),scc_sizes(k),best.get(k,(None,))[0]) for k in sources))
if __name__=='__main__':main(sys.argv[1])
