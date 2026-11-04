import itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
def conv(r):return[(1<<len(r))-1^(1<<i)^x for i,x in enumerate(r)]
def inflate(Q,blocks):
    off=[];n=0
    for B in blocks:off.append(n);n+=len(B)
    rows=[0]*n
    for j,B in enumerate(blocks):
        for v,x in enumerate(B):
            rows[off[j]+v]=x<<off[j]
            for k,C in enumerate(blocks):
                if Q[j]>>k&1:rows[off[j]+v]|=((1<<len(C))-1)<<off[k]
    return rows
def tt(m):return[((1<<m)-1)^((1<<(i+1))-1)for i in range(m)]
data=json.loads((P/'selected.json').read_text());cases=[];queries=[];dedup={}
def put(r,cap,physical,a,c):
    key=(tuple(r),tuple(cap))
    if key not in dedup:
        dedup[key]=len(queries);queries.append(dict(compressed_rows=r,capacities=cap,physical_rows=physical,a_mask=a,c_mask=c))
    return dedup[key]
for index,d in enumerate(data):
    for m,ambient in ([(3,False),(5,False),(3,True)]if index==0 else [(3,False)]):
        physical_blocks=[tt(m)]+[[2,0]for _ in range(8)]
        compressed_blocks=[tt(m)]+[[0]for _ in range(8)]
        rootcap=[1]*m+[2]*8;old=m+16
        g=d['certificate_group'];groupnames=[m+2*(g-1),m+2*(g-1)+1];cg=m+g-1
        case=dict(node=d['node'],m=m,ambient=ambient,group=g,group_names=groupnames,terms=[])
        for r,weight in [(m-1,1)]+([(2,3)]if m==5 else []):
            term=dict(r=r,weight=weight,queries=[])
            for aa in range(4):
                for cc in range(4):
                    a=(aa&1)|((aa>>1&1)<<1);c=((1<<r)if cc&1 else 0)|(sum(1<<x for x in groupnames)if cc&2 else 0)
                    ca=a;cy=((1<<r)if cc&1 else 0)|((1<<cg)if cc&2 else 0)
                    cap=rootcap[:];pold=old;ctypes=m+8
                    if ambient:
                        a|=1<<(old+1);c|=1<<(old+1);ca|=1<<(ctypes+1);cy|=1<<(ctypes+1);cap +=[1,1];pold+=2;ctypes+=2
                    values=[]
                    for target in [False,True]:
                        Q=conv(d['rows'])if target else d['rows'];physical=inflate(Q,physical_blocks);comp=inflate(Q,compressed_blocks)
                        if ambient:
                            physical +=[1<<(old+1),(1<<old)-1]
                            for v in range(old):physical[v]|=1<<old
                            base=m+8;comp +=[1<<(base+1),(1<<base)-1]
                            for v in range(base):comp[v]|=1<<base
                        physical +=[1<<(pold+1),0];comp +=[1<<(ctypes+1),0]
                        for rows,n,x,y in [(physical,pold,a,c),(comp,ctypes,ca,cy)]:
                            for v in range(n):
                                if x>>v&1:rows[v]|=1<<n
                                else:rows[n]|=1<<v
                                if y>>v&1:rows[v]|=1<<(n+1)
                                else:rows[n+1]|=1<<v
                        values.append(put(comp,cap+[1,1],physical,a,c))
                    term['queries'].append(dict(a=a,c=c,sign=(-1)**(4-aa.bit_count()-cc.bit_count()),source=values[0],target=values[1]))
            case['terms'].append(term)
        cases.append(case)
(P/'physical_queries.json').write_text(json.dumps(dict(cases=cases,queries=queries),separators=(',',':'))+'\n')
(P/'physical_input.txt').write_text(str(len(queries))+'\n'+'\n'.join(str(len(x['compressed_rows']))+'\n'+' '.join(map(str,x['compressed_rows']))+'\n'+' '.join(map(str,x['capacities']))for x in queries)+'\n')
print(len(cases),len(queries))
