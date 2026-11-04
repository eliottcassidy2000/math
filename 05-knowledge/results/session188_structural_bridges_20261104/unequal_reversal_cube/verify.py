"""Exact literal controls for unequal cyclic reversal cubes."""
import itertools,json
from functools import lru_cache

def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def opposite(a):
    return tuple(((1 << len(a))-1) ^ row ^ (1 << i) for i,row in enumerate(a))

def tournament(n, bits):
    out = [0]*n
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        u,v = (j,i) if bits >> bit & 1 else (i,j)
        out[u] |= 1 << v
    return tuple(out)

def hamilton_subsets(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1,1 << n):
        for v in range(n):
            count = dp[mask][v]
            if not count:
                continue
            unused = a[v] & ~mask
            while unused:
                bit = unused & -unused
                unused -= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return [1]+[sum(row) for row in dp[1:]]

@lru_cache(None)
def hamilton(a):
    return hamilton_subsets(a)[-1]

@lru_cache(None)
def covers(a):
    n = len(a)
    H = hamilton_subsets(a)
    dp = [[0]*(n+1) for _ in range(1 << n)]
    dp[0][0] = 1
    for mask in range(1,1 << n):
        sub = mask
        while sub:
            rest = mask ^ sub
            for k in range(1,n+1):
                dp[mask][k] += H[sub]*dp[rest][k-1]
            sub = (sub-1)&mask
    return tuple(dp[-1][1:])

def substitute(quotient, blocks):
    offsets = [0]
    for block in blocks:
        offsets.append(offsets[-1]+len(block))
    out = [0]*offsets[-1]
    for i, block in enumerate(blocks):
        for u, row in enumerate(block):
            out[offsets[i]+u] |= row << offsets[i]
            for j in range(len(blocks)):
                if quotient[i] >> j & 1:
                    out[offsets[i]+u] |= ((1 << len(blocks[j]))-1) << offsets[j]
    return tuple(out)

def extend(a, x, c):
    n = len(a)
    out = list(a)+[1 << (n+1),0]
    for v in range(n):
        if x >> v & 1:
            out[v] |= 1 << n
        else:
            out[n] |= 1 << v
        if c >> v & 1:
            out[v] |= 1 << (n+1)
        else:
            out[n+1] |= 1 << v
    return tuple(out)


def bits_of_blocks(blocks):
 offsets=[];masks=[];off=0
 for b in blocks:offsets.append(off);masks.append(((1<<len(b))-1)<<off);off+=len(b)
 return offsets,masks

def cross(blocks):
 qs=[(0,)+covers(b) for b in blocks]
 at=lambda q,k:q[k] if k<len(q) else 0
 return all(at(qs[i],k+2)*at(qs[j],k)<=at(qs[i],k+1)*at(qs[j],k+1)
            for i,j in itertools.permutations(range(3),2) for k in range(1,max(map(len,qs))))

def case(blocks):
 need(cross(blocks),'case lacks theorem hypothesis')
 active=[i for i,b in enumerate(blocks) if len(b)>=2];r=len(active)
 offsets,masks=bits_of_blocks(blocks);Q=((2,4,1),(4,1,2))
 states={}
 for theta in range(2):
  for bits in range(1<<r):
   changed={i for j,i in enumerate(active) if bits>>j&1}
   state=substitute(Q[theta],tuple(opposite(b) if i in changed else b for i,b in enumerate(blocks)))
   states[theta,bits]=state
 need(len(set(states.values()))==2**(1+r),'actual state count')
 root_queries=[(masks[i],masks[j]) for i,j in ((0,1),(1,2),(2,0))]
 queries=list(root_queries);decoders=[]
 for theta in range(2):
  for j,i in enumerate(active):
   T,U=states[theta,0],states[theta,1<<j]
   rep=offsets[i];outside=sum(1<<v for v in range(len(T)) if not masks[i]>>v&1 and T[v]>>rep&1)
   choice=None
   for a,c in itertools.product(range(1<<len(blocks[i])),repeat=2):
    x=outside|(a<<offsets[i]);y=outside|(c<<offsets[i])
    f,g=hamilton(extend(T,x,y)),hamilton(extend(U,x,y))
    if (f-g)%4==2:
     choice=(x,y);break
   need(choice is not None,'literal lifted child separator not found')
   queries.append(choice);decoders.append({'theta':theta,'child':i,'bit_index':j,'query':choice,'values':[f,g]})
 records=[];signatures=set();root_values={0:set(),1:set()}
 for (theta,bits),T in states.items():
  vals=tuple(hamilton(extend(T,a,c)) for a,c in queries)
  need(vals not in signatures,'fixed query response collision');signatures.add(vals)
  root_values[theta].add(vals[:3])
  for d in decoders:
   if d['theta']==theta:
    idx=queries.index(d['query']);bit=(bits>>d['bit_index'])&1
    need(vals[idx]==d['values'][bit],'child decoder affected by other children')
  records.append({'theta':theta,'bits':bits,'rows':T,'fixed_values':vals})
 need(all(len(v)==1 for v in root_values.values()),'root vector depends on child bit')
 root_sums=[sum(next(iter(root_values[t]))) for t in range(2)]
 need(root_sums[0]<root_sums[1],'root cyclic sum wrong order')
 need(len(queries)==3+2*r,'fixed count')
 return {'block_rows':blocks,'nontrivial_children':r,'states':len(states),'adaptive_query_bound':3+r,
         'fixed_query_count':len(queries),'distinct_fixed_queries':len(set(queries)),
         'root_query_sums':root_sums,'queries':queries,'child_decoders':decoders,'records':records}

if __name__=='__main__':
 cases=[(tournament(2,0),tournament(2,0),(2,4,1)),
        (tournament(2,0),tournament(3,0),(2,4,1)),
        ((0,),tournament(2,0),tournament(2,0)),((0,),(0,),(0,))]
 print(json.dumps({'status':'PASS','cases':[case(c) for c in cases],
                   'scope':'Conditional unequal cyclic reversal cubes; no unrestricted injectivity claim.'},indent=2,sort_keys=True))
