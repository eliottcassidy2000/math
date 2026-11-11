#!/usr/bin/env python3
"""Literal independent checks for core/ready-atlas maps; no repository imports."""
import argparse
import itertools
import json
from functools import lru_cache

CORES = {
    'K1': (0,), 'C3': (2, 4, 1), 'S4': (2, 4, 9, 3),
    'S5a': (2, 4, 9, 19, 7), 'S5b': (2, 4, 25, 3, 11),
}
B0 = (16, 1, 3, 7, 14)
B1 = (20, 1, 2, 7, 14)

def need(condition, label):
    if not condition:
        raise RuntimeError(label)

def dp(rows):
    n = len(rows)
    full = (1 << n) - 1
    end = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        end[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v, count in enumerate(end[mask]):
            avail = rows[v] & ~mask & full
            while count and avail:
                bit = avail & -avail
                avail -= bit
                end[mask | bit][bit.bit_length() - 1] += count
    hp = [sum(x) for x in end]
    hp[0] = 1
    return hp, end

def induced(rows, ids):
    return tuple(sum(((rows[i] >> j) & 1) << b for b, j in enumerate(ids)) for i in ids)

def iso(left, right):
    if len(left) != len(right):
        return None
    n = len(left)
    for p in itertools.permutations(range(n)):
        if all(((left[i] >> j) & 1) == ((right[p[i]] >> p[j]) & 1)
               for i in range(n) for j in range(n)):
            return list(p)
    return None

def transitive_dictionary(rows):
    ans = []
    for root in range(len(rows)):
        rest = [v for v in range(len(rows)) if v != root]
        for p in itertools.permutations(rest):
            if all(rows[p[i]] >> p[j] & 1 for i in range(len(p)) for j in range(i+1, len(p))):
                ans.append({'root': root, 'spine': list(p),
                            'word': ''.join(str((rows[root] >> v) & 1) for v in p)})
    return ans

def signature(rows, out):
    hp, _ = dp(rows)
    full = (1 << len(rows)) - 1
    val = [0] * 8
    def walk(mask, first, last):
        w = hp[full ^ mask]
        sign = (-1) ** mask.bit_count()
        yeses = (True, not(out >> last & 1), bool(out >> first & 1),
                 bool(out >> first & 1) and not(out >> last & 1))
        for i, yes in enumerate(yeses):
            if yes:
                val[2*i] += w
                val[2*i+1] += sign * w
        avail = rows[last] & ~mask & full
        while avail:
            bit = avail & -avail
            avail -= bit
            walk(mask | bit, first, bit.bit_length()-1)
    for root in range(len(rows)):
        walk(1 << root, root, root)
    return (hp[-1], *val)

def seven(sig):
    d, kp, km, ip, im, op, om, jp, jm = sig
    need(km == -d, 'signed K reciprocity')
    return (d, d+kp, ip, op, -im, -om, d+jp+jm)

@lru_cache(None)
def g(n):
    if n < 2:
        return 1
    return g(n//2) if n & 1 else g(n//2) + n//2 + ((n//2) & 1)

def compiled(sig, a, b):
    h, q, e, s, u, v, t = seven(sig)
    return q*a*b + e*a + s*b + h*(g(a)+g(b)-2) + t + u*(a&1) + v*(b&1)

def attach(rows, out, a=0, b=0):
    left = list(map(int, bin(a)[2:])) if a else []
    right = [1-int(x) for x in bin(b)[2:][::-1]] if b else []
    size = len(left) + len(rows) + len(right)
    ans = [0] * (size+1)
    for i in range(size):
        for j in range(i+1, size):
            internal = len(left) <= i < j < len(left)+len(rows)
            win = bool(rows[i-len(left)] >> (j-len(left)) & 1) if internal else True
            u, v = (i,j) if win else (j,i)
            ans[u] |= 1 << v
    bits = left + [(out >> i) & 1 for i in range(len(rows))] + right
    for v, bit in enumerate(bits):
        ans[size if bit else v] |= 1 << (v if bit else size)
    return tuple(ans)

def root_profiles(rows):
    hp, end = dp(rows)
    n, full = len(rows), (1 << len(rows))-1
    reverse = tuple(sum((i != j and not(rows[i] >> j & 1)) << j for j in range(n)) for i in range(n))
    _, rev_end = dp(reverse)
    ans = []
    for r in range(n):
        F = sum(end[mask][r] * hp[full ^ mask] for mask in range(1,full+1) if mask >> r & 1)
        s, q = rev_end[-1][r], hp[full ^ (1 << r)]
        ans.append({'root':r, 'h':hp[-1], 's':s, 'q':q, 'F':F,
                    'gap':F-hp[-1], 'delta':F-hp[-1]-q+s})
    return ans

def ready_source(tail, k):
    """Letter bit0 is [u->a], bit1 is [u->b]; ports a->b->r and a->r."""
    word = [3]*k + list(tail)
    m = len(word)
    ans = [0]*(m+3)
    root = m+2
    for i in range(m):
        for j in range(i+1,m):
            ans[i] |= 1 << j
        for port in range(2):
            u,v = (i,m+port) if word[i] >> port & 1 else (m+port,i)
            ans[u] |= 1 << v
        ans[root] |= 1 << i
    ans[m] |= (1 << (m+1)) | (1 << root)
    ans[m+1] |= 1 << root
    return tuple(ans)

def tail_core(tail):
    rows = ready_source(tail,1)
    return induced(rows, range(1,len(rows)-1))

def arc_difference(a,b):
    return [(i,j) for i in range(len(a)) for j in range(i+1,len(a))
            if ((a[i] >> j) & 1) != ((b[i] >> j) & 1)]

def birth_beaten_only_by(rows, root):
    n = len(rows)
    return tuple(row | ((1 << n) if v == root else 0) for v,row in enumerate(rows)) + (((1 << n)-1) ^ (1 << root),)

def two_birth_reset(rows, root):
    return birth_beaten_only_by(birth_beaten_only_by(rows,root),len(rows))

def slot_profile(rows, max_slots=6):
    hp,_ = dp(rows)
    full = (1 << len(rows))-1
    cur = [0]*(full+1)
    cur[0] = 1
    vals=[]
    for slots in range(1,max_slots+1):
        nex=[0]*(full+1)
        for mask in range(full+1):
            sub=mask
            while True:
                nex[mask] += hp[sub]*cur[mask^sub]
                if sub == 0:
                    break
                sub=(sub-1)&mask
        cur=nex
        vals.append(cur[full])
    return vals

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fail',action='store_true')
    args=parser.parse_args()
    if args.fail:
        need(False,'intentional hostile control')
    result={'scope':'Literal finite audits plus exact compiler coefficients; no all-target coverage claim.'}
    result['inner_apex_dictionary']={name:transitive_dictionary(rows) for name,rows in CORES.items()}
    result['label_maps']={'S5a_to_B1':iso(CORES['S5a'],B1),'S5b_to_B0':iso(CORES['S5b'],B0)}
    need(all(result['label_maps'].values()),'inherited seed isomorphisms')
    expected=[r**5+6*r**3+2*r for r in range(1,7)]
    need(slot_profile(B0)==expected==slot_profile(B1),'full degree5 profile')
    result['slot_values']=expected
    pentagon=[]
    for tail in ((),(1,),(2,),(1,2),(2,1)):
        W=tail_core(tail)
        record={'tail':tail,'deletion_tail_rows':W,'h':dp(W)[0][-1]}
        record['is_strong_S4']=iso(W,CORES['S4']) is not None
        pentagon.append(record)
    need([x['h'] for x in pentagon]==[1,3,1,5,5],'pentagon deletion quotient')
    result['pentagon_deletion']=pentagon
    atlas={name:[] for name in ('S5a','S5b')}
    for tail in itertools.product(range(4),repeat=3):
        if tail[0]==3 or tail[-1]==0:
            continue
        W=tail_core(tail)
        for name in atlas:
            mapping=iso(W,CORES[name])
            if mapping is not None:
                state=root_profiles(ready_source(tail,1))[-1]
                atlas[name].append({'tail':tail,'root_state':state,'core_iso':mapping})
    need(all(min(x['root_state']['gap'] for x in values)==4 for values in atlas.values()),'first H9 tail budget')
    result['length3_H9_atlas']=atlas
    families=[]
    for tail,coefficient,constant in (((1,0,2),22,5),((1,2,2),25,5),((2,3,1),34,5),((2,2,1),34,5)):
        W=tail_core(tail)
        sig=signature(W,(1 << len(tail))-1)
        d,kp,km,ip,im,op,om,jp,jm=sig
        need((ip,d-ip-im+jp+jm)==(coefficient,constant),'symbolic all-k coefficient')
        states=[]
        for k in range(1,6):
            T=ready_source(tail,k)
            roots=root_profiles(T)
            chosen=roots[-1]
            need((chosen['h'],chosen['s'],chosen['q'],chosen['gap'])==(coefficient*2**k+constant,9,9,4),'C4 family state')
            need(dp(attach(W,(1 << len(tail))-1,2**k-1,0))[0][-1]==chosen['h'],'independent compiler graph')
            states.append({'k':k,'root':chosen,'ready_roots':[r for r in roots if r['gap'] in (1,2)]})
        families.append({'tail':tail,'coefficient':coefficient,'constant':constant,'signature':sig,'seven':seven(sig),'states':states})
    result['C4_families']=families
    seven_a,seven_b=families[2]['seven'],families[3]['seven']
    need(tuple(a-b for a,b in zip(seven_a,seven_b))==(0,0,0,3,0,1,0),'neutral pair exact observer defect')
    for k in range(1,6):
        diff=arc_difference(ready_source((2,3,1),k),ready_source((2,2,1),k))
        need(diff==[(k+1,k+3)],'neutral pair one edge')
    result['C4_neutral_switch']={'tails':[[2,3,1],[2,2,1]],'seven_a':seven_a,'seven_b':seven_b,
                                 'edge_at_k':['k+1','k+3'],'preserved':'root, ordered ports, tape order, h,s_root,q,F',
                                 'lost':'Full flank observer; count difference after right code B is 3B+(B mod2).'}
    base_obstructions=[]
    one_birth_tests=0
    for tail in ((2,3,1),(2,2,1)):
        T=ready_source(tail,1)
        roots=root_profiles(T)
        indegrees=[len(T)-1-row.bit_count() for row in T]
        need([r for r,d in enumerate(indegrees) if d==1]==[0],'only inherited indegree1 root')
        sole_predecessor=next(v for v in range(len(T)) if T[v]&1)
        parent=induced(T,range(1,len(T)))
        predecessor_starts=root_profiles(parent)[sole_predecessor-1]['s']
        need(predecessor_starts==7,'specified parent start obstruction')
        base_obstructions.append({'tail':tail,'rows':T,'starts':[r['s'] for r in roots],
                                  'indegrees':indegrees,'only_indegree1_root':0,
                                  'sole_predecessor':sole_predecessor,'predecessor_starts_after_delete':predecessor_starts})
        for out in range(1 << len(T)):
            child=attach(T,out)
            need(all(r['gap'] not in (1,2) for r in root_profiles(child)),'one birth cannot restore readiness')
            one_birth_tests+=1
    result['addition_obstruction_base_data']=base_obstructions
    result['exhaustive_one_birth_children_without_ready_root']=one_birth_tests
    order7=[]
    for tail in ((),(1,),(2,),(1,2),(2,1)):
        k=4-len(tail)
        count=dp(ready_source(tail,k))[0][-1]
        order7.append(count)
    need(order7==[47,49,39,55,59],'exceptional ready order7 counts')
    result['C2_ready_order7_counts']=order7
    returns=[]
    for k in range(1,6):
        sources=[ready_source(tail,k) for tail in ((2,3,1),(2,2,1))]
        # Removing the middle tail vertex erases the changed edge and restores BA.
        gate2=[induced(T,[v for v in range(len(T)) if v!=k+1]) for T in sources]
        need(gate2[0]==gate2[1]==ready_source((2,1),k),'literal shared 2-ready gate')
        r2=root_profiles(gate2[0])[-1]
        need((r2['h'],r2['s'],r2['q'],r2['gap'])==(14*2**k+3,5,5,2),'2-ready returned response')
        reset2=two_birth_reset(gate2[0],len(gate2[0])-1)
        need(dp(reset2)[0][-1]==42*2**k+13,'2-ready returned reset')
        # Removing port a leaves a common graph with a global sink; removing
        # that last tail vertex gives the pure 10^(k+2) source.
        ancestor=[induced(T,[v for v in range(len(T)) if v!=k+3]) for T in sources]
        need(ancestor[0]==ancestor[1],'shared one-port ancestor')
        need(dp(ancestor[0])[0][-1]==2**(k+2)+1,'Fermat-shaped ancestor count')
        gate1=[induced(T,[v for v in range(len(T)) if v not in (k+2,k+3)]) for T in sources]
        need(gate1[0]==gate1[1],'literal shared 1-ready gate')
        r1=root_profiles(gate1[0])[0]
        need((r1['h'],r1['s'],r1['q'],r1['gap'])==(2**(k+2)+1,2**(k+1),2**(k+1)+1,1),'1-ready returned response')
        reset1=two_birth_reset(gate1[0],0)
        need(dp(reset1)[0][-1]==12*2**k+5,'1-ready returned reset')
        returns.append({'k':k,'source_count':34*2**k+5,'gate2_rows':gate2[0],'gate2_root':r2,
                        'reset2_rows':reset2,'reset2_h':42*2**k+13,
                        'Fermat_ancestor_rows':ancestor[0],'gate1_rows':gate1[0],'gate1_root':r1,
                        'reset1_rows':reset1,'reset1_h':12*2**k+5})
    result['literal_ready_returns']=returns
    switches=[]
    for out,delta in ((9,0),(22,-2)):
        sig0,sig1=signature(B0,out),signature(B1,out)
        for b in range(16):
            T0,T1=attach(B0,out,0,b),attach(B1,out,0,b)
            h0,h1=dp(T0)[0][-1],dp(T1)[0][-1]
            need(h1-h0==delta,'contextual switch difference')
            need(h0==compiled(sig0,0,b) and h1==compiled(sig1,0,b),'compiler independent endpoint count')
            need(arc_difference(T0,T1)==[(0,2)],'switch one edge')
        switches.append({'outmask':out,'signature0':sig0,'signature1':sig1,'delta_at_A0':delta,
                         'root0':root_profiles(attach(B0,out))[-1],
                         'root1':root_profiles(attach(B1,out))[-1]})
    result['labelled_context_switches']=switches
    # Cheap universal-switch hostile: no pair of cuts has equal two linear slopes
    # and a constant +/-2 difference after fixing the two parities.
    bad=[]
    for o0,o1 in itertools.product(range(32),repeat=2):
        s0,s1=seven(signature(B0,o0)),seven(signature(B1,o1))
        if s0[2:4]!=s1[2:4]:
            continue
        for x,y in itertools.product(range(2),repeat=2):
            delta=s1[6]-s0[6]+(s1[4]-s0[4])*x+(s1[5]-s0[5])*y
            if abs(delta)==2:
                bad.append((o0,o1,x,y))
    need(not bad,'no two-free-flank constant2 switch')
    result['two_free_flank_constant2_switches']=bad
    # Five core-free positive families have a sparse overlap, proved in NOTE.
    collisions=[]
    for m in range(1,101):
        h=(4**(m+1)+6*m+5)//9
        need((4**(m+1)+6*m+5)%9==0,'alternating integrality')
        q=h-1
        if q>0 and q&(q-1)==0:
            collisions.append({'m':m,'r':q.bit_length()-1,'h':h})
    need(collisions==[{'m':1,'r':1,'h':3},{'m':2,'r':3,'h':9}],'primitive-family collision control')
    result['primitive_collisions_through_m100']=collisions
    print(json.dumps(result,sort_keys=True,indent=2))

if __name__=='__main__':
    main()
