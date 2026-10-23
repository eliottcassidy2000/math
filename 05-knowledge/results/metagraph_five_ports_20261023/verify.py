"""Independent literal certificate consumer; does not canonicalize tournaments.

Checks every source isomorphism and target Hamiltonian word. Endpoint-aware
Held-Karp separately proves that no open target presentation was omitted.
The target set is all 201 classes in the 99 complete n9 ambiguity buckets.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter,defaultdict,deque
import argparse,csv,gzip,hashlib,json,math,struct

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];checks=0
def need(x,why):
    global checks
    checks+=1
    if not x:raise RuntimeError(why)
def rows(code,n=9):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>(n*(n-1)//2-1-k)&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return a
def encode(a,p):
    z=0
    for i,j in combinations(p,2):z=2*z+((a[i]>>j)&1)
    return z
def permutation(z):
    p=[(z>>(4*k))&15 for k in range(9)]
    need(sorted(p)==list(range(9)) and z<1<<36,'literal permutation carrier')
    return p
def path_counts(a):
    total=opened=0
    for first in range(9):
        dp=[[0]*9 for _ in range(512)];dp[1<<first][first]=1
        for mask in range(1,512):
            if not mask>>first&1:continue
            for last in range(9):
                value=dp[mask][last]
                if value:
                    avail=a[last]&~mask
                    while avail:
                        bit=avail&-avail;avail-=bit
                        dp[mask|bit][bit.bit_length()-1]+=value
        total+=sum(dp[511]);opened+=sum(dp[511][last] for last in range(9) if a[first]>>last&1)
    return total,opened
def peel(rc,cc,allowed):
    r=dict(rc);c=dict(cc);rn=defaultdict(set);cn=defaultdict(set);q=deque();fixed={}
    for a,b in allowed:rn[a].add(b);cn[b].add(a)
    q.extend((0,a) for a in r if len(rn[a])==1);q.extend((1,b) for b in c if len(cn[b])==1)
    def remove(a,b):
        rn[a].discard(b);cn[b].discard(a)
        if r[a] and len(rn[a])==1:q.append((0,a))
        if c[b] and len(cn[b])==1:q.append((1,b))
    while q:
        side,v=q.popleft()
        if side==0:
            if not r[v] or len(rn[v])!=1:continue
            a=v;b=next(iter(rn[v]));amount=r[a]
        else:
            if not c[v] or len(cn[v])!=1:continue
            b=v;a=next(iter(cn[v]));amount=c[b]
        need(amount<=r[a] and amount<=c[b],'infeasible capacity peel');r[a]-=amount;c[b]-=amount;fixed[a,b]=fixed.get((a,b),0)+amount
        if not r[a]:
            for bb in list(rn[a]):remove(a,bb)
        if not c[b]:
            for aa in list(cn[b]):remove(aa,b)
    need(all(not z or rn[a] for a,z in r.items()) and all(not z or cn[b] for b,z in c.items()),'unsupported residual')
    return fixed,{k:v for k,v in r.items() if v},{k:v for k,v in c.items() if v}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    up=ROOT/'05-knowledge/results/blue_quotient_20261020b/n9_classes.tsv.gz';universe={int(r['code']):r for r in csv.DictReader(gzip.open(up,'rt'),delimiter='\t')}
    need(len(universe)==191536 and sum(math.factorial(9)//int(r['Aut']) for r in universe.values())==1<<36,'inherited class orbit mass')
    menu=json.loads((HERE/'menu_n9.json').read_text());base=next(x for x in menu['menus'] if x['name']=='v4' and x['round']==1);need(base['all_collisions_recorded'],'missing collisions');groups=base['collisions'];gid={c:i for i,g in enumerate(groups) for c in g}
    need(len(gid)==201 and len(groups)==99,'scope');dual={c:int(r['converse_code']) for c,r in universe.items()};need(all(dual[c] in gid for c in gid),'selected targets not converse closed');mats={c:rows(c) for c in gid};tab={int(r['target']):r for r in csv.DictReader((HERE/'focused_targets.tsv').open(),delimiter='\t')};need(set(tab)==set(gid),'target inventory')
    path=HERE/'focused_ports.bin.gz';raw=gzip.decompress(path.read_bytes()) if path.exists() else (HERE/'focused_ports.bin').read_bytes();need(len(raw)%36==0,'binary row width')
    cells=list((i,j) for i in range(9) for j in range(i+2,9));ix={e:k for k,e in enumerate(cells)};seen=set();tc=Counter();bc=Counter();actual=[Counter() for _ in groups];source_set=set();examples={};allmass=0
    for mask,source,target,order,targetpath in struct.iter_unpack('<IQQQQ',raw):
        need(mask<1<<28 and mask not in seen,'repeated/invalid source presentation');seen.add(mask);need(target in gid and source in universe,'class endpoint');p=permutation(order);P=permutation(targetpath)
        a=[1<<(i+1) if i<8 else 0 for i in range(9)]
        for k,(i,j) in enumerate(cells):
            if mask>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
        need(not a[0]>>8&1,'source closure');need(encode(a,p)==source,'source explicit isomorphism')
        b=[x^511^(1<<i) for i,x in enumerate(a)]
        for i in range(8):b[i]^=1<<(i+1);b[i+1]^=1<<i
        inv=[P.index(i) for i in range(9)];need(encode(b,inv)==target,'target explicit path isomorphism');need(all(mats[target][P[i]]>>P[i+1]&1 for i in range(8)),'literal target path')
        ref=sum(((mask>>k)&1)<<ix[8-j,8-i] for k,(i,j) in enumerate(cells));blue=int(ref==mask);sa=a[0].bit_count();sb=a[8].bit_count();need(b[0].bit_count()==9-sa and b[8].bit_count()==7-sb,'endpoint transport')
        row=(source,sa,sb,blue);col=(target,9-sa,7-sb,blue);actual[gid[target]][row,col]+=1;tc[target]+=1;bc[target]+=blue;source_set.add(source);examples.setdefault(target,{'source':source,'mask':mask,'source_isomorphism':p,'target_hamiltonian_word':P,'source_scores':[sa,sb],'blue':blue});allmass+=1
    totals=[0,0]
    for target in sorted(gid):
        H,opened=path_counts(mats[target]);aut=int(universe[target]['Aut']);need(opened%aut==0 and opened//aut==tc[target],'complete target fibre by independent endpoint DP');r=tab[target];need((H,opened,tc[target],bc[target],aut)==tuple(int(r[k]) for k in ('H','open_paths','open_presentations','blue_presentations','aut')),'target summary');totals[0]+=H;totals[1]+=opened
    flat=Counter()
    for edges in actual:
        for (row,col),amount in edges.items():flat[row+(col[0],)]+=amount
    for (source,a,b,color,target),amount in flat.items():need(flat[dual[source],8-b,8-a,color,dual[target]]==amount,'complete reflected native port capacity')
    stats=Counter();residual=[];counterfeit=None
    for g,edges in zip(groups,actual):
        rc=Counter();cc=Counter()
        for (row,col),amount in edges.items():rc[row]+=amount;cc[col]+=amount
        allowed=[(row,(target,9-row[1],7-row[2],row[3])) for row in rc for target in g if (target,9-row[1],7-row[2],row[3]) in cc]
        fixed,rr,cr=peel(rc,cc,allowed);need(all(edges[e]==v for e,v in fixed.items()),'native capacity agreement');stats.update({'source_buckets':len(rc),'actual_port_groups':len(edges),'allowed_entries':len(allowed),'forced_entries':len(fixed),'forced_mass':sum(fixed.values()),'residual_mass':sum(rr.values()),'residual_source_buckets':len(rr),'residual_target_buckets':len(cr)})
        if rr:residual.append({'targets':g,'mass':sum(rr.values()),'source_buckets':len(rr),'target_buckets':len(cr)})
        # A concrete feasible switch inside one residual capacity rectangle.
        if counterfeit is None and len(g)==2:
            by=defaultdict(lambda:defaultdict(list))
            for (row,col),amount in edges.items():
                if row in rr and col in cr:by[row[1:]][col[0]].append(row)
            for typ,columns in sorted(by.items()):
                if set(columns)!=set(g) or typ[2]!=0:continue
                found=False
                for a in columns[g[0]]:
                    for b in columns[g[1]]:
                        if a==b:continue
                        delta=Counter({a+(g[0],):-1,a+(g[1],):1,b+(g[0],):1,b+(g[1],):-1})
                        for (s,x,y,color,t),d in list(delta.items()):delta[dual[s],8-y,8-x,color,dual[t]]+=d
                        if any(flat[e]+d<0 for e,d in delta.items()):continue
                        rd=Counter();cd=Counter();merged=Counter()
                        for (s,x,y,color,t),d in delta.items():
                            rd[s,x,y,color,gid[t]]+=d;cd[t,9-x,7-y,color]+=d;merged[tuple(sorted((min(s,dual[s]),min(t,dual[t]))))]+=d
                        need(not any(rd.values()) and not any(cd.values()),'counterfeit margins');need(all(d%2==0 for d in merged.values()),'counterfeit black twins')
                        if not any(merged.values()):continue
                        counterfeit={'target_pair':g,'edits':[{'source':s,'scores':[x,y],'blue':color,'target':t,'delta':d,'native_mass':flat[s,x,y,color,t]} for (s,x,y,color,t),d in sorted(delta.items()) if d],
                            'merged_black_delta':[[list(e),d] for e,d in sorted(merged.items()) if d],
                            'preserved':'one-round v4 source buckets, target endpoint/color capacities, explicit converse-reflection pairs, even merged black multiplicities'};found=True;break
                    if found:break
                if found:break
    need(stats['forced_mass']+stats['residual_mass']==allmass,'capacity conservation')
    probe=json.loads((HERE/'probe.json').read_text());rec={r['code']:r['rooted5'] for r in probe['orders']['9']['records']};new_signatures={c:tuple(sorted((row[6],row[24]) for row in rec[c])) for c in gid}
    need(all(len({new_signatures[c] for c in g})==len(g) for g in groups),'two motif selector incomplete')
    need(counterfeit is not None,'baseline has no demonstrated feasible switch')
    out={'status':'FINITE-EXACT all open presentations of 201 n9 target classes in the 99 complete v4-one-round ambiguity buckets; not full n9 wiring',
         'input_universe_sha256':hashlib.sha256(up.read_bytes()).hexdigest(),'port_bytes_sha256':hashlib.sha256(raw).hexdigest(),'port_record_width':36,'target_classes':len(gid),'source_classes':len(source_set),'literal_ports':allmass,'blue_ports':sum(bc.values()),'hamiltonian_paths':totals[0],'open_hamiltonian_paths':totals[1],
         'baseline_capacity':dict(stats),'baseline_residuals':residual,'baseline_counterfeit':counterfeit,'two_rooted_motifs':{'words':[66,204],'forced_port_groups':stats['actual_port_groups'],'forced_mass':allmass,'residual_mass':0},'literal_examples':{str(k):v for k,v in sorted(examples.items())},'checks':checks}
    if args.write:(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads((HERE/'verification.json').read_text())==out,'frozen verification receipt')
    print('literal ports',allmass,'source classes',len(source_set),'baseline',dict(stats),'two motifs force all',flush=True)
if __name__=='__main__':main()
