#!/usr/bin/env python3
"""Small exact controls for the conditional reformulation; no spectrum search."""
from itertools import permutations
from pathlib import Path
import hashlib, json, math
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent

def need(b,msg):
    if not b: raise RuntimeError(msg)

def prime(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))

def graph(n,mask):
    out=[0]*n;k=0
    for a in range(n):
        for b in range(a+1,n):
            if mask>>k&1:out[a]|=1<<b
            else:out[b]|=1<<a
            k+=1
    return out

def strong(out):
    n=len(out)
    for v in range(n):
        seen=1<<v
        while True:
            new=seen
            for u in range(n):
                if seen>>u&1:new|=out[u]
            if new==seen:break
            seen=new
        if seen!=(1<<n)-1:return False
    return True

def literal(out):
    n=len(out);h=0;s=[0]*n;F=[0]*n
    for p in permutations(range(n)):
        good=[bool(out[p[i]]>>p[i+1]&1) for i in range(n-1)]
        if all(good):h+=1;s[p[0]]+=1
        for k in range(1,n+1):
            if all(good[:k-1]) and all(good[k:]):F[p[k-1]]+=1
    q=[]
    for v in range(n):
        q.append(sum(all(out[p[i]]>>p[i+1]&1 for i in range(n-2)) for p in permutations([u for u in range(n) if u!=v])))
    return h,s,q,F

def subset(out):
    n=len(out);full=(1<<n)-1
    end=[[0]*n for _ in range(full+1)];H=[0]*(full+1);H[0]=1
    for mask in range(1,full+1):
        for v in range(n):
            if mask>>v&1:
                old=mask^(1<<v)
                end[mask][v]=1 if old==0 else sum(end[old][u] for u in range(n) if old>>u&1 and out[u]>>v&1)
        H[mask]=sum(end[mask])
    F=[sum(end[mask][v]*H[full^mask] for mask in range(1,full+1) if mask>>v&1) for v in range(n)]
    return H[full],F

def all_profiles(out):
    # Independent enumeration of actual directed simple paths, retaining start/end names.
    n=len(out);full=(1<<n)-1;H=[0]*(full+1);H[0]=1
    end=[[0]*n for _ in range(full+1)];start=[0]*n
    def visit(mask,last,first):
        H[mask]+=1;end[mask][last]+=1
        if mask==full:start[first]+=1
        more=out[last]&~mask
        while more:
            bit=more&-more;more-=bit;visit(mask|bit,bit.bit_length()-1,first)
    for v in range(n):visit(1<<v,v,v)
    F=[sum(end[S][v]*H[full^S] for S in range(1,full+1) if S>>v&1) for v in range(n)]
    return H,[{'root':v,'profile':[H[full],start[v],H[full^(1<<v)],F[v]]} for v in range(n)]

def restrict(rows,order):
    return [sum((rows[u]>>v&1)<<i for i,v in enumerate(order)) for u in order]

def positive_both_branches():
    base=[102,124,120,113,33,64,16]
    count,profiles=all_profiles(base);need(count[-1]==35,'positive branch input')
    direct=literal(base)
    need([a['profile'] for a in profiles]==[[direct[0],direct[1][v],direct[2][v],direct[3][v]] for v in range(7)],'literal source profiles')
    out=[]
    inputs=[(107,[486,508,504,497,353,448,384,16,128],191),
            (109,[470,508,504,497,480,321,128,288,64],239)]
    for target,raw,mask in inputs:
        old=[i for i in range(9) if mask>>i&1]
        order=old+[i for i in range(9) if not mask>>i&1]
        rows=restrict(raw,order)
        need(restrict(rows,list(range(7)))==base,'old named source preservation')
        h,p=all_profiles(rows);need(h[-1]==target,'positive target count')
        dh,dF=subset(rows)
        need(dh==target and [a['profile'][3] for a in p]==dF,'independent child subset counts')
        birth=[]
        for v in (7,8):
            old_names=[f'u{i}' for i in range(7)]+(['new_a'] if v==8 else [])
            birth.append({'new_name':'new_a' if v==7 else 'new_b',
                          'beats_old':[old_names[u] for u in range(v) if rows[v]>>u&1],
                          'loses_to_old':[old_names[u] for u in range(v) if not rows[v]>>u&1],
                          'prefix_H':h[(1<<(v+1))-1]})
        balanced=[a['root'] for a in p if a['profile'][1]==a['profile'][2] and a['profile'][0]==a['profile'][3]]
        need(not balanced,'unexpected balanced endpoint')
        root=p[0]['profile'];need(root[1]==root[2],'old marked cut equality')
        out.append({'target':target,'rows':rows,'source_bank_rows':raw,'source_bank_subset':mask,
                    'bank_to_named_order':order,'old_to_new':list(range(7)),
                    'names':[f'u{i}' for i in range(7)]+['new_a','new_b'],
                    'retained_root':'u0','retained_inports':['u3','u4'],'births':birth,
                    'all_root_profiles':p,'balanced_roots':balanced})
    # Cheap complete scan only of the inherited 40 target representatives.
    bank=ROOT/'05-knowledge/results/metagraph_primegame_20261102/low_sources.txt'
    census={107:{'classes':0,'classes_containing_H35':0},109:{'classes':0,'classes_containing_H35':0}}
    for line in bank.read_text().splitlines():
        vals=list(map(int,line.split()));n,H=vals[:2]
        if H not in census:continue
        rows=vals[2:2+n];a,p=all_profiles(rows);need(a[-1]==H,'inherited bank target count')
        census[H]['classes']+=1;census[H]['classes_containing_H35']+=35 in a
    need(census=={107:{'classes':14,'classes_containing_H35':2},109:{'classes':26,'classes_containing_H35':4}},'target bank scan')
    return {'source_rows':base,'source_names':[f'u{i}' for i in range(7)],'source_all_root_profiles':profiles,
            'retained_root':'u0','retained_inports':['u3','u4'],'children':out,'inherited_target_bank_scan':census}

def main():
    levels=[];balanced=set();witnesses={}
    for n in range(3,6):
        ns=0;bc=0
        for mask in range(1<<(n*(n-1)//2)):
            out=graph(n,mask)
            if not strong(out):continue
            ns+=1;h,s,q,F=literal(out);dh,dF=subset(out)
            need((h,F)==(dh,dF),'independent literal/subset discrepancy')
            for r in range(n):
                if s[r]==q[r] and F[r]==h:
                    bc+=1;balanced.add(h)
                    witnesses.setdefault(h,{'rows':out,'root':r,'profile':[h,s[r],q[r],F[r]]})
        levels.append({'n':n,'strong_labelled':ns,'balanced_rooted_labelled':bc})
    need(11 not in balanced and 13 not in balanced,'balanced hostile failed')
    # These values cannot occur above order five: cited sharp strong minimum f(6)=15.
    need(prime(107) and prime(109) and prime(38501),'prime bridge controls')
    need(12833==41*313 and 38501==3*12833+2,'arithmetic bridge identity')
    exceptional=[];lane_predecessors=[]
    for N in range(5,10001,2):
        lane=prime(N) and N!=7 or N%7==0 and prime(N//7) and N//7!=3
        if not lane:continue
        need(N%3!=0,'unexpected lane divisible by3')
        c=2 if N%3==2 else 4;h=(N-c)//3
        need(3*h+c==N and h%2==1 and 0<h<N,'ternary predecessor')
        if h in (7,21):exceptional.append(N)
        lane_predecessors.append((N,h,c))
    need(exceptional==[23,67],'lane predecessor exceptions')
    # Readiness r gives output N=3h+2r and reset pair (N,h).
    gcd_tests=0
    for h in range(1,100,2):
        for r in range(1,100):
            N=3*h+2*r
            if N%3==0:continue
            need(math.gcd(N,h)==math.gcd(N,r),'reset gcd identity')
            gcd_tests+=1
    # Inherited all-order empty ready fibres, used as scoped inputs, not rerun.
    positive=positive_both_branches()
    refs=['05-knowledge/results/ready_return_pentagon_20261103/OBSTRUCTION.md',
          '05-knowledge/results/ready_return_pentagon_20261103/inverse.json',
          '05-knowledge/results/second_feedback_18979_20261103/PROOFS.md',
          '05-knowledge/results/second_feedback_18979_20261103/next_fibres.json',
          '05-knowledge/results/adaptive_metagraph_reformulation_20261103/ADAPTIVE.md',
          '05-knowledge/results/adaptive_metagraph_reformulation_20261103/METAGRAPH.md',
          '05-knowledge/results/coverage_machine_induction_20261103/axiom/AXIOM.md',
          '01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md',
          '05-knowledge/results/metagraph_primegame_20261102/low_sources.txt',
          '05-knowledge/results/metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md']
    inv=json.loads((ROOT/refs[1]).read_text())
    need(35 in inv['neither_ready'],'inherited35 receipt inconsistent')
    nxt=json.loads((ROOT/refs[3]).read_text())
    need(not any(x.get('kind')=='hit' and x.get('target')==12833 for x in nxt['bank']),'inherited12833 receipt inconsistent')
    print(json.dumps({'status':'PASS','scope':'Literal all labelled strong tournaments n3..5; independent subset controls; finite arithmetic checks. Ready fibre exclusions are inherited, not recomputed.',
        'levels':levels,'balanced_H_through_order5':sorted(balanced),
        'positive_named_H35_both_branches':positive,
        'balanced_witnesses':witnesses,'all_order_no_balanced_H_using_f6_15':[11,13],
        'ready_gap_to_prime_targets':[{'source':35,'targets':[107,109]},{'source':12833,'targets':[38501]}],
        'lane_predecessors_checked_through10000':len(lane_predecessors),'forbidden_predecessor_lane_exceptions':exceptional,
        'gcd_identities_checked':gcd_tests,
        'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in refs},
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
