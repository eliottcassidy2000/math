"""Independent even-path cut audit of the named bounded strong frontiers.

No repository implementation is imported. The originating JSON is read only
after each independent frontier has been generated, for exact comparison.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCES=((126,124,88,112,32,196,144,31),
         (126,124,88,80,32,204,144,31))

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def subset_counts(rows):
    # First-vertex recurrence, distinct from the originating last-vertex DP.
    n=len(rows);full=(1<<n)-1
    p=[[0]*n for _ in range(full+1)];h=[0]*(full+1);h[0]=1
    for mask in range(1,full+1):
        vv=mask
        while vv:
            bit=vv&-vv;vv-=bit;v=bit.bit_length()-1;rest=mask^bit
            if not rest:p[mask][v]=1
            else:
                choices=rows[v]&rest
                while choices:
                    b=choices&-choices;choices-=b
                    p[mask][v]+=p[rest][b.bit_length()-1]
        h[mask]=sum(p[mask])
    return h

def even_capacities(rows,h):
    # For each directed path P of even order from i to j, weight by H(S-P).
    n=len(rows);full=(1<<n)-1;C=[[0]*n for _ in rows]
    for start in range(n):
        p=[[0]*n for _ in range(full+1)];p[1<<start][start]=1
        for mask in range(1,full+1):
            if not(mask>>start&1):continue
            for end,value in enumerate(p[mask]):
                if not value:continue
                if not(mask.bit_count()&1):C[start][end]+=value*h[full^mask]
                choices=rows[end]&~mask&full
                while choices:
                    bit=choices&-choices;choices-=bit
                    p[mask|bit][bit.bit_length()-1]+=value
    return C

def attach(rows,O):
    n=len(rows)
    return tuple(r|(0 if O>>v&1 else 1<<n) for v,r in enumerate(rows))+(O,)

def strong(rows):
    full=(1<<len(rows))-1
    for root in range(len(rows)):
        seen=1<<root
        while True:
            enlarged=seen
            for v,r in enumerate(rows):
                if seen>>v&1:enlarged|=r
            if enlarged==seen:break
            seen=enlarged
        if seen!=full:return False
    return True

def cut_value(h,C,O):
    return h+2*sum(x for i,row in enumerate(C) if O>>i&1
                   for j,x in enumerate(row) if not(O>>j&1))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fail',action='store_true')
    if parser.parse_args().fail:need(False,'intentional independent frontier hostile')
    bound=427;reports=[];total_cuts=0;literal_controls=0
    inherited=json.loads((HERE.parent/'alternatives/frontier.json').read_text())
    for si,source in enumerate(SOURCES):
        need(strong(source),'strong source')
        hs=subset_counts(source);h=hs[-1];need(h==141,'source Hamilton count')
        full=(1<<len(source))-1
        q2=sum(hs[mask]*hs[full^mask] for mask in range(full+1))
        ear=attach(attach(source,0),full)
        ear_h=subset_counts(ear)[-1];literal_controls+=1
        need(strong(ear),'cyclic two-ear strongness')
        need(ear_h==h+q2,'two-ear identity')
        need(ear_h>=(len(source)+2)*h>bound,'two-ear threshold')
        layer=[(source,h,())];layers=[];reached=set()
        while layer:
            children=[];cuts=0
            for rows,parent_h,path in layer:
                hs=subset_counts(rows);need(hs[-1]==parent_h,'parent value')
                C=even_capacities(rows,hs);full=(1<<len(rows))-1
                # Check a nonuniform child independently of whether cutoff retains it.
                control=(19*len(rows)+7)&full
                need(cut_value(parent_h,C,control)==subset_counts(attach(rows,control))[-1],
                     'independent high/low cut control')
                literal_controls+=1
                for O in range(full+1):
                    value=cut_value(parent_h,C,O);cuts+=1
                    if O in (0,full):
                        need(value==parent_h,'uniform attachment control')
                        continue
                    if value>bound:continue
                    child=attach(rows,O)
                    need(strong(child),'mixed cut gives strong child')
                    need(subset_counts(child)[-1]==value,'retained literal child')
                    literal_controls+=1
                    need(value>parent_h,'strict growth in strong extension')
                    children.append((child,value,path+(O,)));reached.add(value)
            li=len(layers);origin=inherited['sources'][si]['layers'][li]
            want={(tuple(x['rows']),x['h'],tuple(x['path'])) for x in origin['children']}
            need(set(children)==want,'complete labelled frontier comparison')
            need(len(children)==len(origin['children']),'labelled multiplicity comparison')
            need(cuts==origin['cuts_examined'],'all cuts compared')
            layers.append({'parent_count':len(layer),'cuts_examined':cuts,
                           'children_count':len(children),
                           'count_histogram':dict(sorted(Counter(v for _,v,_ in children).items()))})
            total_cuts+=cuts;layer=children
        need(len(layers)==len(inherited['sources'][si]['layers']),'empty closure layer')
        need(425 not in reached,'425 absent')
        need((427 in reached)==(si==1),'427 distinguishes sources')
        reports.append({'source_rows':source,'h':h,'cyclic_two_ear_h':ear_h,
                        'two_ear_lower_bound':(len(source)+2)*h,'layers':layers,
                        'reached_counts':sorted(reached),
                        'all_labelled_children_match':True})
    return {'status':'PASS','bound':bound,'total_cuts_examined':total_cuts,
            'literal_controls':literal_controls,'sources':reports,
            'scope':'Independent even-directed-path cut polynomial for every labelled cut from all retained strong parents; exact frontier closure. All-order implication additionally uses the written strong-ear threshold proof.'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
