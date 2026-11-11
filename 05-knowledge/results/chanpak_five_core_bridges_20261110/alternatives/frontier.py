"""Complete bounded strong extension frontiers of two labelled H141 sources."""
import argparse
from collections import Counter
import json
import check as a
from two_birth import all_birth_counts

def strong(rows):
    full=(1 << len(rows))-1
    for root in range(len(rows)):
        seen=1 << root
        while True:
            fresh=seen
            for v,row in enumerate(rows):
                if seen >> v & 1:
                    fresh |= row
            if fresh == seen:
                break
            seen=fresh
        if seen != full:
            return False
    return True

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fail',action='store_true')
    args=parser.parse_args()
    if args.fail:
        a.need(False,'intentional frontier hostile')
    bound=427
    report={'bound':bound,'scope':'All strong supertournaments of either specified source with H<=427, using the proved two-ear cost bound.'}
    sources=[]
    total_control=0
    for tail in ((2,3,1),(2,2,1)):
        source=a.ready_source(tail,2)
        hp,_=a.dp(source)
        full=(1 << len(source))-1
        h=hp[-1]
        q2=sum(hp[mask]*hp[full^mask] for mask in range(full+1))
        cyclic=a.attach(a.attach(source,0),full)
        a.need(a.dp(cyclic)[0][-1]==h+q2,'cyclic two-ear count')
        a.need((len(source)+2)*h>bound,'all later two-ear steps exceed bound')
        layer=[{'rows':source,'h':h,'path':[]}]
        layers=[]
        while layer:
            children=[]
            examined=0
            for parent in layer:
                rows=parent['rows']
                counts=all_birth_counts(rows)
                examined+=len(counts)
                a.need(counts[0]==counts[-1]==parent['h'],'uniform birth controls')
                # One nonuniform direct recount per parent, independent of cutoff.
                control=(37*len(rows)+11) % len(counts)
                a.need(a.dp(a.attach(rows,control))[0][-1]==counts[control],'quadratic cutoff direct control')
                total_control+=1
                for out,value in enumerate(counts):
                    if value>bound:
                        continue
                    child=a.attach(rows,out)
                    if not strong(child):
                        continue
                    a.need(value==a.dp(child)[0][-1],'retained frontier independent DP')
                    a.need(value>parent['h'],'strong extension strict count')
                    children.append({'rows':child,'h':value,'path':parent['path']+[out]})
            layers.append({'extension_number':len(layers)+1,'parent_count':len(layer),'cuts_examined':examined,
                           'children_count':len(children),'count_histogram':dict(sorted(Counter(x['h'] for x in children).items())),
                           'children':children})
            layer=children
        expected=[25,4,0] if tail==(2,3,1) else [28,4,0]
        a.need([x['children_count'] for x in layers]==expected,'finite frontier closure sizes')
        reached={node['h'] for layer in layers for node in layer['children']}
        a.need(425 not in reached,'425 actual-source obstruction')
        a.need((427 in reached)==(tail==(2,2,1)),'427 distinguishes the sources')
        sources.append({'tail':tail,'rows':source,'h':h,'ordered_slot_profile_at2':q2,
                        'first_cyclic_two_ear_h':h+q2,'universal_two_ear_lower_bound':(len(source)+2)*h,
                        'layers':layers,'reached_counts':sorted(reached),'targets':{'425':False,'427':427 in reached}})
    report['sources']=sources
    report['independent_parent_controls']=total_control
    print(json.dumps(report,sort_keys=True,indent=2))

if __name__=='__main__':
    main()
