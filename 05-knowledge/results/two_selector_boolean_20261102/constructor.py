"""Shortest actual Boolean inverse in the declared seven-vertex operation graph."""
from collections import deque
import argparse
import json
from experiment import rows, ham, path_cube_counts, remark, R, need

S=(0,3,4,1,2,5,6)


def construct(tile, target=17, selectors=2):
    need(0<=tile<32768,'tile is a 15-bit fixed-path chart')
    need(selectors in (0,1,2),'number of named selectors')
    counts=path_cube_counts(7,True)
    previous={tile:None};todo=deque([tile]);destination=None
    while todo:
        x=todo.popleft()
        if counts[x]==target:destination=x;break
        candidates=[('C',x^32767)]
        g=rows(7,x,True)
        for name,p in list(zip(('R','S'),(R,S)))[:selectors]:
            y=remark(g,p)
            if y>=0:candidates.append((name,y))
        for name,y in candidates:
            if y not in previous:previous[y]=(x,name);todo.append(y)
    result={'source_tile':tile,'source_rows':rows(7,tile,True),'target_H':target,
            'allowed_commands':['C']+list(('R','S')[:selectors]),
            'scope':'exact shortest word in this fixed-path command language only',
            'reachable':destination is not None}
    if destination is None:
        result['exhausted_component_size']=len(previous)
    else:
        trace=[];x=destination
        while previous[x] is not None:
            parent,name=previous[x];trace.append({'command':name,'source':parent,'target':x,
                                                 'H_before':counts[parent],'H_after':counts[x]})
            x=parent
        trace.reverse()
        order=list(range(7))
        for step in trace:
            step['source_vertex_order']=order.copy()
            p=R if step['command']=='R' else S if step['command']=='S' else tuple(range(7))
            step['selected_path_in_current_labels']=list(p)
            if step['command']!='C':order=[order[i] for i in p]
            step['target_vertex_order']=order.copy()
        for x in [tile]+[step['target'] for step in trace]:
            need(ham(rows(7,x,True))==counts[x],'literal endpoint-DP output check')
        result.update(word=[step['command'] for step in trace],trace=trace,
                      target_rows=rows(7,destination,True),target_vertex_order=order,
                      minimum_commands=len(trace))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tile',type=int)
    parser.add_argument('--target',type=int,default=17)
    parser.add_argument('--selectors',type=int,choices=(0,1,2),default=2)
    args=parser.parse_args()
    print(json.dumps(construct(args.tile,args.target,args.selectors),indent=2,sort_keys=True))
