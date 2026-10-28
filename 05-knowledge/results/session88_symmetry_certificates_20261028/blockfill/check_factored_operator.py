#!/usr/bin/env python3
"""Exact application from one immutable named row per block."""
import json
import check_blockfill as b


def factors(result):
    # Each row is sorted and includes its own state. The least support state
    # occurs exactly once, so factor extraction itself is linear in |X|.
    return tuple(row for source,row in result['Q'].items() if source==row[0][0])


def apply(blocks,vector):
    output={}
    for row in blocks:
        value=sum(weight*vector[state] for state,weight in row)
        for state,weight in row:
            output[state]=value
    return output


def main():
    parents=[()]; results=[]
    for n in range(8):
        result=b.build(parents,n,'blockfill'); blocks=factors(result)
        keys=sorted(result['states']); vector={k:((i*37)%101)-50 for i,k in enumerate(keys)}
        direct={i:sum(w*vector[j] for j,w in row) for i,row in result['Q'].items()}
        fast=apply(blocks,vector)
        b.need(fast==direct,'factored operator mismatch')
        twice=apply(blocks,fast)
        b.need(twice=={k:(n+1)*v for k,v in fast.items()},'factored power mismatch')
        b.need(sum(len(row) for row in blocks)==len(keys),'factor support overlap/gap')
        b.need(len({id(row) for row in result['Q'].values()})==len(blocks),'immutable block alias pattern')
        results.append({'order':n+1,'pointed_states':len(keys),'stored_distinct_immutable_rows':len(blocks),
                        'stored_factor_entries':sum(len(row) for row in blocks),
                        'explicit_sparse_Q_entries':sum(len(row) for row in result['Q'].values()),
                        'exact_Qv_and_Q_squared_v_agree':True})
        parents=list(result['bank'])
    return {'status':'FINITE-EXACT','checks':results,
            'complexity':'Qv requires one multiply-add per state plus one broadcast assignment per state. Bit complexity depends on integer sizes. Constructing Q and transporter witnesses has separate costs.'}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
