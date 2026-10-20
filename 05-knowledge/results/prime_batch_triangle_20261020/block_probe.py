#!/usr/bin/env python3
"""Declared finite source scan; no universal-entry inference."""
from collections import Counter
import json
from separate_controls import odd_step, optimal_blocks

sources=list(range(3,2002,2))+[10087,35655,52527,77031,837799,63728127,670617279]
hist=Counter();records=[];max_blocks=0;max_depth=0
for n in sources:
    x=n;word=[]
    while x>=n:
        x,a=odd_step(x);word.append(a)
        if len(word)>1000:raise RuntimeError('declared first-exit route cap')
    rows=optimal_blocks(word,n)
    if rows[-1][0]>=n:raise RuntimeError('finest envelope must equal the paid endpoint')
    count=len(rows);hist[count]+=1
    if count>max_blocks or len(word)>max_depth or n>2001:
        records.append({'source':n,'depth':len(word),'clock':sum(word),
                        'endpoint':x,'minimum_blocks':count,'cuts':rows[-1][1],
                        'upper':str(rows[-1][0]),'word':word})
    max_blocks=max(max_blocks,count);max_depth=max(max_depth,len(word))
print(json.dumps({'status':'FINITE-EXACT; listed sources only',
                  'universe':{'all_odds':[3,2001],'extra':sources[1000:],'route_cap':1000},
                  'sources':len(sources),'block_histogram':dict(sorted(hist.items())),
                  'maximum_blocks':max_blocks,'maximum_first_exit_depth':max_depth,
                  'records':records},indent=2,sort_keys=True))
