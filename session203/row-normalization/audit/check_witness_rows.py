"""Verify every archived base-999 row against explicit physical operations."""
from pathlib import Path
import json

def need(b,s):
    if not b: raise RuntimeError(s)

def pad(a,r,z):
    b=a[:]+[0];n=len(a)
    for i in range(n):
        if z and i==r:b[n]|=1<<i
        else:b[i]|=1<<n
    return b

def amalgam(a,r,b,t):
    result=a[:]+[0]*(len(b)-1)
    mapping={t:r};next_index=len(a)
    for j in range(len(b)):
        if j!=t:mapping[j]=next_index;next_index+=1
    for i in range(len(b)):
        for j in range(len(b)):
            if b[i]>>j&1:result[mapping[i]]|=1<<mapping[j]
    for i in range(len(a)):
        for j in range(len(b)):
            if i!=r and j!=t:result[i]|=1<<mapping[j]
    return result

p=Path(__file__).resolve().parent
cases=json.loads((p/'REVIEWED_BASE999_WITNESSES.json').read_text())['cases']
for c in cases:
    a=c['seed_rows'][:]
    for letter in c['suffix_word']:
        need(letter in 'OZ','word alphabet');a=pad(a,c['root'],letter=='Z')
    need(a==c['source_rows'],'source rows')
    need(amalgam([2,8,3,5],0,a,c['root'])==c['compiled_rows'],'compiled rows')
    need(0<=c['j']<=8 and c['s']==c['r']+3*c['j'],'actual endpoint shift')
    need(c['h']+c['L']-c['s']==c['B']-2*c['j'],'actual base shift')
    need(c['target']==3*c['B']+2*c['r']==c['literal_H'],'target arithmetic')
print(json.dumps({'status':'PASS','targets':[c['target'] for c in cases],
                  'source_orders':[len(c['source_rows']) for c in cases],
                  'child_orders':[len(c['compiled_rows']) for c in cases],
                  'scope':'Exact row construction and archived arithmetic; independent counts are checked by check.cpp.'},indent=2,sort_keys=True))
