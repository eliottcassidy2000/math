#!/usr/bin/env python3
"""Fast witness/audit verifier; exhaustive minimal-order replay is separate."""
from pathlib import Path
from math import isqrt
import json
D=Path(__file__).resolve().parent

def check(p,msg):
 if not p:raise RuntimeError(msg)
R=json.loads((D/'COLOR_FIRST_FAILURE_CENSUS.json').read_text());W=json.loads((D/'birth_witnesses_25_80.json').read_text());base=json.loads((D/'receipt.json').read_text())['rows'];st={r['n']:bool(r['count']) for r in base}
for n in range(33,37):
 p=W[str(n)];check(len(p)==n and set(p)==set(range(1,n+1)) and all(isqrt(a+b)**2==a+b for a,b in zip(p,p[1:])),'birth witness');st[n]=True
check([r['n'] for r in R]==list(range(25,37)),'census range')
check([(r['path_count'],r['word_labels'],r['max_word_fibre']) for r in R[-4:]]==[(669,669,1),(4041,4041,1),(17175,17175,1),(12960,12957,2)],'recorded census')
check(all(r['path_count']==r['word_labels'] for r in R[:-1]),'earlier collision')
H=R[-1]['first_collision'];P=H['paths'];check(len(P)==2 and P[0]!=P[1],'distinct witnesses')
words=[]
for p in P:
 check(len(p)==36 and set(p)==set(range(1,37)),'vertex coverage');check((p[0],p[-1])==(6,32),'endpoints')
 check(all(isqrt(a+b)**2==a+b for a,b in zip(p,p[1:])),'square edges')
 word=[0 if st[max(a,b)] else 1 if isqrt(2*max(a,b)-1)-isqrt(max(a,b))==2 else 2 for a,b in zip(p,p[1:])];words.append(word)
check(words[0]==words[1]==H['word'],'word mismatch')
X=[{tuple(sorted(e)) for e in zip(p,p[1:])} for p in P];diff=sorted(X[0]^X[1]);check(diff,'same edge set')
print(json.dumps(dict(status='Explicit36hostileverified; first-order claim additionally depends on exhaustive census replay/audit',endpoints=[6,32],word=[c+1 for c in words[0]],paths=P,symmetric_difference=diff,recorded_census_last_four=R[-4:]),indent=2))
