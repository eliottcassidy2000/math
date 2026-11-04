import json,pathlib
P=pathlib.Path(__file__).resolve().parent
bank=json.loads((P/'bank.json').read_text())
with (P/'independent_queries.txt').open('w') as f:
    for g in bank:
        n=len(g['rows'])
        for a in g['masks']:
            for c in g['masks']:
                r=g['rows'][:]+[1<<(n+1),0]
                for i in range(n):
                    if a>>i&1:r[i]|=1<<n
                    else:r[n]|=1<<i
                    if c>>i&1:r[i]|=1<<(n+1)
                    else:r[n+1]|=1<<i
                print(n+2,*r,file=f)
