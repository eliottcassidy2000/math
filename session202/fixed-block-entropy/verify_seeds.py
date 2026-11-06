import json
from verify import hamilton, need

def build(word):
    rows=[2,0];root=0
    for letter in word:
        n=len(rows);rows.append(0)
        for v in range(n):
            if letter!='O' and v==root: rows[n]|=1<<v
            else: rows[v]|=1<<n
        if letter=='J':root=n
    return tuple(rows)
results=[]
for n,word in [(49,'OJJOOJ'),(63,'OZZOOJ'),(343,'OOJOZJJOJ')]:
    rows=build(word);value=hamilton(rows);need(value==n,'seed graph count')
    results.append({'target':n,'word':word,'rows':rows,'literal_H':value,'order':len(rows)})
print(json.dumps({'status':'PASS','seeds':results},indent=2,sort_keys=True))
