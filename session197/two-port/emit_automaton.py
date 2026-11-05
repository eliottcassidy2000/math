#!/usr/bin/env python3
from pathlib import Path
import argparse,json
from interface import shapes,transitions,brute
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
k=2;seed=brute([2,0],[0,1])
data={'k':k,'shapes':shapes(k),'seed':seed,'generators':[{'S':S,'drop':d,'transitions':transitions(k,S,d)} for S in range(4) for d in range(3)]}
js=json.dumps(data,indent=2)+'\n'
s='const int hc[]={' + ','.join(str(i) for i,a in enumerate(shapes(k)) if len(a)==1) + '};\n'
s+='const State seed={{' + ','.join(map(str,seed)) + '}};\n'
s+='const vector<vector<pair<int,int>>> tr[12]={\n'
for S in range(4):
 for d in range(3):s+='{'+','.join('{'+','.join('{'+str(j)+','+str(m)+'}' for j,m in row)+'}' for row in transitions(k,S,d))+'},\n'
s+='};\n'
for name,content in [('automaton.json',js),('automaton.inc',s)]:
 if args.check:
  if (ROOT/name).read_text()!=content:raise RuntimeError(name+' does not match')
 else:(ROOT/name).write_text(content)
print('Both serialized automata match the direct combinatorial generator.' if args.check else 'Wrote both serialized automata.')
