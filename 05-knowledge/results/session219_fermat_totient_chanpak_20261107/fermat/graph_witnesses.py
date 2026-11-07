import json
from pathlib import Path
from verify_fermat_grammar import state,F,require
ROOT=Path(__file__).parent

def graph(w):
    # Literal birth model recovered from Session200 independent verifier.
    rows=[2,0];root=0
    for c in w:
        n=len(rows);rows.append(0)
        for v in range(n):
            if c!='O' and v==root:rows[n]|=1<<v
            else:rows[v]|=1<<n
        if c=='J':root=n
    return rows,root
prefix='OOZJOZJJOJJOZ'
words={'short_root_change_escape':prefix+'O'*20,
       'five_run_escape':'O'*31+'Z'*15+'O'*9+'Z'*5+'O'*9}
items=[]
for name,w in words.items():
    rows,root=graph(w+'J');n=len(rows)
    for i in range(n):
        require(not(rows[i]>>i&1),'no diagonal')
        for j in range(i):require(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1,'complete tournament')
    item={'name':name,'source_word':w,'closing_letter':'J','n':n,'root':root,'H':str(F(state(w))),
          'row_encoding':'integer bit j of row i is 1 iff edge i->j',
          'outgoing_rows_decimal':[str(x) for x in rows]}
    items.append(item)
require(int(items[1]['H'])==1182825546590390812645==176530139331685*6700417,'five-run value')
require(int(items[1]['H'])%641==275,'five-run prime separation')
(ROOT/'graph_witnesses.json').write_text(json.dumps(items,indent=2)+'\n')
# Independent subset-DP inputs; compare literal graph path count to profile formulas.
checks=[]
for k in range(7):
    w=prefix[:k];checks.extend([(w,state(w)[0]),(w+'J',F(state(w)))])
checks.extend([(prefix,state(prefix)[0]),(prefix+'J',F(state(prefix))),
               (prefix+'OJ',4909+2*6013),(prefix+'OOJ',4909+4*6013)])
with (ROOT/'literal_dp_inputs.txt').open('w') as f:
    f.write(str(len(checks))+'\n')
    for w,H in checks:
        rows,_=graph(w);f.write(f'{len(rows)} {H}\n'+' '.join(map(str,rows))+'\n')
print('Generated 36- and 72-vertex literal witnesses; 18 small/prefix DP checks through order18.')
