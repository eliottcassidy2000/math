#!/usr/bin/env python3
"""Literal adjacency reconstruction; no transfer or inverse imports."""
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261030'))
from audit_feedback import require,subset_dp


def adjacency(word):
 rows=[2,0];x=1
 for kind,cut in word:
  block=(8,1,3,6) if kind=='D' else (0,)
  require(kind in ('D','S') and 0<=cut<2**len(block),'generator encoding')
  old=len(rows);rows.extend([0]*len(block))
  for i in range(old):
   if i!=x:rows[i]|=sum(1<<(old+j) for j in range(len(block)))
  for i,mask in enumerate(block):
   for j in range(len(block)):
    if mask>>j&1:rows[old+i]|=1<<(old+j)
   if cut>>i&1:rows[x]|=1<<(old+i)
   else:rows[old+i]|=1<<x
 return rows


def main():
 folder=HERE/'inverse_witnesses';folder.mkdir(exist_ok=True);receipts=[]
 for file in ('diamond_inverse.json','diamond_inverse_103061.json'):
  for record in json.loads((HERE/file).read_text())['results']:
   if record['word'] is None:continue
   rows=adjacency(record['word']);N=len(rows);F=record['target']
   require(N==record['vertices'],'retained vertex clock')
   for i in range(N):
    require(not(rows[i]>>i&1),'loop')
    for j in range(i):require((rows[i]>>j&1)+(rows[j]>>i&1)==1,'pair orientation')
   direct=None
   if N<=13:
    direct=subset_dp(rows);require(direct==F,'literal subset-DP count')
   name=f'H{F}.txt'
   lines=[f'{N} {F}']+[''.join(str(row>>j&1) for j in range(N)) for row in rows]
   (folder/name).write_text('\n'.join(lines)+'\n')
   receipts.append({'target':F,'vertices':N,'rows':rows,'adjacency':name,'python_literal_count':direct,
                    'large_case_audit':'separate checked C++ adjacency DP' if direct is None else None})
 print(json.dumps({'status':'passed','witnesses':receipts},indent=2,sort_keys=True))


if __name__=='__main__':main()
