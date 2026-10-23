"""Exhaustive one/two actual vertex marking probe on the explicit p23 pair."""
from pathlib import Path
from collections import Counter
from itertools import permutations,combinations
import argparse,json
HERE=Path(__file__).resolve().parent;CHECKS=0
def need(ok,why):
 global CHECKS
 CHECKS+=1
 if not ok:raise ValueError(why)

def refine(a,marks):
 n=len(a);colors=[len(marks)]*n
 for i,v in enumerate(marks):colors[v]=i
 profiles=[];rounds=0
 while True:
  profiles.append(tuple(sorted(Counter(colors).values())))
  keys=[(colors[v],tuple(sorted(colors[w]for w in range(n)if a[v][w])),
          tuple(sorted(colors[w]for w in range(n)if a[w][v])))for v in range(n)]
  ids={key:i for i,key in enumerate(sorted(set(keys)))};new=[ids[key]for key in keys]
  need(all(colors[i]==colors[j] or new[i]!=new[j] for i in range(n)for j in range(n)),'refinement never merges')
  if len(set(new))==len(set(colors)):break
  colors=new;rounds+=1
  need(rounds<n,'finite strict partition refinement')
 code=None
 if len(set(colors))==n:
  order=sorted(range(n),key=lambda v:colors[v])
  code=hex(sum(int(a[i][j])<<b for b,(i,j)in enumerate(combinations(order,2))))
 return rounds,profiles,code

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args();cases=[]
 for mask in (382046,431774):
  p=23;a=[[i!=j and ((mask>>((j-i)%p))&1)for j in range(p)]for i in range(p)]
  zero=refine(a,());need(zero[1]==[(23,)],'unmarked regular hostile')
  record={'mask':mask,'unmarked_profile':[23]}
  for count in (1,2):
   rounds=Counter();profiles=Counter();codes=Counter();initial_to_code=[]
   for marks in permutations(range(p),count):
    depth,history,code=refine(a,marks)
    need(history[-1]==(1,)*p and code is not None,'all actual marks become discrete')
    rounds[depth]+=1;profiles[str(history)]+=1;codes[code]+=1
    initial_to_code.append({'marks':list(marks),'rounds_to_discrete':depth,'canonical_word':code})
   need(sum(rounds.values())==(23 if count==1 else 506),'complete marked universe')
   record[f'{count}_marks']={'cases':sum(rounds.values()),'rounds_histogram':dict(sorted(rounds.items())),
     'stable_cell_sizes':[1]*23,'profile_history_histogram':dict(sorted(profiles.items())),
     'canonical_word_multiplicities':dict(sorted(codes.items()))}
  cases.append(record)
 need(set(cases[0]['1_marks']['canonical_word_multiplicities']).isdisjoint(cases[1]['1_marks']['canonical_word_multiplicities']),'single-mark incidence reconstructions separate')
 need(set(cases[0]['2_marks']['canonical_word_multiplicities']).isdisjoint(cases[1]['2_marks']['canonical_word_multiplicities']),'ordered-pair incidence reconstructions separate')
 p=31;S={x*x%p for x in range(1,p)};a=[[i!=j and (j-i)%p in S for j in range(p)]for i in range(p)]
 controls=[]
 for root in range(p):
  depth,history,code=refine(a,(root,));need(history[-1]==(1,15,15) and code is None,'Paley31 one-mark hostile');controls.append(depth)
 result={'status':'FINITE-EXACT all23 single marks and506 ordered distinct marks on each p23 host',
  'initial_colors':'ordered marked vertices have distinct colors; all other vertices share common rooted-five record',
  'cases':cases,'Paley31_control':{'roots':31,'stable_cell_sizes':[1,15,15],'rounds':dict(Counter(controls))},'checks':CHECKS}
 data=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
 if args.write:(HERE/'marked_probe.json').write_bytes(data)
 else:need((HERE/'marked_probe.json').read_bytes().replace(b'\r\n',b'\n')==data,'frozen marked receipt')
 print(json.dumps({'checks':CHECKS,'rounds':[{'mask':c['mask'],'one':c['1_marks']['rounds_histogram'],'two':c['2_marks']['rounds_histogram']}for c in cases],'Paley31_stable':[1,15,15]}))
if __name__=='__main__':main()
