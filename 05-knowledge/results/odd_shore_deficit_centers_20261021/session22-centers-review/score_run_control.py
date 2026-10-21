import itertools,json,pathlib
cases=0
for q in range(5,9):
 for c in range(1,5):
  for indices in itertools.combinations(range(q),c):
   cycle=[int(i in indices) for i in range(q)]*5
   for cut in range(len(cycle)):
    word=cycle[cut:]+cycle[:cut];run=maximum=0
    for x in word:
     run=run+1 if x else 0;maximum=max(maximum,run)
    if maximum>c:raise RuntimeError('cyclic-cut score-run bound failed')
    cases+=1
out={'scope':'binary score words only, not tournament enumeration','q_range':[5,8],'shore_multiplicities':[1,4],'repetitions':5,'cut_cases_checked':cases}
p=pathlib.Path(__file__).resolve().parent
(p/'score_run_results.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
