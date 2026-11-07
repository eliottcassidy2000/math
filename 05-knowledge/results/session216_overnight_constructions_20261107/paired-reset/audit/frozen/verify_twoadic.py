from inherited import profile,need
import json
x=list(map(int,open('../session216-carriers/BANK.txt').read().split()));bank=[];i=1
for _ in range(x[0]):n=x[i];bank.append(x[i+1:i+1+n]);i+=n+1
g=bank[125];mask=16;rows=[mask<<1]+[(g[j]<<1)|(not(mask>>j&1)) for j in range(len(g))]
p=profile(rows);need(p==(123,4,33,630,123,33,2760,672,1422),'congruence witness')
strong=(48,5,9,3,14,30);q=profile(strong);need(q[0]==33 and q[6]==360,'strong counterexample')
for v in range(6):
 reached=1<<v
 while True:
  nxt=reached
  for i in range(6):
   if reached>>i&1:nxt|=strong[i]
  if nxt==reached:break
  reached=nxt
 need(reached==63,'strong connectivity')
print(json.dumps({'status':'PASS','necessary_congruence_witness':{'rows':rows,'profile':p},'strong_even_div16_counterexample':{'rows':strong,'profile':q}},indent=2))
