import json
from verify_seeds import circulate,walk_contract,fixed_contract
from verify_pathcovers import need
p=json.load(open(__file__.replace('check_r5_45.py','receipt.json')))
r=next(x for x in p['rows'] if x['n']==9);c,d=r['c'],r['d'];A=circulate(5,{1,2})
H=walk_contract(A,c);D=fixed_contract(A,c,d);a=5*81**5;e=3*81**2
need(H%a==0 and D%e==0,'orbit normalization');h,b=H//a,D//e;need((h-b)%2==0,'pair parity')
print(json.dumps(dict(n=45,H=H,d1=D,aut=a,centralizer=e,h=h,b=b,paired_black_degree=(h-b)//2),indent=2))
