from search_nine_offsets import data_all
from module_tools import *
import json
HERE=Path(__file__).resolve().parent
T=(6,12,24,33,107,7,47);n=len(T);hits=[]
for bits in range(1<<n):
 Z=tuple(T[i]|((1<<n) if bits>>i&1 else 0) for i in range(n))+(sum(1<<i for i in range(n) if not bits>>i&1),)
 h,ws,us=data_all(Z)
 if h<129:continue
 for O in range(1,1<<(n+1)):
  q=us[O];j=(q&-q).bit_length()-1 if q else 0
  if h-ws[O]==129 and q>>j==9:
   z={'rows':Z,'O':O,'D':129,'H':h,'U':q,'j':j,'orientation_mask':bits};hits.append(z);print(z,flush=True)
HERE.joinpath('offset_extension_receipt.json').write_text(json.dumps({'source':T,'target_D':129,'extensions':128,'hits':hits},indent=2))
print('END',len(hits))
