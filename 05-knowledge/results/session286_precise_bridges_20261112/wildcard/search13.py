from scout import run
import json
for mask in range(32):
 S=[j if mask>>(j-1)&1 else 13-j for j in range(1,7)]
 ans=run(13,S);print(json.dumps(ans),flush=True)
 if ans['best'] and ans['best'][0]==2:break
