from screen import *
if __name__=='__main__':
 bank=load_bank()
 profiles={}
 for r in bank:
  if r['n']>8:continue
  h,w=taxmat(np.array(r['rows'],dtype=np.int64));c=int(np.sum(w))//2
  profiles.setdefault((int(h),c),r)
 out=[]
 for (h,c),r in profiles.items():
  for (g,d),s in profiles.items():
   if h*g+c*g-h*d==1 and (h>1 or g>1):out.append({'A':r['rows'],'B':s['rows'],'profile_A':[h,c],'profile_B':[g,d]})
 print('profiles',len(profiles),'nontrivial_pairs',len(out),out[:20])
 (HERE/'module_pairs.json').write_text(json.dumps(out,indent=2))
