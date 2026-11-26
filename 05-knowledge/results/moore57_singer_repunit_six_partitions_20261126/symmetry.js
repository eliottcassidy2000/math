// Fiber-preserving automorphism test for both predecessor prism covers.
const fs=require("node:fs"), vm=require("node:vm"), path=require("node:path");
const code=fs.readFileSync(path.join(__dirname,"../moore57_sharp_prism_boundary_20261126/prism_boundary.js"),"utf8");
const cases=vm.runInNewContext(code.slice(0,code.indexOf("for(const e of [true,false])"))+"\n[splice(true),splice(false)]");
for(const gs of cases){
  const roots=[];
  for(let y=0;y<54;y++){
    const P=Array(54).fill(-1),q=[0];P[0]=y;
    let ok=true;
    for(const x of q){
      for(const g of gs){
        const j=g[x],z=g[P[x]];
        if(P[j]<0){P[j]=z;q.push(j);}
        else if(P[j]!==z){ok=false;break;}
      }
      if(!ok)break;
    }
    if(ok&&q.length===54&&new Set(P).size===54)roots.push(y);
  }
  if(roots.length!==1||roots[0]!==0)throw Error("unexpected centralizer");
  console.log("fiber-preserving sheet actions: identity only");
}
