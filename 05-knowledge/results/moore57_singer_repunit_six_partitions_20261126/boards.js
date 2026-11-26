// Rebuild the input boards from the pinned predecessor; optional: --check to compare copies.
const fs = require("node:fs"), vm = require("node:vm"), path = require("node:path");
const frame = path.join(__dirname, "../moore57_sharp_prism_boundary_20261126/prism_boundary.js");
let src=fs.readFileSync(frame,"utf8");
const marker=src.indexOf("for(const e of [true,false])");
if(marker<0)throw Error("boundary replay entry point not found");
src=src.slice(0,marker);
const tables=vm.runInNewContext(src+`
const ans=[];
for(const e of [true,false])for(const t of [true,false]) {
  const c=prism(e), b=boundary(c,t,e?(t?401:411):(t?601:611)).output;
  const taken=new Set();
  for(const r of b)for(let i=0;i<6;i++)if(r[i]<54)
    for(let j=i+1;j<6;j++)if(r[j]<54)taken.add(c.code(i,r[i],j,r[j]));
  const data=[];
  for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)for(let x=0;x<n;x++){
    let mask=0n;
    for(let y=0;y<n;y++)if(c.allowed(i,x,j,y)&&!taken.has(c.code(i,x,j,y)))
      mask|=1n<<BigInt(y);
    data.push(mask.toString(16));
  }
  ans.push(data.join(" "));
}
ans;
`);
for(const [i,val] of tables.entries()){
  const fn=path.join(__dirname,`b${i}.txt`);
  if(process.argv.includes("--check")){
    if(val!==fs.readFileSync(fn,"utf8").trim())throw Error("board mismatch");
  }else fs.writeFileSync(fn,val+"\n");
}
console.log(tables.length,"exact boundary boards");
