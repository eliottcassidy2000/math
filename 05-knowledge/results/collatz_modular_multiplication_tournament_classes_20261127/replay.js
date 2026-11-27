const check=x=>{if(!x)throw Error("failed replay")};
function gcd(a,b){while(b)[a,b]=[b,a%b];return a}
function mod(a,b){return (a%b+b)%b}
function order(g,m){if(gcd(g,m)!==1)return 0;let r=g,k=1;for(;r!==1;r=r*g%m)k++;return k}
function orb(g,m){const all=[],seen=new Set();for(let v=1;v<m;v++)if(!seen.has(v)){const C=[];for(let x=v;!seen.has(x);x=x*g%m){C.push(x);seen.add(x)}all.push(C)}return all}
function generator(p){for(let g=2;g<p;g++)if(order(g,p)===p-1)return g}
function getClasses(p){
 const inv2=(p+1)/2,a=3*inv2%p,g=a*a%p,h=order(g,p);
 if(!h||h%2===0)return [];
 const u=generator(p),powers=[1];for(let i=1;i<p-1;i++)powers.push(u*powers.at(-1)%p);
 const t=(p-1)/h,r=t/2, seen=new Set(), rows=[];
 for(let mask=0;mask<2**r;mask++){
   const v=Array.from({length:t},(_,i)=>i<r?(mask>>i&1):1-(mask>>(i-r)&1)).join("");
   const key=Array.from({length:t},(_,j)=>v.slice(j)+v.slice(0,j)).sort()[0];
   if(seen.has(key))continue;seen.add(key);
   rows.push(powers.filter((_,i)=>v[i%t]==="1").sort((a,b)=>a-b));
 }
 let burn=0;for(let j=0;j<t;j++){const d=gcd(t,j);if(r%d)burn+=2**(d/2);}
 check(rows.length*t===burn);return rows
}
function stats(p,D,doPath){
 const set=new Set(D),hist={};
 for(let t=1;t<p;t++){let k=D.filter(a=>set.has(mod(t-a,p))).length+Number(set.has(t));hist[k]=(hist[k]||0)+1}
 if(!doPath)return {D,hist};
 const E=Array.from({length:p},(_,x)=>D.map(y=>(x+y)%p).reduce((m,i)=>m|(1<<i),0)), full=(1<<p)-1;
 const F=new Float64Array((1<<p)*p);F[p]=1;
 for(let S=1;S<full;S+=2)for(let ends=S;ends;ends&=ends-1){
   const x=31-Math.clz32(ends&-ends),c=F[S*p+x];if(!c)continue;
   for(let todo=E[x]&~S;todo;todo&=todo-1){const b=todo&-todo,y=31-Math.clz32(b);F[(S|b)*p+y]+=c}
 }
 let h=0;for(let x=0;x<p;x++)h+=F[full*p+x];h*=p;
 check(Number.isSafeInteger(h));return {D,hist,H:h}
}
for(const p of [5,7,11,13,17,19,23,29,31]){
 const a=3*(p+1)/2%p,g=a*a%p, R=getClasses(p);
 console.log("p, ord(a), ord(g), classes",p,order(a,p),order(g,p),R.length);
 if(p<=19)for(const D of R){const row=stats(p,D,true);check(row.H%(p*order(g,p))===0);console.log(row);}
}
for(let m=5;m<=401;m+=2)if(m%3!==0){
 const a=3*(m+1)/2%m,g=a*a%m, O=orb(g,m);
 check(O.every(C=>!C.includes(m-C[0]))===(order(g,m)%2===1));
 if(m===25||m===161)console.log("Composite",m,order(g,m),O.length,
  O.reduce((s,C)=>{const k=String(C.length);s[k]=(s[k]||0)+1;return s},{}))
}

