
const seed=[1,1,1,1];
const mod=(x,n)=>(x%n+n)%n;
function one(v,c,n){let [h,e,u,g]=v;return (c=="P"?[h,0,2*u,g]:c=="K"?[h-e+u,g,2*u-g,g]:[h-e+u,h-e,h+u,h]).map(x=>mod(x,n));}
function go(w,n,v=seed){for(const c of w)v=one(v,c,n);return v.map(x=>mod(x,n));}
function eq(v,w){return v.every((x,i)=>x===w[i]);}
function check(ok,x){if(!ok)throw Error(JSON.stringify(x));}
function wall(h,g,n){
 if(n==2)return "";
 if(n==4)return ["","PK","PKPM","PPKMM"][((g-1)/2)*2+(h-1)/2];
 let z=mod((h-g)/2,n/2), w=z, r=0;
 if(w==0){w=1;r=Math.log2(n)-1;}
 else while(w%2==0){w/=2;r++;}
 let half=n/2,q=n/4;
 let v=K(mod(g-q,half),mod(w+q,half),w,half)+"K";
 if(go(v,n)[0]!==g)v+="K";
 check(go(v+"P".repeat(r)+"PM",n).join()===[h,g,h,g].join(),[v,h,g,n]);
 return v+"P".repeat(r)+"PM";
}
function M(h,e,g,n){
 let u=mod(h+g-e,n),state=[h,e,u,g],t=0;
 if(mod(u-h,2)){state=one(state,"M",n);t++;}
 for(let bit=2;bit<n;bit*=2){
  if(mod(state[2]-state[0],2*bit)>=bit){
   state=go("M".repeat(2*bit),n,state);t+=2*bit;
  }
 }
 check(state[0]===state[2]&&state[1]===state[3],state);
 return wall(state[0],state[3],n)+"M".repeat(mod(-t,2*n));
}
function P(h,u,g,n){return M(h,mod(h+g-u/2,n),g,n)+"P";}
function K(h,u,g,n){
 if(n==2)return "";
 let x=mod(u-g,n),y=mod(h-u,n),r=0;
 if(x===0){x=2;r=Math.log2(n)-1;}
 else while(x%4===0){x/=2;r++;}
 // uOld=(x+2g)/2 modulo n/2 is even.
 let b=mod((x+2*g)/2,n/2),a=mod(y+b-g,n);
 return P(a,b,g,n)+"K".repeat(r+1);
}
function census(n){
 let targets=new Set(),odd=Array.from({length:n/2},(_,i)=>2*i+1);let max=0;
 for(const h of odd)for(const g of odd){
  for(let e=0;e<n;e++){
   let w=M(h,e,g,n),v=[h,e,mod(h+g-e,n),g];
   check(eq(go(w,n),v),[n,v,w,go(w,n)]);targets.add(v.join());max=Math.max(max,w.length);
  }
  for(let u=0;u<n;u++){
   let v=[h,u%2?g:0,u,g],w=u%2?K(h,u,g,n):P(h,u,g,n);
   check(eq(go(w,n),v),[n,v,w,go(w,n)]);targets.add(v.join());max=Math.max(max,w.length);
  }
 }
 check(targets.size===n*n*(n-1)/2,[n,targets.size]);
 return [targets.size,max];
}
function rooted(w){
 let g=[[false]],a=0;
 for(const c of w){let b=g.length;g.push(Array(b+1).fill(false));for(let i=0;i<b;i++){g[i].push(c==="P"||i!==a);g[b][i]=!g[i][b];}if(c==="M")a=b;}
 return [g,a];
}
function iso(prefix){
 let [a,s]=rooted(prefix+"MK"),[b,t]=rooted(prefix+"PM"),n=a.length;
 function p(i){return i===n-1?n-2:i===n-2?n-1:i;}
 check(p(s)===t&&a.every((r,i)=>r.every((x,j)=>x===b[p(i)][p(j)])),prefix);
}
function run(){
 for(let n of [2,4,8,16,32])console.log("2-power",n,...census(n));
 for(let n of [2,4,8,16,32,64]){
  let seen=new Set([seed.join()]), q=[seed];
  for(let i=0;i<q.length;i++)for(let c of "PKM"){let v=one(q[i],c,n),s=v.join();if(!seen.has(s)){seen.add(s);q.push(v);}}
  check(seen.size===n*n*(n-1)/2,[n,seen.size]);
  console.log("independent dyadic BFS",n,seen.size);
 }
 let q=[""];for(let i=0;i<4;i++)q=q.concat(q.filter(s=>s.length===i).flatMap(s=>Array.from("PKM",x=>s+x)));
 for(const s of q)iso(s);
 console.log("rooted swaps",q.length);
 for(let m of [3,5,9,15]){
  let n=2*m,S=new Set([seed.join()]),W=[seed];
  for(let i=0;i<W.length;i++)for(let c of "PKM"){let v=one(W[i],c,n),t=v.join();if(!S.has(t)){S.add(t);W.push(v);}}
  let odd=Array.from({length:n/2},(_,i)=>2*i+1);
  for(let h of odd)for(let g of odd)for(let u=0;u<n;u++){
   check(S.has([h,u%2?g:0,u,g].join()),["PK",n,h,u,g]);
   check(S.has([h,mod(h+g-u,n),u,g].join()),["M",n,h,u,g]);
  }
  check(S.size===n*n*(n-1)/2,[n,S.size]);console.log("mixed even",n,S.size);
 }
}
run();
