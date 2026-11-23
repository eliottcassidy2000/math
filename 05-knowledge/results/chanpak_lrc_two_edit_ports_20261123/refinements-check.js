"use strict";
function eq(a, b) {
  if (JSON.stringify(a) !== JSON.stringify(b))
    throw Error(JSON.stringify({actual: a, expected: b}));
}
function tournament(r) {
  const n = r.length, all = (1 << n) - 1;
  const reach = r.map((row, i) => row | 1 << i);
  const dp = Array.from({length: all+1}, () => new Float64Array(n));
  for (let i=0;i<n;i++) {
    eq((r[i] >> i) & 1, 0);
    for (let j=0;j<i;j++) eq(((r[i]>>j)&1)+((r[j]>>i)&1), 1);
    dp[1<<i][i] = 1;
  }
  for (let k=0;k<n;k++) for (let i=0;i<n;i++)
    if (reach[i] & 1 << k) reach[i] |= reach[k];
  for (let s=1;s<=all;s++) for (let i=0;i<n;i++) if (dp[s][i]) {
    for (let z=r[i]&~s;z;z&=z-1) {
      const bit=z&-z, j=31-Math.clz32(bit);
      dp[s|bit][j]+=dp[s][i];
    }
  }
  return [dp[all].reduce((a,b)=>a+b,0), [...new Set(reach.map((x,i) =>
    reach.reduce((s,z,j) => s | (((x>>j&1)&&(z>>i&1)) ? 1<<j : 0), 0)))]];
}
const rows = [4094,4092,4088,4080,4064,3392,3200,3360,3648,3296,2048,0];
for (let mask=0;mask<4;mask++) {
  const v=rows.slice();
  [[7,11],[0,7]].forEach(([i,j],b) => {
    if (mask>>b&1) { v[i]^=1<<j; v[j]^=1<<i; }
  });
  const [h,scc]=tournament(v);
  eq([h, scc.filter(s=>s&(s-1))],
    [[13,[992]],[43,[4064]],[301,[1023]],[1099,[4095]]][mask]);
}
const prime = n => {
  if (n < 2) return false;
  for(let d=2;d*d<=n;d++) if (n%d===0) return false;
  return true;
};
for (let n=1;n<=9;n++) {
  const h=Array.from({length:1<<n}, (_,w) => {
    let a=1;
    for(let i=0;i<n;i++) if(w>>i&1) for(let j=i+1;j<n;j++) if(!(w>>j&1))
      a+=j-i===1 ? 2 : 1<<(j-i-1);
    return a;
  });
  const hits=[];
  for(let w=0;w<h.length;w++) if(h[w]>=5 && h[w]!==7 && prime(h[w])) {
    for(let i=0;i<n;i++) {
      if(h[w^1<<i]===7*h[w]) hits.push([w,i]);
      for(let j=i+1;j<n;j++)
        if(h[w^1<<i^1<<j]===7*h[w]) hits.push([w,i,j]);
    }
  }
  eq(hits.length, n===9 ? 2 : 0);
  if(n===9) {
    eq([h[152],h[153],h[156],h[157]], [37,217,81,259]);
    eq(hits.some(a=>JSON.stringify(a)===JSON.stringify([152,0,2])),true);
    const row=w=>Array.from({length:10},(_,i)=> i===9 ? w :
      ((511^((1<<(i+1))-1)) | ((w>>i&1) ? 0 : 512)));
    for(const w of [152,153,156,157]) eq(tournament(row(w))[0],h[w]);
  }
}
const affine=[[1,-40],[3,-68],[5,-148],[67,-2368],[19,80],[23,-80],[1,0],[0,1]];
const small=new Set();
for(let i=0;i<8;i++) for(let j=0;j<i;j++) {
  const k=Math.abs(affine[i][0]*affine[j][1]-affine[j][0]*affine[i][1]);
  eq(k>0,true);
  for(let d=2;d<=k;d++) if(k%d===0&&prime(d)) small.add(d);
}
const ex = [2,3,5,7,11,13,17,19,23,31,37,67,73,331,383,751,1049];
eq([...small].sort((a,b)=>a-b),ex);
const c = {5:3,7:4,11:5,13:4,17:6,19:6,23:6,31:6,37:5,67:6,73:6,
  331:6,383:6,751:6,1049:6};
for(let p=5;p<=1100;p++) if(prime(p)) {
  const n=Array.from({length:p},(_,q)=>q).filter(q=>affine.some(([a,b])=>(a+b*q)%p===0)).length;
  eq(n,c[p]??7);
}
for(let r of [1,2,4,5,7,8]) {
  let good=0;
  for(let q=0;q<9;q++) if((q-2*r)%3===0) {
    const d=r-40*q, num=5*r-148*q;
    eq(num%3,0);
    if(num%9===0) continue;
    eq([d,3*r-68*q,num/3,67*r-2368*q,19*r+80*q,23*r-80*q,r,q]
        .every(x=>x%3!==0),true);
    good++;
  }
  eq(good,2);
}
