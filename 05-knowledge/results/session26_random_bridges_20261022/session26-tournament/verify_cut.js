'use strict';
function H(w){let h=1n,a=0n,p=0n;for(const b of w){if(b==='1'){a=2n*a+1n;p=1n;}else {h+=a+p;a*=2n;p=0n;}}return h;}
function cap(w){let c=0n;for(let j=0;j<w.length;j++)if(w[j]==='1')for(let i=j+1;i<w.length;i++)if(w[i]==='0')c+=i-j===1?2n:1n<<BigInt(i-j-1);return c;}
function dp(w){let n=w.length+1,N=1<<n;let a=Array(n).fill(0);for(let i=0;i<n-1;i++){for(let j=i+1;j<n-1;j++)a[i]|=1<<j;if(w[i]==='1')a[n-1]|=1<<i;else a[i]|=1<<(n-1);}let d=Array.from({length:N},()=>Array(n).fill(0));for(let i=0;i<n;i++)d[1<<i][i]=1;for(let mask=1;mask<N;mask++)for(let i=0;i<n;i++)if(d[mask][i])for(let j=0;j<n;j++)if(!(mask>>j&1)&&(a[i]>>j&1))d[mask|(1<<j)][j]+=d[mask][i];return BigInt(d[N-1].reduce((s,x)=>s+x,0));}
let words=0,swaps=0,dps=0,corners=0,extrema=[];
for(let m=0;m<=12;m++){let min=null,max=null,minw=[],maxw=[];for(let k=0;k<2**m;k++){let w=m?k.toString(2).padStart(m,'0'):'';let h=H(w),c=cap(w);if(h!==c+1n)throw Error('cut '+w);words++;if(m<=8){if(dp(w)!==h)throw Error('DP '+w);dps++;}for(let j=0;j+1<m;j++)if(w.slice(j,j+2)==='01'){let v=w.slice(0,j)+'10'+w.slice(j+2);if(cap(v)<c+2n)throw Error('swap '+w);swaps++;}
const ones=[...w].flatMap((b,i)=>b==='1'?[i]:[]),zeros=[...w].flatMap((b,i)=>b==='0'?[i]:[]).reverse(),rows=ones.map(i=>zeros.filter(j=>i<j).length);let f=1n;
for(let r=0;r<rows.length;r++)for(let j=0;j<rows[r];j++){let hook=rows[r]-j+rows.slice(r+1).filter(t=>t>j).length;if(hook!==zeros[j]-ones[r])throw Error('hook '+w);f+=1n<<BigInt(hook-1);if(hook===1){f++;corners++;}}
if(f!==h)throw Error('Ferrers '+w);
if(m>=2&&w[0]==='1'&&w[m-1]==='0'){if(min===null||h<min){min=h;minw=[w];}else if(h===min)minw.push(w);if(max===null||h>max){max=h;maxw=[w];}else if(h===max)maxw.push(w);}}
if(m>=2)extrema.push({m,min:String(min),max:String(max),minw,maxw});}
console.log(JSON.stringify({status:'PASS',words,swaps,dps,corners,extrema},null,2));
