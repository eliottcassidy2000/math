const assert=require('node:assert/strict');
function counts(n,mask=0,flip=null){
 const all=(1<<n)-1, p=new Int32Array(n);
 for(let i=0;i<n;i++) for(let j=i+1;j<n;j++){
  let b=(i===0&&j===n-1)||(n%2&&i>=1&&j===i+(n-1)/2&&((mask>>(i-1))&1));
  if(flip&&i===flip[0]&&j===flip[1]) b=!b;
  if(b)p[i]|=1<<j;else p[j]|=1<<i;
 }
 const h=new Float64Array(all+1),dp=new Float64Array((all+1)*n);h[0]=1;
 for(let s=1;s<=all;s++)for(let I=s;I;I&=I-1){
  const b=I&-I,i=31-Math.clz32(b),t=s^b;let d=t?0:1;
  for(let J=t&p[i];J;J&=J-1)d+=dp[t*n+31-Math.clz32(J&-J)];
  h[s]+=dp[s*n+i]=d;
 }
 const f=(z)=>{let g=h;for(let k=1;k<z;k++){let temp=new Float64Array(all+1);for(let s=0;s<=all;s++){let t=s;for(;;t=(t-1)&s){temp[s]+=h[t]*g[s^t];if(!t)break;}}g=temp;}return g[all];};
 return {H:h[all],f};
}
function choose(n,k){let v=1;for(let i=1;i<=k;i++)v=v*(n-i+1)/i;return v;}
function parity(m,z,p){let v=0;for(let k=p;k<=m;k+=2)v+=choose(m,k)*z**(m-k);return v;}
function formula(n,i,j,z){
 if(i===0&&j===n-1)return -2*z*parity(n-2,z,1);
 if(i===0){i=n-1-j;j=n-1;}
 if(j===n-1){const a=i-1,b=n-i-2;return 2*(z**(a+2)*parity(b,z,1)-z**(b+1)*parity(a,z,0));}
 const b=j-i-1,s=n-b-4,o=(m)=>parity(m,z,1);
 let q=0;for(let t=1;t<=b;t+=2)q+=choose(b,t)*o(s+b-t);
 return z*((z+1)**s*(z+2)**b-(z-1)**s*(z-2)**b)+2*z**(s+3)*o(b)-2*z**(b+1)*o(s)+4*z*z*q;
}
let checks=0;
for(let n=3;n<=10;n++){
 const base=counts(n);const ties=[];
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  const f=counts(n,0,[i,j]);
  for(let z=1;z<=(n<=8?3:2);z++){assert.equal(f.f(z)-base.f(z),formula(n,i,j,z),[n,i,j,z].join(','));checks++;}
  if(f.H===base.H)ties.push([i,j]);
 }
 console.log('n='+n+' ties='+JSON.stringify(ties));
}
console.log('single-flip checks='+checks);
const eu=[1,1,1,2,5,16,61,272,1385,7936];
for(let r=1;r<=6;r++){
 let v=0;
 for(let m=0;m<1<<r;m++)v+=(-1)**m.toString(2).split('1').length*(-counts(2*r+1,m).H);
 assert.equal(v,2*(-1)**(r-1)*eu[r-1]);console.log('r='+r+' contrast='+v);
}

for (const flip of [null,[3,6]]) {
 const m=counts(7,0,flip),a=[];
 for(let j=1;j<=5;j++){
  let t=0;for(let k=1;k<=j;k++)t+=(-1)**(j-k)*choose(j,k)*m.f(k);
  a.push(t);
 }
 const low3=[141,136,38].reduce((t,c,i)=>t+c*a[i],0);
 const low4=[2917,4638,3274,984,96].reduce((t,c,i)=>t+c*a[i],0);
 const gap=m.f(2);
 assert.equal(m.H,33);assert.equal(gap,flip?660:612);
 console.log('module '+JSON.stringify({flip,a,gap,cube3:[low3,low3+gap],cube4:[low4,low4+gap]}));
}

const F=x=>(x**3-x)/6;
for(const [N,g] of [[104548,306],[113284,330],[13100266,306],[13930714,330]]){
 const f=new Set;for(let x=1;F(x)<=N;x++)f.add(F(x));
 let hits=0,placements=0;
 for(let a=1;a**3<=g;a++)if(g%a===0)for(let b=a;b*b<=g/a;b++)if(g%(a*b)===0){
  const c=g/a/b;
  for(let u=1;4*F(u)<=N;u++){
   const low=[u,u+a+b,u+a+c,u+b+c].reduce((x,y)=>x+F(y),0);
   if(low>N)break;
   placements++;
   assert.equal(low+g,[u+a,u+b,u+c,u+a+b+c].reduce((x,y)=>x+F(y),0));
   if(f.has(N-low))hits++;
  }
 }
 assert.equal(hits,0);console.log('common fifth: '+JSON.stringify({N,g,placements,hits}));
}
