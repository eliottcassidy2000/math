function Hdp(r){const n=r.length,F=(1<<n)-1,D=new Int32Array((F+1)*n);
for(let i=0;i<n;i++)D[(1<<i)*n+i]=1;
for(let s=1;s<=F;s++)for(let i=0;i<n;i++)if(D[s*n+i])for(let b=r[i]&(F^s);b;b&=b-1){const t=b&-b;D[(s|t)*n+31-Math.clz32(t)]+=D[s*n+i];}
return D;}
function H(r){const n=r.length,D=Hdp(r);return D.slice(((1<<n)-1)*n).reduce((x,y)=>x+y,0);}
function strong(r){let s=1,z=0;while(z!==s){z=s;for(let i=0;i<r.length;i++)if(r[i]&s)s|=1<<i;}return s===(1<<r.length)-1;}
function rows(n,t){const R=Array(n).fill(0);let k=0;for(let j=1;j<n;j++)for(let i=0;i<j;i++){if(j===i+1||(t>>k++)&1)R[i]|=1<<j;else R[j]|=1<<i;}return R;}
function rooted(r){const n=r.length,F=(1<<n)-1,D=Hdp(r),S=Hdp(r.map((a,i)=>a^F^(1<<i))),hs=new Int32Array(F+1);hs[0]=1;const u=Array(n).fill(0);for(let m=1;m<=F;m++)for(let i=0;i<n;i++)hs[m]+=D[m*n+i];for(let m=1;m<=F;m++)for(let i=0;i<n;i++)u[i]+=hs[F^m]*S[m*n+i];return [hs[F],u,D.slice(F*n)];}
function prime(x){if(x<2)return false;for(let i=2;i*i<=x;i++)if(x%i===0)return false;return true;}



const witnesses=[[5,[94,4,8,18,6,95,30],[[4,5],[0,4]],29],[11,[30,44,8,0,46,77,31],[[3,4],[0,2]],45],[13,[118,4,0,23,38,78,30],[[2,3],[0,2]],49],[17,[238,44,8,0,111,76,142,62],[[4,6],[3,4]],23],[19,[126,76,8,16,38,78,28],[[0,3],[0,2]],77],[23,[10,4,9,18,39,15,191,63],[[5,6],[0,7]],77],[29,[42,12,9,0,47,78,159,63],[[3,4],[0,7]],135],[31,[34,48,43,19,37,8,191,63],[[5,6],[1,2]],177],[37,[42,36,9,18,39,12,191,63],[[5,6],[0,6]],179],[41,[32,173,9,17,39,12,191,61],[[5,6],[0,1]],163],[43,[210,20,17,23,0,95,158,62],[[4,5],[2,3]],245],[47,[66,156,25,17,1,223,158,29],[[5,7],[4,5]],119],[53,[162,172,9,1,47,12,191,60],[[5,6],[3,4]],205],[59,[210,156,137,17,4,223,158,24],[[5,7],[4,5]],229],[61,[414,284,280,272,288,335,415,318,0],[[0,1],[7,8]],93],[67,[254,204,104,16,38,74,24,124],[[6,7],[0,3]],111],[71,[510,508,232,48,36,64,152,312,124],[[0,4],[0,7]],369],[73,[454,412,392,401,421,463,414,256,0],[[0,2],[6,8]],97],[79,[502,428,392,401,422,460,414,256,0],[[1,2],[6,8]],101],[83,[486,412,392,401,421,462,414,256,0],[[0,8],[3,8]],447],[89,[510,508,200,48,36,68,152,312,124],[[0,2],[0,6]],373]];
function toggled(r,arcs){
 r=r.slice();for(const [u,v] of arcs){r[u]^=1<<v;r[v]^=1<<u;}return r;
}
for(const [p,r,e,m] of witnesses){
 if(!prime(p)||H(r)!==p||H(toggled(r,[e[0]]))!==m||
    H(toggled(r,e))!==7*p||!(p<m&&m<7*p))throw Error("bad literal sequence "+p);
 console.log(p,r.length,m,H(toggled(r,[e[1]])),7*p);
}
const E=witnesses.find(([p])=>p===89),r=E[1],e=E[2];
const u=toggled(r,[e[0]]),v=toggled(r,[e[1]]);
if(H(u)!==373||H(v)!==373)throw Error("equal-middle square");
const pop=x=>x?1+pop(x&(x-1)):0;
function signature(r){const d=r.map(pop);return r.map((w,i)=>d[i]+":"+d.filter((_,j)=>w>>j&1).sort().join(",")).sort().join(";");}
if(signature(u)===signature(v))throw Error("need nonisomorphic middle certificates");
const C=[18,20,25,19,0];
function port(a){return [a,4];}
for(const [i,j,end] of [[0,2,13],[1,2,11]]){
 if(H(C)!==5||H(toggled(C,[port(i)]))!==9||H(toggled(C,[port(j)]))!==11||
 H(toggled(C,[port(i),port(j)]))!==end)throw Error("marginal twin");
}
console.log("same-host marginal twin: 5; 9 and 11; then 13 or 11.");
