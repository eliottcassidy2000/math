'use strict';
function perms(a){if(!a.length)return [[]];return a.flatMap((x,i)=>perms(a.filter((_,j)=>i!==j)).map(p=>[x,...p]));}
let receipt=[];
for(let n=3;n<=5;n++){
let pairs=[],ix=Array.from({length:n},()=>Array(n).fill(-1));for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){ix[i][j]=ix[j][i]=pairs.length;pairs.push([i,j]);}
let count=1<<pairs.length,taus=[];
function edge(t,i,j){return i<j?!!(t>>ix[i][j]&1):!(t>>ix[i][j]&1);}
for(let t=0;t<count;t++){let best=0;for(let s=0;s<1<<n;s++){let vs=[];for(let i=0;i<n;i++)if(s>>i&1)vs.push(i);if(vs.length<=best)continue;let good=true;for(let i=0;i<vs.length;i++)for(let j=i+1;j<vs.length;j++)for(let k=j+1;k<vs.length;k++){let a=edge(t,vs[i],vs[j]),b=edge(t,vs[j],vs[k]),c=edge(t,vs[i],vs[k]);if(a===b&&a!==c)good=false;}if(good)best=vs.length;}taus.push(n-best);}
let edges=0,slack={},ps=perms(Array.from({length:n},(_,i)=>i));for(let t=0;t<count;t++)for(const p of ps){let ok=true,path=0;for(let j=0;j<n-1;j++){if(!edge(t,p[j],p[j+1]))ok=false;path|=1<<ix[p[j]][p[j+1]];}if(!ok)continue;let u=t^((count-1)^path),s=taus[t]+taus[u]-Math.floor(n/3);if(s<0)throw Error('wiring');edges++;slack[s]=(slack[s]||0)+1;}let trans=count-1,path=0;for(let i=0;i<n-1;i++)path|=1<<ix[i][i+1];let u=trans^((count-1)^path);if(taus[u]!==Math.floor(n/3))throw Error('sharp');receipt.push({n,labelledTournaments:count,directedPathComplementInstances:edges,slackDistribution:slack,transitiveNeighborTau:taus[u]});
}console.log(JSON.stringify({status:'PASS',receipt},null,2));
