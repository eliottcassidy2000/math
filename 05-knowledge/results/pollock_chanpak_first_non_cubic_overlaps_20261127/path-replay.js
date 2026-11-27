function pathCount(v,mask){
 const word=(2*v).toString(2), n=word.length;
 const edge=Array.from({length:6},()=>Array(6).fill(false));
 let k=0;
 for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){
   edge[i][j]=!(mask>>k++&1);edge[j][i]=!edge[i][j];
 }
 const rows=Array.from({length:64},()=>[]);
 for(let S=0;S<64;S++){
   function visit(T,first,last){
     if(T===S){rows[S].push([first,last]);return}
     for(let i=0;i<6;i++)if((S>>i&1)&&!(T>>i&1)&&(last<0||edge[last][i]))
       visit(T|1<<i,first<0?i:first,i);
   }
   visit(0,-1,-1);
 }
 const good=[1,1,0,0,1,0,1,0], g=Array.from({length:2},()=>Array(2).fill(0));
 for(let S=0;S<256;S++){
   const T=255^S;
   for(const [first,last] of rows[S&63])for(const [f,l] of rows[T&63]){
     const end = S&128?7:S&64?6:last;
     const begin = f>=0?f:T&64?6:T&128?7:-1;
     for(let leftOK=0;leftOK<2;leftOK++)
     for(let rightEmpty=0;rightEmpty<2;rightEmpty++)
       if((end<0?leftOK:!good[end])&&(!rightEmpty||begin<0||good[begin]))g[leftOK][rightEmpty]++;
   }
 }
 let result=0, full=(1<<n)-1;
 for(let S=0;S<=full;S++){
   const T=full^S;
   const leftOK=(!S||word[31-Math.clz32(S)]==="0")?1:0;
   if(T && word[31-Math.clz32(T & -T)]!=="1")continue;
   result+=g[leftOK][T===0?1:0];
 }
 return result;
}
for(const [v,h] of [[5185,8541511],[10259,16818421],[10415,17021085],[10337,17023573]]){
 const a=pathCount(v,177),b=pathCount(v,181); if(a!==h||b!==h+2)throw Error("path count");console.log(v,a,b);
}
