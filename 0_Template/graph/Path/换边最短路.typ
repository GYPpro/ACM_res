#import "@preview/codelst:2.0.1": sourcecode
// Display inline code in a box
#set text(font:("Times New Roman","Source Han Serif SC"))
#show raw.where(block: false): box.with(
  fill: luma(230),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
)
#show raw.where(block: true): block.with(
  fill: luma(240),
  inset: 10pt,
  radius: 4pt,
)
#show raw: set text(
    font: ("consolas", "Source Han Serif SC")
  )
#set page(
  // flipped: true,
  // background: [#image("background.png")]
  paper: "a4",
)
#set text(
    font:("Times New Roman","Source Han Serif SC"),
    style:"normal",
    weight: "regular",
    size: 13pt,
)
#show math.equation:set text(font:("New Computer Modern Math","Source Han Serif SC"))
#let nxtIdx(name) = box[ #counter(name).step()#counter(name).display()]
#set math.equation(numbering: "(1)")

#set page(
  paper:"a4",
  number-align: right,
  margin: (x:2cm,y:2.5cm),
  header: [
    #box(baseline:5pt)[#set text(
      size: 11pt,
    )
    #align(
      left+bottom, 
      [
        #smallcaps[ ]
        #h(1fr)#text(" ",fill:rgb("#898989"));
      ]
    )]
    #line(start: (0pt,-10pt),end:(483pt,-10pt))
  ],
  numbering: "1/1"
)
#set math.mat(delim: "[")
#set math.vec(delim: "[")

#set page(
  paper:"a4",
  number-align: right,
  margin: (x:2cm,y:2.5cm),
  header: [
    #box(baseline:5pt)[#set text(
      size: 11pt,
    )
    #align(
      left+bottom,
      [
        #smallcaps[Templetes]
        #h(1fr)#text("Github GYPpro/Acm_res",fill:rgb("#898989"));
      ]
    )]
    #line(start: (0pt,-10pt),end:(483pt,-10pt))
  ],
  numbering: "1/1"
)
///MAIN---MAIN///

给你一个n个点，m条边的无向图，每条边连接点u、v，并且有个长度w。
有q次询问，每次询问给你一对t、x,表示仅当前询问下，将t这条边的长度修改为x,请你输出当前1到n的最短路长度。

数据范围

2 ≤ n ≤ 2e5
1 ≤ m, q ≤ 2e5 
1 ≤ wi,xi ≤ 1e9

我们先不考虑修改，考虑给出指定的边，求出经过这条边的最短路
设这条边连接u,v,很容易想到这样的话只需要从1与n分别跑一遍Dijkstra,最短路长度就是 min(1到u的距离+边长+v到n的距离,1到v的距离+边长+u到n的距离)。

我们可以把修改分为以下几类：
   - 1.修改的边在1到n的最短路上，边的长度变大了。
   - 2.修改的边在1到n的最短路上，边的长度变小了。
   - 3.修改的边不在1到n的最短路上，边的长度变大了。
   - 4.修改的边不在1到n的最短路上，边的长度变小了。

很容易知道：
	对于 2， 原最短路长度-原边长+新边长 就是答案。
	对于 3， 原最短路长度 就是答案。
	对于 4， 由前面的思考得 min(原最短路长度，min(1到u的距离+新边长+v到n的距离,1到v的距离+新边长+u到n的距离) 就是答案。

都是O(1)得出答案


于是只剩下1了。
对于原问题，我们可以得到一些简单的结论。
令原最短路为E,其路径为E1,E2,E3……Ek.
对于每个不是E上的点u,1到u的最短路必定会使用E的一段前缀（可以为空）。
令这个前缀为Lu,1到u的最短路经过E1,E2,……E(Lu)。

同理可得u到n的最短路必定会使用E的一段后缀（可以为空）。
令这个后缀为Ru,u到n的最短路经过E1,E2,……E(Ru)。

特别地，对于E上的点u，令其L为E上连接自己的边的编号，R为自己连到下一个E上点的边的编号



关于如何求出每个点的L与R,我们可以先从n到1跑一遍Dijkstra,求出E的路径。
然后分别从1到n,n到1各跑一遍Dijkstra,过程中分别更新每个非E上点的L,R。

有了每个点的L,R后，考虑如何解决1。
1就是求min(E的长度-原边长+新边长，不经过修改的这条边的最短路长度)

问题从而转化成如何快速求出 不经过E上某条边的最短路长度

考虑使用线段树，树上的每段区间 l,r 的值表示 整个图不经过E上l到r这段的最短路长度

有了每个点的L,R我们很容易用图中非E上的边更新线段树

例如：一条边连接点u,v,经过这条边的最短路长度为len，

     我们可以把树上 Lu+1,Rv 的区间用len更新,比个min

     同样地，可以把 Lv+1,Ru 的区间用len更新

  由于代码中点1的l,r为0，且l=r,所以l要加1

  如果图方便，可以把l[1]=1,r[1]=0,每个点的l比r大1

  不懂的话可以自己画图比划比划

从而我们可以在O(logn)的时间内回答每个问题

总的复杂度为O((m+n+q)logn)
==== solution:
#sourcecode[
```cpp
#include<bits/stdc++.h>
using namespace std;
#define ll long long
const ll inf=1e18;
#define N 200505
#define mod 1000000007
int main()
{
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int t=1;
    void solve();
    // cin>>t;
    while(t--)solve();
}
//换边最短路
//题目：有一张n个点，m条带权无向边的图，现在有qu次询问，每次询问将边ti的权值修改为ch后的最短路，
//      询问对原图无影响

//      可用于删边、换边的最短路

//思想是先跑一遍dj，找到原始的最短路，用线段维护路上每条边删除后的最短路的距离
//此模板建议原封不动地抄
int n,m,qu;   //n点数，m边数，qu询问次数
struct Original_edge{int u,v;ll len;}OE[N];  //存原边
struct tabl{int to,id;ll len;};   //用于邻接表的结构体、存边
struct node{int id; ll len; bool operator < (const node &x) const {return x.len<len;}};//用于优先队列的结构体，存点和距离
vector<tabl> edge[N<<1];   //邻接表
priority_queue<node> q;    //用于dj的优先队列
int lstvis[N]={0};
int l[N]={0};
int r[N]={0};    //
int ind[N]={0};   //标记，ind[i]==j，表示边i在原始最短路上，且是路上的第j条边
int vis[N];       //用于dj的vis数组

ll t[N<<4]={0};       //这个是线段树

ll disT[N],disN[N];   //disT[i]从起点到点i的最短距离，disN[i]从终点到点i的最短距离
bool on_path[N];     //记录原始最短路，on_path[i]==true,表示点i在最短路上
int path_cnt=0;     //表示原始最短路上的边数

void dj(int p,ll dis[],int f)//起点， 对应数组， 操作编号
{
    for(int i=1;i<=n;i++)//初始化
    {
        vis[i]=0;
        dis[i]=inf;
    }
    dis[p]=0;
    q.push((node){p,0});  //加入起点
    while(!q.empty())
    {
        node temp=q.top();
        q.pop();
        int u=temp.id;
        ll w=temp.len;
        if(vis[u])continue;//跳过处理过的点
        vis[u]++;
        dis[u]=w;//记录距离
        for(auto re:edge[u])
        {
            int v=re.to;
            int id=re.id;
            ll tw=re.len;

            if(dis[v]<=w+tw)continue;
            dis[v]=w+tw;  //松弛
            lstvis[v]=id;  //记前驱边，表示点v从边id到达
            q.push((node){v,w+tw});

            if(f==1&&!on_path[v])l[v]=l[u];  //操作1，此时是从起点出发，需要记住前缀
            if(f==2&&!on_path[v])r[v]=r[u];  //操作2，此时是从终点出发，需要记住后缀
        }
    }
}
void trace()
{
    int u=1;    //u初始为起点
    on_path[u]=true;  //做好起点的初始化
    l[u]=r[u]=0;   

    for(int i=1;u!=n;i++) //第一次dj是从终点开始，所以这里要从起点开始，向终点找路
    {
        int e_id=lstvis[u];  //取前驱边
        ind[e_id]=i;      //给前驱标记

        u^=OE[e_id].u^OE[e_id].v;  //取前驱点
        on_path[u]=true;       //前驱点做标记
        l[u]=r[u]=i;      //做标记
        path_cnt++;      //路长增加
    }
}
void build(int le,int ri,int p)
{
    t[p]=inf;  //初始化为无穷大
    if(le==ri)return;
    int mid=(le+ri)>>1,lef=(p<<1),rig=lef|1;
    build(le,mid,lef);
    build(mid+1,ri,rig);
}
void update(int L,int R,int x,int y,int p,ll k)
{
    if(x>y)return;
    if(x<=L&&R<=y)  //区间修改
    {
        t[p]=min(t[p],k);
        return;
    }
    int mid=(L+R)>>1,lef=(p<<1),rig=lef|1;
    if(x<=mid)update(L,mid,x,y,lef,k);
    if(y>mid)update(mid+1,R,x,y,rig,k);
}
ll query(int L,int R,int x,int y,int p)
{
    ll ans=t[p];
    if(L==R)    //单点查询
    {
        return ans;
    }
    int mid=(L+R)>>1,lef=(p<<1),rig=lef|1;
    if(y<=mid)ans=min(ans,query(L,mid,x,y,lef));   //全程取最小值
    else ans=min(ans,query(mid+1,R,x,y,rig));
    return ans;
}
void solve()
{
    memset(on_path,false,sizeof(on_path));  //初始化

    cin>>n>>m>>qu;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        ll w;
        cin>>u>>v>>w;

        OE[i].u=u;   //记录原边，后续要用下标查询边
        OE[i].v=v;
        OE[i].len=w;

        edge[u].push_back({v,i,w});  //邻接表
        edge[v].push_back({u,i,w});
    }   

    dj(n,disN,0);   //先dj一次，找到原始最短路
    trace();
    dj(1,disT,1);   //分别对起点和终点进行dj，得到每条边的
    dj(n,disN,2);

    build(1,path_cnt,1);  //建立线段树，维护不经过某些边的最短路

    for(int i=1;i<=m;i++)
    {
        int u=OE[i].u,v=OE[i].v;
        
        ll w=OE[i].len;

        if(ind[i])continue;  //在路上的边不做更新

        //为线段树加入元素
        update(1,path_cnt,l[u]+1,r[v],1,disT[u]+w+disN[v]);//注意这里左边界要加一
        update(1,path_cnt,l[v]+1,r[u],1,disT[v]+w+disN[u]);
    }

    while(qu--)
    {
        int ti;
        ll ch;
        ll ans;
        cin>>ti>>ch;
        //询问会有两种大情况，换了最短路上的边，以及，换了最短路外的边

        if(ind[ti])   //若换了路上的边
        {
            ans=disT[n]-OE[ti].len+ch;//最理想的情况就是将原最短路的边给换一下
                                        //若新的边权更小，则不用比较了

            if(ch>OE[ti].len)    //若新的边权更大，则需要比较
            {
                ans=min(ans,query(1,path_cnt,ind[ti],ind[ti],1));   //查询不经过ti边的最短路
            }
        }
        else   //若换了最短路外的边
        {
            ans=disT[n];   //最理想是原最短路
                            //换的边比原来的边大，则无需考虑

            if(OE[ti].len>ch)  //若换的边更小，则需要比较
            {
                int u=OE[ti].u,v=OE[ti].v;

                ans=min(ans,min(disT[u]+ch+disN[v],disT[v]+ch+disN[u]));
                //新的最短路可能是   disT[u]+ch+disN[v] 表示从起点沿最优路径到达u，再经过边ti，再从v沿着最优路径到达终点
                //       也可能是   disT[v]+ch+disN[u]  表示从起点沿着最优路径到达v，再经过边ti，再从u沿着最优路径到达终点           
            }
        }
        cout<<ans<<"\n";
    }
}
```]