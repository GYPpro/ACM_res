#include<bits/stdc++.h>
using namespace std;
#define ll long long
/////////////////////////////////////
struct p4779SPFA  //链式前向星
{
    int to;
    int nex;
    ll w;
}star[550000];
int head[100005]={0},cnt=0;
void add(int u,int v,ll w)
{
    star[++cnt].to=v;
    star[cnt].nex=head[u];
    star[cnt].w=w;
    head[u]=cnt;
}
///////////////////////////////

//SPFA,一款暴力的最短路算法
//常用于网络流、判负环、带负边的最短路
//暴力版的dj，每次松弛都让未入队的点入队，直到无法松弛为止
int n,m,s;
ll dis[100005]={0};   //距离
int times[100005]={0}; //遍历次数
bool vis[100005];  //入队情况
void SPFA()
{
    memset(vis, 0 ,sizeof(vis)); //清空vis
    for(int i=1;i<=n;i++)dis[i]=2147483647;  //距离初始化为无穷大
    dis[s]=0;  //起点归零
    queue<int> q;  //队列
    q.push(s);
    while(!q.empty())
    {
        int temp=q.front();
        q.pop();
        times[temp]++;
        vis[temp]=false;  //出队判定

        if(times[temp]>2+n)  //判负环  ，遍历次数过多说明出现了负环
        {
            times[s]=n+10;  //做标记
            break;
        }

        for(int i=head[temp];i;i=star[i].nex)  //遍历儿子
        {
            int v=star[i].to;
            if(dis[v]>dis[temp]+star[i].w)  //可松弛
            {
                dis[v]=dis[temp]+star[i].w;//松弛
                if(!vis[v])  //未入队则入队
                {
                    q.push(v);
                    vis[v]=true;
                }
            }
        }
    }

}
void solve()
{
    //spfa用于求最短路
    cin>>n>>m>>s;
    for(int i=1;i<=m;i++)
    {
        int u,v;
        ll w;
        cin>>u>>v>>w;
        add(u,v,w);
    }

    SPFA();//建边后调用即可

    if(times[s]>n)  //找负环
    {
        cout<<"minus circle!\n";
        return;
    }
    for(int i=1;i<=n;i++)cout<<dis[i]<<" ";
    cout<<"\n";
}