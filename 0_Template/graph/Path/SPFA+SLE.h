#include<bits/stdc++.h>
using namespace std;
#define ll long long 
const int inf =1e9+7;
////////////////////
struct p3008  //链式前向星
{
    int to;
    int nex;
    int v;
}a[300000];
int head[25006]={0},cnt=0;

int dis[25006]={0};  //距离

void add(int u,int v,int w)
{
    a[++cnt].nex=head[u];
    a[cnt].to=v;
    a[cnt].v=w;
    head[u]=cnt;
}
/////////////////////////////////
//SPFA，但是双端队列优化
//小的放队头，大的放队尾
void spfa(int t,int s)
{
    for(int i=1;i<=t;i++)dis[i]=inf;  //初始化距离
    deque<int> q;
    vector<int> vis(t+1),cou(t+1);

    q.push_front(s);  //起点入队

    dis[s]=0;  //起点初始化
    vis[s]++;

    while(!q.empty())
    {
        int u=q.front();  //取出队头
        q.pop_front();

        cou[u]++;

        if(cou[u]>t)return;  //遍历次数过大，出现负环

        for(int i=head[u];i;i=a[i].nex)  //遍历儿子
        {
            int v=a[i].to;
            int w=a[i].v;
            if(dis[v]>dis[u]+w)  //可松弛
            {
                dis[v]=dis[u]+w;  //松弛操作

                if(!vis[v])
                {
                    if(q.empty())q.push_back(v);       //特判！队空则随便入队
                    else if(dis[q.front()]>=dis[v])q.push_front(v);  //小的放队头
                    else q.push_back(v);  //大的放队尾
                    vis[v]++;
                }
            }
        }
        vis[u]--;  //取消入队标记
    }
}
void solve()
{
    int t,r,p,s;
    cin>>t>>r>>p>>s;
    for(int i=1;i<=r;i++)
    {
        int u,v;
        int w;
        cin>>u>>v>>w;
        add(u,v,w);
        add(v,u,w);
    }
    for(int i=1;i<=p;i++)
    {
        int u,v;
        int w;
        cin>>u>>v>>w;
        add(u,v,w);
    }

    spfa(t,s);  //建好边调用即可，t是点数，s是起点
    
    for(int i=1;i<=t;i++)
    {
        if(dis[i]>=inf)cout<<"NO PATH\n";
        else cout<<dis[i]<<"\n";
    }
}