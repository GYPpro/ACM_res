
//比ISAP和DINIC快一点，但不稳定 O(根(m)*n^2)
//思想是从源点开始，给每条边都推入满的流量，并设定每个点都可以存储一定的流量ei，称超额量
//且节点会伺机将超额量推向深度低的节点
int n,m,s,t;
int vis[1500]={0};
int head[1500]={0},cnt=1;
int h[1500]={0},gap[1500]={0};
int sum=0;
int e[1500]={0};
struct p4722
{
    int to,nex;
    int w;
}a[540005];
struct node
{
    bool operator() (int x,int y)const{return h[x]<h[y];}
};
priority_queue<int,vector<int>,node> q; 
queue<int> que;
void add(int u,int v,int w)   //链式前向星
{
    a[++cnt].to=v;
    a[cnt].nex=head[u];
    a[cnt].w=w;
    head[u]=cnt;
    a[++cnt].to=u;
    a[cnt].nex=head[v];
    a[cnt].w=0;
    head[v]=cnt;
}
void bfs_rebel()  //深度标签初始化，从汇点广搜
{
    que.push(t);
    h[t]=0;
    while(!que.empty())
    {
        int temp=que.front();
        que.pop();
        gap[h[temp]]++;
        for(int i=head[temp];i;i=a[i].nex)
        {
            int v=a[i].to;
            if(h[v]==0&&v!=t&&a[i^1].w)
            {
                h[v]=h[temp]+1;
                que.push(v);
            }
        }
    }
}
void push(int u)   //推流
{
    int i;
    for(i=head[u];i;i=a[i].nex)
    {
        int v=a[i].to;
        if(a[i].w<=0||h[u]!=h[v]+1)continue;//跳过0容边和深度不递减的点

        int f=min(e[u],a[i].w); //取残留容量和超额量中的小者

        e[v]+=f;    //推流，更新
        e[u]-=f;
        a[i].w-=f;
        a[i^1].w+=f;

        if(!vis[v]&&v!=t&&v!=s)  //子节点返回队列
        {
            vis[v]=1;
            q.push(v);
        }

        if(e[u]==0)break;   //超额量分发完毕退出函数
    }
}
void rebel(int u)   //重贴标签
{
    int i;
    h[u]=inf;   //标签初始化
    for(i=head[u];i;i=a[i].nex)
    {
        int v=a[i].to;
        if(a[i].w&&h[v]+1<h[u])h[u]=h[v]+1;  //找到儿子中的最小者，使u下一次恰好可以给它推流
    }
}
void hlpp()   //核心函数
{   
    h[s]=n;e[s]=inf; //源点s深度为n，超额量为无穷
    for(int i=head[s];i;i=a[i].nex)  //源点不入队，因此在外处理
    {
        int v=a[i].to;
        if(int f=a[i].w)
        {
            a[i].w-=f;     //更新网络
            a[i^1].w+=f;
            e[s]-=f;
            e[v]+=f;
            if(v!=s&&v!=t&&!vis[v])  //注意汇点和重复点不可入队
            {
                q.push(v);
                vis[v]=1;
            }
        }
    }
    while(!q.empty()) //计算最大流
    {
        int u=q.top();
        q.pop();
        vis[u]=0;  //出队取消标记
        push(u);//推流
        if(e[u]<=0)continue;  //无超额流就跳过，不入队
        gap[h[u]]--;//准备更新深度
        if(!gap[h[u]])  //本深度无其他点，则比u高的点都不可能再给u推流，可以把它们放在n+1处
        {
            for(int v=1;v<=n;v++)
            {
                if(v!=s&&v!=t&&h[v]>h[u]&&h[v]<n+1)
                {
                    h[v]=n+1;
                }
            }
        }
        rebel(u);  //重贴标签
        gap[h[u]]++;
        q.push(u);//入队
        vis[u]=1;
    }
    sum=e[t];
}