
const int INF = 0x3f3f3f3f

class PD//AC
{
public:
    class edge
    {
    public:
        int v, f, c, next;
        edge(int _v,int _f,int _c,int _next)
        {
            v = _v,f = _f,c = _c,next = _next;
        }
        edge() { }
    } ;

    void vecset(int value,vector<int> &arr)
    {
        for(int i = 0;i < arr.size();i ++) arr[i] = value;
        return;
    }

    class node
    {
    public:
        int v, e;
    } ;

    class mypair
    {
    public:
        int dis, id;

        bool operator<(const mypair &a) const { return dis > a.dis; }

        mypair(int d, int x)
        {
            dis = d;
            id = x;
        }
    };

    vector<int> head,dis,vis,h;
    vector<edge> e;
    vector<node> p;
    int n, m, s, t, cnt = 1, maxf, minc;

    PD(int _n,int _m,int _s,int _t)
    {
        n = _n, m = _m,s = _s,t = _t;
        maxf = 0, minc = 0;
        head.resize(n+2), dis.resize(n+2),vis.resize(n+2);
        e.resize(2);
        h.resize(n+2), p.resize(m+2);
    }

    void addedge(int u, int v, int f, int c)
    {
        e.push_back(edge(v,f,c,head[u]));
        head[u] = e.size()-1;
        e.push_back(edge(u,0,-c,head[v]));
        head[v] = e.size()-1;
    }

    bool dijkstra()
    {
        priority_queue<mypair> q;
        vecset(INF,dis);
        vecset(0,vis);
        dis[s] = 0;
        q.push(mypair(0, s));
        while (!q.empty())
        {
            int u = q.top().id;
            q.pop();
            if (vis[u])
                continue;
            vis[u] = 1;
            for (int i = head[u]; i; i = e[i].next)
            {
                int v = e[i].v, nc = e[i].c + h[u] - h[v];
                if (e[i].f && dis[v] > dis[u] + nc)
                {
                    dis[v] = dis[u] + nc;
                    p[v].v = u;
                    p[v].e = i;
                    if (!vis[v])
                        q.push(mypair(dis[v], v));
                }
            }
        }
        return dis[t] != INF;
    }

    void spfa()
    {
        queue<int> q;
        vecset(63,h);
        h[s] = 0, vis[s] = 1;
        q.push(s);
        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            vis[u] = 0;
            for (int i = head[u]; i; i = e[i].next)
            {
                int v = e[i].v;
                if (e[i].f && h[v] > h[u] + e[i].c)
                {
                    h[v] = h[u] + e[i].c;
                    if (!vis[v])
                    {
                        vis[v] = 1;
                        q.push(v);
                    }
                }
            }
        }
    }

    int pd()
    {
        spfa();
        while (dijkstra())
        {
            int minf = INF;
            for (int i = 1; i <= n; i++)
                h[i] += dis[i];
            for (int i = t; i != s; i = p[i].v)
                minf = min(minf, e[p[i].e].f);
            for (int i = t; i != s; i = p[i].v)
            {
                e[p[i].e].f -= minf;
                e[p[i].e ^ 1].f += minf;
            }
            maxf += minf;
            minc += minf * h[t];
        }
        return 0;
    }

    pair<int,int> get()
    {
        return {maxf,minc};
    }
};
