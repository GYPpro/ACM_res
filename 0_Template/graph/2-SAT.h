#include <bits/stdc++.h>
using namespace std;


using pii = pair<int,int>;
class graph {
public:
    vector<vector<pii>> cnj;
    vector<int> dfn,scc,sz;
    int sc;
    int n;

    graph() { }

    graph(int _n) {
        n = _n;
        cnj.resize(n+1);
    };

    graph(vector<vector<int>> _cnj){

    }
private:
    reset(auto &t){t.clear();t.resize(n+1);};
public:
    void addOrdEdge(int u,int v,int w = 1) {

    }

    void addUnordEdge(int u,int v,int w = 1) {

    }

    //tarjin缩点
    void runSCC(){
        reset(dfn),reset(scc),reset(sz);
        vector<int> low(n+1),s,vst(n+1),sz(n+1);
        sc = 0;
        int tp = 0,dfncnt = 0;
        using itf = function<void(int)>;
        itf dfs = [&](int u) -> void {
            low[u] = dfn[u] = ++dfncnt;
            s.pb(u),vst[u] = 1;
            for(auto [v,_]:cnj[u])
            {
                if(!dfn[v])
                {
                    dfs(v);
                    low[u] = min(low[u],low[v]);
                } else if(vst[v]){
                    low[u] = min(low[u],dfn[v]);
                }
            }
            if(dfn[u] == low[u]){
                sc ++;
                while(s.back()!=u) {
                    scc[s.back()] = sc;
                    sz[sc]++;
                    vst[s.back()] = 0;
                    s.pop_back();
                }
                scc[s.back()] = sc;
                sz[sc]++;
                vst[s.back()] = 0;
                s.pop_back();
            }
        };
        for(int i = 1;i <= n;i ++)
            if(!dfn[i]) dfs(i);
    }

    /**
     * @param method:string 
     * Kruskal:O(m log m) for 稠密图
     * Prim:O((n+m) log n) for 稀疏图
     * @return vector<vector<int>> 生成树邻接表
     */
    vector<vector<pii>> runMST(string method = "Kruskal"){
        if(method == "Kruskal") {

        } else if(method == "Prim"){

        }
    }


    vector<pii> findBridge(string method = "tarjin"){
        
    }

    vector<int> findCut(string method = "tarjin"){

    }

    void eraseEdge(vector<pii> lst) {

    }

    void eraseVertice(vector<int> lst) {

    }
};

class two_SAT{
public:
    graph g;
    int n;
    two_SAT(int _n) {
        assert(_n%2 == 0);
        g = graph(_n);
        n = _n;
    }

    int tr(int p){
		if(p%2) return p+1;
		else return p-1;
    }

    void addDistinckPair(int u,int v){
        g.addOrdEdge(u,tr(v));
        g.addOrdEdge(v,tr(u));
    }
     
    /**
     * @brief 
     * 图的2-SAT：选择一组节点，对于所有选取的节点满足：
     * 任意x，节点编号2x和2x+1中必选一个
     * 所有的节点，其后缀一定在选取集合中。
     * @return pair<bool,vector<int>>
     * 无解返回{0,{}}
     * 有解返回{1,{无序的节点集合}}
     */
    pair<bool,vector<int>> run2SAT(){
        g.runSCC();
        for(int i = 1;i <= n;i += 2) if(g.scc[i] == g.scc[tr(i)]){
            return{0,{}};
        }
        vector<int> cel(g.sc+1);
        vector<int> res;
        for(int i = 1;i <= n/2;i ++)
        {
            READD:;
            if(cel[g.scc[i*2-1]]) res.pb(i*2-1);
            else if(cel[g.scc[i*2]]) res.pb(i*2);
            else {
                cel[min(g.scc[i*2-1],g.scc[i*2])] = 1;
                goto READD;
            }
        }
        return {1,res};
    }

};