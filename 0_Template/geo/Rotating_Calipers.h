
//Rotating_Calipers
#ifndef pii
#define pii pair<int, int>
#endif
class Rotating_Calipers
{
public:
    using vec_pii = vector<pair<int, int>>;
    vec_pii p;

    static int cross(pair<int, int> p1, pair<int, int> p2, pair<int, int> p0)
    {
        pii t1 = {p1.fi - p0.fi, p1.se - p0.se},
            t2 = {p2.fi - p0.fi, p2.se - p0.se};
        return t1.fi * t2.se - t1.se * t2.fi;
    }

    static int dis(const pii &p1,const pii &p2){
        return (p1.fi - p2.fi) * (p1.fi - p2.fi) + (p1.se - p2.se) * (p1.se - p2.se);
    };

public:
    
    Rotating_Calipers() {}

    Rotating_Calipers(vec_pii _A) {
        build(_A);
    }

    void build(const vec_pii & _A) {
        p = ConvexHull(_A);
    }

    static vec_pii ConvexHull(vec_pii A, int flag = 1)
    {
        int n = A.size();
        if (n <= 2) return A; 
        vec_pii ans(n * 2);
        sort(A.begin(), A.end());
        int now = -1;
        for (int i = 0; i < n; i++)
        { // 维护下凸包
            while (now > 0 && cross(A[i], ans[now], ans[now - 1]) < flag) now--;
            ans[++now] = A[i];
        }
        int pre = now;
        for (int i = n - 2; i >= 0; i--)
        { // 维护上凸包
            while (now > pre && cross(A[i], ans[now], ans[now - 1]) < flag) now--;
            ans[++now] = A[i];
        }
        ans.resize(now);
        return ans;
    }

    int getDiameter() {
        int j = 1,ans = 0;
        int m = p.size();
        p.push_back(p[0]);
        for(int i = 0;i < m;i ++)
        {
            while( cross(p[i+1],p[j],p[i]) > cross(p[i+1],p[j+1],p[i]) ) j = (j+1)%m;
            ans = max(ans, max( dis(p[i],p[j]) , dis(p[i+1],p[j]) ) );
        }
        p.pop_back();
        return ans;
    }

};
