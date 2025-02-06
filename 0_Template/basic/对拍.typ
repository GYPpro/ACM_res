==== 数据生成器：
#sourcecode[```cpp
#include <bits/stdc++.h>
using namespace std;

class generator{
public:
    // using std::mt19937 
    std::mt19937 mt;
    generator(){ mt.seed(std::random_device()()); };
    generator(int n) { mt.seed(n); };

    int randi(int l,int r) { return std::uniform_int_distribution<int>(l,r)(mt); }

    vector<int> randi(int n,int l,int r) {
        vector<int> rt;
        while(n--) rt.push_back(randi(l,r));
    }

    double randf(double l,double r) { return std::uniform_real_distribution<double>(l,r)(mt); }

    vector<double> randf(int n,double l,double r) {
        vector<double> rt;
        while(n--) rt.push_back(randf(l,r));
    }

    string rands(int l,bool ifa = 1,bool ifA = 0,bool ifd = 0) {
        string rt;
        while(l--) {
            int t = randi(0,61);
            if(t < 26) rt.push_back('a'+t);
            else if(t < 52) rt.push_back('A'+t-26);
            else rt.push_back('0'+t-52);
        }
        return rt;
    }

    vector<int> randp(int n) {
        vector<int> rt;
        for(int i = 1;i <= n;i ++) rt.push_back(i);
        std::shuffle(rt.begin(),rt.end(),mt);
        return rt;
    }

    vector<int> randt(int n) {
        vector<int> rt;
        for(int i = 2;i <= n;i ++) rt.push_back(randi(1,i-1));
        return rt;
    }

    vector<vector<int>> randg(int n,int m,bool forceconnected = 0) {
        vector<vector<int>> rt(n+1);
        vector<int> p = randp(n);
        for(int i = 2;i <= n;i ++) {
            int t = randi(1,i-1);
            rt[p[i]].push_back(p[t]);
            rt[p[t]].push_back(p[i]);
        }
        for(int i = n+1;i <= m;i ++) {
            int t = randi(1,n);
            rt[p[t]].push_back(p[i]);
            rt[p[i]].push_back(p[t]);
        }
        if(forceconnected) {
            vector<int> vis(n+1);
            std::queue<int> q;
            q.push(1);
            vis[1] = 1;
            while(q.size()) {
                int x = q.front();q.pop();
                for(auto y:rt[x]) if(!vis[y]) {
                    vis[y] = 1;
                    q.push(y);
                }
            }
            for(int i = 1;i <= n;i ++) if(!vis[i]) {
                int t = randi(1,n);
                rt[i].push_back(t);
                rt[t].push_back(i);
            }
        }
        return rt;
    }
} gc;

int main()
{
    freopen("G.A.in","w",stdout);
    
}

```]
==== 检查器：
#sourcecode[```cpp
#include <bits/stdc++.h>
using namespace std;

const string QUES_NAME = "E";

int system(string s) {
    return system(s.c_str());
}

using namespace std;
int main()
{
    cout << "compiling...\n";
    system("g++ " + QUES_NAME + ".ptc.cpp -o " + QUES_NAME + ".ptc.exe -std=c++2a -DFC");
    system("g++ " + QUES_NAME + ".std.cpp -o " + QUES_NAME + ".std.exe -std=c++2a -DFC");
    system("g++ " + QUES_NAME + ".gnc.cpp -o " + QUES_NAME + ".gnc.exe -std=c++2a -DFC");
    cout << "compile complete\n";
    int t = 1;
    while(++t){
        system(".\\" + QUES_NAME + ".gnc.exe");
        system(".\\" + QUES_NAME + ".ptc.exe");
        system(".\\" + QUES_NAME + ".std.exe");
        system("cls");
        if (system("fc " + QUES_NAME + ".A.std " + QUES_NAME + ".A.ptc")) {
            cout << "WA\n";
            system("pause");
            return 0;
        } else cout << "AC at test:" << t-1 << "\n";
    }
}


```]