#include <stdio.h>
#include <vector>
#include <string>
#include <string.h>
#include <iostream>
// #include <myLib\myLib.h>
using namespace std;

#define ll long long
const ll int maxn = (ll)(1e7 + 1);

vector<int> pri;
vector<bool> vis(maxn, false);

int main()
{
    // 
    freopen("D:\\Desktop\\Document\\Coding\\C++\\ACM\\0_list\\Prime_Number_List\\list.md", "w", stdout);
    ll step = (ll)maxn/100;
    for (ll int i = 2; i < maxn; i++)
    {
        // if ((i % step) == 0)
        // {
        //     UpdatePrc(maxn, i, 50, string("generating"));
        //     if(step <= i/100) step += step;
            
        // }
        if(i % step == 0) cout << i << " " << pri.size() << endl;
        if (!vis[i])
        {
            pri.push_back(i);
        }
        for (ll int j = 0; j < pri.size(); j++)
        {
            if ((ll)i * pri[j] >= maxn)
                break;
            vis[(ll)i * pri[j]] = true;
            if (i % pri[j] == 0)
                break;
        }
    }
    
    // cout << endl << "Generating Accomplished\n    Writing File\n";
    // step = 100;
    // int pow = 2;
    // ll int pop = 0;
    // ll int idx = 0;
    // ll int hcu = 0;
    // freopen("D:\\Desktop\\Document\\Coding\\C++\\ACM\\0_list\\Prime_Number_List\\list.md", "w", stdout);
    // cout << "# prime number list\n|rank|range|number|\n|:-:|:-:|:-|\n";
    // cout << "| 1e" << pow << " | " << step+step*(pop-1) << "~" << step+step*pop-1 << " | ";
    // for(;idx < pri.size();idx ++)
    // {
    //     hcu ++;
    //     if(pri[idx] > step+step*pop) {
    //         pop ++;
    //         hcu = step / 800;
    //         if(pop == 10){
    //             step *= 10;
    //             pow ++;
    //             pop = 1;
    //         }
    //         cout << " | ";
    //         cout << endl;
    //         cout << "| 1e" << pow << " | " << step+step*(pop-1) << "~" << step+step*pop-1 << " | ";
    //     }
    //     if(idx <= 1e3) cout << pri[idx] << " ";
    //     else if(step % hcu == 0) cout << pri[idx] << " ";
    // }
    // cout << "|\n";
    fclose(stdout);
}