#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>   
using namespace std;
#define _USE_MATH_DEFINES
//Вывод вектора
void printVector(vector<double>a)
{
    for(int i = 0; i<a.size(); i++)
    {
        cout<< a[i]<< endl;
    }
}
//Вывод 2д вектора
void printMatrix(vector<vector<double> > a)
{
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < a[i].size(); j++)
        {
            cout << a[i][j] << "  ";
        }
        cout<<endl;
    }
}
//Проверка равенства векторов
void EqualityOfVectors(vector<vector<double> >g1, vector<vector<double> >g)
{
    double ag1;
    double ag;
    for(int i=0;i<g1.size();i++)
    {
        ag=sqrt(g[i][0]*g[i][0]+g[i][1]*g[i][1]+g[i][2]*g[i][2]);
        ag1=sqrt(g1[i][0]*g1[i][0]+g1[i][1]*g1[i][1]+g1[i][2]*g1[i][2]);
        if((ag1-ag)*(ag1-ag)>1.0e-8)
        {
            cout<<"MISTAKE";
        }
    }
}
//Увеличение 2д вектора
vector<vector<double> > Nlrg(vector<vector<double> > v,int v_index, int v_coordinates)
{ 
    v.resize(v_index);
    for (int j = 0; j < v.size(); j++)
    {
            v[j].resize(v_coordinates);
    }
    return v;
}
//Модуль вектора
double ModVector(vector<double> v, int n,int m)//Первая-Номер 1-й координаты (не с нуля) и последней
{
    double g=0.0;
    for(int i=n-1;i<m;i++)
    {
        g+=v[i]*v[i];
    }
    return sqrt(g);
}