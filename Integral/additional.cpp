#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>   
using namespace std;
#define _USE_MATH_DEFINES
void ZSI(vector<vector<double> >a, vector<vector<double> >b)
{
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            if((a[i][j]+a[i][j+3]-b[i][j]-b[i][j+3])*(a[i][j]+a[i][j+3]-b[i][j]-b[i][j+3])>1.0e-8)
            {
                cout<<"MISTAKE  ";
            }
        }
    }
}
void ZSE(vector<vector<double> > a, vector<vector<double> > E1)
{
    double E1_2, E_2;
    for (int i=0;i<E1.size();i++)
    {
        E1_2=E1[i][0]*E1[i][0] + E1[i][1]*E1[i][1] + E1[i][2]*E1[i][2] + E1[i][3]*E1[i][3] + E1[i][4]*E1[i][4] + E1[i][5]*E1[i][5];//Верно
        E_2=a[i][0]*a[i][0] + a[i][1]*a[i][1] + a[i][2]*a[i][2] + a[i][3]*a[i][3] + a[i][4]*a[i][4] + a[i][5]*a[i][5];//Верно
        if((E1_2-E_2)*(E1_2-E_2)>1.0e-8)
        {
            cout<<"IT DOES NOT WORK";
            cout<<(E1_2-E_2)<<endl;
        }
    }
}
//Бинарный поиск для приближения к сетке
double getClosest(double,double,double); 
double findClosest(vector<double> VelNet, double target) 
{ 
    int n=VelNet.size();
	if (target <= VelNet[0]) 
    {
		return VelNet[0]; 
    }
	if (target >= VelNet[n - 1]) 
    {
		return VelNet[n - 1]; 
    }
	int i = 0, j = n, mid = 0; 
	while (i < j) { 
		mid = (i + j) / 2; 

		if (VelNet[mid] == target) 
			return VelNet[mid]; 
        if (target < VelNet[mid]) {  
			if (mid > 0 && target > VelNet[mid - 1]) 
				return getClosest(VelNet[mid - 1], 
								VelNet[mid], target); 
			j = mid; 
		} 

		else { 
			if (mid < n - 1 && target < VelNet[mid + 1]) 
				return getClosest(VelNet[mid], 
								VelNet[mid + 1], target);  
			i = mid + 1; 
		} 
	} 
    return VelNet[mid];
     
}  
double getClosest(double val1, double val2, double target) 
{ 
	if ((target - val1) >= (val2 - target)) 
		return val2; 
	else
		return val1; 
} 