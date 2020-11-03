#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>   
#include "additional.h"
#include "workwVectors.h"
using namespace std;
#define _USE_MATH_DEFINES
const int N=20;//Размер скоростной сетки
const double Bltzmn=1.380649*1.0e-23;
const double mass=2.6*1.0e-26;
const double Temp=300.0;
const double stretch =2;//Для изменения максимальной скорости входящей в область
const double E_cut=5;//Параметр для обрезания сферы
//Метод для быстрого удаления элемента из матрицы
template <class ForwardIt, class SortUniqIndsFwdIt>
inline ForwardIt remove_at(
    ForwardIt first,
    ForwardIt last,
    SortUniqIndsFwdIt ii_first,
    SortUniqIndsFwdIt ii_last)
{
    if(ii_first == ii_last)
        return last;
    typedef typename std::iterator_traits<ForwardIt>::difference_type diff_t;
    typedef typename std::iterator_traits<SortUniqIndsFwdIt>::value_type ind_t;
    ForwardIt destination = first + static_cast<diff_t>(*ii_first);
    while(ii_first != ii_last)
    {
        for(ind_t cur = *ii_first++; ii_first != ii_last; ++ii_first)
        {
            const ind_t nxt = *ii_first;
            if(nxt - cur > 1)
                break;
            cur = nxt;
        }
        const ForwardIt source_first =
            first + static_cast<diff_t>(*(ii_first - 1)) + 1;
        const ForwardIt source_last =
            ii_first != ii_last ? first + static_cast<diff_t>(*ii_first) : last;
        std::move(source_first, source_last, destination);
        destination += source_last - source_first;
    }
    return destination;
}
int modpow(int x,int n, int m)//Вычисление x^n mod m рекурсией
{
    if(n==0) return 1%m;
    long long u=modpow(x,n/2,m);
    u=(u*u)%m;
    if(n%2==1) u=(u*x)%m;
    return u;
}
double frac(double a)//Вычисление дробного остатка
{
    return a-(int)a;
}
int I(double a)//Возвращает индекс скорости по сетке скоростей 
{
    return (E_cut+a)/(2*E_cut/N);
}
//MAIN
int main()
{
    //Переменнные
    int s = 0;//Для проверки каких-то параметров
    int N_x=N,N_y=N,N_z=N;//Размеры сетки
    double delta_E=E_cut*2.0/N;//Расстояние между скоростями
    vector<vector<double> > a;//Основной массив
    vector<vector<double> > g;//Относительные скорости
    vector<double> Def_ang;//Угол отклонения
    vector<vector<double> > E1;//Матрица скоростей после соударения
    vector<double> VelNet;//Скоростная сетка
    vector<vector<double> > E_near;//Ближайшие точки к разлетным скоростям
    vector<vector<double> > Eta_cube;//Куб скоростей вокруг одной из разлетных
    vector<double> Energy_cube;//Значение энергий вокруг 
    double E_zero;//Энергия точной разлетной скорости
    vector<double> E_cm;//Энергия центра масс
    double minDist;
    double E_check;
    
    
    long double Dlm, Dlms,Dab;
    vector<vector<double> >lam;//Скорости после проецирования
    vector<vector<double> >lam_s;//Скорости после проецирования
    vector<vector<double> >mu;//Скорости после проецирования
    vector<vector<double> >mu_s;//Скорости после проецирования
    double E_near_energy=0.0;
    double dXi_x,dXi_y,dXi_z;//Для вычисления интеграла
    long double C;//Константа при интегрировании
    double ag;//Модуль вектора g
    double gxy;//gxy^2=gx^2+gy^2
    double E1_2, E_2;//Проверка 
    double E_0v, E_1v,E_2v;//Для вычисления r
    vector<int> Elements_to_erase;//Массив индексов элементов, которые будут удалены
    double nc=0.0;// Численная плотность
    double R[8];//Вектор случайного сдвига сетки коробова
    vector<vector<double> > g1;//Относительные скорости после соударения
    long int n; //Число пар соударяющихся частиц1
    int Appr_dot;
    double N_to_fin;
    double konst_to_ln=0.0;//Для проверки аргумента логарифма
    int p=100003;//Размер сетки коробова
    int K_b[8];
    int uy=0;
    int b=20285;//Первый коэфициент сетки
    int t_max=1;//Число шагов по времени
    cout<<"Enter max t-> ";
    cin>>t_max;
    double Nu;//Число точек сетки скоротей, попадающих во вписанную сферу
    double tau=1.0/(sqrt(Bltzmn*Temp/mass)*M_PI*sqrt(2));//Шаг интегрирования по времени(длина свободного пробега на скорость)
    double rel[3];//Лишний массив для относительных скоростей
    long double Delta;//Параметр для интегрирования
    double f[N_x][N_y][N_z];//Функции распределения
    VelNet.resize(N); 
    vector<double> T_long;//Продольная температура
    vector<double> H;//Н-функции
    vector<double> T;//Температура
    H.resize(t_max);
    T.resize(t_max);
    double T_Zero=0.0;//Начальная температура(вычисляется по начальному значению сетки скоростей)
    T_long.resize(t_max); 
    for(int i=0;i<VelNet.size();i++)//Заполняем сетку скоростей
    {
        VelNet[i]=-E_cut+(i+0.5)*2*E_cut/N;
    }
    dXi_x=4.0/N*sqrt(3*Bltzmn*Temp/mass);//Шаг по рельной скорости вдоль оси
    dXi_y=4.0/N*sqrt(3*Bltzmn*Temp/mass);
    dXi_z=4.0/N*sqrt(3*Bltzmn*Temp/mass);
    double msr=dXi_x*dXi_y*dXi_z;
//Функция распределния вида f(x,y,z)=A{exp(-(Vx-U)^2-Vy^2-Vz^2)+exp(-(Vx+U)^2-Vy^2-Vz^2)} Она симметрична, но я это не использую
    for(int i=0;i<N_x;i++)//Заполняем массив функци распределения н.у.
    {
        for(int j=0;j<N_y;j++)
        {
            for(int k=0;k<N_z;k++)
            {
                if(sqrt((VelNet[i])*(VelNet[i])+VelNet[j]*VelNet[j]+VelNet[k]*VelNet[k])<E_cut)//Попадает ли скорость в сферу
                {
                f[i][j][k]=2.5*1.0e25*0.5*pow(1.0/(M_PI),1.5)*
                (
                    exp(
                    -((VelNet[i]*stretch*sqrt(3.0)/E_cut-sqrt(1.5))*(VelNet[i]*stretch*sqrt(3.0)/E_cut-sqrt(1.5))+VelNet[j]*stretch*sqrt(3.0)/E_cut*VelNet[j]*stretch*sqrt(3.0)/E_cut+VelNet[k]*stretch*sqrt(3.0)/E_cut*VelNet[k]*stretch*sqrt(3.0)/E_cut)
                    )+
                    exp(
                    -((VelNet[i]*stretch*sqrt(3.0)/E_cut+sqrt(1.5))*(VelNet[i]*stretch*sqrt(3.0)/E_cut+sqrt(1.5))+VelNet[j]*stretch*sqrt(3.0)/E_cut*VelNet[j]*stretch*sqrt(3.0)/E_cut+VelNet[k]*stretch*sqrt(3.0)/E_cut*VelNet[k]*stretch*sqrt(3.0)/E_cut)
                    )
                );
                Nu++;//Сколько точек попадает в сферу
                nc+=f[i][j][k]*msr;
                }
                else
                {
                f[i][j][k]=0.0;    
                nc+=f[i][j][k]*msr;
                }
                T_Zero+=(f[i][j][k])*((VelNet[i]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch))*(VelNet[i]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch))+VelNet[j]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)*VelNet[j]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)+VelNet[k]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)*VelNet[k]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch))*msr;                
            }
        }
    }
     for(int i=0;i<N_x;i++)//Нормируем, чтобы плотность была равна 1
    {
        for(int j=0;j<N_y;j++)
        {
            for(int k=0;k<N_z;k++)
            {
                f[i][j][k]=f[i][j][k]/(nc);
            }
        }
    }

    C=tau*M_PI*M_PI/3.0*Nu*(2*sqrt(3.0*Bltzmn*Temp/mass))*(2*sqrt(3.0*Bltzmn*Temp/mass))*(2*sqrt(3.0*Bltzmn*Temp/mass))/p;//Параметр для вычисления интеграла
    T_Zero=T_Zero*mass/(3*Bltzmn*nc);//Итоговое значение начальной температуры
    cout<<"T[0]="<<T_Zero<<endl;
    cout<<"N in sphere="<<Nu<<endl;
    cout<<"C="<<C<<endl;
    cout<<"tau="<<tau<<endl;
    vector<double> r;//Индексы для интерполяции
//Печать в файл
    ofstream out;
    out.open("resultfv");
    for(int i=0;i<N_x;i++)
        {
            for(int j=0;j<N_y;j++)
            {
                for(int k=0;k<N_z;k++)
                {
                    out<<VelNet[i]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)<<" "<<f[i][j][k]<<endl;
                }
            }
        }
    out.close();
    K_b[0]=1;//Коэфиициенты сетки коробова
    K_b[1]=b;
    for(int i=2;i<8;i++)
    {
        K_b[i]=modpow(b,i,p);
    }
cout<<"n before"<<nc<<endl;
//Цикл по времени
for(int t=0;t<t_max;t++)
{
    a=Nlrg(a,p,8);//Число соударений
    Eta_cube=Nlrg(Eta_cube,8,3);//Число ближайших точек
    E_cm.resize(3);
    //Получим массив вида n строк и 8 столбцов, первое число индексирует номер соударяющейся пары, а второй элемент точки
    //Заполняем массив сеткой Коробова размера p
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<a[i].size();j++)
        {
            a[i][j]=frac(1.0*K_b[j]*(i+1)/p);
        }
    }
    //Сдвигаем на случайный вектор каждый шаг по времени
    srand(static_cast<unsigned int>(clock()));
    for(int i=0;i<8;i++)
    {
        R[i] = double(rand()) / (double(RAND_MAX) + 1.0);
    }
    //Сдвигаем по модулю 1
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<a[i].size();j++)
        {
            a[i][j]=a[i][j]+R[j];
            if(a[i][j]>=1) a[i][j]-=1;
        }
    }
    //Теперь будем преобразововать массив
    for(int j=0;j<a.size();j++)
    {
        for(int i=0;i<6;i++)
        {
            a[j][i]=2.0*E_cut*a[j][i]-E_cut;//Преобразование координат
        }
        a[j][7]=a[j][7]*2.0*M_PI;
        a[j][6]=sqrt(a[j][6]);
    }
    //Теперь приближаем к сетке скоростей
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<6;j++)
        {
            a[i][j]=findClosest(VelNet,a[i][j]);
        }
    }
    Elements_to_erase.clear();
    //Обрезаем в сферу тоже с большой вероятностью верно
    for(int i=0;i<a.size();i++)
    {
        if((ModVector(a[i],1,3)>=E_cut)||(ModVector(a[i],4,6)>=E_cut))
        {
            Elements_to_erase.push_back(i);
        }
    }
    //Удаляем элементы, скорость которых не попадает в сферу
    a.erase(remove_at(a.begin(), a.end(), Elements_to_erase.begin(), Elements_to_erase.end()), a.end());
    g=Nlrg(g,a.size(),3);
    //Теперь находим относительные скорости
    for(int j=0;j<a.size();j++)
    {
        for(int i=0;i<3;i++)
        {
            g[j][i]=a[j][i+3]-a[j][i];
        }
    }
    Elements_to_erase.clear();
    //Отбрасываем нулевые относительные скорости
    for(int i=0;i<a.size();i++)
    {
        if(ModVector(g[i],1,3)<1.0e-8)//Проверка на равенство нулю
        {
            Elements_to_erase.push_back(i);
        }
    }
    a.erase(remove_at(a.begin(), a.end(), Elements_to_erase.begin(), Elements_to_erase.end()), a.end());
    g.erase(remove_at(g.begin(), g.end(), Elements_to_erase.begin(), Elements_to_erase.end()), g.end());
    cout<<a.size()<<endl;//Остается приблизительно 1/4 от всех точек
    //Тут уже отброшены все наподходящие а и g
    g1=Nlrg(g1,a.size(),3);
    E1=Nlrg(E1,a.size(),6);
    Def_ang.resize(a.size());
    E_near=Nlrg(E_near,a.size(),6);
    r.resize(a.size());
    //Находим углы отклонения
    for(int i=0;i<Def_ang.size();i++)
    {
        Def_ang[i]=2.0*acos(a[i][6]);
    }
    //находим относительные скорости после соударения
    for(int i=0;i<a.size();i++)
    {
        gxy=sqrt(g[i][0]*g[i][0]+g[i][1]*g[i][1]);//sqrt(g_x^2+g_y^2)
        ag=sqrt(g[i][0]*g[i][0]+g[i][1]*g[i][1]+g[i][2]*g[i][2]);//Модуль вектора g
        if(gxy>1.0e-7)//Проверка на простой случай
        {
            g1[i][0]=g[i][0]*cos(Def_ang[i])-(g[i][0]*g[i][2]/gxy)*cos(a[i][7])*sin(Def_ang[i])+(ag*g[i][1]/gxy)*sin(Def_ang[i])*sin(a[i][7]);//Верно
            g1[i][1]=g[i][1]*cos(Def_ang[i])-(g[i][1]*g[i][2]/gxy)*cos(a[i][7])*sin(Def_ang[i])-(ag*g[i][0]/gxy)*sin(Def_ang[i])*sin(a[i][7]);//Верно
            g1[i][2]=g[i][2]*cos(Def_ang[i])+gxy*cos(a[i][7])*sin(Def_ang[i]);//Верно
        }
        else//Если это простой случай
        {
           g1[i][0]=ag*sin(a[i][7])*sin(Def_ang[i]);//Верно
           g1[i][1]=ag*cos(a[i][7])*sin(Def_ang[i]);//Верно
           g1[i][2]=ag*cos(Def_ang[i]); //Верно
        }
    }
    //Тут имея относительные скорости получим новые скорсти частиц после соударения
    for(int i=0;i<E1.size();i++)
    {
        for(int k=0;k<3;k++)
        {
            E1[i][k]=(a[i][k]+a[i][k+3])*0.5-g1[i][k]*0.5;         
            E1[i][k+3]=(a[i][k]+a[i][k+3])*0.5+g1[i][k]*0.5;
        }
    }
    Elements_to_erase.clear();
    for(int i=0;i<a.size();i++)
    {
        if((ModVector(E1[i],1,3)>=E_cut)||(ModVector(E1[i],4,6)>=E_cut))
        {
            Elements_to_erase.push_back(i);
        }
    }
    a.erase(remove_at(a.begin(), a.end(), Elements_to_erase.begin(), Elements_to_erase.end()), a.end());
    E1.erase(remove_at(E1.begin(), E1.end(), Elements_to_erase.begin(), Elements_to_erase.end()), E1.end());
    E_near.erase(remove_at(E_near.begin(), E_near.end(), Elements_to_erase.begin(), Elements_to_erase.end()), E_near.end());
//Тут будет реализован проекционный метод
    //Ищем ближайшие точки к  E1
    for(int i=0;i<E1.size();i++)
    {
        for(int j=0;j<6;j++)
        {
            E_near[i][j]=findClosest(VelNet,E1[i][j]);
        }
    }
    Elements_to_erase.clear();
    //Проверяем, лежит ли ближайшая точка в сфере (если не лежит то мы сразу отбрасываем соударение)
    for(int i=0;i<a.size();i++)
    {
        if((ModVector(E_near[i],1,3)>=E_cut)||(ModVector(E_near[i],4,6)>=E_cut))//Проверка на равенство нулю
        {
            Elements_to_erase.push_back(i);
        }
    }
    a.erase(remove_at(a.begin(), a.end(), Elements_to_erase.begin(), Elements_to_erase.end()), a.end());
    E1.erase(remove_at(E1.begin(), E1.end(), Elements_to_erase.begin(), Elements_to_erase.end()), E1.end());
    E_near.erase(remove_at(E_near.begin(), E_near.end(), Elements_to_erase.begin(), Elements_to_erase.end()), E_near.end());
    //Создаем матрицы для скоростей в проекционном методе
    lam=Nlrg(lam,a.size(),3);
    lam_s=Nlrg(lam_s,a.size(),3);
    mu=Nlrg(mu,a.size(),3);
    mu_s=Nlrg(mu_s,a.size(),3);
    Elements_to_erase.clear();
    //Получили Новые массивы и набор точек в прве скоростей и нужного размера
    //Далее идет проекционный метод
    for(int i=0;i<a.size();i++)
    {
        Energy_cube.resize(8);
        Eta_cube=Nlrg(Eta_cube,8,3);//Аналогично
        //Находим покомпонентно ближайшие точки к разлетной скорости Е1[i][0..2]
        //Т.е. берем точку, затем смотрим как она расположена относительно точной и строим куб
        //Все углы проиндексированы и в зависимости от их расположения добавляется или вычитается шаг сетки
        if(E1[i][0]<E_near[i][0])
        {
            Eta_cube[0][0]=E_near[i][0]-delta_E;
            Eta_cube[4][0]=E_near[i][0]-delta_E;
            Eta_cube[3][0]=E_near[i][0]-delta_E;
            Eta_cube[7][0]=E_near[i][0]-delta_E;
            Eta_cube[1][0]=E_near[i][0];
            Eta_cube[2][0]=E_near[i][0];
            Eta_cube[5][0]=E_near[i][0];
            Eta_cube[6][0]=E_near[i][0];
        }
        else
        {
            Eta_cube[1][0]=E_near[i][0]+delta_E;
            Eta_cube[2][0]=E_near[i][0]+delta_E;
            Eta_cube[5][0]=E_near[i][0]+delta_E;
            Eta_cube[6][0]=E_near[i][0]+delta_E;
            Eta_cube[0][0]=E_near[i][0];
            Eta_cube[4][0]=E_near[i][0];
            Eta_cube[3][0]=E_near[i][0];
            Eta_cube[7][0]=E_near[i][0];

        }
        if(E1[i][1]<E_near[i][1])
        {
            Eta_cube[0][1]=E_near[i][1]-delta_E;
            Eta_cube[1][1]=E_near[i][1]-delta_E;
            Eta_cube[2][1]=E_near[i][1]-delta_E;
            Eta_cube[3][1]=E_near[i][1]-delta_E;
            Eta_cube[4][1]=E_near[i][1];
            Eta_cube[5][1]=E_near[i][1];
            Eta_cube[6][1]=E_near[i][1];
            Eta_cube[7][1]=E_near[i][1];
        }
        else
        {
            Eta_cube[4][1]=E_near[i][1]+delta_E;
            Eta_cube[5][1]=E_near[i][1]+delta_E;
            Eta_cube[6][1]=E_near[i][1]+delta_E;
            Eta_cube[7][1]=E_near[i][1]+delta_E;
            Eta_cube[0][1]=E_near[i][1];
            Eta_cube[1][1]=E_near[i][1];
            Eta_cube[2][1]=E_near[i][1];
            Eta_cube[3][1]=E_near[i][1];
        }
        if(E1[i][2]<E_near[i][2])
        {
            Eta_cube[0][2]=E_near[i][2]-delta_E;
            Eta_cube[1][2]=E_near[i][2]-delta_E;
            Eta_cube[4][2]=E_near[i][2]-delta_E;
            Eta_cube[5][2]=E_near[i][2]-delta_E;
            Eta_cube[2][2]=E_near[i][2];
            Eta_cube[3][2]=E_near[i][2];
            Eta_cube[6][2]=E_near[i][2];
            Eta_cube[7][2]=E_near[i][2];
        }
        else
        {
            Eta_cube[2][2]=E_near[i][2]+delta_E;
            Eta_cube[3][2]=E_near[i][2]+delta_E;
            Eta_cube[6][2]=E_near[i][2]+delta_E;
            Eta_cube[7][2]=E_near[i][2]+delta_E;
            Eta_cube[0][2]=E_near[i][2];
            Eta_cube[1][2]=E_near[i][2];
            Eta_cube[4][2]=E_near[i][2];
            Eta_cube[5][2]=E_near[i][2];
        }
        //Тут получаем массив состоящий из ближайших элементов к скоростям
        //Теперь находим энергиии 
        for(int k=0;k<3;k++)
        {
            E_cm[k]=0.5*(a[i][k]+a[i][k+3]);//Скорость центра масс
        }
        for(int j=0;j<8;j++)
        {
            Energy_cube[j]=(Eta_cube[j][0]-E_cm[0])*(Eta_cube[j][0]-E_cm[0])+
            (Eta_cube[j][1]-E_cm[1])*(Eta_cube[j][1]-E_cm[1])+
            (Eta_cube[j][2]-E_cm[2])*(Eta_cube[j][2]-E_cm[2]);
        }
        E_zero=0.25*((a[i][0]-a[i][3])*(a[i][0]-a[i][3])+
        (a[i][1]-a[i][4])*(a[i][1]-a[i][4])+
        (a[i][2]-a[i][5])*(a[i][2]-a[i][5]));
        E_near_energy=(E_near[i][0]-E_cm[0])*(E_near[i][0]-E_cm[0])+
        (E_near[i][1]-E_cm[1])*(E_near[i][1]-E_cm[1])+
        (E_near[i][2]-E_cm[2])*(E_near[i][2]-E_cm[2]);
        //Анализируем три случая
/*
1-Энергия точного = энергии ближайшего узла
2-Энергия точного < энергии ближайшего узла
3-Энергия точного > энергии ближайшего узла
*/

        if((E_zero-E_near_energy)*(E_zero-E_near_energy)<1.0e-10)
        {
            for(int k=0;k<3;k++)
            {
                lam[i][k]=E_near[i][k];
                lam_s[i][k]=E_near[i][k];
                mu[i][k]=E_near[i][k+3];
                mu_s[i][k]=E_near[i][k+3];
            }
            r[i]=1.0;
        }
//Второй случай
        else if(E_zero<E_near_energy)
        {
           for(int k=0;k<3;k++)
            {
                lam_s[i][k]=E_near[i][k];
                mu_s[i][k]=E_near[i][k+3];
            }
            //Удаляем элементы с энергией превышающей условную
            for(int j=0;j<Energy_cube.size();j++)
            {
                if(E_zero<Energy_cube[j])
                {
                    Energy_cube.erase(Energy_cube.begin()+j);
                    Eta_cube.erase(Eta_cube.begin()+j);
                    j--;
                }
            }
            for(int j=0;j<Energy_cube.size();j++)
            {
                //Тут будет условие на попадание скорости точки и парной к ней в сетку
                E_check=sqrt((2.0*E_cm[0]-Eta_cube[j][0])*(2.0*E_cm[0]-Eta_cube[j][0])+
                (2.0*E_cm[1]-Eta_cube[j][1])*(2.0*E_cm[1]-Eta_cube[j][1])+
                (2.0*E_cm[2]-Eta_cube[j][2])*(2.0*E_cm[2]-Eta_cube[j][2]));
                //Вычисляем скорость узла и парную к ней
                if(ModVector(Eta_cube[j],1,3)>=E_cut||E_check>=E_cut)
                {
                    Energy_cube.erase(Energy_cube.begin()+j);
                    Eta_cube.erase(Eta_cube.begin()+j);
                    j--;
                }
            }
            if(Energy_cube.size()!=0)
            {
                N_to_fin=sqrt((Eta_cube[0][0]-E1[i][0])*(Eta_cube[0][0]-E1[i][0])+
                (Eta_cube[0][1]-E1[i][1])*(Eta_cube[0][1]-E1[i][1])+
                (Eta_cube[0][2]-E1[i][2])*(Eta_cube[0][2]-E1[i][2]));
                minDist=N_to_fin;
                Appr_dot=0;
            }
            else 
            {
                //Если вообще никаких точек для проецирования нет, то записываем точку как подлежащую удалению
                minDist=-1.0;
                Elements_to_erase.push_back(i);
            }
            for(int j=0;j<Energy_cube.size();j++)
            {
                N_to_fin=sqrt((Eta_cube[j][0]-E1[i][0])*(Eta_cube[j][0]-E1[i][0])+
                (Eta_cube[j][1]-E1[i][1])*(Eta_cube[j][1]-E1[i][1])+
                (Eta_cube[j][2]-E1[i][2])*(Eta_cube[j][2]-E1[i][2]));
                //Находим точку с минимальным расстоянием до разлетного
                if(N_to_fin<minDist)
                {
                    minDist=N_to_fin;
                    Appr_dot=j;
                }
            }
            if(Energy_cube.size()!=0)
            {
                for(int k=0;k<3;k++)
                {
                    lam[i][k]=Eta_cube[Appr_dot][k];
                    mu[i][k]=2.0*E_cm[k]-lam[i][k];
                    E_0v=a[i][0]*a[i][0]+a[i][1]*a[i][1]+a[i][2]*a[i][2]+
            a[i][3]*a[i][3]+a[i][4]*a[i][4]+a[i][5]*a[i][5];
        E_1v=lam[i][0]*lam[i][0]+lam[i][1]*lam[i][1]+lam[i][2]*lam[i][2]+
            mu[i][0]*mu[i][0]+mu[i][1]*mu[i][1]+mu[i][2]*mu[i][2];
        E_2v=lam_s[i][0]*lam_s[i][0]+lam_s[i][1]*lam_s[i][1]+lam_s[i][2]*lam_s[i][2]+
            mu_s[i][0]*mu_s[i][0]+mu_s[i][1]*mu_s[i][1]+mu_s[i][2]*mu_s[i][2];
         r[i]=(E_0v-E_1v)/(E_2v-E_1v);
         cout<<r[i]<<" "<<E_0v<<" "<<E_1v<<" "<<E_2v<<endl;
                }
            }
            else
            {
               for(int k=0;k<3;k++)
                {
                    lam_s[i][k]=0;
                    mu_s[i][k]=0;
                    r[i]=M_PI;
                } 
            }
            
        }
//Третий случай
        else if(E_zero>E_near_energy)
        {
            for(int k=0;k<3;k++)
            {
                lam[i][k]=E_near[i][k];
                mu[i][k]=E_near[i][k+3];
            }
        //Удаляем элементы с энергией  не превышающей условную
            for(int j=0;j<Energy_cube.size();j++)
            {
                if(E_zero>Energy_cube[j])
                {
                    Energy_cube.erase(Energy_cube.begin()+j);
                    Eta_cube.erase(Eta_cube.begin()+j);
                    j--;
                }
            }
            for(int j=0;j<Energy_cube.size();j++)
            {
                //Тут будет условие на попадание скорости точки и парной к ней в сетку
                E_check=sqrt((2.0*E_cm[0]-Eta_cube[j][0])*(2.0*E_cm[0]-Eta_cube[j][0])+
                (2.0*E_cm[1]-Eta_cube[j][1])*(2.0*E_cm[1]-Eta_cube[j][1])+
                (2.0*E_cm[2]-Eta_cube[j][2])*(2.0*E_cm[2]-Eta_cube[j][2]));
                //Вычисляем скорость узла и парную к ней
                if(ModVector(Eta_cube[j],1,3)>=E_cut||E_check>=E_cut)
                {
                    Energy_cube.erase(Energy_cube.begin()+j);
                    Eta_cube.erase(Eta_cube.begin()+j);
                    j--;
                }
            }
            if(Energy_cube.size()!=0)
            {
                N_to_fin=sqrt((Eta_cube[0][0]-E1[i][0])*(Eta_cube[0][0]-E1[i][0])+
                (Eta_cube[0][1]-E1[i][1])*(Eta_cube[0][1]-E1[i][1])+
                (Eta_cube[0][2]-E1[i][2])*(Eta_cube[0][2]-E1[i][2]));
                minDist=N_to_fin;
                Appr_dot=0;
            }
            else 
            {
                //Если вообще никаких точек для проецирования нет, то записываем точку как подлежащую удалению
                minDist=-1.0;
                Elements_to_erase.push_back(i);
            }
            for(int j=0;j<Energy_cube.size();j++)
            {
                N_to_fin=sqrt((Eta_cube[j][0]-E1[i][0])*(Eta_cube[j][0]-E1[i][0])+
                (Eta_cube[j][1]-E1[i][1])*(Eta_cube[j][1]-E1[i][1])+
                (Eta_cube[j][2]-E1[i][2])*(Eta_cube[j][2]-E1[i][2]));
                //Находим точку с минимальным расстоянием до разлетного
                if(N_to_fin<minDist)
                {
                    minDist=N_to_fin;
                    Appr_dot=j;
                }
            }
            if(Energy_cube.size()!=0)
            {
                for(int k=0;k<3;k++)
                {
                    lam_s[i][k]=Eta_cube[Appr_dot][k];
                    mu_s[i][k]=2.0*E_cm[k]-lam_s[i][k];
                }
                E_0v=a[i][0]*a[i][0]+a[i][1]*a[i][1]+a[i][2]*a[i][2]+
            a[i][3]*a[i][3]+a[i][4]*a[i][4]+a[i][5]*a[i][5];
        E_1v=lam[i][0]*lam[i][0]+lam[i][1]*lam[i][1]+lam[i][2]*lam[i][2]+
            mu[i][0]*mu[i][0]+mu[i][1]*mu[i][1]+mu[i][2]*mu[i][2];
        E_2v=lam_s[i][0]*lam_s[i][0]+lam_s[i][1]*lam_s[i][1]+lam_s[i][2]*lam_s[i][2]+
            mu_s[i][0]*mu_s[i][0]+mu_s[i][1]*mu_s[i][1]+mu_s[i][2]*mu_s[i][2];
         r[i]=(E_0v-E_1v)/(E_2v-E_1v);
         cout<<r[i]<<" "<<E_0v<<" "<<E_1v<<" "<<E_2v<<endl;
            }
            else
            {
               for(int k=0;k<3;k++)
                {
                    lam_s[i][k]=0;
                    mu_s[i][k]=0;
                } 
            }
        }
        
    }
    a.erase(remove_at(a.begin(), a.end(), Elements_to_erase.begin(), Elements_to_erase.end()), a.end());
    E1.erase(remove_at(E1.begin(), E1.end(), Elements_to_erase.begin(), Elements_to_erase.end()), E1.end());
    E_near.erase(remove_at(E_near.begin(), E_near.end(), Elements_to_erase.begin(), Elements_to_erase.end()), E_near.end());
    lam.erase(remove_at(lam.begin(), lam.end(), Elements_to_erase.begin(), Elements_to_erase.end()), lam.end());
    lam_s.erase(remove_at(lam_s.begin(), lam_s.end(), Elements_to_erase.begin(), Elements_to_erase.end()), lam_s.end());
    mu.erase(remove_at(mu.begin(), mu.end(), Elements_to_erase.begin(), Elements_to_erase.end()), mu.end());
    mu_s.erase(remove_at(mu_s.begin(), mu_s.end(), Elements_to_erase.begin(), Elements_to_erase.end()), mu_s.end());
    r.erase(remove_at(r.begin(), r.end(), Elements_to_erase.begin(), Elements_to_erase.end()), r.end());
//Теперь вычисляем интеграл
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            rel[j]=a[i][j]-a[i][j+3];//Находим относительную скорость
        }
        if((f[I(lam[i][0])][I(lam[i][1])][I(lam[i][2])]>1.0e-15)&&(f[I(mu[i][0])][I(mu[i][1])][I(mu[i][2])]
        >1.0e-15))//Проверяем деление на ноль
        {
            Delta=(pow(f[I(lam_s[i][0])][I(lam_s[i][1])][I(lam_s[i][2])]*f[I(mu_s[i][0])][I(mu_s[i][1])][I(mu_s[i][2])]/(f[I(lam[i][0])][I(lam[i][1])][I(lam[i][2])]*f[I(mu[i][0])][I(mu[i][1])][I(mu[i][2])]),r[i])*f[I(lam[i][0])][I(lam[i][1])][I(lam[i][2])]*f[I(mu[i][0])][I(mu[i][1])][I(mu[i][2])]-
            f[I(a[i][0])][I(a[i][1])][I(a[i][2])]*f[I(a[i][3])][I(a[i][4])][I(a[i][5])])*
            sqrt(rel[0]*rel[0]+rel[1]*rel[1]+rel[2]*rel[2])*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch);
        }
        else 
        {
            Delta=-f[I(a[i][0])][I(a[i][1])][I(a[i][2])]*f[I(a[i][3])][I(a[i][4])][I(a[i][5])]*
            sqrt(rel[0]*rel[0]+rel[1]*rel[1]+rel[2]*rel[2])*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch);
        }
        //Проверяем неотрицательность
        if(f[I(a[i][0])][I(a[i][1])][I(a[i][2])]+C*Delta<0||
        f[I(a[i][3])][I(a[i][4])][I(a[i][5])]+C*Delta<0||
        f[I(lam[i][0])][I(lam[i][1])][I(lam[i][2])]-(1-r[i])*C*Delta<0||
        f[I(mu[i][0])][I(mu[i][1])][I(mu[i][2])]-(1-r[i])*C*Delta<0||
        f[I(lam_s[i][0])][I(lam_s[i][1])][I(lam_s[i][2])]-r[i]*C*Delta<0||
        f[I(mu_s[i][0])][I(mu_s[i][1])][I(mu_s[i][2])]-r[i]*C*Delta<0)
        {
            continue;
        }
        else
        {
            f[I(a[i][0])][I(a[i][1])][I(a[i][2])]+=C*Delta;
            f[I(a[i][3])][I(a[i][4])][I(a[i][5])]+=C*Delta;
            f[I(lam[i][0])][I(lam[i][1])][I(lam[i][2])]-=(1-r[i])*C*Delta;
            f[I(mu[i][0])][I(mu[i][1])][I(mu[i][2])]-=(1-r[i])*C*Delta;
            f[I(lam_s[i][0])][I(lam_s[i][1])][I(lam_s[i][2])]-=r[i]*C*Delta;
            f[I(mu_s[i][0])][I(mu_s[i][1])][I(mu_s[i][2])]-=r[i]*C*Delta;
        }
    }

//Проверки
    ZSE(a,E1);
    EqualityOfVectors(g,g1);
    ZSI(E1,a);
    ZSI(E1,E_near);
    ZSI(a,E_near);
    cout<<t<<endl;
    T[t]=0;
    T_long[t]=0;
    H[t]=0;
    //Находим Температуры и Н-функцию в момент времени t
    for(int i=0;i<N_x;i++)
    {
        for(int j=0;j<N_y;j++)
        {
            for(int k=1;k<N_z;k++)
            {
                konst_to_ln=(f[i][j][k]*f[i][j][k]);
                if(konst_to_ln>0)
                {
                    H[t]+=f[i][j][k]*0.5*log(konst_to_ln)*msr;
                }  
                T[t]+=msr*(f[i][j][k])*((VelNet[i]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch))*(VelNet[i]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch))+VelNet[j]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)*VelNet[j]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)+VelNet[k]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)*VelNet[k]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch));

                T_long[t]+=f[i][j][k]*VelNet[i]*VelNet[i]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)*msr;
            }
        }
    }
    //Проверяем, что сохраняется плотность
    nc=0.0;
    for(int i=0;i<N_x;i++)
    {
        for(int j=0;j<N_y;j++)
        {
            for(int k=0;k<N_z;k++)
            {
                nc+=f[i][j][k]*msr;
            }
        }
    }    
    int h=0;
    cout<<"n after"<<nc<<endl;
}    
//Конец цикла
nc=0.0;
    for(int i=0;i<N_x;i++)
    {
        for(int j=0;j<N_y;j++)
        {
            for(int k=0;k<N_z;k++)
            {
                nc+=f[i][j][k]*msr;
            }
        }
    }    
    cout<<"n after"<<nc<<endl;
    //Результаты выводятся
    out.open("resultfvfin");
    for(int i=0;i<N_x;i++)
    {      
        for(int j=0;j<N_y;j++)
        {
            for(int k=0;k<N_z;k++)
            {
                out<<VelNet[i]*sqrt(3*Bltzmn*Temp/mass)/(E_cut/stretch)<<" "<<f[i][j][k]<<endl;
            }
        }
    }
    out.close();
    out.open("Templong");
    for(int i=0;i<t_max;i++)
    {
        out<<i<<" "<<T_long[i]*(mass/(3*Bltzmn*nc))<<endl;
    }
    out.close();
    out.open("Temp");
    for(int i=0;i<t_max;i++)
    {
        out<<i<<" "<<(T[i]*(mass/(3*Bltzmn*nc)))/*/T_Zero-1)*100*/<<endl;
    }
    out.close();
    out.open("TempR");
    for(int i=0;i<t_max;i++)
    {
        out<<i<<" "<<(T[i]*(mass/(3*Bltzmn*nc))/T_Zero-1)*100<<endl;
    }
    out.close();
    out.open("Hfunc");
    for(int i=0;i<t_max;i++)
    {
        out<<i<<" "<<H[i]<<endl;
    }
    out.close();
}
