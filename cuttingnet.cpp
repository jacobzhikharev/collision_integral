#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;
#define _USE_MATH_DEFINES
//Метод для быстрого удаления элемента из матрицы
template <class ForwardIt, class SortUniqIndsFwdIt>
inline ForwardIt remove_at(
    ForwardIt first,
    ForwardIt last,
    SortUniqIndsFwdIt ii_first,
    SortUniqIndsFwdIt ii_last)
{
    if(ii_first == ii_last) // no indices-to-remove are given
        return last;
    typedef typename std::iterator_traits<ForwardIt>::difference_type diff_t;
    typedef typename std::iterator_traits<SortUniqIndsFwdIt>::value_type ind_t;
    ForwardIt destination = first + static_cast<diff_t>(*ii_first);
    while(ii_first != ii_last)
    {
        // advance to an index after a chunk of elements-to-keep
        for(ind_t cur = *ii_first++; ii_first != ii_last; ++ii_first)
        {
            const ind_t nxt = *ii_first;
            if(nxt - cur > 1)
                break;
            cur = nxt;
        }
        // move the chunk of elements-to-keep to new destination
        const ForwardIt source_first =
            first + static_cast<diff_t>(*(ii_first - 1)) + 1;
        const ForwardIt source_last =
            ii_first != ii_last ? first + static_cast<diff_t>(*ii_first) : last;
        std::move(source_first, source_last, destination);
        // std::copy(source_first, source_last, destination) // c++98 version
        destination += source_last - source_first;
    }
    return destination;
}
//Вывод вектора Работает тоже верно
void printVector(vector<double>a)
{
    for(int i = 0; i<a.size(); i++)
    {
        cout<< a[i]<< endl;
    }
}
//Вывод матрицы любого размера
//Этот войд работает точно верно
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
//Проверка равенства векторов с какой-то точностью. ТОже работает верно
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
//Проверка закона сохранения энергии
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
//Увеличение размера массива тоже вроде верно
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
//Проверка зси
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
int modpow(int x,int n, int m)
{
    if(n==0) return 1%m;
    long long u=modpow(x,n/2,m);
    u=(u*u)%m;
    if(n%2==1) u=(u*x)%m;
    return u;
}
double frac(double a)
{
    return a-(int)a;
}
//MAIN
int main()
{
    //Переменнные
    int s = 0;//Для проверки каких-то параметров
    int N=20;//Число точек в пр-ве скоростей вдоль оси
    double E_cut=4.8;//Параметр для обрезания сферы
    vector<vector<double> > a;//Основной массив
    vector<vector<double> > g;//Относительные скорости
    vector<double> Def_ang;//Угол отклонения
    vector<vector<double> > E1;//Матрица скоростей после соударения
    vector<double> VelNet;//Скоростная сетка
    vector<vector<double> > E_near;//Ближайшие точки к разлетным скоростям
    vector<vector<double> > Eta_cube;//Куб скоростей вокруг одной из разлетных
    vector<double> Energy_cube;//Значение энергий вокруг 
    double E_zero;//Энергия точной разлетной скорости
    vector<double> E_cm;
    double minDist;
    double E_check;
    vector<vector<double> >lam;
    vector<vector<double> >lam_s;
    vector<vector<double> >mu;
    vector<vector<double> >mu_s;
    double E_near_energy=0.0;
    double ag;//Модуль вектора g
    double gxy;//gxy^2=gx^2+gy^2
    double E1_2, E_2;//Проверка 
    vector<int> Elements_to_erase;
    double R[8];//Вектор случайного сдвига сетки коробова
    vector<vector<double> > g1;//Относительные скорости после соударения
    long int n; //Число пар соударяющихся частиц1
    int Appr_dot;
    double N_to_fin;
    int p=50021;//Размер сетки коробова
    cout<<"In p: ";
    cin>>p;
    int K_b[8];
    int b=11281;
    cout<<"In b: ";
    cin>>b;
    //Какая-то обработка
    n=p;
    //Изменение размера массивов
    a=Nlrg(a,n,8);
    VelNet.resize(N); 
    Eta_cube=Nlrg(Eta_cube,8,3);
    E_cm.resize(3);
    //Получим массив вида n строк и 8 столбцов, первое число индексирует номер соударяющейся пары, а второй элемент точки
    //Первые 6 - скорости, оставшиеся углы и прицельные парамет
    
    srand(static_cast<unsigned int>(clock()));//Что-то для случайного числа
    for(int i=0;i<8;i++)
    {
        R[i] = double(rand()) / (double(RAND_MAX) + 1.0);//Случайный элемент массива
        cout<<R[i]<<endl;
    }
    //Тут имеем массив a[n][8] со случайными числами из(0,1)(или из [0,1))
    //Заполняем массив сеткой Коробова размера p
    //Находим коэффициенты
    K_b[0]=1;
    K_b[1]=b;
    for(int i=2;i<8;i++)
    {
        K_b[i]=modpow(b,i,p);
        cout<<K_b[i]<<endl;
    }
    //Заполняем массив
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<a[i].size();j++)
        {
            a[i][j]=frac(1.0*K_b[j]*(i+1)/p);
        }
    }
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
    }
    //Приближение к дискретной сетке скоростей
    //Создаем саму сетку скоростей
    for(int i=0;i<20;i++)
    {
        VelNet[i]=-E_cut+(i+0.5)*2*E_cut/N;
    }
    //Теперь приближаем к ней 
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<6;j++)
        {
            a[i][j]=findClosest(VelNet,a[i][j]);
        }
    }
    //Обрезаем в сферу тоже с большой вероятностью верно
    for(int i=0;i<a.size();i++)
    {
        if((ModVector(a[i],1,3)>=E_cut)||(ModVector(a[i],4,6)>=E_cut))
        {
            Elements_to_erase.push_back(i);
        }
    }
    a.erase(remove_at(a.begin(), a.end(), Elements_to_erase.begin(), Elements_to_erase.end()), a.end());
    n=a.size();
    g=Nlrg(g,n,3);
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
    cout<<endl<<Elements_to_erase.size()<<" 0 rel."<<endl;
    a.erase(remove_at(a.begin(), a.end(), Elements_to_erase.begin(), Elements_to_erase.end()), a.end());
    g.erase(remove_at(g.begin(), g.end(), Elements_to_erase.begin(), Elements_to_erase.end()), g.end());
    cout<<a.size()<<endl;
    //Тут уже отброшены все наподходящие а и g
    n=a.size();
    g1=Nlrg(g1,n,3);
    E1=Nlrg(E1,n,6);
    Def_ang.resize(n);
    E_near=Nlrg(E_near,n,6);
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
        //Верно
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
    g.erase(remove_at(g.begin(), g.end(), Elements_to_erase.begin(), Elements_to_erase.end()), g.end());
    g1.erase(remove_at(g1.begin(), g1.end(), Elements_to_erase.begin(), Elements_to_erase.end()), g1.end());
    Def_ang.erase(remove_at(Def_ang.begin(), Def_ang.end(), Elements_to_erase.begin(), Elements_to_erase.end()), Def_ang.end());
    n=a.size();
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
    //Проверяем, лежит ли ближайшая в точке(если не лежит то мы сразу отбрасываем соударение)
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
    g.erase(remove_at(g.begin(), g.end(), Elements_to_erase.begin(), Elements_to_erase.end()), g.end());
    g1.erase(remove_at(g1.begin(), g1.end(), Elements_to_erase.begin(), Elements_to_erase.end()), g1.end());
    Def_ang.erase(remove_at(Def_ang.begin(), Def_ang.end(), Elements_to_erase.begin(), Elements_to_erase.end()), Def_ang.end());
    n=a.size();
    //Создаем матрицы для скоростей в проекционном методе
    lam=Nlrg(lam,n,3);
    lam_s=Nlrg(lam_s,n,3);
    mu=Nlrg(mu,n,3);
    mu_s=Nlrg(mu_s,n,3);
    cout<<n<<" Before proection"<<endl;
    Elements_to_erase.clear();
    //Получили Новые массивы и набор точек в прве скоростей и нужного размера
//Далее идет проекционный метод
    //Начало самого проекционного метода
    for(int i=0;i<a.size();i++)
    {
        Energy_cube.resize(8);//Тк алгоритм подразумевает отбрасывание лишних элементов массива
        Eta_cube=Nlrg(Eta_cube,8,3);//Аналогично
        //Находим покомпонентно ближайшие точки к разлетной скорости Е1[i][0..2]
        if(E1[i][0]<E_near[i][0])
        {
            Eta_cube[0][0]=E_near[i][0]-0.48;
            Eta_cube[4][0]=E_near[i][0]-0.48;
            Eta_cube[3][0]=E_near[i][0]-0.48;
            Eta_cube[7][0]=E_near[i][0]-0.48;
            Eta_cube[1][0]=E_near[i][0];
            Eta_cube[2][0]=E_near[i][0];
            Eta_cube[5][0]=E_near[i][0];
            Eta_cube[6][0]=E_near[i][0];
        }
        else
        {
            Eta_cube[1][0]=E_near[i][0]+0.48;
            Eta_cube[2][0]=E_near[i][0]+0.48;
            Eta_cube[5][0]=E_near[i][0]+0.48;
            Eta_cube[6][0]=E_near[i][0]+0.48;
            Eta_cube[0][0]=E_near[i][0];
            Eta_cube[4][0]=E_near[i][0];
            Eta_cube[3][0]=E_near[i][0];
            Eta_cube[7][0]=E_near[i][0];

        }
        if(E1[i][1]<E_near[i][1])
        {
            Eta_cube[0][1]=E_near[i][1]-0.48;
            Eta_cube[1][1]=E_near[i][1]-0.48;
            Eta_cube[2][1]=E_near[i][1]-0.48;
            Eta_cube[3][1]=E_near[i][1]-0.48;
            Eta_cube[4][1]=E_near[i][1];
            Eta_cube[5][1]=E_near[i][1];
            Eta_cube[6][1]=E_near[i][1];
            Eta_cube[7][1]=E_near[i][1];
        }
        else
        {
            Eta_cube[4][1]=E_near[i][1]+0.48;
            Eta_cube[5][1]=E_near[i][1]+0.48;
            Eta_cube[6][1]=E_near[i][1]+0.48;
            Eta_cube[7][1]=E_near[i][1]+0.48;
            Eta_cube[0][1]=E_near[i][1];
            Eta_cube[1][1]=E_near[i][1];
            Eta_cube[2][1]=E_near[i][1];
            Eta_cube[3][1]=E_near[i][1];
        }
        if(E1[i][2]<E_near[i][2])
        {
            Eta_cube[0][2]=E_near[i][2]-0.48;
            Eta_cube[1][2]=E_near[i][2]-0.48;
            Eta_cube[4][2]=E_near[i][2]-0.48;
            Eta_cube[5][2]=E_near[i][2]-0.48;
            Eta_cube[2][2]=E_near[i][2];
            Eta_cube[3][2]=E_near[i][2];
            Eta_cube[6][2]=E_near[i][2];
            Eta_cube[7][2]=E_near[i][2];
        }
        else
        {
            Eta_cube[2][2]=E_near[i][2]+0.48;
            Eta_cube[3][2]=E_near[i][2]+0.48;
            Eta_cube[6][2]=E_near[i][2]+0.48;
            Eta_cube[7][2]=E_near[i][2]+0.48;
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
        if((E_zero-E_near_energy)*(E_zero-E_near_energy)<1.0e-7)
        {
            for(int k=0;k<3;k++)
            {
                lam[i][k]=E_near[i][k];
                lam_s[i][k]=E_near[i][k];
                mu[i][k]=E_near[i][k+3];
                mu_s[i][k]=E_near[i][k+3];
            }
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
                }
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

//Проверки
    ZSE(a,E1);
    EqualityOfVectors(g,g1);
    ZSI(E1,a);
    ZSI(E1,E_near);
    ZSI(a,E_near);
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            if((lam[i][j]+mu[i][j]-a[i][j]-a[i][j+3])*(lam[i][j]+mu[i][j]-a[i][j]-a[i][j+3])>1.0e-8)
            {
                cout<<"MISTAKE  ";
            }
        }
    }
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            if((lam_s[i][j]+mu_s[i][j]-a[i][j]-a[i][j+3])*(lam_s[i][j]+mu_s[i][j]-a[i][j]-a[i][j+3])>1.0e-8)
            {
                cout<<"MISTAKE s-a ";
            }
        }
    }
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            if((lam_s[i][j]+mu_s[i][j]-lam[i][j]-mu[i][j])*(lam_s[i][j]+mu_s[i][j]-lam[i][j]-mu[i][j])>1.0e-8)
            {
                cout<<"MISTAKE  lam-mu ";//<<lam_s[][]+mu_s[][]<<" "<<lam[];
            }
        }
    }
    for(int i=0;i<a.size();i++)
    {
        if(ModVector(lam_s[i],1,3)>=E_cut||ModVector(lam[i],1,3)>=E_cut||ModVector(mu[i],1,3)>=E_cut||ModVector(mu_s[i],1,3)>=E_cut)
        {
            cout<<"Nu ty i Mudak ";
        }
    }
    //printMatrix(E_near);//Вывод массива
    cout <<endl
        << "ty  "
        << M_PI
        <<"door"
        <<endl
        <<a.size()
        <<endl;
    
        
}
