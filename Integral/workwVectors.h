#ifndef V_H
#define V_H
void printVector(std::vector<double>a);
void printMatrix(std::vector<std::vector<double> > a);
void EqualityOfVectors(std::vector<std::vector<double> >g1, std::vector<std::vector<double> >g);
std::vector<std::vector<double> > Nlrg(std::vector<std::vector<double> > v,int v_index, int v_coordinates);
double ModVector(std::vector<double> v, int n,int m);
#endif