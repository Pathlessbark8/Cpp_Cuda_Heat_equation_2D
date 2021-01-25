#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include <chrono> 
#include <iomanip>
using namespace std::chrono; 
using namespace std;



int main(){
    int n,ntime;
    double delta, sigma,nu,dom_len,dt,r;



    fstream fin;

    fin.open("input.dat",ios::in);

    fin>>n>>sigma>>nu>>dom_len>>ntime;

    fin.close();


    delta=dom_len/(n-1);

    dt=(sigma*pow(delta,2))/nu;



    double *T = (double *)malloc((n+1) * (n+1 )* sizeof(double));
    double *x = (double *)malloc((n+1) *( n+1) * sizeof(double));
    double *y = (double *)malloc((n+1) * (n+1) * sizeof(double));
    double *T_old = (double *)malloc((n+1) * (n+1) * sizeof(double));


auto start = high_resolution_clock::now(); 

    for(int i=1;i<n+1;++i)
    {
        *(x+n+i)=0.0;
       
        *(x+n*n+i)=2.0;
      
        *(y+n*i+1)=0.0;
      
        *(y+n+n*i)=2.0;
    }
   

    for(int i=2;i<n;++i)
    {
        for(int j=1;j<n+1;++j)
        {    
            *(x+i*n+j)=*(x+(i-1)*n+j)+delta;
            *(y+j*n+i)=*(y+(j)*n+i-1)+delta;
        }
    }

   
    
    for(int i=1;i<n+1;++i)
    {
        for(int j=1;j<n+1;++j)
        {
            if(*(x+i*n+j)<=1.5 && *(x+i*n+j)>=0.5 && *(y+i*n+j)<=1.5 && *(y+i*n+j) >=0.5)
            {
                *(T+i*n+j)=2.0;
            }           
            else{
                *(T+i*n+j)=1.0;
            }
        }
    }
    
    fstream foutw;
    foutw.open("int.dat",ios::out);
    
    for(int i=1;i<n+1;++i)
    {
        for(int j=1;j<n+1;++j)
        {
            foutw<<*(x+i*n+j)<<" "<<*(y+i*n+j)<<" "<<*(T+i*n+j)<<"\n";
        }
    }
    foutw.close();
  
    r=(nu*dt)/pow(delta,2);


    for(int i=1;i<=ntime;++i)
    {
        cout<<"time_it:"<<i<<endl;
        for(int j=1;j<n+1;++j)
        {
            for(int k=1;k<n+1;++k)
            {
                *(T_old+j*n+k)=*(T+j*n+k);
            }
        }

        for(int j=2;j<n;j++)
        {
            for(int k=2;k<n;k++)
            {
                *(T+j*n+k)=*(T_old+j*n+k)+r*( *(T_old+(j+1)*n+k)+*(T_old+j*n+k+1)+*(T_old+(j-1)*n+k)+*(T_old+j*n+k-1)- 4* *(T_old+j*n+k));
            }
        }
    }


    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<milliseconds>(stop - start); 

    cout<<"simulation completed"<<endl;
    cout<<"Time_Taken : "<<duration.count()<<endl;

    fstream fout;
    fout.open("soln.dat",ios::out);
    
    for(int i=1;i<n+1;++i)
    {
        for(int j=1;j<n+1;++j)
        {
            fout<<std::scientific<<*(x+i*n+j)<<" "<<*(y+i*n+j)<<" "<<*(T+i*n+j)<<"\n";
        }
    }
    fout.close();
    return 0;
}