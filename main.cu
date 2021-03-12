// #include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include <chrono> 
#include <iomanip>
using namespace std::chrono; 
using namespace std;


__global__ void Calculate(double *T,double *T_old,int r,int n)
{
    int j = blockDim.x*(blockIdx.x) + threadIdx.x;
    int k = blockDim.y*(blockIdx.y) + threadIdx.y;

    if(j*n+k>n && j*n+k<=n*n)
    {
        
        *(T+j*n+k)=*(T_old+j*n+k)+r*( *(T_old+(j+1)*n+k)+*(T_old+j*n+k+1)+*(T_old+(j-1)*n+k)+*(T_old+j*n+k-1)- 4* *(T_old+j*n+k));
    }
    __syncthreads();
    
}


int main(){
    int n,ntime;
    double delta, sigma,nu,dom_len,dt,r;



    fstream fin;

    fin.open("input.dat",ios::in);

    fin>>n>>sigma>>nu>>dom_len>>ntime;

    fin.close();


    delta=dom_len/(n-1);

    dt=(sigma*pow(delta,2))/nu;



auto start = high_resolution_clock::now(); 

    double *T = (double *)malloc((n+1) * (n+1 )* sizeof(double));
    double *x = (double *)malloc((n+1) *( n+1) * sizeof(double));
    double *y = (double *)malloc((n+1) * (n+1) * sizeof(double));
    double *T_old = (double *)malloc((n+1) * (n+1) * sizeof(double));




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

    



    double *temp;
    double *dev_T;
    double *dev_t_old;


    
    cudaMalloc(&dev_T,(n+1)*(n+1)*sizeof(double));
    cudaMalloc(&dev_t_old,(n+1)*(n+1)*sizeof(double));
    cudaMalloc(&temp,(n+1)*(n+1)*sizeof(double));

    // for(int j=1;j<n+1;++j)
    //     {
    //         for(int k=1;k<n+1;++k)
    //         {
    //             *(T_old+j*n+k)=*(T+j*n+k);
    //         }
    //     }
    
    cudaMemcpy(dev_T, T, (n+1)*(n+1)*sizeof(double), cudaMemcpyHostToDevice);

    for(int i=1;i<=ntime;++i)
    {
        //cout<<"time_it:"<<i<<endl;
        
       
        temp=dev_T;
        dev_T=dev_t_old;
        dev_t_old=temp;
        dim3 a(32,8);
        dim3 b(n/a.x,n/a.y);
        // cudaMemcpy(dev_t_old, T_old, (n+1)*(n+1)*sizeof(double), cudaMemcpyHostToDevice);

        Calculate<<<b,a>>>(dev_T,dev_t_old,r,n);

        cudaDeviceSynchronize();
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
        fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
        }

      
        // cudaMemcpy(T_old, dev_t_old, (n+1)*(n+1)*sizeof(double), cudaMemcpyDeviceToHost);
        
    }
    cudaMemcpy(T, dev_t_old, (n+1)*(n+1)*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(&dev_T);
    cudaFree(&dev_t_old);
    cudaFree(&temp);

    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<milliseconds>(stop - start); 

    printf("simulation completed\n");
    printf("Time_Taken :  %d\n",duration.count());

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
    delete x,y,T,T_old,temp;
    return 0;
}