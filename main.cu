#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include <chrono> 
#include <iomanip>
using namespace std::chrono; 
using namespace std;


// __global__ void Calculate(double *T,double *T_old,int r,int n)
// {
//     int j = blockDim.x*(blockIdx.x) + threadIdx.x;
//     int k = blockDim.y*(blockIdx.y) + threadIdx.y;
    
//     // printf("%d %d Hi\n",j,k);
//     if(j*n+k>n && j*n<=n*n)
//     {
        
//         *(T+j*n+k)=*(T_old+j*n+k)+r*( *(T_old+(j+1)*n+k)+*(T_old+j*n+k+1)+*(T_old+(j-1)*n+k)+*(T_old+j*n+k-1)- 4* *(T_old+j*n+k));
//         // T[j][k]= T_old[j][k]+r*(T_old[j+1][k]+ T_old[j][k+1]+ T_old[j-1][k]+ T_old[j][k-1]-4*T_old[j][k]);
//     }

// }


int main(){
    int n,ntime;
    double delta, sigma,nu,dom_len,dt,r;



    fstream fin;

    fin.open("input.dat",ios::in);

    fin>>n>>sigma>>nu>>dom_len>>ntime;

    fin.close();

    // cout<<n<<" "<<sigma<<" "<<nu<<" "<<dom_len<<" "<<ntime<<endl;;

    delta=dom_len/(n-1);

    dt=(sigma*pow(delta,2))/nu;



    double *T = (double *)malloc((n+1) * (n+1 )* sizeof(double));
    double *x = (double *)malloc((n+1) *( n+1) * sizeof(double));
    double *y = (double *)malloc((n+1) * (n+1) * sizeof(double));
    double *T_old = (double *)malloc((n+1) * (n+1) * sizeof(double));


auto start = high_resolution_clock::now(); 

// cout<<delta<<endl;
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
            // cout<<i<<" "<<j<<" "<<x[i][j]<<endl;
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

    // for(int i=1;i<2;++i)
    // {
    //     for(int j=1;j<n+1;++j)
    //     {    
    //         cout<<*(y+i*n+j+1)<<" ";
    //     }
    //     cout<<endl;
    // }
    
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

    

    // double*devPtr_T;
    // double*devPtr_T_old;
    // size_t pitch_T;
    // size_t pitch_T_old;
    // size_t host_pitch_T = n * sizeof(double);
    // size_t host_pitch_T_old = n * sizeof(double);


    


    // // cout << "sizeof(float): " << sizeof(double)<< endl;
    // // cout << "width: " << n << endl;
    // // cout << "height: " << n << endl;
    // // cout << "pitch:  " << pitch_T << endl;
    // // cout << "pitch:  " << pitch_T_old << endl;

  

    // cudaMallocPitch(&devPtr_T,&pitch_T,n * sizeof(double),n);
    // cudaMallocPitch(&devPtr_T_old,&pitch_T_old,n * sizeof(double),n);

    double *dev_T;
    double *dev_t_old;
  

    for(int i=1;i<=ntime;++i)
    {
        // cout<<"time_it:"<<i<<endl;
        for(int j=1;j<n+1;++j)
        {
            for(int k=1;k<n+1;++k)
            {
                *(T_old+j*n+k)=*(T+j*n+k);
            }
        }

        // cudaMemcpy2D(devPtr_T, pitch_T, &T, host_pitch_T, n*sizeof(double) , n, cudaMemcpyHostToDevice);
        // cudaMemcpy2D(devPtr_T_old, pitch_T_old,&T_old, host_pitch_T_old, n*sizeof(double) , n, cudaMemcpyHostToDevice);

        // dim3 a(32,32);
        // dim3 b(n/a.x,n/a.y);
        // Calculate<<<a,b>>>(devPtr_T,devPtr_T_old,r,n,pitch_T,pitch_T_old);
        // cudaDeviceSynchronize();

        // cudaMalloc(&dev_T,(n+1)*(n+1)*sizeof(double));
        // cudaMalloc(&dev_t_old,(n+1)*(n+1)*sizeof(double));
    
        // dim3 a(32,32);
        // dim3 b(n/a.x,n/a.y);
        // cudaMemcpy(dev_T, T, (n+1)*(n+1)*sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(dev_t_old, T_old, (n+1)*(n+1)*sizeof(double), cudaMemcpyHostToDevice);

        // Calculate<<<a,b>>>(dev_T,dev_t_old,r,n);

        // cudaDeviceSynchronize();
        // cudaError_t error = cudaGetLastError();
        // if (error != cudaSuccess) {
        // fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
        // }


        // cudaMemcpy2D(T,host_pitch_T, devPtr_T, pitch_T, n*sizeof(double), n, cudaMemcpyDeviceToHost);
        // cudaMemcpy2D(T_old, host_pitch_T_old, devPtr_T_old, pitch_T_old, n*sizeof(double), n, cudaMemcpyDeviceToHost);


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