/*
 * Felipe Gimenez
 * 02 - 13 - 2019
 * Algorithm to:
 * 		create a Queue of methods to compare
 * 		calc the average of time elapsed to solve linear problems
 */

#ifndef TIME_H
#include <sys/time.h>
#define TIME_H
#endif

#ifndef VECTOR
#define VECTOR
#include <vector>
#endif

#ifndef OMP_H
#include "omp.h"
#define OMP_H
#endif


#ifndef MATRIX_H
#include "matrix.h"
#define MATRIX_H
#endif

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#ifndef ASSERT_H
#define ASSERT_H
#include <assert.h>
#endif

#ifndef METHODS_H
#define METHODS_H

/* 
 * this section contains all configurations of linear system
 */
#define ITER 10000 // max iteration to solve
#define ERRORMAX 0.0001 // max error

class Methods
{

	// used only by specifications
	public:
	virtual int solve(Matrix *matrix, Result *result)
	{}
	
	// apply little and especific changes on matrix
	public:
	virtual void changes(Matrix *matrix)
	{}
	
	public:
	virtual void save(int n,double time, int k)
	{}
	
};

class CompareMethods
{
	// quantity of tests to the same size
	private:
	int repeat;
	
	// initial value of n (order to matrix n x n)
	private:
	int ini;
	
	// next ini = ini + step
	private:
	int step;
	
	// testList
	std::vector<Methods*> testList;
	
	// testAverage
	std::vector<double> average;
	
	// the function put the initial values
	public:
	CompareMethods(int ini, int step, int repeat)
	{
		this->repeat = repeat;
		this->ini = ini;
		this->step = step;
	}
	
	// algorithm to create a list of tests
	public:
	void add(Methods *method)
	{
		testList.push_back(method);
		average.push_back(0.0);
	}
	
	public:
	void start()
	{
		
		Matrix *matrix;
		Result *result;
		int k;

		// this loop reSize the problem
		for(int i = ini; 1 ; i += step )
		{

			matrix = new Matrix(i);
			result = new Result(i);
			
			// this loop repeat the problem with the same size
			for(int j = 0; j < repeat ; j++)
			{
				// this loop solve all testList
				for(int p = 0; p < testList.size(); p++)
				{
					// make changes if necessary
					testList[p]->changes(matrix);
					
					// start chronometer
					double ti,tf,dt;
					dt = ti = tf  = 0;
					struct timeval tii,tff;
					gettimeofday(&tii,NULL);

					// start to solve
					k = testList[p]->solve(matrix,result);
	
					// stop chronometer
					gettimeofday(&tff,NULL);
					tf = (double)tff.tv_usec +((double)tff.tv_sec*(1000000.0));
					ti = (double)tii.tv_usec +((double)tii.tv_sec*(1000000.0));
					dt=(tf-ti)/1000;
					
					// saving test : order , time , k_iteration
					testList[p]->save(i,dt,k);
					
					// add the time elapsed (to calc the average)
					average[p]+=dt;
					
					// reset all the vectors of result class
					result->reset();
					
					// undo matrix changes if necessary
					testList[p]->changes(matrix);
				}
				matrix->values();
			}
			//print size of matrix
			printf("%i, ",i);
			for(int j = 0 ; j < testList.size(); j++)
			{
				// calc the average and show
				average[j] /= repeat;
				printf("%lf%s",average[j],j==testList.size()-1?"\n":", ");
				average[j] = 0.0;
			}
			// deleting objects to recreate  with new size
			delete matrix;
			delete result;
		}
	}
			
};

class JacobiSequential:public Methods
{	
	// solve the linear system using Jacobi Sequencial
	public:
	int solve(Matrix *matrix, Result *result)
	{
		// initial values
		double **a = matrix->getA();
		double *b = matrix->getB();
		int n = matrix->getN();
		double *x0 = result->getX0();
		double *xk = result->getXk();
		int k = 0;
		bool again;
		// end initial values
		
  		while( k < ITER )
  		{
    		again = false;
    		for( int i = 0 ; i < n ; i++ )
    		{
      		xk[i] = b[i] / a[i][i] ;
      		for( int j=0;j<n; j++)
        			if( i != j ) xk[i] -= a[i][j] * x0[j] / a[i][i];
    		}
    		for( int i = 0 ; i < n ; i++)
    		{
    			double err = abs((xk[i]-x0[i])/xk[i]);
      		if( err  > ERRORMAX)
        			again = true;
      		x0[i] = xk[i];
    		}
    		if(!again) break;
    		k++;
  		}
		return k;
	}	
	
	public:
	void save(int n,double time, int k)
	{
		FILE *file = fopen("jsequ.csv","a");
		fprintf(file,"%i, %lf, %i\n",n,time,k);
		fclose(file);
	}
};

class JacobiParallel:public Methods
{

	// solve the linear system using Jacobi Parallel
	public:
	int solve(Matrix *matrix, Result *result)
	{
		// initial values
		double **a = matrix->getA();
		double *b = matrix->getB();
		int n = matrix->getN();
		double *x0 = result->getX0();
		double *xk = result->getXk();
		int k = 0;
		bool again;
		// end initial values
		
  		while( k < ITER )
  		{
    		again = false;
    		#pragma omp parallel for
    		#pragma omp nowait
    		for( int i = 0 ; i < n ; i++ )
    		{
      		xk[i] = b[i] / a[i][i] ;
      		for( int j=0;j<n; j++)
        			if( i != j ) xk[i] -= a[i][j] * x0[j] / a[i][i];
    		}
			#pragma omp parallel for
			#pragma omp nowait
    		for( int i = 0 ; i < n ; i++)
    		{
    			double err = abs((xk[i]-x0[i])/xk[i]);
      		if( err > ERRORMAX)
        			again = true;
      		x0[i] = xk[i];
    		}
    		if(!again) break;
    		k++;
  		}
		return k;
	}
	
	public:
	void save(int n,double time, int k)
	{
		FILE *file = fopen("jpara.csv","a");
		fprintf(file,"%i, %lf, %i\n",n,time,k);
		fclose(file);
	}
	
};

class SeidelSequential:public Methods
{

	// solve the linear system using Seidel Sequential
	public:
	int solve(Matrix *matrix, Result *result)
	{
		// initial values
		double **a = matrix->getA();
		double *b = matrix->getB();
		int n = matrix->getN();
		double *x0 = result->getX0();
		double *xk = result->getXk();
		int k = 0;
		bool again;
		// end initial values
		

  		while( k < ITER )
  		{
    		again = false;
    		for( int i = 0 ; i < n ; i++ )
    		{
      		xk[i] = 0;
      		for( int j=0;j<n; j++)
        			if( i != j ) xk[i] -= a[i][j] * xk[j] / a[i][i];
    			xk[i] += b[i] / a[i][i];
    		
    			double err = abs((xk[i]-x0[i])/xk[i]);
      		if( err > ERRORMAX)
        			again = true;
      		x0[i] = xk[i];
    		}
    		if(!again) break;
    		k++;
  		}
		return k;
	}
	
	public:
	void save(int n,double time, int k)
	{
		FILE *file = fopen("ssequ.csv","a");
		fprintf(file,"%i, %lf, %i\n",n,time,k);
		fclose(file);
	}
	
};

class SeidelParallel:public Methods
{
	public:
	void changes(Matrix *matrix)
	{
		int n = matrix->getN() - 1;
		double **a = matrix->getA();
		double *diag = matrix->getDiag();
		
		// algorithm to change element position
		// i -> n - j - 1
		// j -> n - i - 1
		int x = n;
		for(int j = 0 ; j < x; j++ , x-- )
		{
			for(int i = j+1 ; i < x ; i++)
			{
				double temp = a[i][j];
				a[i][j] = a[n-j][n-i];
				a[n-j][n-i] = temp;
			}
		}
		
		for(int i = 0; i<n+1 ; i++)
		{
			diag[i] = a[i][i]; 
		}
	}

	// solve the linear system using Seidel Parallel with changes
	public:
	int solve(Matrix *matrix, Result *result)
	{
		// initial values
		double *diag = matrix->getDiag();
		double **a = matrix->getA();
		double *b = matrix->getB();
		int n = matrix->getN();
		double *x0 = result->getX0();
		double *xk = result->getXk();
		int k = 0;
		bool again;
		// end initial values
		

  		while( k < ITER )
  		{
    		again = false;
    		
    		#pragma omp parallel for
    		#pragma omp nowait
    		for( int i = 0 ; i < n ; i++ )
    		{
      		xk[i] = b[i] / a[i][i] ;
      		for( int j = i+1 ; j < n ; j++)
        			xk[i] -= (a[i][j] * x0[j]) / a[i][i];
    		}
    		
    		for( int i = n-1 ; i > 0 ; i-- )
    		{
    			#pragma omp parallel for
    			#pragma omp nowait
      		for( int j = i-1 ; j >= 0 ; j--)
      		{ 
      			// changedI and changedJ
      			int cI = n - j - 1 , cJ = n - i - 1;
        			xk[cI] -= (a[i][j] * xk[cJ]) / diag[cI];
      		}
    		}  		
    		
			#pragma omp parallel for
			#pragma omp nowait
    		for( int i = 0 ; i < n ; i++)
    		{
    			double err = abs((xk[i]-x0[i])/xk[i]);
      		if( err > ERRORMAX)
        			again = true;
      		x0[i] = xk[i];
    		}
    		if(!again) break;
    		k++;
  		}
		return k;
	}
	
	public:
	void save(int n,double time, int k)
	{
		FILE *file = fopen("schan.csv","a");
		fprintf(file,"%i, %lf, %i\n",n,time,k);
		fclose(file);
	}
	
};

class SeidelParallel2:public Methods
{	
	// solve the linear system using Seidel Parallel with changes
	public:
	int solve(Matrix *matrix, Result *result)
	{
		// initial values
		double *diag = matrix->getDiag();
		double **a = matrix->getA();
		double *b = matrix->getB();
		int n = matrix->getN();
		double *x0 = result->getX0();
		double *xk = result->getXk();
		int k = 0;
		bool again;
		// end initial values
		

  		while( k < ITER )
  		{
    		again = false;
    		
    		#pragma omp parallel for
    		#pragma omp nowait
    		for( int i = 0 ; i < n ; i++ )
    		{
      		xk[i] = b[i] / a[i][i] ;
      		for( int j = i+1 ; j < n ; j++)
        			xk[i] -= (a[i][j] * x0[j]) / a[i][i];
    		}
    		
    		for( int j=0 ; j<n ; j++ )
    		{
    			#pragma omp parallel for
    			#pragma omp nowait
      		for(int i=j+1; i<n ; i++)
      		{ 
      			xk[i] -= (a[i][j] * xk[j]) / a[i][i];
      		}
    		}  		
    		
			#pragma omp parallel for
			#pragma omp nowait
    		for( int i = 0 ; i < n ; i++)
    		{
    			double err = abs((xk[i]-x0[i])/xk[i]);
      		if( err > ERRORMAX)
        			again = true;
      		x0[i] = xk[i];
    		}
    		if(!again) break;
    		k++;
  		}
		return k;
	}
	
	public:
	void save(int n,double time, int k)
	{
		FILE *file = fopen("snoch.csv","a");
		fprintf(file,"%i, %lf, %i\n",n,time,k);
		fclose(file);
	}
	
};


class SeidelSemiParallel:public Methods
{

	// solve the linear system using Seidel Parallel
	public:
	int solve(Matrix *matrix, Result *result)
	{
		// initial values
		double **a = matrix->getA();
		double *b = matrix->getB();
		int n = matrix->getN();
		double *x0 = result->getX0();
		double *xk = result->getXk();
		int k = 0;
		bool again;
		// end initial values
		

  		while( k < ITER )
  		{
    		again = false;
    		#pragma omp parallel for
      	#pragma omp nowait
      	// U matrix -> parallel
    		for( int i = 0 ; i < n ; i++ )
    		{
      		xk[i] = b[i] / a[i][i];
      		for( int j=i+1;j<n; j++)
        			xk[i] -= a[i][j] * x0[j] / a[i][i];
    		}
    		// L matrix -> sequential
	 		for(int i=1;i<n;i++)
	 			for(int j=0;j<i;j++)
	 				xk[i] -= a[i][j] * xk[j] / a[i][i];
    		
    		// Absolute error -> parallel
    		#pragma omp parallel for
    		#pragma om nowait
    		for(int i=0;i<n;i++)
    		{
    			double err = abs((xk[i]-x0[i])/xk[i]);
      		if( err > ERRORMAX)
        			again = true;
      		x0[i] = xk[i];
    		}
    		if(!again) break;
    		k++;
  		}
		return k;
	}
	
	public:
	void save(int n,double time, int k)
	{
		FILE *file = fopen("ssemi.csv","a");
		fprintf(file,"%i, %lf, %i\n",n,time,k);
		fclose(file);
	}
	
};

class SeidelUnstable:public Methods
{
	// solve the linear system using Seidel Sequential
	public:
	int solve(Matrix *matrix, Result *result)
	{
		// initial values
		double *diag = matrix->getDiag();
		double **a = matrix->getA();
		double *b = matrix->getB();
		int n = matrix->getN();
		double *x0 = result->getX0();
		double *xk = result->getXk();
		int k = 0;
		bool again;
		// end initial values
		

  		while( k < ITER )
  		{
    		again = false;
    		#pragma omp parallel for
    		for( int i = 0 ; i < n ; i++ )
    		{
      		xk[i] = 0;
      		for( int j=0;j<n; j++)
        			if( i != j ) xk[i] -= a[i][j] * x0[j] / a[i][i];
    			xk[i] += b[i] / a[i][i];
    		
    			double err = abs((xk[i]-x0[i])/xk[i]);
      		if( err > ERRORMAX)
        			again = true;
      		x0[i] = xk[i];
    		}
    		if(!again) break;
    		k++;
  		}
		return k;
	}
	
	public:
	void save(int n,double time, int k)
	{
		FILE *file = fopen("sunst.csv","a");
		fprintf(file,"%i, %lf, %i\n",n,time,k);
		fclose(file);
	}
	
};

#endif
