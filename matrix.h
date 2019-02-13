#ifndef STDIO_H
#include <stdio.h>
#define STDIO_H
#endif

#ifndef STDLIB_H
#include <stdlib.h>
#define STDLIB_H
#endif

#ifndef ASSERT_H
#include <assert.h>
#define ASSERT_H
#endif

#ifndef RANDOM_H
#include <random>
#define RANDOM_H
#endif

#ifndef MATRIX_H
#define MATRIX_H
#define abs(x) x<0 ? -x : x // abs function
#define lim 1000 // limit to random (-lim,+lim)

class Result
{
	private:
	int n;
	private:
	double *x0;
	private:
	double *xk;
	
	public:
	Result(int n)
	{
		this->n=n;
		x0 = (double*)malloc(sizeof(double)*n);
		xk = (double*)malloc(sizeof(double)*n);
		assert(x0);
		assert(xk);	
	}
	
	public:
	~Result()
	{
		free(x0);
		free(xk);
	}
	
	public:
	void reset()
	{
		#pragma omp parallel for
		for(int i=0; i<n ;i++)
		{
			xk[i] = 0;
			x0[i] = 0;
		}
	}
	
	public:
	void show()
	{
		for(int i=0; i<n ;i++)
		{
			printf("%f\n",x0[i]);
		}	
	}
	
	public:
	double* getX0()
	{
		return x0;
	}
	
	public:
	double * getXk()
	{
		return xk;
	}
};

class Matrix
{

	// matrix
	private:
	double **a;
	
	// vector of independent terms
	private:
	double *b;
	
	// order of matrix 'a' (n x n)
	private:
	int n;
	
	// vector of diagonal values
	private:
	double *diag;
	

	/*
	 * Algorithm to create a diagonally dominant matrix
	 */	 
	public:
	Matrix(int n)
	{
		this->n = n;
		alloc();
		values();
	}

	/*
	 * algorithm to alloc a (n x n) matrix
	 */
	private:
	void alloc()
	{
		// creating diagonal vector
		diag = (double*)malloc(sizeof(double)*n);
		assert(diag);
	
		// creating all lines
		a = (double**) malloc(sizeof(double*)*n);
		assert(a);
		
		b = (double*) malloc(sizeof(double)*n);
		assert(b);
		
		// creating all collumns
		for(int i = 0 ; i < n ; i++)
		{ 
			a[i] = (double*) malloc(sizeof(double)*n);
			assert(a[i]);
		}
	}
	
	/*
	 * algorithm to fill the matrix with random numbers
	 * and ajust the result to get a diagonally dominant matrix
	 */
	private:
	void values()
	{
		// choose the type of distribution (uniform end real)
		std::random_device generator;
    	//pcg rand(gerenator);
		
		//std::default_random_engine generator;
  		std::uniform_real_distribution<double> distribution(-lim,lim);
  		#pragma omp parallel for
  		for(int i = 0 ; i < n ; i++)
  		{
  			
  			double s = 0;
  			// randomizing numbers
  			for(int j = 0 ; j < n; j++)
  			{
  				a[i][j] = distribution(generator);
  				if(i!=j)
  				{
  					// add all line's elements if i!=j
  					s += abs(a[i][j]);
  				}
  			}
  			
  			// ajust to get diagonally dominant matrix
  			a[i][i] = abs(a[i][i]);
  			if(a[i][i]<s)
  			{
  				a[i][i] /= 100;
  				a[i][i] += s;  
  			}
  			// randomizing a result to line i	
  			b[i] = distribution(generator);
  		}
	}
	
	/*
	 * algorithmm to return the matrix's order
	 */
	public:
	int getN()
	{
		return n;
	}
	
	/*
	 * algorithm to return the matrix's address 'a'
	 */
	 public:
	 double **getA()
	 {
	 	return a;
	 }
	 
	 /*
	  * algorithm to return the independent terms 'b'
	  */
	  public:
	  double *getB()
	  {
	  		return b;
	  }
	  
	 /*
	  * algorithm to return the Diagonal Vector
	  */
	  public:
	  double *getDiag()
	  {
	  		return diag;
	  }
	  
	  /*
	   * algorithm to print a matrix
	   * and write a matrix on HD
	   */
	  public:
	  void show()
	  {	
	  		FILE *arq = fopen("matrix","w");
	  		fprintf(arq,"%i\n",n);
	  		for(int i = 0 ; i  < n ; i++)
	  		{
	  			for(int j = 0 ; j < n ; j++)
	  			{
	  				printf("%f ",a[i][j]);
	  				fprintf(arq,"%.16f\n",a[i][j]);
	  			}
	  			printf(" = %f\n",b[i]);
	  		}
	  		
	  		for(int i=0;i<n;i++)
	  			fprintf(arq,"%.16f\n",b[i]);	
	  		fclose(arq);
	  }
	  
	  /*
	   * algorithm to destruct the matrix
	   */
	  public:
	  ~Matrix()
	  {
	  		// cleaning diag
	  		free(diag);
	  
	  		// cleaning B
	  		free(b);
	  		
	  		// cleaning A
	  		for(int i=0;i<n;i++)
	 			free(a[i]);
	 		free(a);
	 		
	 		// cleaning N
	 		n = 0;
	  }
	  
	  public:
	  void mult(Result *r)
	  {
	  		double *x0 = r->getX0();
	  		for(int i=0;i<n;i++)
	  		{
	  			double s=0;
	  			for(int j=0;j<n;j++)
	  			{
	  				s+= a[i][j]*x0[j];
	  			}
	  			printf("%f\n",s);
	  		}
	  		
	  }	  	
};
#endif

