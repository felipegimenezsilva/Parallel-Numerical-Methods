/*
 * Felipe Gimenez
 * 02 - 13 - 2019
 * Algorithm to manage matrix and vectors. 
 */

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
#define lim 1 // limit to random (-lim,+lim)

/*
 * Result Class is used to store the answer (x0),
 * xk is used temporarily, n is the vector length.
 */
class Result
{
	private:
	int n;
	private:
	double *x0;
	private:
	double *xk;
	
	/*
	 * Constructor method
	 * alloc all the vectors when the object is created
	 * # to 'change' the length, need to create another Result
	 */
	public:
	Result(int n)
	{
		this->n=n;
		x0 = (double*)malloc(sizeof(double)*n);
		xk = (double*)malloc(sizeof(double)*n);
		assert(x0);
		assert(xk);	
		reset();
	}
	
	/*
	 * Destructor Method
	 * free all the vectors
	 */
	public:
	~Result()
	{
		free(x0);
		free(xk);
	}
	
	/*
	 * Reset all the vectors 
	 * # N value not change
	 */
	public:
	void reset()
	{
		#pragma omp parallel for
		#pragma omp nowait
		for(int i=0; i<n ;i++)
		{
			xk[i] = 0;
			x0[i] = 0;
		}
	}
	
	/* display x0 vector */
	public:
	void show()
	{
		for(int i=0; i<n ;i++)
		{
			printf("%f\n",x0[i]);
		}	
	}
	
	// returns x0 pointer
	public:
	double* getX0()
	{
		return x0;
	}
	
	// return xk pointer
	public:
	double * getXk()
	{
		return xk;
	}
};


/*
 * Matrix class is used to:
 *	alloc:
 *		Matrix A (matrix)
 *		Vector B (independent terms)
 *	Randomize:
 *		Diagonally Dominant matrix (A)
 *		Vector (B)
 *	Delete (A,B,...)
 */
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
	 * algorithm to fill the matrix with pseudo random numbers
	 * and ajust the result to get a diagonally dominant matrix
	 */
	public:
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
  				a[i][j] = (double)distribution(generator);
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
  				a[i][i] /= i%5 + i/n + 0.1;
  				a[i][i] += s;  
  			}
  			
  			// randomizing a result to line i	
  			b[i] = distribution(generator);
  		}
	}
	/*
	 * I think is better use random values, but
	 * randomize completely not sure of convergence
	 */
	 /*
	public:
	void values2()
	{
	
		std::random_device generator;
  		std::uniform_real_distribution<double> dist1(-lim,lim);
  		//std::bernoulli_distribution dist1(lim);
	
		double *x0 = (double*)malloc(sizeof(double)*n);
		for(int i = 0; i<n ; i++)
		{
			x0[i] = dist1(generator);
		}
		for(int i = 0; i<n; i++)
		{
			for(int j=0; j<n; j++)
			{
				if(i!=j)
				{
					a[i][j] = dist1(generator);
					b[i] += x0[j] * a[i][j];
				}
			}
			a[i][i] = 1 + i;
			b[i] += a[i][i] * x0[i];
		}
		free(x0);
	}*/

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
	  		for(int i = 0 ; i  < n ; i++)
	  		{
	  			for(int j = 0 ; j < n ; j++)
	  			{
	  				printf("%f ",a[i][j]);
	  			}
	  			printf(" = %f\n",b[i]);
	  		}
	  }
	  
	  /*
	   * algorithm to save a matrix (.txt)
	   */
	  public:
	  void save()
	  {	
	  		FILE *file = fopen("txt_mat.txt","w");
	  		for(int i = 0 ; i  < n ; i++)
	  		{
	  			for(int j = 0 ; j < n ; j++)
	  			{
	  				fprintf(file,"%f ",a[i][j]);
	  			}
	  			fprintf(file,"\n");
	  		}
	  		for(int i=0; i<n ; i++)
	  			fprintf(file,"%f\n",b[i]);
	  		fclose(file);
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
	  
	  /*
	   * Algorithm to Calc and show:
	   *	A * X0
	   * ---- note : A * X0 = B
	   * ---- A and B is member of MATRIX
	   * ---- X0 is member of RESULT
	   */
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
