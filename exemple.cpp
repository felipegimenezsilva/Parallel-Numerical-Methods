/*
 * Felipe Gimenez da Silva
 * 02 - 13 - 2019
 * algorithm: How to use CompareMethods object
 * to compile: g++ exemple.cpp -fopenmp
 */

#include "methods.h"

int main()
{
	/* creating object pointer
	 * parameters:
	 * 		initial matrix size 	( 1 )
	 *		increment size		(200)
	 *		average of		( 5 )
	 */
	//CompareMethods *compare = new CompareMethods(1,1,1);
	CompareMethods compare(1,200,5);
	
	/*
	 * adding list of (objects) Methods to compare
	 *		# to create new object, put your class 
	 *		# in methods.h
	 */
	
	compare.add(new JacobiSequential());	// gauss jacobi sequential
	compare.add(new JacobiParallel());	// gauss jacobi parallel
	compare.add(new SeidelSequential()); 	// gauss seidel sequential
	compare.add(new SeidelParallel());	// gauss seidel parallel 	(change matrix)
	compare.add(new SeidelSemiParallel());	// gauss seidel semi parallel
	compare.add(new SeidelParallel2());	// gauss seidel parallel 2 	(no change matrix)
	compare.add(new SeidelUnstable());	// gauss seidel unistable	(fast and random)
	
	
	/*
	 * starting comparation
	 * 	# this is a loop, whitout END.
	 */
	compare.start();
}
