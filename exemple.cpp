/*
 * Felipe Gimenez da Silva
 * 02 - 13 - 2019
 * algorithm: How to use CompareMethods object
 */

#include "methods.h"

int main()
{
	/* creating object pointer
	 * parameters:
	 * 		initial matrix size ( 1 )
	 *		increment size		(200)
	 *		average of			( 5 )
	 */
	CompareMethods *compare = new CompareMethods(1,200,5);
	
	/*
	 * adding list of (objects) Methods to compare
	 *		# to create new object, put your class 
	 *		# in methods.h
	 */
	compare->add(new JacobiSequential());
	compare->add(new JacobiParallel());
	compare->add(new SeidelSequential());
	compare->add(new SeidelParallel());
	compare->add(new SeidelUnstable());
	
	/*
	 * starting comparation
	 * 		# this is a loop, whitout END.
	 */
	compare->start();

}
