#include "methods.h"

int main()
{
	CompareMethods *compare = new CompareMethods(1,200,10);
	compare->add(new JacobiSequential());
	compare->add(new JacobiParallel());
	compare->add(new SeidelSequential());
	compare->add(new SeidelParallel());
	compare->add(new SeidelUnstable());
	compare->start();

}
