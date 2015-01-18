#include"Euler_timestepper.h"
int main()
{
	Euler_timestepper *obj;
	obj = new Euler_timestepper;

	obj->init("inputs-DG");
	obj->Explicit_timestepping();

	delete(obj);

	return(0);
}
