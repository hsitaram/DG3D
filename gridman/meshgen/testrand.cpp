#include<cstdlib>
#include<iostream>
#include<ctime>

int main()
{
	double randno;

	srand(time(NULL));

	randno = ( (double) rand()/RAND_MAX);

	std::cout<<"rand:"<<randno<<"\n";
}
