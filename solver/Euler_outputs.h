#include<fstream>
#include<cstdlib>
#include"Fieldvar.h"
#include"Gridman.h"

class Euler_outputs
{
	private:


	public:
		Gridman *meshptr;
		std::vector<Fieldvar *> fieldvarptrs;
		int porder;
		void writevtufile(std::string fname);
		void writeonedfile(std::string fname);


};
