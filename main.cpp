


//#include "Vector3.hpp"
//#include "Cell.hpp"
#include "World.hpp"


#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>



int main(int argc, char* argv[])
{

    /*
    vector3<double> a(0.0, 0.0, 1.0);
    vector3<double> b(0.0, 0.0, 0.0);
    b = a;
    b += a;
    std::cout << "a=" << a << " b=" << b << std::endl;
    */

    Cell cell(2.900, 0.040, 800.0);

    World world;
    world.DoEvent();
   

    return 0;
}
