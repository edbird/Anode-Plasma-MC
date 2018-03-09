


//#include "Vector3.hpp"
//#include "Cell.hpp"
#include "World.hpp"

//#include "Generator.hpp"
//#include "Wire.hpp"
//#include "EndCap.hpp"
//#include "IonizationEvent.hpp"
//#include "Cell.hpp"

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

    //Generator generator;
    //Wire wire;
    //EndCap endcap;
    //IonizationEvent ionizationevent(vector3<double>(0.0, 0.0, 0.0));

    //Cell cell(2.900, 0.040, 800.0);

    World world;
    for(int i = 0; i < 10000; ++ i)
        world.DoEvent();
   

    return 0;
}
