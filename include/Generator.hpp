#ifndef GENERATOR_HPP
#define GENERATOR_HPP


#include <random>
#include <cmath>


class Generator
{


    public:

    Generator()
        : _dis_0_(0.0, 1.0)
    {
    }

    vector3<double> GetRandomPosition(const vector3<double>& volume) /*const*/
    {
        double r_x{_dis_0_(_gen_)}; // generate random number r x
        double r_y{_dis_0_(_gen_)}; // generate random number r y
        double r_z{_dis_0_(_gen_)}; // generate random number r z
        //return vector3<double>(r_x * volume.GetX(), r_y * volume.GetY(), r_z * volume.GetZ());
        vector3<double> position(r_x, r_y, r_z);
        position.Scale(volume);
        return position;
    }



    private:

    const double PI{4.0 * std::atan(1.0)};
    std::random_device _random_device_;
    std::mt19937_64 _gen_;
    
    // distribution range 0.0, 1.0
    std::uniform_real_distribution<double> _dis_0_;

};


#endif
