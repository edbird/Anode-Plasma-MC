#ifndef GENERATOR_HPP
#define GENERATOR_HPP


#include "Vector3.hpp"


#include <random>
#include <cmath>


class Generator
{


    public:

    Generator()
        : _random_device_()
        , _gen_(_random_device_())
        , _dis_0_(0.0, 1.0)
    {
    }

    // TODO
    void RandomSeed()
    {
        _gen_.seed(_random_device_());
    }

    void DefaultSeed()
    {
        _gen_.seed(0);
    }

    vector3<double> GetRandomPosition(const vector3<double>& volume) /*const*/
    {
        double r_x{GetRandomUniform()}; // generate random number r x
        double r_y{GetRandomUniform()}; // generate random number r y
        double r_z{GetRandomUniform()}; // generate random number r z
        //return vector3<double>(r_x * volume.GetX(), r_y * volume.GetY(), r_z * volume.GetZ());
        vector3<double> position(r_x, r_y, r_z);
        position.Scale(volume);
        return position;
    }

    // TODO
    double GetRandomExponential(const double lambda)
    {
        std::exponential_distribution<double> dis(lambda);
        return dis(_gen_);
    }

    // TODO
    double GetRandomUniform()
    {
        return _dis_0_(_gen_);
    }

    double GetRandomUniform(const double low, const double high)
    {
        if(high < low) throw std::string("Error in ") + std::string(__func__) + std::string(" high < low");
        return GetRandomUniform() * (high - low) + low;
    }


    private:

    const double PI{4.0 * std::atan(1.0)};
    std::random_device _random_device_;
    std::mt19937_64 _gen_;
    
    // distribution range 0.0, 1.0
    std::uniform_real_distribution<double> _dis_0_;

};


#endif
