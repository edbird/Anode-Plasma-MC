#ifndef ENDCAP_HPP
#define ENDCAP_HPP


class EndCap
{

    public:

    
    EndCap()
    {
    }
    

    EndCap(const vector3<double> &position, const vector3<double> &direction, const double length, const double radius, const double voltage)
        : _cylinder_(position, direction, length, radius)
        , _voltage_{voltage}
    {
    }

    // init function, works like constructor
    void Init(const vector3<double> &position, const vector3<double> &direction, const double length, const double radius, const double voltage)
    {
        _cylinder_.Init(position, direction, length, radius);
        _voltage_ = voltage;
    }

    void SetVoltage(const double voltage)
    {
        _voltage_ = voltage;
    }

    double GetVoltage() const
    {
        return _voltage_;
    }


    private:

    Geometry::Cylinder _cylinder_;
    double _voltage_;

};


#endif
