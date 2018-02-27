#ifndef CELL_HPP
#define CELL_HPP


#include "Vector3.hpp"
#include "Generator.hpp"


class TrackerGas
{


    private:

        
    double _ionization_energy_; // first ionization energy of gas atom / molecule
    double _gamma_mean_free_path_; // assume all energies of photons have same mean free path
    double _gas_fraction_ionizable_;
    double _gas_fraction_quenching_;

};


class IonizationEvent
{

    public:

    IonizationEvent(const vector3<double>& position)
        : _position_{position}
    {

    }

    private:

    vector3<double> _position_;

};


namespace Geometry
{
    
    class Cylinder
    {

        public:

        Cylinder()
            : _length_{0.0}
            , _radius_{0.0}
        {
        }

        Cylinder(const vector3<double> &position, const vector3<double> &direction, const double length, const double radius)
            : _position_{position}
            , _direction_{direction}
            , _length_{length}
            , _radius_{radius}
        {
        }

        void Init(const vector3<double> &position, const vector3<double> &direction, const double length, const double radius)
        {
            _position_ = position;
            _direction_ = direction;
            _length_ = length;
            _radius_ = radius;
        }

        private:

        vector3<double> _position_;
        vector3<double> _direction_;
        double _length_;
        double _radius_;

    };


} // namespace Geometry

class Wire
{

    public:

    
    Wire()
    {
    }
     

    Wire(const vector3<double> &position, const vector3<double> &direction, const double length, const double radius, const double voltage)
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



const double CELL_ENDCAP_LENGTH{0.1}; // 100 mm


class Cell
{


    public:

    Cell(const double length, const double radius, const double anode_voltage)
        : _position_(0.0, 0.0, 0.0)
        , _direction_(0.0, 0.0, 1.0) // anode wire points along z
        , _length_{length}
        , _radius_{radius}
    {
        init_wire(anode_voltage);
        init_endcap();
    }

    void SetPosition(const vector3<double> position)
    {
        _position_ = position;
    }

    IonizationEvent GenerateIonizationEvent(Generator& generator, const vector3<double>& volume) const
    {
        vector3<double> event_position{generator.GetRandomPosition(volume)};
    }


    private:

    void init_wire(const double anode_voltage)
    {
        // initialized anode wire
        _wire_anode_.Init(_position_, _direction_, _length_, _radius_, anode_voltage);
        
        // initialize ground wire vector
        _wire_ground_.clear();

        for(int y{-1}; y <= 1; ++ y)
        {
            for(int x{-1}; x <= 1; ++ x)
            {
                if(x == 0 && y == 0) continue;

                vector3<double> position{_position_};
                const vector3<double> x_delta(0.5 * _radius_, 0.0, 0.0);
                const vector3<double> y_delta(0.0, 0.5 * _radius_, 0.0);
                position += x * x_delta;
                position += y * y_delta;
                _wire_ground_.push_back(Wire(position, _direction_, _length_, _radius_, 0.0));
            }
        }
    }

    void init_endcap()
    {
        _endcap_0_.Init(_position_, _direction_, CELL_ENDCAP_LENGTH, _radius_);
        _endcap_1_.Init(_position_, -_direction_, CELL_ENDCAP_LENGTH, _radius_); // TODO radius should be slightly smaller
    }


    vector3<double> _position_; // vector to anode wire start position
    vector3<double> _direction_; // unit vector in direction along anode wire
    double _length_;
    double _radius_;
    const double _ACTIVE_RADIUS_FRACTION_{0.9}; // the active radius is 

    Wire _wire_anode_;
    /*
    Wire _wire_ground_0_;
    Wire _wire_ground_1_;
    Wire _wire_ground_2_;
    Wire _wire_ground_3_;
    Wire _wire_ground_4_;
    Wire _wire_ground_5_;
    Wire _wire_ground_6_;
    Wire _wire_ground_7_;
    */
    std::vector<Wire> _wire_ground_;


    EndCap _endcap_0_;
    EndCap _endcap_1_;

};


// global definitions of tracker cell size
const double TRACKER_CELL_LENGTH{2.900};
const double TRACKER_CELL_X{0.040}; // x and y size, also radius
const double TRACKER_CELL_Y{TRACKER_CELL_X};

const double TRACKER_CELL_ANODE_WIRE_VOLTAGE{800.0};


class World
{


    public:

    World()
        : _cell_(TRACKER_CELL_LENGTH, TRACKER_CELL_X, TRACKER_CELL_ANODE_WIRE_VOLTAGE) // 2.9 m by 40 mm cell
        , _volume_(3.0 * TRACKER_CELL_X, 3.0 * TRACKER_CELL_Y, 3.0 * TRACKER_CELL_LENGTH) // world volume
    {
        // set cell position to be in center of world volume
        _cell_.SetPosition(vector3<double>(TRACKER_CELL_X, TRACKER_CELL_Y, TRACKER_CELL_LENGTH));
    }


    void DoEvent()
    {
        IonizationEvent event{_cell_.GenerateIonizationEvent(_generator_, _volume_)};
    }


    private:

    Cell _cell_; // world contains a single cell
    Generator _generator_; // event generator

    // world volume
    vector3<double> _volume_; // events are generated in this volume
};


#endif // CELL_HPP
