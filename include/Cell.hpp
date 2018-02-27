#ifndef CELL_HPP
#define CELL_HPP


#include "Vector3.hpp"
#include "Generator.hpp"
#include "Wire.hpp"
#include "EndCap.hpp"
#include "IonizationEvent.hpp"













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
        for(;;)
        {
            // generate event within volume
            vector3<double> event_position{generator.GetRandomPosition(volume)};
            
            // check if event is within cell volume
            // TODO: don't need this because no "track"
            // when generating muon tracks, this was required
            // now it is not

            // check x
            if(_position_.GetX() <= event_position.GetX() && event_position.GetX() < _position_.GetX() + _radius_)
            {
                // check y
                if(_position_.GetY() <= event_position.GetY() && event_position.GetY() < _position_.GetY() + _radius_)
                {
                    // check z
                    if(_position_.GetZ() <= event_position.GetZ() && event_position.GetZ() < _position_.GetZ() + _length_)
                    {
                        return event_position;
                    }
                }
            }
        }
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
        _endcap_0_.Init(_position_, _direction_, CELL_ENDCAP_LENGTH, _radius_, 0.0);
        _endcap_1_.Init(_position_, -_direction_, CELL_ENDCAP_LENGTH, _radius_, 0.0); // TODO radius should be slightly smaller
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



#endif // CELL_HPP
