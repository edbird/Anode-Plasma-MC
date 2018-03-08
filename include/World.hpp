#ifndef WORLD_HPP
#define WORLD_HPP

#include "Vector3.hpp"
#include "Cell.hpp"
#include "Generator.hpp"
#include "IonizationEvent.hpp"


// global definitions of tracker cell size
const double TRACKER_CELL_LENGTH{2.900};
const double TRACKER_CELL_X{0.040}; // x and y size, also radius
const double TRACKER_CELL_Y{TRACKER_CELL_X};

const double TRACKER_CELL_ANODE_WIRE_VOLTAGE{800.0};


class World
{


    public:

    World()
        : _cell_(TRACKER_CELL_LENGTH, 0.5 * TRACKER_CELL_X, TRACKER_CELL_ANODE_WIRE_VOLTAGE) // 2.9 m by 40 mm cell
        , _volume_(3.0 * TRACKER_CELL_X, 3.0 * TRACKER_CELL_Y, 3.0 * TRACKER_CELL_LENGTH) // world volume
        //: _volume_(3.0 * TRACKER_CELL_X, 3.0 * TRACKER_CELL_Y, 3.0 * TRACKER_CELL_LENGTH) // world volume
    {
        // set cell position to be in center of world volume
        _cell_.SetPosition(vector3<double>(TRACKER_CELL_X, TRACKER_CELL_Y, TRACKER_CELL_LENGTH));
    }


    void DoEvent()
    {
        IonizationEvent event{_cell_.GenerateIonizationEvent(_generator_, _volume_)};
        std::cout << event.GetPosition() << std::endl;
    }


    private:

    Cell _cell_; // world contains a single cell
    Generator _generator_; // event generator

    // world volume
    vector3<double> _volume_; // events are generated in this volume
};


#endif
