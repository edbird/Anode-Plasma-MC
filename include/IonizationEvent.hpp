#ifndef IONIZATIONEVENT_HPP
#define IONIZATIONEVENT_HPP


#include "Vector3.hpp"


class IonizationEvent
{

    public:

    IonizationEvent(const vector3<double>& position)
        : _position_{position}
    {

    }

    const vector3<double>& GetPosition() const
    {
        return _position_;
    }


    private:

    vector3<double> _position_;

};

#endif // IONIZATIONEVENT_HPP
