#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP


#include "Vector3.hpp"


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

        // position is position of one end
        // direction is direction vector from the position vector
        // length is the length along the direction vector the cylinder
        // extends
        // radius is the radius of the cylinder
        vector3<double> _position_;
        vector3<double> _direction_;
        double _length_;
        double _radius_;

    };


    class Cube
    {
        
        public:

        Cube()
        {
        }


        void Translate(const vector3<double>& translation)
        {
            _position_ += translation;
        }

        //void Rotate

        private:

        // position is the position of the "origin" corner of the cube
        // _direction_x_ is the unit vector pointing along the x axis of the
        // cube from the position vector (before any rotation / translation
        // this vector also points along the x axis of the "world")
        // _direction_y_ is the same but for the y axis
        vector3<double> _position_;
        vector3<double> _direction_x_;
        vector3<double> _direction_y_;
    };


} // namespace Geometry

#endif // GEOMETRY_HPP
