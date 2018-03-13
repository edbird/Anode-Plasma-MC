#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP


#include "vector3.hpp"


namespace Geometry
{

    // abstract base class for 3d shapes
    class Volume
    {

        public:

        Volume()
        {
        }

        Volume(const vector3<double>& position)
            : _position_{position}
        {
        }

        virtual
        ~Volume() = 0;

        void Init(const vector3<double>& position)
        {
            _position_ = position;
        }
        
        void Translate(const vector3<double>& translation)
        {
            _position_ += translation;
        }

        virtual
        bool PointIntersectionTest(vector3<double> point) const = 0;

        inline
        const vector3<double>& Position() const
        {
            return _position_;
        }

        void SetPosition(const vector3<double>& position)
        {
            _position_ = position;
        }

        private:

        vector3<double> _position_;

    };

    Volume::~Volume() {}

    
    // cylinder class not finished
    class Cylinder : public Volume
    {

        public:

        Cylinder()
            : Volume()
            , _length_{0.0}
            , _radius_{0.0}
        {
        }

        Cylinder(const vector3<double> &position, const vector3<double> &direction, const double length, const double radius)
            : Volume(position)
            , _direction_{direction}
            , _length_{length}
            , _radius_{radius}
        {
        }

        void Init(const vector3<double> &position, const vector3<double> &direction, const double length, const double radius)
        {
            Volume::Init(position);
            _direction_ = direction;
            _length_ = length;
            _radius_ = radius;
        }

        /*
        const vector3<double>& GetPosition() const
        {
            return Volume::GetPosition();
        }
        */

        const vector3<double>& Direction() const
        {
            return _direction_;
        }

        double Radius() const
        {
            return _radius_;
        }

        double Length() const
        {
            return _length_;
        }

        virtual
        bool PointIntersectionTest(vector3<double> point) const
        {
            // TODO: check this function

            point -= Position();
            // TODO: handle rotation

            // check z
            if(0.0 <= point.GetZ() && point.GetZ() < _length_)
            {
                // check radial
                double x{Position().GetX()};
                double y{Position().GetY()};
                double rr{x * x + y + y};
                if(rr <= _radius_ * _radius_)
                {
                    std::cout << __func__ << " return true; point=" << point << std::endl;
                    return true;
                }
            }

            std::cout << __func__ << " return false; point=" << point << std::endl;
            return false;
        }

        private:

        // position is position of one end
        // direction is direction vector from the position vector
        // length is the length along the direction vector the cylinder
        // extends
        // radius is the radius of the cylinder
        vector3<double> _direction_;
        double _length_;
        double _radius_;

    };


    class OpenCylinder : public Cylinder
    {

        public:

        OpenCylinder()
            : Cylinder()
            , _bore_radius_{0.0}
        {
        }

        OpenCylinder(const vector3<double> &position, const vector3<double> &direction, const double length, const double radius, const double bore_radius)
            : Cylinder(position, direction, length, radius)
            , _bore_radius_{bore_radius}
        {
        }

        double GetThickness() const
        {
            return Cylinder::Radius() - _bore_radius_;
        }

        double GetBoreRadius() const
        {
            return _bore_radius_;
        }


        private:

        // radius of bore hole
        double _bore_radius_;

    };


    // cube class not finished
    class Cube : public Volume
    {
        
        public:

        Cube()
        {
        }

        Cube(const vector3<double>& size)
            : _size_{size}
            , Volume(vector3<double>(0.0, 0.0, 0.0))
            , _direction_x_{vector3<double>(1.0, 0.0, 0.0)}
            , _direction_y_{vector3<double>(0.0, 1.0, 0.0)}
        {
        }

        Cube(const vector3<double> &size, const vector3<double>& position, const vector3<double>& direction_x, const vector3<double>& direction_y)
            : _size_{size}
            , Volume(position)
            , _direction_x_{direction_x}
            , _direction_y_{direction_y}
        {
        }

        void Init(const vector3<double>& size, const vector3<double>& position, const vector3<double>& direction_x, const vector3<double> &direction_y)
        {
            _size_ = size;
            Volume::Init(position);
            _direction_x_ = direction_x;
            _direction_y_ = direction_y;
        }

        const vector3<double>& Size() const
        {
            return _size_;
        }

        // TODO: use cross product
        const vector3<double> Direction() const
        {
            return vector3<double>(0.0, 0.0, 1.0);
        }

        bool PointIntersectionTest(vector3<double> point) const
        {

            point -= Volume::Position();
            // TODO rotation

            // check x
            if(0.0 <= point.GetX() && point.GetX() < _size_.GetX())
            {
                // check y
                if(0.0 <= point.GetY() && point.GetY() < _size_.GetY())
                {
                    // check z
                    if(0.0 <= point.GetZ() && point.GetZ() < _size_.GetZ())
                    {
                        std::cout << __func__ << " return true; point=" << point << std::endl;
                        return true;
                    }
                    else
                    {
                        std::cerr << "PointIntersectionTest: test Z failed: z=" << point.GetZ() << " range=[" << 0.0 << "," << _size_.GetZ() << "]";
                        if(!(0.0 <= point.GetZ()))
                            std::cerr << " !0.0 <= " << point.GetZ();
                        if(!(point.GetZ() < _size_.GetZ()))
                            std::cerr << " !" << point.GetZ() << " < " << _size_.GetZ();
                        std::cerr << std::endl;
                    }
                }
                else
                {
                    std::cerr << "PointIntersectionTest: test Y failed: y=" << point.GetY() << " range=[" << 0.0 << "," << _size_.GetY() << "]";
                    if(!(0.0 <= point.GetY()))
                        std::cerr << " !0.0 <= " << point.GetY();
                    if(!(point.GetY() < _size_.GetY()))
                        std::cerr << " !" << point.GetY() << " < " << _size_.GetY();
                    std::cerr << std::endl;
                }
            }
            else
            {
                std::cerr << "PointIntersectionTest: test X failed: x=" << point.GetX() << " range=[" << 0.0 << "," << _size_.GetX() << "]";
                if(!(0.0 <= point.GetX()))
                    std::cerr << " !0.0 <= " << point.GetX();
                if(!(point.GetX() < _size_.GetX()))
                    std::cerr << " !" << point.GetX() << " < " << _size_.GetX();
                std::cerr << std::endl;
            }

            std::cout << __func__ << " return false; point=" << point << std::endl;
            return false;

        }
        
        //void Rotate

        private:

        // position is the position of the "origin" corner of the cube
        // _direction_x_ is the unit vector pointing along the x axis of the
        // cube from the position vector (before any rotation / translation
        // this vector also points along the x axis of the "world")
        // _direction_y_ is the same but for the y axis
        vector3<double> _direction_x_;
        vector3<double> _direction_y_;
        vector3<double> _size_;
    };


} // namespace Geometry

#endif // GEOMETRY_HPP
