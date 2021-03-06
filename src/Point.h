#ifndef Point_h
#define Point_h
#include<assert.h>
class MyPoint : public std::array<double, 2>
{
public:

    // dimension of space (or "k" of k-d tree)
    // KDTree class accesses this member
    static const int DIM = 2;

    // the constructors
    MyPoint() {}
    MyPoint(double x, double y)
    {
        (*this)[0] = x;
        (*this)[1] = y;
    }
};

struct Point {
  double x, y, z, z_r;
  int hl;
  
  Point() {
    x = 0.0;
    y = 0.0;
    z = 0.0;
    z_r = 0.0;
    hl = 0;
  }
  
  Point(const double &x, const double &y, const double &z,const double &z_r, const int &hl) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->z_r = z_r;
    this->hl = hl;
  }
  

  
  float operator[](const int &coordinate) const {
    if (coordinate == 0) return x;
    else if (coordinate == 1) return y;
    else if (coordinate == 2) return z;
    else assert(false);
  }

    bool operator<(const Point &o)  const {
        //return x < o.x || y < o.y || z < o.z;

      return x < o.x || (x == o.x && y < o.y) || (x == o.x && y == o.y && z < o.z);}
    bool const operator==(const Point &o) const{
        return x == o.x && y == o.y && z == o.z;
    }
  /*
  const Point operator+(const Point &other) const {
    return Point(x+other.x, y+other.y, z+other.z);
  }
  
  const Point operator-(const Point &other) const {
    return Point(x-other.x, y-other.y, z-other.z);
  }
  
  const Point operator*(const float &other) const {
    return Point(x*other, y*other, z*other);
  }
  
  const Point operator/(const float &other) const {
    return Point(x/other, y/other, z/other);
  }
  
  float dot(const Point &other) const {
    return x*other.x + y*other.y + z*other.z;
  }
  
  const Point cross(const Point &other) const {
    return Point(y*other.z-z*other.y, -(x*other.z-z*other.x), x*other.y-y*other.x);
  }
   */
};

struct Face{
    std::vector<unsigned int> Exterior;
    std::vector<std::vector<unsigned int>> Interior;
    int hole_no;
    Face() {
        Exterior = {};
        Interior = {};
        hole_no = 0;
    }
    Face(const std::vector<unsigned int> &Exterior, const std::vector<std::vector<unsigned int>> &Interior, const int& hole_no) {
        this->Exterior = Exterior;
        this->Interior = Interior;
        this->hole_no = hole_no;
    }
};

std::ostream& operator<<(std::ostream& os, const Point& p) {
  os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
  return os;
}

#endif /* Point_h */
