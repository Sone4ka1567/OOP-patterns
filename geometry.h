#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo>

const double accuracy = 0.0001;
using std::abs;
#define pi M_PI
struct Point {
    double x = 0.0;
    double y = 0.0;
    explicit Point(double _x = 0.0, double _y = 0.0) : x(_x), y(_y) {};
    double distance(const Point& another) const; //redo
    void rotate(const Point& center, double angle);
    void scale(const Point& p, double k);
    bool Between(Point first, Point second) const;
};

bool operator==(const Point& a, const Point& b) {
    if(abs(a.x - b.x) < accuracy && abs(a.y - b.y) < accuracy) return true;
    return false;
}

bool operator!=(const Point& a, const Point& b) {
    return !(a == b);
}

double Point::distance(const Point &another) const {
    return sqrt((x - another.x) * (x - another.x) + (y - another.y) * (y - another.y));
}

void Point::rotate(const Point &center, double angle) {
    Point copy = Point(this->x - center.x, this->y - center.y);
    double new_x = copy.x * cos(angle) - copy.y * sin(angle);
    double new_y = copy.x * sin(angle) + copy.y * cos(angle);
    copy.x = new_x;
    copy.y = new_y;
    copy.x += center.x;
    copy.y += center.y;
    (*this) = copy;
}

void Point::scale(const Point &point, double k) {
    double new_x = x - point.x;
    double new_y = y - point.y;
    x = point.x + new_x * k;
    y = point.y + new_y * k;
}


class Line {
public:
    double a = 0.0 , b = 0.0, c = 0.0; // ax + by + c = 0
    Line() = default;
    Line(double a, double b, double c) : a(a), b(b), c(c) {}
    Line(const Point& p1, const Point& p2); //задать двумя точками.
    Line(const double& angle, const double& move); // y = kx + b
    Line(const Point& p, const double& angle); //задать точкой и тангенсом угла наклона
    friend bool operator==(const Line& l1, const Line& l2);
    Line perp(const Point& p) const;
    bool containsPoint (Point point) const ;
};


Line::Line(const Point& p1,const Point& p2) {
    a = p2.y - p1.y;
    b = p1.x - p2.x;
    c = p1.y * (p2.x - p1.x) - p1.x * (p2.y - p1.y);
}

Line::Line(const double& angle, const double& move) {
    a = angle;
    b = -1;
    c = move;
}

Line::Line(const Point& p, const double& angle) {
    a = angle;
    b = -1;
    c = p.y - p.x * angle;
}

bool operator==(const Line& l1, const Line& l2) {
    if(abs(l1.a * l2.b - l2.a * l1.b) < accuracy && abs(l1.b * l2.c - l2.b * l1.c) < accuracy && abs(l1.a * l2.c - l2.a * l1.c) < accuracy) return true;
    return false;
}

bool operator!=(const Line& l1, const Line& l2) {
    return !(l1 == l2);
}

Point cross(const Line& l1, const Line& l2) {
    double new_x = 0.0, new_y = 0.0;
    if((l1.b == 0.0 && l2.b == 0.0) || abs(l1.a / l1.b - l2.a / l2.b) < accuracy)
        std::cout << "parallel";
    else if(l1.b == 0) {
        new_x = -l1.c / l1.a;
        new_y = (-l2.a) * (new_x) / l2.b - l2.c / l2.b;
    }
    else if(l2.b == 0) {
        new_x = -l2.c / l2.a;
        new_y = (-l1.a) * (new_x) / l1.b - l1.c / l1.b;
    }
    else {
        new_x = (l1.c * l2.b / l1.b - l2.c) / (l2.a - l1.a * l2.b / l1.b);
        new_y = (-l1.a) * (new_x) / l1.b - l1.c / l1.b;
    }
    return Point(new_x, new_y);
}

Line Line::perp(const Point& p) const { //AX + BY + C = 0;
    Line result = Line(-b, a, -a * p.y + b * p.x);
    return result;
}

bool Line::containsPoint (Point point) const {
    return (abs(a * point.x + b * point.y + c) < accuracy);
}

double AngleBetweenVectors(double x1, double y1, double x2, double y2) {
    return (x1 * x2 + y1 * y2) / (sqrt(x1 * x1 + y1 * y1) * sqrt(x2 * x2 + y2 * y2));
}

Point ReflexPoint(const Point& point, const Line& l) {
    Line p = l.perp(point);
    Point h = cross(l, p);

    double new_x = h.x + (h.x - point.x);
    double new_y = h.y + (h.y - point.y);

    Point new_p(new_x, new_y);
    return new_p;
}

bool Point::Between(Point first, Point second) const {
    Line l(first, second);
    if (!l.containsPoint((*this))) return false;
    if (first.x >= second.x + accuracy) {
        Point temp = first;
        first = second;
        second = temp;
    }
    if ((x >= second.x - accuracy) || (x <= first.x + accuracy)) return false;
    if ((y <= second.y + accuracy) && (y >= first.y - accuracy)) return true;
    if ((y <= first.y + accuracy) && (y >= second.y - accuracy)) return true;
    return false;
}



class  Shape {
public:
    Shape () {};
    virtual ~Shape () {};
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool operator!=(const Shape& another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual  bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    virtual void rotate(Point center, double angle) = 0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
};


class Polygon: public Shape {
protected:
    std::vector<Point> vertices;
public:
    Polygon() = default;
    Polygon(std::vector<Point> _vertices): vertices(std::move(_vertices)) {};
    Polygon(std::initializer_list<Point> _vertices): vertices(_vertices) {};
    long long verticesCount() const; //сколько вершин
    std::vector<Point> getVertices() const; //вывести вершины
    bool isConvex() const; //выпуклый ли


    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape& another) const override;
    bool operator!=(const Shape& another) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;
    void rotate(Point center, double angle) override;
    void reflex(Point center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;
};

long long Polygon::verticesCount() const {
    return vertices.size();
}

std::vector<Point> Polygon::getVertices() const {
    return vertices;
}

bool Polygon::isConvex() const {
    int sz = vertices.size();
    for(int i = 0; i < sz; ++i) {
        int next_i = (i + 1) % sz;
        int which_side = 1;
        Line checker = Line(vertices[i], vertices[next_i]);
        if(checker.a * vertices[(i + 2) % sz].x + checker.b * vertices[(i + 2) % sz].y + checker.c < 0) which_side = -1;
        for(int j = i + 2; j <= sz + i - 1; ++j) {
            int jj = j % sz;
            int side = 1;
            if(checker.a * vertices[jj].x + checker.b * vertices[jj].y + checker.c < 0) which_side = -1;
            if(side != which_side) return false;
        }
    }
    return true;
}


double Polygon::perimeter() const {
    double ans = 0.0;
    int sz = vertices.size();
    for(int i = 0; i < sz; ++i) {
        int next_i = (i + 1) % sz;
        ans += vertices[i].distance(vertices[next_i]);
    }
    return ans;
}

double Polygon::area() const {
    double ans = 0.0;
    int sz = vertices.size();
    for(int i = 0; i < sz; ++i) {
        int next_i = (i + 1) % sz;
        ans += (vertices[i].x + vertices[next_i].x) * (vertices[i].y - vertices[next_i].y) / 2;
    }
    return abs(ans);
}


bool Polygon::operator==(const Shape &another) const {
    Polygon first = (*this);
    try {
        Polygon second = dynamic_cast<const Polygon&>(another);
        if (second.verticesCount() != first.verticesCount()) return false;

        size_t start_ind = 0;
        while (second.vertices[start_ind++] != first.vertices[0])
            if (start_ind == second.vertices.size()) return false;

        --start_ind;
        bool to_right = true, to_left = true;
        int sz = first.vertices.size();
        for (int i = 0; i < sz; i++)
            if (first.vertices[i] != second.vertices[(i + start_ind) % sz]) {
                to_right = false;
                break;
            }

        for (int i = 0; i < sz; i++)
            if (first.vertices[i] != second.vertices[(sz - i + start_ind) % sz]) {
                to_left = false;
                break;
            }

        return (to_left || to_right);
    }
    catch (...) {
        return false;
    }
}


bool Polygon::operator!=(const Shape &another) const {
    return !(this->operator==(another));
}

bool Polygon::isCongruentTo(const Shape &another) const {
    Polygon first = (*this);
    try {
        Polygon second = dynamic_cast<const Polygon&> (another);
        if(first.verticesCount() != second.verticesCount()) return false;
        size_t sz = first.verticesCount();
        std::vector<std::pair<double, double>> sides_and_angles_first(sz);
        std::vector<std::pair<double, double>> sides_and_angles_second(sz);
        for(size_t i = 0; i < sz; ++i) {
            size_t next_i = (i + 1) % sz;
            size_t prev_i = (i + sz - 1) % sz;
            double dist1 = first.vertices[i].distance(first.vertices[next_i]);
            double dist2 = second.vertices[i].distance(second.vertices[next_i]);
            double angle1 = AngleBetweenVectors(first.vertices[next_i].x - first.vertices[i].x, first.vertices[next_i].y - first.vertices[i].y, first.vertices[prev_i].x - first.vertices[i].x, first.vertices[prev_i].y - first.vertices[i].y);
            double angle2 = AngleBetweenVectors(second.vertices[next_i].x - second.vertices[i].x, second.vertices[next_i].y - second.vertices[i].y, second.vertices[prev_i].x - second.vertices[i].x, second.vertices[prev_i].y - second.vertices[i].y);
            sides_and_angles_first[i] = {dist1, angle1};
            sides_and_angles_second[i] = {dist2, angle2};
        }

        for(size_t start_index = 0; start_index < sz; ++start_index) {
                bool to_left = true, to_right = true;
                for(size_t j = 0; j < sz; ++j) { //to_right
                    if(abs(sides_and_angles_first[j].first - sides_and_angles_second[(start_index + j) % sz].first) >= accuracy || abs(sides_and_angles_first[j].second - sides_and_angles_second[(start_index + j) % sz].second) >= accuracy)
                        to_right = false;
                }
                for(size_t j = 0; j < sz; ++j) { //to_left
                    if(abs(sides_and_angles_first[j].first - sides_and_angles_second[(start_index - 1 - j + sz) % sz].first) >= accuracy || abs(sides_and_angles_first[j].second - sides_and_angles_second[(start_index - j + sz) % sz].second) >= accuracy)
                        to_left = false;
                }
                if(to_left || to_right)
                    return true;
        }
        return false;
    }
    catch(...) {
        return false;
    }
}

bool Polygon::isSimilarTo(const Shape &another) const {
    Polygon first = (*this);
    try {
        Polygon second = dynamic_cast<const Polygon&> (another);
        if(first.verticesCount() != second.verticesCount()) return false;
        size_t sz = first.verticesCount();
        std::vector<double> check1, check2; //delete
        std::vector<std::pair<double, double>> sides_and_angles_first(sz);
        std::vector<std::pair<double, double>> sides_and_angles_second(sz);
        for(size_t i = 0; i < sz; ++i) {
            size_t next_i = (i + 1) % sz;
            size_t prev_i = (i + sz - 1) % sz;
            double dist1 = first.vertices[i].distance(first.vertices[next_i]);
            double dist2 = second.vertices[i].distance(second.vertices[next_i]);
            double angle1 = AngleBetweenVectors(first.vertices[next_i].x - first.vertices[i].x, first.vertices[next_i].y - first.vertices[i].y, first.vertices[prev_i].x - first.vertices[i].x, first.vertices[prev_i].y - first.vertices[i].y);
            double angle2 = AngleBetweenVectors(second.vertices[next_i].x - second.vertices[i].x, second.vertices[next_i].y - second.vertices[i].y, second.vertices[prev_i].x - second.vertices[i].x, second.vertices[prev_i].y - second.vertices[i].y);
            sides_and_angles_first[i] = {dist1, angle1};
            sides_and_angles_second[i] = {dist2, angle2};
        }

        for(size_t start_index = 0; start_index < sz; ++start_index) {
            double coefficient = sides_and_angles_first[0].first / sides_and_angles_second[start_index].first;
            bool to_left = true, to_right = true;
            for(size_t j = 0; j < sz; ++j) { //to_right
                if(abs(sides_and_angles_first[j].first / sides_and_angles_second[(start_index + j) % sz].first - coefficient) >= accuracy || abs(sides_and_angles_first[j].second - sides_and_angles_second[(start_index + j) % sz].second) >= accuracy)
                    to_right = false;
            }
            coefficient = sides_and_angles_first[0].first / sides_and_angles_second[(start_index - 1) % sz].first;
            for(size_t j = 0; j < sz; ++j) { //to_left
                if(abs(sides_and_angles_first[j].first / sides_and_angles_second[(start_index - 1 - j + sz) % sz].first - coefficient) >= accuracy || abs(sides_and_angles_first[j].second - sides_and_angles_second[(start_index - j + sz) % sz].second) >= accuracy)
                    to_left = false;
            }
            if(to_left || to_right) {
                return true;
            }
        }
        return false;
    }
    catch(...) {
        return false;
    }
}

bool Polygon::containsPoint(Point point) const {
    double angle = 0.0;
    int sz = vertices.size();
    for (int i = 0; i < sz; i++) {
        Point next_vertice = vertices[(i + 1) % sz];
        Point this_vertice = vertices[i];
        if (point == this_vertice || point.Between(next_vertice, this_vertice)) return true;
        double vect = (this_vertice.x - point.x) * (next_vertice.y - point.y) - (this_vertice.y - point.y) * (next_vertice.x - point.x);
        double scal = (this_vertice.x - point.x) * (next_vertice.x - point.x) + (this_vertice.y - point.y) * (next_vertice.y - point.y);
        angle += atan2 (vect, scal);
    }
    return (abs(abs(angle) - 2 * pi) < accuracy);
}

void Polygon::rotate(Point center, double angle) {
    angle *= pi / 180;
    for(size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].rotate(center, angle);
    }
}

void Polygon::reflex(Point center) {
    for(size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].rotate(center, pi);
    }
}

void Polygon::reflex(Line axis) {
    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i] = ReflexPoint(vertices[i], axis);
    }
}

void Polygon::scale(Point center, double coefficient) {
    for(size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].scale(center, coefficient);
    }
}




class Ellipse: public Shape {
protected:
    Point f1, f2;
    double distance = 0.0;// сумма расстояний от точки эллипса до них
public:
    Ellipse() = default;
    Ellipse(const Point& _f1, const Point& _f2, double _distance) : f1(_f1), f2(_f2),  distance(_distance) {};
    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    virtual Point center() const;


    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape& another) const override;
    bool operator!=(const Shape& another) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;
    void rotate(Point center, double angle) override;
    void reflex(Point center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;
};

std::pair<Point, Point> Ellipse::focuses() const{
    return std::make_pair(f1, f2);
}

std::pair<Line, Line> Ellipse::directrices() const{
    double new_x = distance / (2 * eccentricity()) ;
    Point cent = this->center();
    Point vec1 = Point((f1.x - cent.x) * new_x / f1.distance(cent) + cent.x, (f1.y - cent.y) * new_x / f1.distance(cent) + cent.y);
    Point vec2 = Point((f2.x - cent.x) * new_x / f2.distance(cent) + cent.x, (f2.y - cent.y) * new_x / f2.distance(cent) + cent.y);
    return std::make_pair(Line(f1,f2).perp(vec1), Line(f1,f2).perp(vec2));
}

double Ellipse::eccentricity() const{
    return f1.distance(f2) / distance;
}

Point Ellipse::center() const {
    Point c;
    c.x = (f1.x + f2.x) / 2;
    c.y = (f1.y + f2.y) / 2;
    return c;
}

double Ellipse::perimeter() const { // pi(3(a+b) - sqrt((3a + b)(a + 3b)))
    double a = distance / 2;
    double c = f1.distance(f2) / 2;
    double b = sqrt(a * a - c * c);
    double ans = pi * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b)));
    return ans;
}


double Ellipse::area() const { // a * b * pi
    double a = distance / 2;
    double c = f1.distance(f2) / 2;
    double b = sqrt(a * a - c * c);
    double ans = a * b * pi;
    return ans;
}


bool Ellipse::operator==(const Shape &another) const {
    Ellipse first = (*this);
    try {
        Ellipse second = dynamic_cast<const Ellipse&>(another);
        bool flag1 = (first.f1 == second.f1) && (first.f2 == second.f2);
        bool flag2 = (first.f2 == second.f1) && (first.f1 == second.f2);
        return ((flag1 || flag2) && (first.distance - second.distance) < accuracy);
    }
    catch(...) {
        return false;
    }
}


bool Ellipse::operator!=(const Shape &another) const {
    return !(this->operator==(another));
}


bool Ellipse::isCongruentTo(const Shape &another) const { //сравнить distance и расстояние между фокусами
    try {
        Ellipse second = dynamic_cast<const Ellipse&> (another);
        double c_first = f1.distance(f2);
        double c_second = second.f1.distance(second.f2);
        return (abs(distance - second.distance) < accuracy) &&
               (abs(c_first - c_second) < accuracy);
    }
    catch (...) {
        return false;
    }
}

bool Ellipse::isSimilarTo(const Shape &another) const {
    try {
        Ellipse second = dynamic_cast<const Ellipse&>(another);
        double k = distance / second.distance;
        double c_first = f1.distance(f2);
        double c_second = second.f1.distance(second.f2);
        return abs( k * c_second - c_first) < accuracy ? true : false;
    }
    catch (...) {
        return false;
    }
}


bool Ellipse::containsPoint(Point point) const {
    double supposed_dist = 0.0;
    supposed_dist += point.distance(f1);
    supposed_dist += point.distance(f2);
    return supposed_dist < distance + accuracy;
}


void Ellipse::rotate(Point center, double angle) {
    angle *= pi / 180;
    f1.rotate(center, angle);
    f2.rotate(center, angle);
}

void Ellipse::reflex(Point center) {
    f1.rotate(center, pi);
    f2.rotate(center, pi);
}

void Ellipse::reflex(Line axis) {
    f1 = ReflexPoint(f1, axis);
    f2 = ReflexPoint(f2, axis);
}

void Ellipse::scale(Point center, double coefficient) {
    f1.scale(center, coefficient);
    f2.scale(center, coefficient);
    distance = distance * coefficient;
}



class Circle: public Ellipse {
public:
    Point center_;
    double rad = 0.0;

    Circle() = default;
    Circle(Point _center, double _rad): Ellipse(), center_(_center), rad(_rad) {
        f1 = f2 = _center;
        distance = _rad * 2;
    };

    double radius() const {
        return distance / 2; //return rad
    }

    Point center() const override {
        return center_;
    };
};


class Rectangle: public Polygon {
protected:
    Point a, b, c, d;
public:
    Rectangle() = default;
    Rectangle(Point a, Point d, double k);
    Point center() const;
    std::pair<Line, Line> diagonals();
};

Rectangle::Rectangle(Point _a, Point _c, double k):Polygon() {
    a = _a; //ab - меньшая стороны
    c = _c;
    double diag = a.distance(c); //AB^2 + BC^2 = DIAG^2
    if(k < 1) k = 1 / k; //bc = k * ab - > ab^2 + ab^2 * k ^ 2 = diag^2-> ab^2 = diag^2 / (1 + k^2)
    double ab = sqrt(diag * diag / (1 + k * k));
    Point vec_ac = Point(c.x - a.x, c.y - a.y);
    Point a_new_b;
    a_new_b.x = vec_ac.x * ab / diag + a.x;
    a_new_b.y = vec_ac.y * ab / diag + a.y;
    k = atan(k);
    a_new_b.rotate(a, k);
    b = a_new_b;
    Point vec_ca = Point( a.x - c.x, a.y - c.y);
    Point c_new_d;
    c_new_d.x = vec_ca.x * ab / diag + c.x;
    c_new_d.y = vec_ca.y * ab / diag + c.y;
    c_new_d.rotate(c, k);
    d = c_new_d;
    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(c);
    vertices.push_back(d);
}

Point Rectangle::center() const {
    Point cen;
    cen.x = (a.x + c.x) / 2;
    cen.y = (a.y + c.y) / 2;
    return cen;
}

std::pair<Line, Line> Rectangle::diagonals() {
    Line d1 = Line(a, c);
    Line d2 = Line(b, d);
    return std::make_pair(d1, d2);
}


class Square: public Rectangle {
private:
    Point a, b, c, d;
public:
    Square() = default;
    Square(Point _a, Point _c);
    Circle circumscribedCircle();
    Circle inscribedCircle();
};


Square::Square(Point _a, Point _c) {
    a = _a;
    c = _c;
    double diag = a.distance(c); //AB^2 + BC^2 = DIAG^2
    double ab = diag / sqrt(2);
    Point vec_ac = Point(c.x - a.x, c.y - a.y);
    Point a_new_b;
    a_new_b.x = vec_ac.x * ab / diag + a.x;
    a_new_b.y = vec_ac.y * ab / diag + a.y;
    a_new_b.rotate(a, pi / 4);
    b = a_new_b;
    Point vec_ca = Point( a.x - c.x, a.y - c.y);
    Point c_new_d;
    c_new_d.x = vec_ca.x * ab / diag + c.x;
    c_new_d.y = vec_ca.y * ab / diag + c.y;
    c_new_d.rotate(c, pi / 4);
    d = c_new_d;
    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(c);
    vertices.push_back(d);
}

Circle Square::circumscribedCircle() {
    Circle circ;
    circ.center_ = Point((a.x + c.x) / 2, (a.y + c.y) / 2);
    circ.rad = a.distance(circ.center_);
    return circ;
}

Circle Square::inscribedCircle() {
    Circle circ;
    circ.center_ = Point((a.x + c.x) / 2, (a.y + c.y) / 2);
    circ.rad = a.distance(b);
    return circ;
}



class Triangle: public Polygon {
private:
    Point a, b , c ;
public:
    Triangle() = default;
    Triangle(Point _a, Point _b, Point _c): a(_a), b(_b), c(_c) {
        vertices.push_back(a);
        vertices.push_back(b);
        vertices.push_back(c);
    };
    Point circumscribedCircleCenter() const;
    Circle circumscribedCircle(); // Описанная
    Circle inscribedCircle() const;
    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;
};

Point Triangle::circumscribedCircleCenter() const {
    double d = (a.x - b.x) * (c.y - a.y) - (a.y - b.y)*(c.x - a.x);
    double _x = (a.y - b.y) * (c.x * c.x + c.y * c.y) + (b.y - c.y) * (a.x * a.x  + a.y * a.y) + (c.y - a.y) * (b.x * b.x + b.y * b.y);
    double _y = (a.x - b.x) * (c.x * c.x + c.y * c.y) + (b.x - c.x) * (a.x * a.x  + a.y * a.y) + (c.x - a.x) * (b.x * b.x + b.y * b.y);
    double new_x = - _x / (2 * d);
    double new_y = _y / (2 * d);
    Point center = Point(new_x, new_y);
    return center;
}

Circle Triangle::circumscribedCircle() {
    double r = a.distance(b) * b.distance(c) * c.distance(a) / (4 * this->area());
    Circle ans = Circle(this->circumscribedCircleCenter(), r);
    return ans;
}

Circle Triangle::inscribedCircle() const {
    double r = 2 * this->area() / this->perimeter();
    double ab = a.distance(b);
    double ac = a.distance(c);
    double bc = b.distance(c);
    double new_x = (bc * a.x + ac * b.x + ab * c.x) / (ab + ac + bc);
    double new_y = (bc * a.y + ac * b.y + ab * c.y) / (ab + ac + bc);
    Point cen = Point(new_x, new_y);
    Circle ans = Circle(cen, r);
    return ans;
}

Point Triangle::centroid() const {
    Point _centroid = Point ((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3);
    return _centroid;
}

double det(double a, double b, double c, double d, double e, double f, double g, double h, double i) { // a b c
    double first_part = a * e * i +  b * f * g + c * d * h;                                           //  d e f
    double second_part = c * e * g + b * d * i + a * f * h;                                           //  g h i
    return first_part - second_part;
}


Point Triangle::orthocenter() const { //через определители матриц
    double det_x  = det(a.y, b.x * c.x + a.y * a.y, 1.0,
                        b.y, c.x * a.x + b.y * b.y, 1.0,
                        c.y, a.x * b.x + c.y * c.y, 1.0);

    double det_xy = det(a.x, a.y, 1.0,
                        b.x, b.y, 1.0,
                        c.x, c.y, 1.0);

    double det_y  = det(b.y * c.y + a.x * a.x, a.x, 1.0,
                        c.y * a.y + b.x * b.x, b.x, 1.0,
                        a.y * b.y + c.x * c.x, c.x, 1.0);

    return Point (det_x / det_xy, det_y / det_xy);
}

Line Triangle::EulerLine() const {
    Line euler = Line(this->orthocenter(), this->circumscribedCircleCenter());
    return euler;
}

Circle Triangle::ninePointsCircle() const {
    Point ortho = this->orthocenter();
    Point circum_center = this->circumscribedCircleCenter();
    double  new_x = (ortho.x + circum_center.x) / 2;
    double  new_y = (ortho.y + circum_center.y) / 2;
    Point cent = Point(new_x, new_y);
    double r = a.distance(b) * b.distance(c) * c.distance(a) / (4 * this->area()) / 2;
    Circle circ = Circle(cent, r);
    return circ;
}
