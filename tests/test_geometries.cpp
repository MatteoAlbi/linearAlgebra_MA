#include <gtest/gtest.h>
#include "linear_algebra_ma/geometries/point.hpp"
#include "linear_algebra_ma/geometries/segment.hpp"

using namespace MA::geometries;

using std::cout, std::endl, std::string, std::pair;
using std::invalid_argument, std::runtime_error, std::out_of_range;
using pdd = pair<double, double>;
using ppp = pair<Point, Point>;

// custom equal operator
namespace std{
    template<typename T>
    bool operator==(pair<const T &, const T &> p1, pair<T,T> p2){return p1.first==p2.first && p1.second==p2.second;}
}


TEST(infinity, expected_behavior){
    /**
     * if one of this test fails, following 
     * functions will have undefined behavior:
     * Point::distance
     * Point::slope
     * distance(Point,Point)
     * Segment::length
     * Segment::slope
    */

    EXPECT_TRUE(INFINITY > 0);
    EXPECT_FALSE(INFINITY < 0);
    EXPECT_TRUE(-INFINITY < 0);
    EXPECT_FALSE(-INFINITY > 0);
    EXPECT_TRUE((-INFINITY)*INFINITY < 0);
    EXPECT_TRUE(INFINITY*INFINITY > 0);
    EXPECT_TRUE((-INFINITY)*(-INFINITY) > 0);
}


TEST(Point, constructor_get) {
    EXPECT_NO_THROW(Point p);
    Point p;
    EXPECT_TRUE(std::isnan(p.x()));
    EXPECT_TRUE(std::isnan(p.y()));

    EXPECT_NO_THROW(p = Point(2,3));
    EXPECT_EQ(p.x(), 2);
    EXPECT_EQ(p.y(), 3);
    EXPECT_NO_THROW(p = Point(9,9));
    EXPECT_EQ(p.get(), (pdd{9,9}));
}

TEST(Point, set) {
    Point p;

    p.x(5);
    p.y(7);
    EXPECT_EQ(p.x(), 5);
    EXPECT_EQ(p.y(), 7);
    p.set(1,1);
    EXPECT_EQ(p.get(), (pdd{1,1}));
}

TEST(Point, isnan){
    Point p;
    EXPECT_TRUE(p.isnan());
    p.x(3);
    EXPECT_TRUE(p.isnan());
    p.y(INFINITY);
    EXPECT_FALSE(p.isnan());
    p.x(NAN);
    EXPECT_TRUE(p.isnan());
}

TEST(Point, isinf){
    Point p(2,2);
    EXPECT_FALSE(p.isinf());
    p.x(INFINITY);
    EXPECT_TRUE(p.isinf());
    p.y(INFINITY);
    EXPECT_TRUE(p.isinf());
    p.x(NAN);
    EXPECT_TRUE(p.isinf());
}

TEST(Point, equal_operator){
    Point p1(2,2), p2(2,2);
    EXPECT_TRUE(p1==p2);
    p2.x(3);
    EXPECT_FALSE(p1==p2);

    // NAN == NAN is evaluated as false
    p2.x(NAN);
    p1.x(NAN);
    EXPECT_FALSE(p1==p2);

    // INF == INF is evaluated as true
    p2.x(INFINITY);
    p1.x(INFINITY);
    EXPECT_TRUE(p1==p2);
    p2.x(-INFINITY);
    p1.x(-INFINITY);
    EXPECT_TRUE(p1==p2);
    // there is distintion between +/- INF
    p2.x(INFINITY);
    p1.x(-INFINITY);
    EXPECT_FALSE(p1==p2);
}

TEST(Point, unequal_operator){
    Point p1(2,2), p2(2,2);
    EXPECT_FALSE(p1!=p2);
    p2.x(3);
    EXPECT_TRUE(p1!=p2);

    // NAN != NAN is evaluated as true
    p2.x(NAN);
    p1.x(NAN);
    EXPECT_TRUE(p1!=p2);

    // INF != INF is evaluated as false
    p2.x(INFINITY);
    p1.x(INFINITY);
    EXPECT_FALSE(p1!=p2);
    p2.x(-INFINITY);
    p1.x(-INFINITY);
    EXPECT_FALSE(p1!=p2);
    // there is distintion between +/- INF
    p2.x(INFINITY);
    p1.x(-INFINITY);
    EXPECT_TRUE(p1!=p2);
}

TEST(Point, distance){
    Point p1(2,2), p2(2,2);
    EXPECT_EQ(p1.distance(p2), 0);
    p2.x(5);
    EXPECT_EQ(p2.distance(p1), 3);
    p2.y(6);
    EXPECT_EQ(p1.distance(p2), 5);
    p2.y(-2);
    EXPECT_EQ(distance(p2,p1), 5);

    // only one point has infinite coordinate
    p2.x(INFINITY);
    EXPECT_EQ(distance(p2,p1), INFINITY);
    // points have different infinite coordinates
    p1.y(INFINITY);
    EXPECT_EQ(p1.distance(p2), INFINITY);
    // points have same infinite coordinate with same sign
    p1.y(7);
    p1.x(INFINITY);
    EXPECT_TRUE(std::isnan(distance(p2,p1)));
    // points have same infinite coordinate with opposite sign
    p2.x(-INFINITY);
    EXPECT_EQ(distance(p1,p2), INFINITY);

    p2.x(NAN);
    EXPECT_THROW(p1.distance(p2), invalid_argument);
}

TEST(Point, slope){
    Point p1(1,1), p2(1,1);

    // same points
    EXPECT_TRUE(std::isnan(p1.slope(p2)));
    p1.x() = NAN;
    // one of the point is NAN
    EXPECT_THROW(p2.slope(p1), invalid_argument);

    // no coord = inf
    p1.x(1);
    p2 = Point(2,2);
    EXPECT_EQ(p1.slope(p2), 1);
    EXPECT_EQ(p1.slope(p2), p2.slope(p1));
    p2.y() = 1;
    EXPECT_EQ(slope(p2,p1), 0);
    EXPECT_EQ(slope(p2,p1), slope(p1,p2));
    // vertical
    p2 = Point(1,2);
    EXPECT_EQ(p1.slope(p2), INFINITY);

    // both inf
        // one as both coordinate = inf
    p1 = Point(2, INFINITY);
    p2 = Point(INFINITY, INFINITY);
    EXPECT_TRUE(std::isnan(p1.slope(p2)));
        // diff coordinate = inf
    p2 = Point(INFINITY, 2);
    EXPECT_TRUE(std::isnan(p2.slope(p1)));
        // y coordinate, same sign
    p2 = Point(1, INFINITY);
    EXPECT_TRUE(std::isnan(p2.slope(p1)));
        // y coordinate, diff sign
    p1 = Point(2, -INFINITY);
    EXPECT_EQ(slope(p1,p2), INFINITY);
        // x coordinate, same sign
    p1 = Point(-INFINITY, 3);
    p2 = Point(-INFINITY, 1);
    EXPECT_TRUE(std::isnan(p2.slope(p1)));
        // x coordinate, diff sign
    p1 = Point(INFINITY, 3);
    EXPECT_EQ(slope(p2,p1), 0);

    // only first is inf
    p2 = Point(4,5);
        // both coordinate = inf
    p1 = Point(INFINITY, -INFINITY);
    EXPECT_TRUE(std::isnan(slope(p1,p2)));
        // only x
    p1 = Point(INFINITY, 3);
    EXPECT_EQ(slope(p1,p2), 0);
        // only y
    p1 = Point(6, -INFINITY);
    EXPECT_EQ(slope(p1,p2), INFINITY);

    // only second is inf
    p1 = Point(4,5);
        // both coordinate = inf
    p2 = Point(INFINITY, -INFINITY);
    EXPECT_TRUE(std::isnan(slope(p1,p2)));
        // only x
    p2 = Point(INFINITY, 3);
    EXPECT_EQ(slope(p1,p2), 0);
        // only y
    p2 = Point(6, -INFINITY);
    EXPECT_EQ(slope(p1,p2), INFINITY);

}
// test to_vec
// test sum and diff operator




TEST(Segment, constructor_get) {
    EXPECT_NO_THROW(Segment s);
    Segment s;
    EXPECT_TRUE(s.p1().isnan());
    EXPECT_TRUE(s.p2().isnan());

    EXPECT_NO_THROW(s = Segment({2,3},{4,5}));
    EXPECT_EQ(s.p1(), Point(2,3));
    EXPECT_EQ(s.p2(), Point(4,5));
    EXPECT_NO_THROW(s = Segment({6,7},{8,9}));
    EXPECT_EQ(s.get(), (ppp{{6,7},{8,9}}));
}

TEST(Segment, set) {
    Segment s;

    s.p1({1,1});
    s.p2({2,2});
    EXPECT_EQ(s.get(), (ppp{{1,1},{2,2}}));
    s.set({3,3},{4,4});
    EXPECT_EQ(s.get(), (ppp{{3,3},{4,4}}));
}

TEST(Segment, equal_operator) {
    Segment s1, s2;
    // points are NAN
    EXPECT_FALSE(s1==s2);

    s1 = Segment({1,1},{2,2});
    s2 = Segment({1,1},{2,2});
    EXPECT_TRUE(s1==s2);
    s2.p1().x(3);
    EXPECT_FALSE(s1==s2);
    // one coordinate is NAN
    s2.p1().x(NAN);
    EXPECT_FALSE(s1==s2);
    // same coordinates are NAN
    s1.p1().x(NAN);
    EXPECT_FALSE(s1==s2);
    // one coordinate is NAN
    s2.p1().x(INFINITY);
    EXPECT_FALSE(s1==s2);
    // same coordinates are +INF
    s1.p1().x(INFINITY);
    EXPECT_TRUE(s1==s2);
    // same coordinates are +/-INF
    s1.p1().x(-INFINITY);
    EXPECT_FALSE(s1==s2);
}

TEST(Segment, unequal_operator) {
    Segment s1, s2;
    // points are NAN
    EXPECT_TRUE(s1!=s2);

    s1 = Segment({1,1},{2,2});
    s2 = Segment({1,1},{2,2});
    EXPECT_FALSE(s1!=s2);
    s2.p1().x(3);
    EXPECT_TRUE(s1!=s2);
    // one coordinate is NAN
    s2.p1().x(NAN);
    EXPECT_TRUE(s1!=s2);
    // same coordinates are NAN
    s1.p1().x(NAN);
    EXPECT_TRUE(s1!=s2);
    // one coordinate is NAN
    s2.p1().x(INFINITY);
    EXPECT_TRUE(s1!=s2);
    // same coordinates are +INF
    s1.p1().x(INFINITY);
    EXPECT_FALSE(s1!=s2);
    // same coordinates are +/-INF
    s1.p1().x(-INFINITY);
    EXPECT_TRUE(s1!=s2);
}

TEST(Segment, length){
    Segment s({2,2}, {2,2});

    EXPECT_EQ(s.length(), 0);
    s.p2().x(5);
    EXPECT_EQ(s.length(), 3);
    s.p2().y(6);
    EXPECT_EQ(s.length(), 5);
    s.p2().y(-2);
    EXPECT_EQ(s.length(), 5);

    s.p2().x(NAN);
    EXPECT_THROW(s.length(), invalid_argument);
}

TEST(Segment, slope){
    Segment s({1,1},{1,1});

    // same points
    EXPECT_TRUE(std::isnan(s.slope()));
    s.p1().x() = NAN;
    // one of the point is NAN
    EXPECT_THROW(s.slope(), invalid_argument);

    s = Segment({1,1},{2,2});
    EXPECT_EQ(s.slope(), 1);
    s.p2().y() = 1;
    EXPECT_EQ(s.slope(), 0);
    s.p2() = {1,2};
    EXPECT_EQ(s.slope(), INFINITY);
}

TEST(Segment, intersection){
    // intersection at (0.95, 1.25)
    Segment p({1.8,2.1},{0.8,1.1});
    Segment q({1,1.25},{0,1.25});
    // cout << p.intersection(q) << endl;
    EXPECT_EQ(p.intersection(q), Point(0.95, 1.25));

    // no intersection
    p = Segment({-1,0.5},{1,0.5});
    q = Segment({0,1},{0,2});
    // cout << p.intersection(q) << endl;
    EXPECT_TRUE(p.intersection(q).isnan());

    // overlapping
    p = Segment({-1,1},{1,1});
    q = Segment({0,1},{2,1});
    EXPECT_EQ(p.intersection(q), Point(0, 1));
}