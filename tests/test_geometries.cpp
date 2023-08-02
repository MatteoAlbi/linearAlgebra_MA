#include <gtest/gtest.h>
#include "linear_algebra_ma/geometries/point.hpp"
#include "linear_algebra_ma/geometries/segment.hpp"

using namespace MA::geometries;

using std::cout, std::endl, std::string;
using std::invalid_argument, std::runtime_error, std::out_of_range;
using pdd = std::pair<double, double>;

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
    EXPECT_EQ(p1.distance(p2), 3);
    p2.y(6);
    EXPECT_EQ(p1.distance(p2), 5);
    p2.y(-2);
    EXPECT_EQ(p1.distance(p2), 5);

    p2.x(NAN);
    EXPECT_THROW(p1.distance(p2), invalid_argument);
    p2.x(INFINITY);
    EXPECT_THROW(p1.distance(p2), invalid_argument);
}


TEST(Segment, constructor_get) {
    EXPECT_NO_THROW(Segment s);
    Segment s;
    EXPECT_TRUE(s.p1().isnan());
    EXPECT_TRUE(s.p2().isnan());

    EXPECT_NO_THROW(s = Point(2,3));
    EXPECT_EQ(s.x(), 2);
    EXPECT_EQ(s.y(), 3);
    EXPECT_NO_THROW(s = Point(9,9));
    EXPECT_EQ(s.get(), (pdd{9,9}));
}