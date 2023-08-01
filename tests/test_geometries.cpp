#include <gtest/gtest.h>
#include "linear_algebra_ma/geometries/point.hpp"
#include "linear_algebra_ma/geometries/line.hpp"

TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}