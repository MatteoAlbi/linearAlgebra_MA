#include <gtest/gtest.h>
#include "linear_algebra_ma/linearAlgebra.hpp"

TEST(HelloTest2, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}