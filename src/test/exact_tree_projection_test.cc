#include "../core/exact_tree_projection.h"

#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "test_helpers.h"

using std::vector;
using tree_sparsity::exact_tree_projection;

void RunAlgo(const vector<double>& x,
             size_t degree,
             size_t sparsity,
             const vector<bool>& expected_result) {
  std::unique_ptr<bool[]> result(new bool[x.size()]);
  //ASSERT_TRUE(exact_tree_projection(x.data(), x.size(), degree, sparsity, 
  exact_tree_projection(x.data(), x.size(), degree, sparsity, result.get());
  //    result.data());
  vector<bool> res(x.size());
  std::copy(result.get(), result.get() + x.size(), res.begin());
  CheckResult(expected_result, res);
}

// d = 2

TEST(ExactTreeProjectionTest, SimpleBinaryTest) {
  const double x2[] = {1, 1, 0, 1, 1, 0, 0};
  const bool res2[] = {1, 1, 0, 1, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 4, res);
}

TEST(ExactTreeProjectionTest, SimpleEmptyNodeTest) {
  const double x2[] = {1, 0, 0, 1, 1, 0, 0};
  const bool res2[] = {1, 1, 0, 1, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 4, res);
}

TEST(ExactTreeProjectionTest, SimpleDecisionTest) {
  const double x2[] = {1, 0, 0, 2, 3, 0, 0};
  const bool res2[] = {1, 1, 0, 0, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 3, res);
}

TEST(ExactTreeProjectionTest, SimpleNoDecisionTest) {
  const double x2[] = {1, 0, 0, 0, 0, 2, 3};
  const bool res2[] = {1, 0, 1, 0, 0, 1, 1};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 4, res);
}

TEST(ExactTreeProjectionTest, SimpleAlmostFullTest) {
  const double x2[] = {1, 1, 1, 1, 1, 1, 0};
  const bool res2[] = {1, 1, 1, 1, 1, 1, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 6, res);
}

TEST(ExactTreeProjectionTest, NotFullTree) {
  // Layout the memory so that the non-existing child has a large value.
  const double x2[] = {1, 0, 1, 0, 0, 1, 100};
  vector<double> x(begin(x2), end(x2));
  x.resize(6);
  const bool res2[] = {1, 0, 1, 0, 0, 1};
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 3, res);
}

TEST(ExactTreeProjectionTest, MediumSizeTest) {
  const double x2[] = {50, 97, 50, 43, 42, 45, 75, 73, 48, 50, 78, 57, 89,
                       69, 38, 17, 52, 31, 81, 68,  8, 96, 41, 79, 62, 69, 54,
                       51, 82, 94, 38};
  const bool res2[] = {1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0,
                       1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                       0, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 9, res);
}


// d = 3

TEST(ExactTreeProjectionTest, SimpleTest3) {
  const double x2[] = {1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0};
  const bool res2[] = {1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 3, 4, res);
}

TEST(ExactTreeProjectionTest, SimpleDecisionTest3) {
  const double x2[] = {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0};
  const bool res2[] = {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 3, 4, res);
}
