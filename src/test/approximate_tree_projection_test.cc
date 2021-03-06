#include "../core/approximate_tree_projection.h"

#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "test_helpers.h"

using std::vector;
using tree_sparsity::binsearch_options;
using tree_sparsity::approximate_tree_projection;

void RunAlgo(const vector<double>& x,
             int degree,
             int k_low,
             int k_high,
             const vector<bool>& expected_result) {
  binsearch_options opts;
  opts.layout = binsearch_options::kCompleteTree;
  opts.verbose = true;
  opts.output_function = WriteToStderr;

  std::unique_ptr<bool[]> result(new bool[x.size()]);
  double llow;
  double lhigh;
  int numiter;
  
  ASSERT_TRUE(approximate_tree_projection(x.data(), x.size(), degree, k_low,
        k_high, opts, result.get(), &llow, &lhigh, &numiter));
  vector<bool> res(x.size());
  std::copy(result.get(), result.get() + x.size(), res.begin());
  CheckResult(expected_result, res);
}

void RunAlgoWavelet(const vector<double>& x,
                    int degree,
                    int k_low,
                    int k_high,
                    const vector<bool>& expected_result) {
  binsearch_options opts;
  opts.layout = binsearch_options::kWaveletTree;
  opts.verbose = true;
  opts.output_function = WriteToStderr;

  std::unique_ptr<bool[]> result(new bool[x.size()]);
  double llow;
  double lhigh;
  int numiter;
  
  ASSERT_TRUE(approximate_tree_projection(x.data(), x.size(), degree, k_low,
        k_high, opts, result.get(), &llow, &lhigh, &numiter));
  vector<bool> res(x.size());
  std::copy(result.get(), result.get() + x.size(), res.begin());
  CheckResult(expected_result, res);
}

// d = 2

TEST(TreeapproxBinsearchTest, SimpleBinaryTest) {
  const double x2[] = {1, 1, 0, 1, 1, 0, 0};
  const bool res2[] = {1, 1, 0, 1, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, SimpleBinaryTest2) {
  const double x2[] = {10, 10, 7, 10, 10, 7, 7};
  const bool res2[] = {1, 1, 0, 1, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, EmptyParentTest) {
  const double x2[] = {1, 0, 0, 0, 0, 1, 1};
  const bool res2[] = {1, 0, 1, 0, 0, 1, 1};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, NotConvexTest) {
  const double x2[] = {0.5, 0, 0, 2, 3, 0, 0};
  const bool res2[] = {0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 3, 3, res);
}

TEST(TreeapproxBinsearchTest, NotConvexTest2) {
  const double x2[] = {0.5, 0, 0, 2, 3, 0, 0};
  const bool res2[] = {1, 1, 0, 1, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 3, 4, res);
}

TEST(TreeapproxBinsearchTest, NotFullTree) {
  // Layout the memory so that the non-existing child has a large value.
  const double x2[] = {1, 0, 1, 0, 0, 1, 100};
  vector<double> x(begin(x2), end(x2));
  x.resize(6);
  const bool res2[] = {1, 0, 1, 0, 0, 1};
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 2, 3, 4, res);
}

// d = 3

TEST(TreeapproxBinsearchTest, SimpleBinaryTestD3) {
  const double x2[] = {1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0};
  const bool res2[] = {1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 3, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, SimpleBinaryTest2D3) {
  const double x2[] = {10, 10, 7, 7, 10, 10, 7, 7, 7, 7, 7, 7, 7};
  const bool res2[] = {1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 3, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, EmptyParentTestD3) {
  const double x2[] = {1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0};
  const bool res2[] = {1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 3, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, NotConvexTestD3) {
  const double x2[] = {0.5, 0, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0};
  const bool res2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 3, 3, 3, res);
}

TEST(TreeapproxBinsearchTest, NotConvexTest2D3) {
  const double x2[] = {0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3};
  const bool res2[] = {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 3, 3, 4, res);
}

TEST(TreeapproxBinsearchTest, NotFullTreeD3) {
  // Layout the memory so that the non-existing child has a large value.
  const double x2[] = {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 100};
  vector<double> x(begin(x2), end(x2));
  x.resize(12);
  const bool res2[] = {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1};
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 3, 3, 4, res);
}

// d = 4

TEST(TreeapproxBinsearchTest, SimpleBinaryTestD4) {
  const double x2[] = {1, 1, 0, 0, 0,
                       1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const bool res2[] = {1, 1, 0, 0, 0,
                       1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 4, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, SimpleBinaryTest2D4) {
  const double x2[] = {10, 7, 10, 7, 7,
                       7, 7, 7, 7, 10, 10, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};
  const bool res2[] = {1, 0, 1, 0, 0,
                       0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 4, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, EmptyParentTestD4) {
  const double x2[] = {1, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  const bool res2[] = {1, 0, 0, 1, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 4, 4, 5, res);
}

TEST(TreeapproxBinsearchTest, NotConvexTestD4) {
  const double x2[] = {0.5, 0, 0, 0, 0,
                       2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const bool res2[] = {0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 4, 3, 3, res);
}

TEST(TreeapproxBinsearchTest, NotConvexTest2D4) {
  const double x2[] = {0.5, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0};
  const bool res2[] = {1, 0, 1, 0, 0,
                       0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 4, 3, 4, res);
}

TEST(TreeapproxBinsearchTest, NotFullTreeD4) {
  // Layout the memory so that the non-existing child has a large value.
  const double x2[] = {1, 0, 0, 0, 1,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 100};
  vector<double> x(begin(x2), end(x2));
  x.resize(20);
  const bool res2[] = {1, 0, 0, 0, 1,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  vector<bool> res(begin(res2), end(res2));
  RunAlgo(x, 4, 3, 4, res);
}



// d = 2

TEST(TreeapproxBinsearchWaveletTest, SimpleBinaryTest) {
  const double x2[] = {0, 1, 1, 0, 1, 1, 0, 0};
  const bool res2[] = {1, 1, 1, 0, 1, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 2, 5, 6, res);
}

TEST(TreeapproxBinsearchWaveletTest, SimpleBinaryTest2) {
  const double x2[] = {0, 10, 10, 7, 10, 10, 7, 7};
  const bool res2[] = {1, 1, 1, 0, 1, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 2, 5, 6, res);
}

TEST(TreeapproxBinsearchWaveletTest, EmptyParentTest) {
  const double x2[] = {1, 1, 0, 0, 0, 0, 1, 1};
  const bool res2[] = {1, 1, 0, 1, 0, 0, 1, 1};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 2, 5, 6, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotConvexTest) {
  const double x2[] = {0.5, 1, 0, 0, 2, 3, 0, 0};
  const bool res2[] = {0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 2, 3, 4, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotConvexTest2) {
  const double x2[] = {0.5, 1, 0, 0, 2, 3, 0, 0};
  const bool res2[] = {1, 1, 1, 0, 1, 1, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 2, 4, 5, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotFullTree) {
  // Layout the memory so that the non-existing child has a large value.
  const double x2[] = {1, 1, 0, 1, 0, 0, 1, 100};
  vector<double> x(begin(x2), end(x2));
  x.resize(7);
  const bool res2[] = {1, 1, 0, 1, 0, 0, 1};
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 2, 4, 5, res);
}

// d = 3

TEST(TreeapproxBinsearchWaveletTest, SimpleBinaryTestD3) {
  const double x2[] = {1, 1, 0, 0, 1, 1, 0, 0, 0};
  const bool res2[] = {1, 1, 0, 0, 1, 1, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 3, 4, 5, res);
}

TEST(TreeapproxBinsearchWaveletTest, SimpleBinaryTest2D3) {
  const double x2[] = {10, 10, 7, 10, 10, 7, 7, 7};
  const bool res2[] = {1, 1, 0, 1, 1, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 3, 4, 5, res);
}

TEST(TreeapproxBinsearchWaveletTest, EmptyParentTestD3) {
  const double x2[] = {1, 0, 0, 0, 0, 0, 0, 1, 1};
  const bool res2[] = {1, 0, 1, 0, 0, 0, 0, 1, 1};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 3, 4, 5, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotConvexTestD3) {
  const double x2[] = {0.5, 0, 0, 0, 2, 3, 0, 0, 0};
  const bool res2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 3, 3, 3, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotConvexTest2D3) {
  const double x2[] = {0.5, 0, 0, 0, 2, 3, 0, 0, 0};
  const bool res2[] = {1, 1, 0, 0, 1, 1, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 3, 3, 4, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotFullTreeD3) {
  // Layout the memory so that the non-existing child has a large value.
  const double x2[] = {1, 0, 1, 0, 0, 0, 0, 1, 100};
  vector<double> x(begin(x2), end(x2));
  x.resize(8);
  const bool res2[] = {1, 0, 1, 0, 0, 0, 0, 1};
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 3, 3, 4, res);
}

// d = 4

TEST(TreeapproxBinsearchWaveletTest, SimpleBinaryTestD4) {
  const double x2[] = {1, 1, 0, 0,
                       1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const bool res2[] = {1, 1, 0, 0,
                       1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 4, 4, 5, res);
}

TEST(TreeapproxBinsearchWaveletTest, SimpleBinaryTest2D4) {
  const double x2[] = {10, 7, 10, 7,
                       7, 7, 7, 7, 10, 10, 7, 7, 7, 7, 7, 7};
  const bool res2[] = {1, 0, 1, 0,
                       0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 4, 4, 5, res);
}

TEST(TreeapproxBinsearchWaveletTest, EmptyParentTestD4) {
  const double x2[] = {1, 0, 0, 0,
                       0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  const bool res2[] = {1, 0, 1, 0,
                       0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 4, 4, 5, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotConvexTestD4) {
  const double x2[] = {0.5, 0, 0, 0,
                       2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const bool res2[] = {0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 4, 3, 3, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotConvexTest2D4) {
  const double x2[] = {0.5, 0, 0, 0,
                       0, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0};
  const bool res2[] = {1, 0, 1, 0,
                       0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0};
  vector<double> x(begin(x2), end(x2));
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 4, 3, 4, res);
}

TEST(TreeapproxBinsearchWaveletTest, NotFullTreeD4) {
  // Layout the memory so that the non-existing child has a large value.
  const double x2[] = {1, 0, 0, 1,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 100};
  vector<double> x(begin(x2), end(x2));
  x.resize(15);
  const bool res2[] = {1, 0, 0, 1,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  vector<bool> res(begin(res2), end(res2));
  RunAlgoWavelet(x, 4, 3, 4, res);
}
