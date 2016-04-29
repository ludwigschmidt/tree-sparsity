#ifndef __APPROXIMATE_TREE_PROJECTION_H__
#define __APPROXIMATE_TREE_PROJECTION_H__

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace tree_sparsity {

struct binsearch_options {
  enum TreeLayout {
    kCompleteTree = 0,
    kWaveletTree,
  };

  int max_num_iterations;
//  double lambda_low;
//  double lambda_high;
  TreeLayout layout;
  bool verbose;
  void (*output_function)(const char*);

  binsearch_options() : max_num_iterations(-1),
                     // lambda_low(-1), lambda_high(-1),
                        layout(kCompleteTree),
                        verbose(false),
                        output_function(NULL) {}
};


bool approximate_tree_projection(const double* x,
                                 int size,
                                 int degree,
                                 int sparsity_low,
                                 int sparsity_high,
                                 const binsearch_options& options,
                                 bool* support,
                                 double* final_lambda_low,
                                 double* final_lambda_high,
                                 int* num_iterations);


namespace internal {

int compute_tree_d(const double* const x,
                   int size,
                   int degree,
                   double lambda,
                   int last_parent,
                   bool* const support,
                   std::vector<double>* const subtree_weights,
                   std::vector<int>* const bfs_queue) {
  std::vector<double>& w = *subtree_weights;
  std::vector<int>& q = *bfs_queue;

  // compute subtree weights
  for (int ii = size - 1; ii > last_parent; --ii) {
    w[ii] = x[ii] - lambda;
  }

  int child_index;
  
  // last parent might not be full
  w[last_parent] = x[last_parent] - lambda;
  child_index = last_parent * degree;
  for (int jj = 1; jj <= degree; ++jj) {
    child_index += 1;
    if (child_index >= size) {
      break;
    }
    if (w[child_index] > 0.0) {
      w[last_parent] += w[child_index];
    }
  }

  // other nodes are full
  if (last_parent > 0) {
    for (int ii = last_parent - 1; ; --ii) {
      w[ii] = x[ii] - lambda;

      child_index = ii * degree;
      for (int jj = 1; jj <= degree; ++jj) {
        child_index += 1;
        if (w[child_index] > 0.0) {
          w[ii] += w[child_index];
        }
      }
      if (ii == 0) {
        break;
      }
    }
  }

  // compute supports with a BFS starting at the root
  std::fill(support, support + size, false);
  int support_size = 0;
  int q_next = 0;
  int q_end = 0;
  if (w[0] > 0.0) {
    q[q_end] = 0;
    q_end += 1;
  }

  int cur;
  while (q_next < q_end) {
    cur = q[q_next];
    q_next += 1;

    support_size += 1;
    support[cur] = true;

    if (cur > last_parent) {
      continue;
    }

    child_index = cur * degree;

    if (cur == last_parent) {
      for (int ii = 1; ii <= degree; ++ii) {
        child_index += 1;  
        if (child_index >= size) {
          break;
        }
        if (w[child_index] > 0.0) {
          q[q_end] = child_index;
          q_end += 1;
        }
      }
    } else {
      for (int ii = 1; ii <= degree; ++ii) {
        child_index += 1;  
        if (w[child_index] > 0.0) {
          q[q_end] = child_index;
          q_end += 1;
        }
      }
    }
  }

  return support_size;
}


int compute_tree_wavelet_d(const double* x,
                           int size,
                           int degree,
                           double lambda,
                           int last_parent,
                           bool* support,
                           std::vector<double>* subtree_weights,
                           std::vector<int>* bfs_queue) {
  std::vector<double>& w = *subtree_weights;
  std::vector<int>& q = *bfs_queue;

  // compute subtree weights
  for (int ii = size - 1; ii > last_parent; --ii) {
    w[ii] = x[ii] - lambda;
  }

  int child_index;
  int num_children;
  // last parent might not be full
  w[last_parent] = x[last_parent] - lambda;
  if (last_parent != 0) {
    child_index = last_parent * degree;
    num_children = degree;
  } else {
    child_index = last_parent * degree + 1;
    num_children = degree - 1;
  }
  for (int jj = 1; jj <= num_children; ++jj) {
    if (child_index >= size) {
      break;
    }
    if (w[child_index] > 0.0) {
      w[last_parent] += w[child_index];
    }
    child_index += 1;
  }

  // other nodes are full
  if (last_parent > 0) {
    for (int ii = last_parent - 1; ii > 0; --ii) {
      w[ii] = x[ii] - lambda;

      child_index = ii * degree;
      for (int jj = 1; jj <= degree; ++jj) {
        if (w[child_index] > 0.0) {
          w[ii] += w[child_index];
        }
        child_index += 1;
      }
    }
  }

  // handle root separately (has only d-1 children)
  w[0] = x[0] - lambda;
  for (int ii = 1; ii < degree; ++ii) {
    if (w[ii] > 0.0) {
      w[0] += w[ii];
    }
  }

  // compute supports with a BFS starting at the root
  std::fill(support, support + size, false);
  int support_size = 0;
  int q_next = 0;
  int q_end = 0;
  if (w[0] > 0.0) {
    q[q_end] = 0;
    q_end += 1;
  }

  int cur;
  while (q_next < q_end) {
    cur = q[q_next];
    q_next += 1;

    support_size += 1;
    support[cur] = true;

    if (cur > last_parent) {
      continue;
    }

    if (cur == 0) {
      num_children = degree - 1;
      child_index = cur * degree + 1;
    } else {
      num_children = degree;
      child_index = cur * degree;
    }

    if (cur == last_parent) {
      for (int ii = 1; ii <= num_children; ++ii) {
        if (w[child_index] > 0.0) {
          q[q_end] = child_index;
          q_end += 1;
        }
        child_index += 1;
        if (child_index >= size) {
          break;
        }
      }
    } else {
      for (int ii = 1; ii <= num_children; ++ii) {
        if (w[child_index] > 0.0) {
          q[q_end] = child_index;
          q_end += 1;
        }
        child_index += 1;  
      }
    }
  }

  return support_size;
}


int compute_tree_2(const double* x,
                   int size,
                   double lambda,
                   int last_parent,
                   bool* support,
                   std::vector<double>* subtree_weights,
                   std::vector<int>* bfs_queue) {
  std::vector<double>& w = *subtree_weights;
  std::vector<int>& q = *bfs_queue;

  // compute subtree weights
  for (int ii = size - 1; ii > last_parent; --ii) {
    w[ii] = x[ii] - lambda;
  }

  int child_index;
  // last parent might not be full
  w[last_parent] = x[last_parent] - lambda;
  child_index = last_parent * 2 + 1;
  if (w[child_index] > 0.0) {
    w[last_parent] += w[child_index];
  }
  child_index += 1;
  if (child_index < size && w[child_index] > 0.0) {
    w[last_parent] += w[child_index];
  }

  if (last_parent > 0) {
      for (int ii = last_parent - 1; ; --ii) {
      w[ii] = x[ii] - lambda;

      child_index = 2 * ii + 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      if (ii == 0) {
        break;
      }
    }
  }

  // compute supports with a BFS starting at the root
  std::fill(support, support + size, false);
  int support_size = 0;
  int q_next = 0;
  int q_end = 0;
  if (w[0] > 0.0) {
    q[q_end] = 0;
    q_end += 1;
  }

  int cur;
  while (q_next < q_end) {
    cur = q[q_next];
    q_next += 1;

    support_size += 1;
    support[cur] = true;

    if (cur > last_parent) {
      continue;
    }

    child_index = 2 * cur + 1;
    if (w[child_index] > 0.0) {
      q[q_end] = child_index;
      q_end += 1;
    }

    child_index += 1;
    if (cur != last_parent || child_index < size) {
      if (w[child_index] > 0.0) {
        q[q_end] = child_index;
        q_end += 1;
      }
    }
  }

  return support_size;
}


int compute_tree_wavelet_2(const double* x,
                           int size,
                           double lambda,
                           int last_parent,
                           bool* support,
                           std::vector<double>* subtree_weights,
                           std::vector<int>* bfs_queue) {
  std::vector<double>& w = *subtree_weights;
  std::vector<int>& q = *bfs_queue;

  // compute subtree weights
  for (int ii = size - 1; ii > last_parent; --ii) {
    w[ii] = x[ii] - lambda;
  }

  int child_index;
  // last parent might not be full
  w[last_parent] = x[last_parent] - lambda;
  if (last_parent == 0) {
    if (w[1] > 0) {
      w[0] += w[1];
    }
  } else {
    child_index = last_parent * 2;
    if (w[child_index] > 0.0) {
      w[last_parent] += w[child_index];
    }
    child_index += 1;
    if (child_index < size && w[child_index] > 0.0) {
      w[last_parent] += w[child_index];
    }

    for (int ii = last_parent - 1; ii > 0; --ii) {
      w[ii] = x[ii] - lambda;

      child_index = 2 * ii;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }
    }

    // handle root separately
    w[0] = x[0] - lambda;
    if (w[1] > 0) {
      w[0] += w[1];
    }
  }

  // compute supports with a BFS starting at the root
  std::fill(support, support + size, false);
  int support_size = 0;
  int q_next = 0;
  int q_end = 0;
  if (w[0] > 0.0) {
    support_size = 1;
    support[0] = true;
    
    if (1 < size && w[1] > 0.0) {
      q[q_end] = 1;
      q_end += 1;
    }
  } else {
    return 0;
  }

  int cur;
  while (q_next < q_end) {
    cur = q[q_next];
    q_next += 1;

    support_size += 1;
    support[cur] = true;

    if (cur > last_parent) {
      continue;
    }

    child_index = 2 * cur;
    if (w[child_index] > 0.0) {
      q[q_end] = child_index;
      q_end += 1;
    }

    child_index += 1;
    if (cur != last_parent || child_index < size) {
      if (w[child_index] > 0.0) {
        q[q_end] = child_index;
        q_end += 1;
      }
    }
  }

  return support_size;
}


int compute_tree_4(const double* x,
                   int size,
                   double lambda,
                   int last_parent,
                   bool* support,
                   std::vector<double>* subtree_weights,
                   std::vector<int>* bfs_queue) {
  std::vector<double>& w = *subtree_weights;
  std::vector<int>& q = *bfs_queue;

  // compute subtree weights
  for (int ii = size - 1; ii > last_parent; --ii) {
    w[ii] = x[ii] - lambda;
  }

  int child_index;

  w[last_parent] = x[last_parent] - lambda;
  child_index = last_parent * 4 + 1;
  if (w[child_index] > 0.0) {
    w[last_parent] += w[child_index];
  }
  child_index += 1;
  if (child_index < size) {
    if (w[child_index] > 0.0) {
      w[last_parent] += w[child_index];
    }

    child_index += 1;
    if (child_index < size) {
      if (w[child_index] > 0.0) {
        w[last_parent] += w[child_index];
      }
      
      child_index += 1;
      if (child_index < size) {
        if (w[child_index] > 0.0) {
          w[last_parent] += w[child_index];
        }
      }
    }
  }

  if (last_parent > 0) {
    for (int ii = last_parent - 1; ; --ii) {
      w[ii] = x[ii] - lambda;

      child_index = 4 * ii + 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      if (ii == 0) {
        break;
      }
    }
  }

  // compute supports with a BFS starting at the root
  std::fill(support, support + size, false);
  int support_size = 0;
  int q_next = 0;
  int q_end = 0;
  if (w[0] > 0.0) {
    q[q_end] = 0;
    q_end += 1;
  }

  int cur;
  while (q_next < q_end) {
    cur = q[q_next];
    q_next += 1;

    support_size += 1;
    support[cur] = true;

    if (cur > last_parent) {
      continue;
    }

    child_index = 4 * cur + 1;
    if (w[child_index] > 0.0) {
      q[q_end] = child_index;
      q_end += 1;
    }
    child_index += 1;

    if (cur == last_parent) {
      if (child_index < size) {
        if (w[child_index] > 0.0) {
          q[q_end] = child_index;
          q_end += 1;
        }

        child_index += 1;
        if (child_index < size) {
          if (w[child_index] > 0.0) {
            q[q_end] = child_index;
            q_end += 1;
          }

          child_index += 1;
          if (child_index < size) {
            if (w[child_index] > 0.0) {
              q[q_end] = child_index;
              q_end += 1;
            }
          }
        }
      }
    } else {
      if (w[child_index] > 0.0) {
        q[q_end] = child_index;
        q_end += 1;
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        q[q_end] = child_index;
        q_end += 1;
      }

      child_index += 1;
      if ( w[child_index] > 0.0) {
        q[q_end] = child_index;
        q_end += 1;
      }
    }
  }

  return support_size;
}


int compute_tree_wavelet_4(const double* x,
                           int size,
                           double lambda,
                           int last_parent,
                           bool* support,
                           std::vector<double>* subtree_weights,
                           std::vector<int>* bfs_queue) {
  std::vector<double>& w = *subtree_weights;
  std::vector<int>& q = *bfs_queue;

  // compute subtree weights
  for (int ii = size - 1; ii > last_parent; --ii) {
    w[ii] = x[ii] - lambda;
  }

  int child_index;

  w[last_parent] = x[last_parent] - lambda;
  if (last_parent == 0) {
    child_index = last_parent * 4 + 1;
    if (w[child_index] > 0.0) {
      w[last_parent] += w[child_index];
    }
    child_index += 1;
    if (child_index < size) {
      if (w[child_index] > 0.0) {
        w[last_parent] += w[child_index];
      }

      child_index += 1;
      if (child_index < size) {
        if (w[child_index] > 0.0) {
          w[last_parent] += w[child_index];
        }
      }
    }
  } else {
    child_index = last_parent * 4;
    if (w[child_index] > 0.0) {
      w[last_parent] += w[child_index];
    }
    child_index += 1;
    if (child_index < size) {
      if (w[child_index] > 0.0) {
        w[last_parent] += w[child_index];
      }

      child_index += 1;
      if (child_index < size) {
        if (w[child_index] > 0.0) {
          w[last_parent] += w[child_index];
        }
        
        child_index += 1;
        if (child_index < size) {
          if (w[child_index] > 0.0) {
            w[last_parent] += w[child_index];
          }
        }
      }
    }
    for (int ii = last_parent - 1; ii > 0; --ii) {
      w[ii] = x[ii] - lambda;

      child_index = 4 * ii;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        w[ii] += w[child_index];
      }

    }

    // handle root separately
    w[0] = x[0] - lambda;
    if (w[1] > 0.0) {
      w[0] += w[1];
    }
    if (w[2] > 0.0) {
      w[0] += w[2];
    }
    if (w[3] > 0.0) {
      w[0] += w[3];
    }
  }

  // compute supports with a BFS starting at the root
  std::fill(support, support + size, false);
  int support_size = 0;
  int q_next = 0;
  int q_end = 0;
  if (w[0] > 0.0) {
    support_size = 1;
    support[0] = true;

    if (1 < size) {
      if (w[1] > 0.0) {
        q[q_end] = 1;
        q_end += 1;
      }

      if (2 < size) {
        if (w[2] > 0.0) {
          q[q_end] = 2;
          q_end += 1;
        }

        if (3 < size) {
          if (w[3] > 0.0) {
            q[q_end] = 3;
            q_end += 1;
          }
        }
      }
    }
  }

  int cur;
  while (q_next < q_end) {
    cur = q[q_next];
    q_next += 1;

    support_size += 1;
    support[cur] = true;

    if (cur > last_parent) {
      continue;
    }

    child_index = 4 * cur;
    if (w[child_index] > 0.0) {
      q[q_end] = child_index;
      q_end += 1;
    }
    child_index += 1;

    if (cur == last_parent) {
      if (child_index < size) {
        if (w[child_index] > 0.0) {
          q[q_end] = child_index;
          q_end += 1;
        }

        child_index += 1;
        if (child_index < size) {
          if (w[child_index] > 0.0) {
            q[q_end] = child_index;
            q_end += 1;
          }

          child_index += 1;
          if (child_index < size) {
            if (w[child_index] > 0.0) {
              q[q_end] = child_index;
              q_end += 1;
            }
          }
        }
      }
    } else {
      if (w[child_index] > 0.0) {
        q[q_end] = child_index;
        q_end += 1;
      }

      child_index += 1;
      if (w[child_index] > 0.0) {
        q[q_end] = child_index;
        q_end += 1;
      }

      child_index += 1;
      if ( w[child_index] > 0.0) {
        q[q_end] = child_index;
        q_end += 1;
      }
    }
  }

  return support_size;
}


int compute_tree(const double* x,
                 int size,
                 int degree,
                 double lambda,
                 binsearch_options::TreeLayout layout,
                 int last_parent,
                 bool* support,
                 std::vector<double>* subtree_weights,
                 std::vector<int>* bfs_queue) {
  if (layout == binsearch_options::kCompleteTree) {
    if (degree == 2) {
      return compute_tree_2(x, size, lambda, last_parent, support,
          subtree_weights, bfs_queue);
    } else if (degree == 4) {
      return compute_tree_4(x, size, lambda, last_parent, support,
          subtree_weights, bfs_queue);
    } else {
      return compute_tree_d(x, size, degree, lambda, last_parent, support,
          subtree_weights, bfs_queue);
    }
  } else {
    if (degree == 2) {
      return compute_tree_wavelet_2(x, size, lambda, last_parent, support,
          subtree_weights, bfs_queue);
    } else if (degree == 4) {
      return compute_tree_wavelet_4(x, size, lambda, last_parent, support,
          subtree_weights, bfs_queue);
    } else {
      return compute_tree_wavelet_d(x, size, degree, lambda, last_parent,
          support, subtree_weights, bfs_queue);
    }
  }
}

}  // namespace internal


bool approximate_tree_projection(const double * const x,
                                 int size,
                                 int degree,
                                 int k_low,
                                 int k_high,
                                 const binsearch_options& options,
                                 bool* support,
                                 double* final_lambda_low,
                                 double* final_lambda_high,
                                 int* num_iterations) {
  const int kOutputBufferSize = 10000;
  char output_buffer[kOutputBufferSize];

  if (size < degree) {
    snprintf(output_buffer, kOutputBufferSize, "Currently, the tree must "
        "have at least d nodes.\n");
    options.output_function(output_buffer);
    return false;
  }

  int max_num_iterations = 0;
  if (options.max_num_iterations < 0) {
    max_num_iterations = 2 * static_cast<int>(std::ceil(log2(size)));
  } else {
    max_num_iterations = options.max_num_iterations;
  }
  // TODO: uncomment
  /*if (options.lambda_low < 0) {
    new_options.lambda_low = 0.0;
  }
  if (options.lambda_high < 0) {
    new_options.lambda_high = 0.0;
    for (size_t ii = 0; ii < size; ++ii) {
      new_options.lambda_high = max(new_options.lambda_high, x[ii]);
    }
  }*/
  if (options.verbose) {
    if (options.output_function == NULL) {
      return false;
    }
  }
  if (k_low > k_high) {
    snprintf(output_buffer, kOutputBufferSize, "Error: k_low > k_high.\n");
    options.output_function(output_buffer);
    return false;
  }
  if (k_low > size) {
    snprintf(output_buffer, kOutputBufferSize, "Error: k_low > n.\n");
    options.output_function(output_buffer);
    return false;
  }
  if (k_high > size) {
    snprintf(output_buffer, kOutputBufferSize, "Error: k_high > n.\n");
    options.output_function(output_buffer);
    return false;
  }

  std::vector<double> subtree_weights(size);
  std::vector<int> bfs_queue(size);

  // TODO: use options
  double lambda_low = 0.0;
  double lambda_high = 0.0;
  for (int ii = 0; ii < size; ++ii) {
    subtree_weights[ii] = x[ii];
  }
  std::nth_element(subtree_weights.begin(),
                   subtree_weights.begin() + size - k_low / 2,
                   subtree_weights.end());
  lambda_high = subtree_weights[size - k_low / 2] / 2.0;

  if (options.verbose) {
    snprintf(output_buffer, kOutputBufferSize, "n: %d  degree: %d  k_low: %d"
        " k_high: %d  l_low: %e  l_high: %e  max_num_iter: %d\n", size,
        degree, k_low, k_high, lambda_low, lambda_high, max_num_iterations);
    options.output_function(output_buffer);
  }

  int last_parent = 0;

  // last parent is the parent of the last leaf
  if (options.layout == binsearch_options::kCompleteTree) {
    last_parent = (size - 2) / degree;
  } else if (options.layout == binsearch_options::kWaveletTree) {
    last_parent = (size - 1) / degree;
  }

  int num_iter = 0;
  double lambda_mid = 0.0;
  int cur_k;

  do {
    num_iter += 1;
    lambda_high = lambda_high * 2.0;
    cur_k = internal::compute_tree(x, size, degree, lambda_high, options.layout,
        last_parent, support, &subtree_weights, &bfs_queue);
    if (options.verbose) {
      snprintf(output_buffer, kOutputBufferSize, "high: l_cur: %e  (l_low: %e, "
          "l_high: %e)  k: %d\n", lambda_high, lambda_low, lambda_high, cur_k);
      options.output_function(output_buffer);
    }
  } while (cur_k > k_high && num_iter < max_num_iterations);
  
  if (cur_k >= k_low) {
    *final_lambda_low = lambda_low;
    *final_lambda_high = lambda_high;
    *num_iterations = num_iter;
    return true;
  }

  while (num_iter < max_num_iterations) {
    num_iter += 1;
    lambda_mid = (lambda_low + lambda_high) / 2.0;

    cur_k = internal::compute_tree(x, size, degree, lambda_mid, options.layout,
        last_parent, support, &subtree_weights, &bfs_queue);

    if (options.verbose) {
      snprintf(output_buffer, kOutputBufferSize, "mid:  l_cur: %e  (l_low: %e, "
          "l_high: %e)  k: %d\n", lambda_mid, lambda_low, lambda_high, cur_k);
      options.output_function(output_buffer);
    }

    if (cur_k <= k_high && cur_k >= k_low) {
      *final_lambda_low = lambda_low;
      *final_lambda_high = lambda_high;
      *num_iterations = num_iter;
      return true;
    }

    if (cur_k > k_high) {
      lambda_low = lambda_mid;
    } else {
      lambda_high = lambda_mid;
    }
  }

  internal::compute_tree(x, size, degree, lambda_high, options.layout,
      last_parent, support, &subtree_weights, &bfs_queue);

  *final_lambda_low = lambda_low;
  *final_lambda_high = lambda_high;
  *num_iterations = num_iter;

  return true;
}


};  // namespace treeapprox

#endif
