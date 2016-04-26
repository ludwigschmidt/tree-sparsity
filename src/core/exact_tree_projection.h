#ifndef __EXACT_TREE_PROJECTION_H__
#define __EXACT_TREE_PROJECTION_H__

#include <cstdio>
#include <queue>
#include <vector>

namespace tree_sparsity {

bool exact_tree_projection(const double* x,
                           int size,
                           int degree,
                           int sparsity,
                           bool* support) {
  if (sparsity > size) {
    return false;
  }

  std::fill(support, support + size, false);

  int last_parent = (size - 2) / degree;

  std::vector<std::vector<double> > table(size);
  std::vector<std::vector<std::vector<size_t> > > num_allocated(size);

  for (int ii = last_parent + 1; ii < size; ++ii) {
    table[ii].resize(2);
    table[ii][0] = 0.0;
    table[ii][1] = x[ii];
  }

  // bottom-up pass: compute best sum for each (node, tree size)

  for (int ii = last_parent; ; --ii) {
    int to_alloc = 1;
    int child_index = ii * degree;
    for (int jj = 1; jj <= degree; ++jj) {
      child_index += 1;
      if (child_index >= size) {
        break;
      }
      to_alloc += (table[child_index].size() - 1);
    }
    to_alloc = std::min(to_alloc, sparsity);

    table[ii].resize(to_alloc + 1);
    std::fill(table[ii].begin(), table[ii].end(), 0.0);
    num_allocated[ii].resize(to_alloc + 1);
    for (int jj = 0; jj < static_cast<int>(num_allocated[ii].size()); ++jj) {
      num_allocated[ii][jj].resize(degree);
    }

    int max_num = 0;
    int prev_maxnum = 0;

    child_index = ii * degree;
    for (int jj = 1; jj <= degree; ++jj) {
      child_index += 1;
      if (child_index >= size) {
        break;
      }

      prev_maxnum = max_num;
      max_num = std::min(sparsity - 1, static_cast<int>(
          max_num + table[child_index].size() - 1));

      for (int cur_num = max_num; ; --cur_num) {
        for (int in_child = std::min(static_cast<int>(
              table[child_index].size() - 1), cur_num); ;
            --in_child) {
          if (table[ii][cur_num] < table[ii][cur_num - in_child]
                                      + table[child_index][in_child]) {
            table[ii][cur_num] = table[ii][cur_num - in_child]
                                    + table[child_index][in_child];
            num_allocated[ii][cur_num][jj - 1] = in_child;
            //printf("new best entry for (%lu, %lu): allocating %lu to %lu "
            //    "value: %lf\n", ii, cur_num, in_child, child_index,
            //    table[ii][cur_num]);
          }

          if (in_child == 0) {
            break;
          }
          if (cur_num - in_child >= prev_maxnum) {
            break;
          }
        }

        if (cur_num <= 1) {
          break;
        }
      }
    }

    for (int jj = table[ii].size() - 1; jj >= 1; --jj) {
      table[ii][jj] = table[ii][jj - 1] + x[ii];
    }
    if (ii == 0) {
      break;
    }
  }

  /*for (size_t ii = 0; ii < x.size(); ++ii) {
    printf("ii = %lu\n", ii);
    for (size_t jj = 0; jj < table[ii].size(); ++jj) {
      printf("  k = %lu: %lf\n", jj, table[ii][jj]);
    }
    printf("\n");
  }*/

  // top-down pass: identify support (BFS)

  std::queue<std::pair<int , int> > q;
  if (table[0][sparsity] > 0.0) {
    q.push(std::make_pair(0, sparsity));
  }

  while (!q.empty()) {
    int cur_node = q.front().first;
    int cur_sparsity = q.front().second;
    q.pop();
    
    //printf("Allocating %lu to node %lu.\n", cur_k, cur_node);

    support[cur_node] = true;
    cur_sparsity -= 1;

    if (cur_node > last_parent) {
      continue;
    }
    
    int child_index = std::min(cur_node * degree + degree, size - 1);
    int start_index = degree - 1 - (cur_node * degree + degree - child_index);
    for (int jj = start_index; ; --jj) {
      int allocated = num_allocated[cur_node][cur_sparsity][jj];
      if (allocated > 0) {
        q.push(std::make_pair(child_index, allocated));
        cur_sparsity -= allocated;
      }

      if (jj == 0) {
        break;
      }
      child_index -= 1;
    }
  }

  return true;
}

};  // namespace tree_sparsity

#endif
