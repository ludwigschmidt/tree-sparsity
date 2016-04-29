# Projections for tree sparsity

This repository contains efficient algorithms for tree sparsity.

## Usage instructions

### Python

#### Installation

Go to the directory `src/python` and run `python setup.py install`.
It is probably a good idea to do this in a virtual environment, e.g., by using virtualenvwrapper.

#### How to use

After a successful installation, you should be able to import a package `treesparsity`.

The exact dynamic program for a tree-sparse projection can then be called via `treesparsity.exact_tree_projection`.
The function expects three parameters:
- A NumPy array containing the node weights.
- The node degree (all nodes except the rightmost leaf are expected to have the same degree).
- The target sparsity.

The function then returns a bool array indicating for each index whether it is in the optimal support.

For instance, `exact_tree_projection` can be used as follows:

```python
x = np.array([1.0, 0, 0, 2, 3, 0, 0])
support = treesparsity.exact_tree_projection(x, 2, 3)
```

The variable `support` then contains

```python
array([ True,  True, False, False,  True, False, False], dtype=bool)
```

Approximate tree-sparse projections are available via `treesparsity.approximate_tree_projection`.
The function expects five parameters:
- A NumPy array containing the node weights.
- The node degree (all nodes except the rightmost leaf are expected to have the same degree).
- A lower bound on the output sparsity (the lower bound might not be achieved in the solutions).
- An upper bound on the output sparsity.
- A boolean value indicating whether verbose output should be printed to stdout.

The function then returns a bool array indicating for each index whether it is in the returned support.
