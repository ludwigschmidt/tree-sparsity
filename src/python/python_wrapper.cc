#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "../core/approximate_tree_projection.h"
#include "../core/exact_tree_projection.h"

PyObject* exact_tree_projection(PyObject*, PyObject* args);
PyObject* approximate_tree_projection(PyObject*, PyObject* args);

void write_to_python_stdout(const char* string) {
  PySys_WriteStdout("%s", string);
}

PyMethodDef module_methods[] = {
  {"exact_tree_projection", exact_tree_projection, METH_VARARGS, ""},
  {"approximate_tree_projection", approximate_tree_projection, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC inittreesparsity() {
  Py_InitModule("treesparsity", module_methods);
  import_array();
}

PyObject* exact_tree_projection(PyObject*, PyObject* args) {
  PyObject* vec;
  int degree;
  int sparsity;

  if (!PyArg_ParseTuple(args, "O!ii", &PyArray_Type, &vec, &degree,
      &sparsity)) {
    PyErr_SetString(PyExc_ValueError, "Could not parse input arguments.");
  }

  PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(PyArray_FROM_OTF(vec,
      NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
  if (!arr) {
    PyErr_SetString(PyExc_ValueError, "Could not convert first argument to "
        "a contiguous double-precision NumPy array.");
    return NULL;
  }

  int ndim = PyArray_NDIM(arr);
  if (ndim != 1) {
    PyErr_SetString(PyExc_ValueError, "Input array must be one-dimensional.");
    Py_DECREF(arr);
    return NULL;
  }
  npy_intp* dims = PyArray_DIMS(arr);

  const double* data = static_cast<double*>(PyArray_DATA(arr));

  PyArrayObject* out_arr = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNew(1,
      dims, NPY_BOOL));
  bool* out_data = static_cast<bool*>(PyArray_DATA(out_arr));

  if (!tree_sparsity::exact_tree_projection(data, dims[0], degree, sparsity,
      out_data)) {
    PyErr_SetString(PyExc_ValueError, "Error in exact tree projection.");
    Py_DECREF(arr);
    Py_DECREF(out_arr);
    return NULL;
  }
  
  Py_DECREF(arr);
  //Py_DECREF(out_arr);

  return reinterpret_cast<PyObject*>(out_arr);
}


PyObject* approximate_tree_projection(PyObject*, PyObject* args) {
  PyObject* vec;
  int degree;
  int sparsity_low;
  int sparsity_high;
  int verbose;

  if (!PyArg_ParseTuple(args, "O!iiii", &PyArray_Type, &vec, &degree,
      &sparsity_low, &sparsity_high, &verbose)) {
    PyErr_SetString(PyExc_ValueError, "Could not parse input arguments.");
  }

  PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(PyArray_FROM_OTF(vec,
      NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
  if (!arr) {
    PyErr_SetString(PyExc_ValueError, "Could not convert first argument to "
        "a contiguous double-precision NumPy array.");
    return NULL;
  }

  int ndim = PyArray_NDIM(arr);
  if (ndim != 1) {
    PyErr_SetString(PyExc_ValueError, "Input array must be one-dimensional.");
    Py_DECREF(arr);
    return NULL;
  }
  npy_intp* dims = PyArray_DIMS(arr);

  const double* data = static_cast<double*>(PyArray_DATA(arr));

  PyArrayObject* out_arr = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNew(1,
      dims, NPY_BOOL));
  bool* out_data = static_cast<bool*>(PyArray_DATA(out_arr));

  double final_lambda_low;
  double final_lambda_high;
  int num_iterations;
  tree_sparsity::binsearch_options opts;
  opts.verbose = verbose > 0;
  opts.output_function = write_to_python_stdout;
  opts.max_num_iterations = 30;

  if (!tree_sparsity::approximate_tree_projection(data, dims[0], degree,
      sparsity_low, sparsity_high, opts, out_data, &final_lambda_low,
      &final_lambda_high, &num_iterations)) {
    PyErr_SetString(PyExc_ValueError, "Error in exact tree projection.");
    Py_DECREF(arr);
    Py_DECREF(out_arr);
    return NULL;
  }
  
  Py_DECREF(arr);
  //Py_DECREF(out_arr);

  return reinterpret_cast<PyObject*>(out_arr);
}
