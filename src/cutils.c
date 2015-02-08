// Implements fast C versions of some of the functions for dms_tools
// Written by Jesse Bloom.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


static char CheckReadQuality_docs[] = "Fast C implementation of *dms_tools.utils.CheckReadQuality*.";

static PyObject *CheckReadQuality(PyObject *self, PyObject *args) {
    // Calling variables: r1, r2, q1, q2, minq, maxlowqfrac, barcodelength
    char *r1, *r2, *q1, *q2, ri;
    int minq, barcodelength;
    double maxlowqfrac, maxlow;
    size_t r1_len, r2_len, nlow, i;
    PyObject *py_newr1, *py_newr2, *py_returntuple;
    // Parse the arguments
    if (! PyArg_ParseTuple(args, "ssssidi", &r1, &r2, &q1, &q2, &minq, &maxlowqfrac, &barcodelength)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to CheckReadQuality");
        return NULL;
    }
    r1_len = strlen(r1);
    r2_len = strlen(r2);
    if ((r1_len != strlen(q1)) || (r2_len != strlen(q2))) {
        PyErr_SetString(PyExc_ValueError, "r1 and q1 must have the same length, as must r2 and q2");
        return NULL;
    }
    char *newr1 = PyMem_New(char, r1_len + 1);
    char *newr2 = PyMem_New(char, r2_len + 1);
    if ((newr1 == NULL) || (newr2 == NULL)) {
        return PyErr_NoMemory();
    }
    // first check r1
    maxlow = maxlowqfrac * r1_len;
    nlow = 0;
    for (i = 0; i < r1_len; i++) {
        ri = toupper(r1[i]);
        if ((ri == 'N') || ((int) q1[i] - 33 < minq)) {
            nlow++;
            if ((i < barcodelength) || (nlow > maxlow)) {
                PyMem_Del(newr1);
                PyMem_Del(newr2);
                Py_RETURN_FALSE;
            }
            newr1[i] = 'N';
        } else {
            newr1[i] = ri;
        }
    }
    newr1[r1_len] = '\0'; // string termination character
    // first check r2
    maxlow = maxlowqfrac * r2_len;
    nlow = 0;
    for (i = 0; i < r2_len; i++) {
        ri = toupper(r2[i]);
        if ((ri == 'N') || ((int) q2[i] - 33 < minq)) {
            nlow++;
            if ((i < barcodelength) || (nlow > maxlow)) {
                PyMem_Del(newr1);
                PyMem_Del(newr2);
                Py_RETURN_FALSE;
            }
            newr2[i] = 'N';
        } else {
            newr2[i] = ri;
        }
    }
    newr2[r2_len] = '\0'; // string termination character
    // Construct return tuple
    py_newr1 = PyString_FromString(newr1);
    py_newr2 = PyString_FromString(newr2);
    PyMem_Del(newr1);
    PyMem_Del(newr2);
    py_returntuple = PyTuple_New(2);
    if (py_returntuple == NULL) {
        return PyErr_NoMemory();
    }
    PyTuple_SetItem(py_returntuple, 0, py_newr1);
    PyTuple_SetItem(py_returntuple, 1, py_newr2);
    return py_returntuple;
}


static char ReverseComplement_docs[] = "Fast C function for DNA sequence reverse complements.\n\nTakes a single calling argument, which must be a Python string.\nThis string should be composed exclusively of DNA nucleotide codes:\nA, T, C, G, a, t, c, g, N, n.\nThe returned variable is a new string of the same length in which the\nnucleotides are reverse-complemented. Case is preserved\nand N/n reverse complements to N/n.";

static PyObject *ReverseComplement(PyObject *self, PyObject *args) {
    // Calling variables: seq
    char *seq;
    size_t seqlen, i, seqlen1;
    PyObject *pyrc;
    // Parse the arguments.  
    if (! PyArg_ParseTuple(args, "s", &seq)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to ReverseComplement.");
        return NULL;
    }
    seqlen = strlen(seq);
    seqlen1 = seqlen - 1;
    char *rc = PyMem_New(char, seqlen + 1);
    if (rc == NULL) {
        return PyErr_NoMemory();
    }
    for (i = 0; i < seqlen; i++) {
        switch (seq[i]) {
            case 'A' : rc[seqlen1 - i] = 'T';
                       break;
            case 'a' : rc[seqlen1 - i] = 't';
                       break;
            case 'C' : rc[seqlen1 - i] = 'G';
                       break;
            case 'c' : rc[seqlen1 - i] = 'g';
                       break;
            case 'G' : rc[seqlen1 - i] = 'C';
                       break;
            case 'g' : rc[seqlen1 - i] = 'c';
                       break;
            case 'T' : rc[seqlen1 - i] = 'A';
                       break;
            case 't' : rc[seqlen1 - i] = 'a';
                       break;
            case 'N' : rc[seqlen1 - i] = 'N';
                       break;
            case 'n' : rc[seqlen1 - i] = 'n';
                       break;
            default : PyErr_SetString(PyExc_ValueError, "Invalid nucleotide code.");
                      return NULL;
        }
    }
    rc[seqlen] = '\0'; // string termination character
    pyrc = PyString_FromString(rc);
    PyMem_Del(rc);
    return pyrc;
}


static char cutils_docs[] = "Fast implementations of some functions in *dms_tools*.\n\n*ReverseComplement* in this module mimics the same function from *dms_tools.utils*.\n\n*CheckReadQuality* in this module mimics the same function from *dms_tools.utils*.";

static PyMethodDef cutils_funcs[] = {
    {"ReverseComplement", (PyCFunction) ReverseComplement, METH_VARARGS, ReverseComplement_docs},
    {"CheckReadQuality", (PyCFunction) CheckReadQuality, METH_VARARGS, CheckReadQuality_docs},
    {NULL}
};

// Initialize the module
void initcutils(void) {
    Py_InitModule3("cutils", cutils_funcs, cutils_docs);
}
