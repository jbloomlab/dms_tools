// Implements fast C versions of some of the functions for dms_tools
// Written by Jesse Bloom.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static char AlignRead_docs[] = "Fast C implementation of *dms_tools.utils.AlignRead*.";

static PyObject *AlignRead(PyObject *self, PyObject *args) {
    PyErr_SetString(PyExc_ValueError, "Although this function appears to work correctly, it has a memory leak and does not improve performance. Use the pure Python implementation instead");
    return NULL;
    // Calling variables: refseq, read, refseqstart, maxmuts, counts, chartype
    char *refseq, *read, *chartype, readcodon[4];
    long refseqstart, maxmuts, nmuts, codonstart, readshift, ncodons, i, arrayindex, icodon, oldcount;
    int has_N;
    PyObject *counts;
    PyObject *charlist = PyList_New(0);
    PyObject *icharlist = PyList_New(0);
    // Parse arguments
    if (! PyArg_ParseTuple(args, "ssllO!s", &refseq, &read, &refseqstart, &maxmuts, &PyDict_Type, &counts, &chartype)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to AlignRead");
        return NULL;
    }
    if (strcmp(chartype, "codon") == 0) {
        nmuts = 0;
        codonstart = (long) (ceil((refseqstart - 1) / 3.0));
        if (((refseqstart - 1) % 3) == 1) {
            readshift = 2;
        } else if (((refseqstart - 1) % 3) == 2) {
            readshift = 1;
        } else {
            readshift = 0;
        }
        if (((strlen(read) - readshift) / 3) < (strlen(refseq) / 3 - codonstart)) {
            ncodons = (strlen(read) - readshift) / 3;
        } else {
            ncodons = strlen(refseq) / 3 - codonstart;
        }
        arrayindex = 0;
        for (icodon = 0; icodon < ncodons; icodon++) {
            has_N = 0;
            for (i = 0; i < 3; i++) {
                readcodon[i] = read[readshift + 3 * icodon + i];
                if (readcodon[i] == 'N') {
                    has_N = 1;
                }
            }
            readcodon[3] = '\0';
            if (! has_N) {
                if ((readcodon[0] != refseq[3 * (codonstart + icodon)]) || (readcodon[1] != refseq[3 * (codonstart + icodon) + 1]) || (readcodon[2] != refseq[3 * (codonstart + icodon) + 2])) {
                    nmuts++;
                    if (nmuts > maxmuts) {
                        Py_RETURN_FALSE;
                    }
                }
                PyList_Append(icharlist, PyInt_FromLong(icodon));
                PyList_Append(charlist, PyString_FromString(readcodon));
                arrayindex++;
            }
        }
        for (i = 0; i < arrayindex; i++) {
            PyObject *codon_number = PyInt_FromLong(PyInt_AsLong(PyList_GetItem(icharlist, i)) + codonstart + 1);
            PyObject *idict = PyDict_GetItem(counts, codon_number);
            PyObject *ipycodon = PyList_GetItem(charlist, i);
            if (PyDict_Contains(idict, ipycodon)) {
                oldcount = PyInt_AsLong(PyDict_GetItem(idict, ipycodon));
                PyDict_SetItem(idict, ipycodon, PyInt_FromLong(1 + oldcount));
            } else {
                PyDict_SetItem(idict, ipycodon, PyInt_FromLong(1));
            }
        }
        Py_RETURN_TRUE;
    } else {
        PyErr_SetString(PyExc_ValueError, "Invalid chartype");
        return NULL;
    }
}


static char AlignSubamplicon_docs[] = "Fast C implementation of *dms_tools.utils.AlignSubamplicon*.";

static PyObject *AlignSubamplicon(PyObject *self, PyObject *args) {
    // Calling variables: refseq, r1, r2, refseqstart, refseqend, maxmuts, maxN, chartype, counts
    char *refseq, *r1, *r2, *chartype;
    char icodon_str[33]; // will hold string representation of icodon
    char mutcodon[4]; // will hold string representation of mutcodon
    long len_refseq, len_r1, len_r2;
    double maxmuts, maxN;
    long refseqstart, refseqend, len_subamplicon, i, n_N, startcodon, codonshift, nmuts, icodon, n, j;
    PyObject *counts, *n_new, *py_returntuple;
    // Parse the arguments
    if (! PyArg_ParseTuple(args, "sssllddsO!", &refseq, &r1, &r2, &refseqstart, &refseqend, &maxmuts, &maxN, &chartype, &PyDict_Type, &counts)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to AlignSubamplicon");
        return NULL;
    }
    len_refseq = (long) strlen(refseq);
    len_r1 = (long) strlen(r1);
    len_r2 = (long) strlen(r2);
    if (refseqstart + len_r1 - 1 > len_refseq) {
        PyErr_SetString(PyExc_ValueError, "R1 extends outside refseq");
        return NULL;
    }
    if (refseqend - len_r2 < 0) {
        PyErr_SetString(PyExc_ValueError, "R2 extends outside refseq");
        return NULL;
    }
    // Build subamplicon of two reads
    len_subamplicon = refseqend - refseqstart + 1;
    char *subamplicon = PyMem_New(char, (size_t) len_subamplicon + 1);
    if (subamplicon == NULL) {
        return PyErr_NoMemory();
    }
    subamplicon[len_subamplicon] = '\0'; // string termination character
    for (i = 0; i < len_subamplicon; i++) {
        if ((i < len_r1) && (i < len_subamplicon - len_r2)) {
            subamplicon[i] = r1[i];
        } else if ((i >= len_r1) && (i < len_subamplicon - len_r2)) {
            subamplicon[i] = 'N';
        } else if ((i < len_r1) && (i >= len_subamplicon - len_r2)) {
            if (r1[i] == r2[i - len_subamplicon + len_r2]) {
                subamplicon[i] = r1[i];
            } else {
                subamplicon[i] = 'N';
            }
        } else {
            subamplicon[i] = r2[i - len_subamplicon + len_r2];
        }
    }
    n_N = 0;
    for (i = 0; i < len_subamplicon; i++) {
        if (subamplicon[i] == 'N') {
            n_N++;
        }
    }
    if (n_N > maxN) {
        PyMem_Del(subamplicon);
        Py_RETURN_FALSE;
    }
    // Count mutations
    if (strcmp(chartype, "codon") == 0) {
        if (refseqstart % 3 == 1) {
            startcodon = (refseqstart + 2) / 3;
            codonshift = 0;
        } else if (refseqstart % 3 == 2) {
            startcodon = (refseqstart + 1) / 3 + 1;
            codonshift = 2;
        } else {
            startcodon = refseqstart / 3 + 1;
            codonshift = 1;
        }
        nmuts = 0;
        for (icodon = startcodon; icodon < refseqend / 3 + 1; icodon++) {
            i = 3 * (icodon - startcodon) + codonshift;
            if (((subamplicon[i] != 'N') && (subamplicon[i + 1] != 'N') && (subamplicon[i + 2] != 'N')) && ((subamplicon[i] != refseq[3 * icodon - 3]) || (subamplicon[i + 1] != refseq[3 * icodon - 2]) || (subamplicon[i + 2] != refseq[3 * icodon - 1]))) {
                nmuts++;
                if (nmuts > maxmuts) {
                    PyMem_Del(subamplicon);
                    Py_RETURN_FALSE;
                }
            }
        }
        for (icodon = startcodon; icodon < refseqend / 3 + 1; icodon++) {
            i = 3 * (icodon - startcodon) + codonshift;
            if (((subamplicon[i] != 'N') && (subamplicon[i + 1] != 'N') && (subamplicon[i + 2] != 'N'))) {
                for (j = 0; j < 3; j++) {
                    mutcodon[j] = subamplicon[i + j];
                }
                mutcodon[3] = '\0'; // string termination character
                sprintf(icodon_str, "%ld", icodon); // icodon_str is now string for icodon
                n = PyInt_AsLong(PyDict_GetItemString(PyDict_GetItemString(counts, icodon_str), mutcodon));
                n_new = PyInt_FromLong(n + 1);
                PyDict_SetItemString(PyDict_GetItemString(counts, icodon_str), mutcodon, n_new);
                Py_DECREF(n_new);
            }
        }
        PyMem_Del(subamplicon);
        py_returntuple = PyTuple_New(2);
        if (py_returntuple == NULL) {
            return PyErr_NoMemory();
        }
        PyTuple_SetItem(py_returntuple, 0, PyBool_FromLong((long) 1));
        PyTuple_SetItem(py_returntuple, 1, PyInt_FromLong(nmuts));
        return py_returntuple;
    } else {
        PyMem_Del(subamplicon);
        PyErr_SetString(PyExc_ValueError, "Invalid chartype");
        return NULL;
    }
    PyMem_Del(subamplicon);
    PyErr_SetString(PyExc_RuntimeError, "Function unexpectedly reached end.");
    return NULL;
}


static char BuildReadConsensus_docs[] = "Fast C implementation of *dms_tools.utils.BuildReadConsensus*.";

static PyObject *BuildReadConsensus(PyObject *self, PyObject *args) {
    // Calling variables: reads, minreadidentity, minreadconcurrence, maxreadtrim
    double minreadidentity, minreadconcurrence, maxnonidentical, min_nt_counts;
    long maxreadtrim, n_nonidentical;
    size_t len_r1, len_r2, max_len_r1, min_len_r1, max_len_r2, min_len_r2, i_len_r1, i_len_r2;
    char *r1_1, *r2_1, *r1_i, *r2_i;
    Py_ssize_t nreads, i, iread;
    PyObject *reads, *tup, *py_r1_consensus, *py_r2_consensus, *py_returntuple;
    const int max_read_len = 2000; // allocate memory for reads up to this length
    long counts_r1_A[max_read_len], counts_r1_C[max_read_len], counts_r1_G[max_read_len], counts_r1_T[max_read_len];
    long counts_r2_A[max_read_len], counts_r2_C[max_read_len], counts_r2_G[max_read_len], counts_r2_T[max_read_len];
    // Parse the arguments
    if (! PyArg_ParseTuple(args, "O!ddl", &PyList_Type, &reads, &minreadidentity, &minreadconcurrence, &maxreadtrim)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to BuildReadConsensus");
        return NULL;
    }
    nreads = PyList_GET_SIZE(reads);
    if (nreads < 2) {
        PyErr_SetString(PyExc_ValueError, "nreads is not >= 2");
        return NULL;
    }
    if ((minreadconcurrence <= 0.5) || (minreadconcurrence > 1.0)) {
        PyErr_SetString(PyExc_ValueError, "minreadconcurrence is not > 0.5 and <= 1.0");
        return NULL;
    }
    tup = PyList_GET_ITEM(reads, 0);
    if (! ((PyTuple_Check(tup)) && (2 == PyTuple_GET_SIZE(tup)) && (PyString_Check(PyTuple_GET_ITEM(tup, 0))) && (PyString_Check(PyTuple_GET_ITEM(tup, 1))))) {
        PyErr_SetString(PyExc_ValueError, "first entry in reads is not a tuple of two strings");
        return NULL;
    }
    r1_1 = PyString_AS_STRING(PyTuple_GET_ITEM(tup, 0));
    r2_1 = PyString_AS_STRING(PyTuple_GET_ITEM(tup, 1));
    max_len_r1 = strlen(r1_1);
    min_len_r1 = strlen(r1_1);
    max_len_r2 = strlen(r2_1);
    min_len_r2 = strlen(r2_1);
    // After this loop, we will have checked that all entries in reads are tuples of two strings
    for (i = 1; i < nreads; i++) {
        tup = PyList_GET_ITEM(reads, i);
        if (! ((PyTuple_Check(tup)) && (2 == PyTuple_GET_SIZE(tup)) && (PyString_Check(PyTuple_GET_ITEM(tup, 0))) && (PyString_Check(PyTuple_GET_ITEM(tup, 1))))) {
            PyErr_SetString(PyExc_ValueError, "an entry in reads is not a tuple of two strings");
            return NULL;
        }
        i_len_r1 = strlen(PyString_AS_STRING(PyTuple_GET_ITEM(tup, 0)));
        i_len_r2 = strlen(PyString_AS_STRING(PyTuple_GET_ITEM(tup, 1)));
        if (i_len_r1 > max_len_r1) {
            max_len_r1 = i_len_r1;
        }
        if (i_len_r2 > max_len_r2) {
            max_len_r2 = i_len_r2;
        }
        if (i_len_r1 < min_len_r1) {
            min_len_r1 = i_len_r1;
        }
        if (i_len_r2 < min_len_r2) {
            min_len_r2 = i_len_r2;
        }
    }
    if (((max_len_r1 - min_len_r1) > maxreadtrim) || ((max_len_r2 - min_len_r2) > maxreadtrim)) {
        // reads cannot be trimmed to compatible lengths
        Py_RETURN_FALSE;
    }
    len_r1 = min_len_r1;
    len_r2 = min_len_r2;
    if ((len_r1 > max_read_len) || (len_r2 > max_read_len)) {
        PyErr_SetString(PyExc_ValueError, "read lengths exceed hardcoded value of max_read_len in BuildReadConsensus");
        return NULL;
    }
    maxnonidentical = (1.0 - minreadidentity) * (len_r1 + len_r2);
    for (i = 0; i < len_r1; i++) {
        counts_r1_A[i] = 0;
        counts_r1_C[i] = 0;
        counts_r1_G[i] = 0;
        counts_r1_T[i] = 0;
        switch (r1_1[i]) {
            case 'A' : counts_r1_A[i]++;
                       break;
            case 'C' : counts_r1_C[i]++;
                       break;
            case 'G' : counts_r1_G[i]++;
                       break;
            case 'T' : counts_r1_T[i]++;
                       break;
            case 'N' : break;
            default : PyErr_SetString(PyExc_ValueError, "Invalid nucleotide");
                      return NULL;
        }
    }
    for (i = 0; i < len_r2; i++) {
        counts_r2_A[i] = 0;
        counts_r2_C[i] = 0;
        counts_r2_G[i] = 0;
        counts_r2_T[i] = 0;
        switch (r2_1[i]) {
            case 'A' : counts_r2_A[i]++;
                       break;
            case 'C' : counts_r2_C[i]++;
                       break;
            case 'G' : counts_r2_G[i]++;
                       break;
            case 'T' : counts_r2_T[i]++;
                       break;
            case 'N' : break;
            default : PyErr_SetString(PyExc_ValueError, "Invalid nucleotide");
                      return NULL;
        }
    }
    for (iread = 1; iread < nreads; iread++) {
        n_nonidentical = 0;
        tup = PyList_GET_ITEM(reads, iread);
        r1_i = PyString_AS_STRING(PyTuple_GET_ITEM(tup, 0));
        r2_i = PyString_AS_STRING(PyTuple_GET_ITEM(tup, 1));
        for (i = 0; i < len_r1; i++) {
            if ((r1_1[i] == 'N') || (r1_i[i] == 'N') || (r1_1[i] != r1_i[i])) {
                n_nonidentical++;
                if (n_nonidentical > maxnonidentical) {
                    Py_RETURN_FALSE;
                }
            }
            switch (r1_i[i]) {
                case 'A' : counts_r1_A[i]++;
                           break;
                case 'C' : counts_r1_C[i]++;
                           break;
                case 'G' : counts_r1_G[i]++;
                           break;
                case 'T' : counts_r1_T[i]++;
                           break;
                case 'N' : break;
                default : PyErr_SetString(PyExc_ValueError, "Invalid nucleotide");
                          return NULL;
            }
        }
        for (i = 0; i < len_r2; i++) {
            if ((r2_1[i] == 'N') || (r2_i[i] == 'N') || (r2_1[i] != r2_i[i])) {
                n_nonidentical++;
                if (n_nonidentical > maxnonidentical) {
                    Py_RETURN_FALSE;
                }
            }
            switch (r2_i[i]) {
                case 'A' : counts_r2_A[i]++;
                           break;
                case 'C' : counts_r2_C[i]++;
                           break;
                case 'G' : counts_r2_G[i]++;
                           break;
                case 'T' : counts_r2_T[i]++;
                           break;
                case 'N' : break;
                default : PyErr_SetString(PyExc_ValueError, "Invalid nucleotide");
                          return NULL;
            }
        }
    }
    // build consensus reads
    min_nt_counts = minreadconcurrence * nreads; // nt must be found >= this many times for consensus
    char *r1_consensus = PyMem_New(char, len_r1 + 1);
    char *r2_consensus = PyMem_New(char, len_r2 + 1);
    if ((r1_consensus == NULL) || (r2_consensus == NULL)) {
        return PyErr_NoMemory();
    }
    r1_consensus[len_r1] = '\0'; // string termination character
    r2_consensus[len_r2] = '\0'; // string termination character
    for (i = 0; i < len_r1; i++) {
        if (counts_r1_A[i] >= min_nt_counts) {
            r1_consensus[i] = 'A';
        } else if (counts_r1_C[i] >= min_nt_counts) {
            r1_consensus[i] = 'C';
        } else if (counts_r1_G[i] >= min_nt_counts) {
            r1_consensus[i] = 'G';
        } else if (counts_r1_T[i] >= min_nt_counts) {
            r1_consensus[i] = 'T';
        } else {
            r1_consensus[i] = 'N';
        }
    }
    for (i = 0; i < len_r2; i++) {
        if (counts_r2_A[i] >= min_nt_counts) {
            r2_consensus[i] = 'A';
        } else if (counts_r2_C[i] >= min_nt_counts) {
            r2_consensus[i] = 'C';
        } else if (counts_r2_G[i] >= min_nt_counts) {
            r2_consensus[i] = 'G';
        } else if (counts_r2_T[i] >= min_nt_counts) {
            r2_consensus[i] = 'T';
        } else {
            r2_consensus[i] = 'N';
        }
    }
    // Construct return tuple
    py_r1_consensus = PyString_FromString(r1_consensus);
    py_r2_consensus = PyString_FromString(r2_consensus);
    PyMem_Del(r1_consensus);
    PyMem_Del(r2_consensus);
    py_returntuple = PyTuple_New(2);
    if (py_returntuple == NULL) {
        return PyErr_NoMemory();
    }
    PyTuple_SetItem(py_returntuple, 0, py_r1_consensus);
    PyTuple_SetItem(py_returntuple, 1, py_r2_consensus);
    return py_returntuple;
}



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


static char cutils_docs[] = "Fast implementations of some functions in *dms_tools*.\n\n*ReverseComplement* in this module mimics the same function from *dms_tools.utils*.\n\n*CheckReadQuality* in this module mimics the same function from *dms_tools.utils*.\n\n*BuildReadConsensus* in this module mimics the same function from *dms_tools.utils*.\n\n*AlignSubamplicon* in this module mimics the same function from *dms_tools.utils*.";

static PyMethodDef cutils_funcs[] = {
    {"ReverseComplement", (PyCFunction) ReverseComplement, METH_VARARGS, ReverseComplement_docs},
    {"CheckReadQuality", (PyCFunction) CheckReadQuality, METH_VARARGS, CheckReadQuality_docs},
    {"BuildReadConsensus", (PyCFunction) BuildReadConsensus, METH_VARARGS, BuildReadConsensus_docs},
    {"AlignSubamplicon", (PyCFunction) AlignSubamplicon, METH_VARARGS, AlignSubamplicon_docs},
    {"AlignRead", (PyCFunction) AlignRead, METH_VARARGS, AlignRead_docs},
    {NULL}
};

// Initialize the module
void initcutils(void) {
    Py_InitModule3("cutils", cutils_funcs, cutils_docs);
}
