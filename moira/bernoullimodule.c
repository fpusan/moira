/*
author = 'Fernando Puente-Sánchez'
email = 'fpusan@gmail.com'
version = '1.0.2'
date = '05-Jan-2016'
license = 'BSD-3'
copyright = 'Copyright 2013-2016 Fernando Puente-Sánchez'

BSD3_LICENSE = """

    Copyright (c) 2015, Fernando Puente Sánchez
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    * Neither the name of moira, Poisson binomial filtering and its Poisson
      approximation, nor the names of its contributors may be used to endorse or
      promote products derived from this software without specific prior
      written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

CITATION = """
Puente-Sánchez F, Aguirre J, Parro V (2015),
A novel conceptual approach to read-filtering in high-throughput amplicon sequencing studies.
Nucleic acids research, gkv1113.
"""
*/

#include <Python.h>
#include <math.h>
#include <string.h>


/// Documentation
static char module_docstring[] = "This module provides an interface for calculating the expected errors of a given sequence using a sum of Bernoulli random variables.";
static char calculate_errors_PB_docstring[] = "This function returns the expected errors of a given sequence with a given confidence value using a sum of Bernoulli random variables.";
/// Compile with: gcc -fpic -shared -I /usr/include/python2.7/ -o bernoulli.so bernoullimodule.c

/// C declarations ////////////////////////////////////////////////////////////////////////////////
double prob_j_errors(double p, int j, int n);
double sum_of_binomials(int j, int k, int n, int(qual_length), double error_probs[], double per_position_accum_probs[][qual_length]);
struct tuple{double expected_errors; int Ns;};
struct tuple test(char contig[], int contig_quals[], double alpha);

/// Python-C interface ////////////////////////////////////////////////////////////////////////////
static PyObject *calculate_errors_PB(PyObject* self, PyObject* args)

{
	char *contig;
	PyObject * contig_quals_listObj; //Python int list.
	PyObject * intObj;	       //A python int inside the list.
	double alpha;
	/// Parse python input tuple into C variables.
	if (!PyArg_ParseTuple(args, "sO!d", &contig, &PyList_Type, &contig_quals_listObj, &alpha)) //The 0! parses for a Python Object (contig_quals_listObj)
	{											  //Checked to be of type PyListType.
		return NULL;
	}
	///Check that alpha is between 0 and 1.
	if (alpha <= 0 || alpha >= 1)
	{
		PyErr_SetString(PyExc_ValueError, "Alpha must be between 0 and 1");
        	return NULL;
	}		
	///Check that sequence and quality are of the same length.
	int contig_quals_size = PyList_Size(contig_quals_listObj);
	if (contig_quals_size != strlen(contig))
	{
		PyErr_SetString(PyExc_ValueError, "contig and contig_quals must have the same length");
        	return NULL;
	}
	///Parse the python list into a C array.
	int contig_quals[contig_quals_size];
	int i;
	for (i = 0; i < contig_quals_size; i++)
	{
		intObj = PyList_GetItem(contig_quals_listObj, i);
		contig_quals [i] = (int)PyInt_AsLong(intObj); //Since python integers are actually implemented as C long integers.
		///If this raised an exception, return it.
		if (PyErr_Occurred())
		{
			return NULL;
		}
		///Having qvalues of 0 (happened in artificial datasets) will lead to divisions by zero in the PB error calculation.
		if (contig_quals [i] == 0)
		{
			contig_quals [i] = 1;
		}
	}
	/// Get results.
	struct tuple results = test(contig, contig_quals, alpha);
	/// Build the output tuple.
	PyObject *result = Py_BuildValue("di", results.expected_errors, results.Ns);
	return result;
}


static PyMethodDef module_methods[] = {
    {"calculate_errors_PB", *calculate_errors_PB, METH_VARARGS, calculate_errors_PB_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initbernoulli(void)
{
	(void) Py_InitModule("bernoulli", module_methods);
}
///////////////////////////////////////////////////////////////////////////////////////////////////



/// C functions ///////////////////////////////////////////////////////////////////////////////////
double prob_j_errors(double p, int j, int n) //Where p is the error probability, j is the number of errors and n the number of observations.
{
	double per_position_accum_probs;
	if (j > n)
	{
		return 0.0; //The formula below would also return 0.
	}
	else 
	{
		per_position_accum_probs = pow((1 - p), n);	//For j == 0.
		int i;
		for(i = 1; i <= j; i++)	//For j > 0.
		{
			per_position_accum_probs = ((n - i + 1) / (1.0*i)) * (p / (1 - p)) * per_position_accum_probs;
		}
		return per_position_accum_probs;
	}		
}	



double sum_of_binomials(int j, int k, int n, int(qual_length), double error_probs[], double per_position_accum_probs[][qual_length]) 
	//#Where j is the number of errors and k is the position in the sequence.
{
	double probability = 0;
	int i;
	
	for(i = 0; i <= j; i++)
	{
		probability += prob_j_errors(error_probs[k], i, n) * per_position_accum_probs[j-i][k-1];
		//Where error_probs[k] is the error probability of the k-th position.
		//Where per_position_accum_probs[j-i][k-1] is the probability that all the bases from position 0 to k-1 had a total of j-i errors.
	}

	return probability;
}



double interpolate(int errors1, double prob1, int errors2, double prob2, double alpha)
{
    	double result = errors1 + ((errors2 - errors1) * ((1 - alpha) - prob1) / (prob2 - prob1));
	if (result < 0) //Happens only for very-short high qual sequences in which the probability of having 0 errors is higher than 1 - alpha.
	{
		result = 0;
	}
	return result;
}



struct tuple test(char contig[], int contig_quals[], double alpha)
{

	///Initialize some variables.
	int i;
	int n = 1; //Since we have a Bernoulli distribution.
	int contig_length = strlen(contig);
	///

	///Translate quality scores into error probabilities.
	int Ns = 0;
	double error_probs [contig_length];
	for(i = 0; i < contig_length; i++)
	{
		if(contig[i] == 78 || contig[i] == 110)	//If position is an "N" or an "n".
		{
			Ns++;
		}
		else	//Will leave error prob at 0 if an N was found at that position.
		{
			error_probs[i - Ns] = pow(10, (contig_quals[i] / -10.0)); //Since we want a continuous list of non-N error_probs.
		}
	}
	///

	///Actually get the job done.
	int max_expected_errors = contig_length + 1;
	int expected_errors = 0;
	double probability;
	double accumulated_probs[max_expected_errors];
	int j;
	int k;
	double per_position_accum_probs[max_expected_errors][contig_length - Ns];
	struct tuple result;
	result.Ns = Ns;
	if((contig_length - Ns) > 0)
	{
		while (1)
		{	
			for(k = 0; k < contig_length - Ns; k++)
			{
				if(k == 0)
				{
					per_position_accum_probs[expected_errors][k] = prob_j_errors(error_probs[k], expected_errors, n);
				}
				else
				{	
					per_position_accum_probs[expected_errors][k] = sum_of_binomials(expected_errors, k, n, contig_length - Ns, error_probs, per_position_accum_probs);
				}
			}	

			probability = per_position_accum_probs[expected_errors][contig_length - Ns - 1];
	
			if(expected_errors == 0)
			{
				accumulated_probs[expected_errors] = probability;
			}
			else
			{
				accumulated_probs[expected_errors] = accumulated_probs[expected_errors - 1] + probability;
			}

			if(accumulated_probs[expected_errors] > (1 - alpha))
			{
				break;
			}
			else
			{
				expected_errors ++;
			}
		}

		result.expected_errors = interpolate(expected_errors - 1, accumulated_probs[expected_errors - 1], expected_errors, accumulated_probs[expected_errors], alpha); 
	}

	else
	{
		result.expected_errors = 0; 
	}

	return result;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
