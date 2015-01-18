/*
author = 'Fernando Puente-Sánchez'
email = 'fpusan@gmail.com'
version = '0.5'
date = '11-Aug-2014'
license = 'GPLv3'
copyright = 'Copyright 2013-2014 Fernando Puente-Sánchez'

GNU_LICENSE = """

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

CITATION = """
Puente-Sánchez F, Aguirre J, Parro V (2014),
A read-filtering algorithm that greatly improves OTU accuracy.
Unpublished.
"""
*/

#include <Python.h>
#include <math.h>
#include <string.h>


/// Documentation
static char module_docstring[] = "This module provides an interface for calculating the expected errors of a given sequence using a sum of Bernoulli random variables.";
static char calculate_errors_PB_docstring[] = "This function returns the expected errors of a given sequence with a given confidence value using a sum of Bernoulli random variables.";
/// Compile with: gcc -fpic -shared -I /usr/include/python2.6/ -o bernoulli.so bernoullimodule.c

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
        PyObject * intObj;               //A python int inside the list.
        double alpha;
        /// Parse python input tuple into C variables.
        if (!PyArg_ParseTuple(args, "sO!d", &contig, &PyList_Type, &contig_quals_listObj, &alpha)) //The 0! parses for a Python Object (contig_quals_listObj)
        {                                                                                          //Checked to be of type PyListType.
                return NULL;
        }
        ///Parse the python list into a C array.
        int contig_quals_size = PyList_Size(contig_quals_listObj);
        int contig_quals[contig_quals_size];
        int i;
        if (contig_quals_size < 0)
	{
		return NULL; //Not a list.
	}
	for (i = 0; i < contig_quals_size; i++)
	{
            intObj = PyList_GetItem(contig_quals_listObj, i);
            contig_quals [i] = (int)PyInt_AsLong(intObj); //Since python integers are actually implemented as C long integers. 
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

	struct tuple result;
	result.expected_errors = interpolate(expected_errors - 1, accumulated_probs[expected_errors - 1], expected_errors, accumulated_probs[expected_errors], alpha); 
	result.Ns = Ns;

	return result;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
