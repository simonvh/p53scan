#include <Python.h>
#include <math.h>

void get_g(double *g, char *seq, int seqlen, int l) {
	// Return uniform background distribution
	int i;

	if (l == 4) {
		for (i = 0; i < l; i++) {
			g[i] = 0.25;
		}
	}
	else {
		for (i = 0; i < l; i++) {
			g[i] = 0.25 * 0.25;
		}
	}
}

void fill_pwm(double pwm[][4], PyObject * pwm_o) {
	
	int pwm_len = PyList_Size(pwm_o);
	int i,j;
	
	for (i = 0; i < pwm_len; i++) {
		
		PyObject *row = PyList_GetItem(pwm_o, i);
		if (!PyList_Check(row)) {
			PyErr_SetString( PyExc_TypeError, "wrong pwm structure");
			return;
		}	
		if (PyList_Size(row) != 4) {
			PyErr_SetString( PyExc_TypeError, "4 nucleotides!");
			return;
		}		
		for (j = 0; j < 4; j++) {
			pwm[i][j] = PyFloat_AsDouble(PyList_GetItem(row, j));
		}
	}
}

double calc_score(char n, int i, double pwm[][4]) {
	double g = 0.25;
	double z = 0.01;
	
	switch(n) {
		case 'A': 
			return log(pwm[i][0] / g + z); 
		case 'C': 
			return log(pwm[i][1] / g + z); 
		case 'G': 
			return log(pwm[i][2] / g + z); 
		case 'T': 
			return log(pwm[i][3] / g + z); 
	}
	return 0;
}

static PyObject * core_maxscore(PyObject *self, PyObject * args)
{
	PyObject *pwm_o;
	int pwm_len;

	if (!PyArg_ParseTuple(args, "O", &pwm_o))
		return NULL;
	
	if (!PyList_Check(pwm_o))
		return NULL;
	
	pwm_len = PyList_Size(pwm_o);
	double pwm[pwm_len][4];
	fill_pwm(pwm, pwm_o);
}



static PyObject * core_pwmscan(PyObject *self, PyObject * args)
{
	PyObject *seq_o;
	PyObject *pwm_o;
	PyObject *cutoff_o;
	char *seq;
	int seq_len;
	int n_report;
	int pwm_len;
	int i, j;
	double zeta = 0.01;

	if (!PyArg_ParseTuple(args, "OOOi", &seq_o, &pwm_o, &cutoff_o, &n_report))
		return NULL;

	// Sequence and length
	if (!PyString_Check(seq_o))
		return NULL;
	seq = PyString_AsString(seq_o);
	seq_len = PyString_Size(seq_o);

	// Retrieve frequency matrix
	if (!PyList_Check(pwm_o))
		return NULL;

	// Weight matrices
	pwm_len = PyList_Size(pwm_o);
	double pwm[pwm_len][4];
	fill_pwm(pwm, pwm_o);

	// Cutoff for every spacer length
	double cutoff;
	cutoff = PyFloat_AsDouble(cutoff_o);
	
	// Scan sequence
	int j_max = seq_len - pwm_len;
	double score_matrix[j_max];
	double rc_score_matrix[j_max];
	double score, rc_score;
	
	if (j_max < 0) { j_max = 0;}
	int m;
	double g = 0.25;
	double z = 0.01;
	for (j = 0; j < j_max; j++) {
		score = 0;
		rc_score = 0;
		for (m = 0; m < pwm_len; m++) {
			switch(seq[j + m]) {
				case 'A': 
					score += log(pwm[m][0] / g + z);
					rc_score += log(pwm[pwm_len - m - 1][3] / g + z); 
					break;
				case 'C': 
					score += log(pwm[m][1] / g + z);
					rc_score += log(pwm[pwm_len - m - 1][2] / g + z); 
					break;
				case 'G': 
					score += log(pwm[m][2] / g + z);
					rc_score += log(pwm[pwm_len - m - 1][1] / g + z); 
					break;
				case 'T': 
					score += log(pwm[m][3] / g + z);
					rc_score += log(pwm[pwm_len - m - 1][0] / g + z); 
					break;
			}
		
		}
		score_matrix[j] = score;
		rc_score_matrix[j] = rc_score;
	}

	// Initialize matrices of n_report highest scores and corresponding positions + strands
	double maxScores[n_report];
	double maxPos[n_report];
	int maxStrand[n_report];
	for (j = 0; j < n_report; j++) {
		maxScores[j] = -100;
		maxPos[j] = -1;
		maxStrand[j] = 1;
	}
	PyObject*  return_list = PyList_New(0);
	
	int p,q;
	for (j = 0; j < j_max; j++) {
		score = score_matrix[j];
		if (n_report > 0) {
			if (score >= cutoff) {
				p = n_report - 1;
				while ((p >= 0) && (score > maxScores[p])) {
					p--;
				}
				if (p < (n_report-1)) {
					for (q = n_report - 1; q > (p + 1); q--) {
						maxScores[q] = maxScores[q - 1];
						maxPos[q] = maxPos[q - 1];
						maxStrand[q] = maxStrand[q - 1];
					}
					maxScores[p + 1] = score;
					maxPos[p + 1] = j;
					maxStrand[p + 1] = 1;
				}
			}
		}
		else {
			if (score >= cutoff) {
				PyObject* row = PyList_New(0);
				PyList_Append(row, PyFloat_FromDouble(score));
				PyList_Append(row, PyInt_FromLong((long) j));
				PyList_Append(row, PyInt_FromLong((long) 1));
				PyList_Append(return_list, row);
			}
		}
	}

	for (j = 0; j < j_max; j++) {
		score = rc_score_matrix[j];
		if (n_report > 0) {
			if (score >= cutoff) {
				p = n_report - 1;
				while ((p >= 0) && (score > maxScores[p])) {
					p--;
				}
				if (p < (n_report-1)) {
					for (q = n_report - 1; q > (p + 1); q--) {
						maxScores[q] = maxScores[q - 1];
						maxPos[q] = maxPos[q - 1];
						maxStrand[q] = maxStrand[q - 1];
					}
					maxScores[p + 1] = score;
					maxPos[p + 1] = j;
					maxStrand[p + 1] = -1;
				}
			}
		}
		else {
			if (score >= cutoff) {
				PyObject* row = PyList_New(0);
				PyList_Append(row, PyFloat_FromDouble(score));
				PyList_Append(row, PyInt_FromLong((long) j));
				PyList_Append(row, PyInt_FromLong((long) 1));
				PyList_Append(return_list, row);
			}
		}
	}

	for (i = 0; i < n_report; i++) {
		if (maxPos[i] > - 1) {
			PyObject* row = PyList_New(0);
			PyList_Append(row, PyFloat_FromDouble(maxScores[i]));
			PyList_Append(row, PyInt_FromLong((long)maxPos[i]));
			PyList_Append(row, PyInt_FromLong((long)maxStrand[i]));
			PyList_Append(return_list, row);
		}
	}
	return return_list;
}


static PyObject * core_scan_triplet(PyObject *self, PyObject * args)
{
	PyObject *seq_o;
	PyObject *pwm_left_o;
	PyObject *pwm_middle_o;
	PyObject *pwm_right_o;
	PyObject *mult_o;
	double cutoff;
	char *seq;
	int seq_len;
	int n_report;
	int pwm_left_len, pwm_middle_len, pwm_right_len;
	int spacer_min, spacer_max;
	int i, j;
	double zeta = 0.01;

	if (!PyArg_ParseTuple(args, "OOOOiidi", &seq_o, &pwm_left_o, &pwm_middle_o, &pwm_right_o, &spacer_min, &spacer_max, &cutoff, &n_report))
		return NULL;

	// Sequence and length
	if (!PyString_Check(seq_o))
		return NULL;
	seq = PyString_AsString(seq_o);
	seq_len = PyString_Size(seq_o);

	// Retrieve frequency matrix
	if (!PyList_Check(pwm_left_o) || !PyList_Check(pwm_middle_o) || !PyList_Check(pwm_right_o))
		return NULL;

	// Weight matrices
	pwm_left_len = PyList_Size(pwm_left_o);
	pwm_middle_len = PyList_Size(pwm_middle_o);
	pwm_right_len = PyList_Size(pwm_right_o);
	double pwm_left[pwm_left_len][4];
	double pwm_middle[pwm_middle_len][4];
	double pwm_right[pwm_right_len][4];
	fill_pwm(pwm_left, pwm_left_o);
	fill_pwm(pwm_middle, pwm_middle_o);
	fill_pwm(pwm_right, pwm_right_o);

	// Scan sequence
	int j_max = seq_len - pwm_right_len - pwm_left_len - pwm_middle_len;
	double score_matrix_left[j_max];
	double score_matrix_middle[j_max];
	double score_matrix_right[j_max];
	double score, rc_score;
	
	if (j_max < 0) { j_max = 0;}
	int m;
	for (j = 0; j < j_max; j++) {
		// left matrix
		score = 0;
		for (m = 0; m < pwm_left_len; m++) {
			score += calc_score(seq[j + m], m, pwm_left);
		}
		score_matrix_left[j] = score;
		
		// middle matrix
		score = 0;
		for (m = 0; m < pwm_middle_len; m++) {
			score += calc_score(seq[j + m], m, pwm_middle);
		}
		score_matrix_middle[j] = score;

		// right matrix
		score = 0;
		for (m = 0; m < pwm_right_len; m++) {
			score += calc_score(seq[j + m], m, pwm_right);
		}
		score_matrix_right[j] = score;
		
	}

	// Initialize matrices of n_report highest scores and corresponding positions + strands
	double maxScores[n_report];
	double maxPos[n_report];
	int maxStrand[n_report];
	int maxSpacer1[n_report];
	int maxSpacer2[n_report];
	for (j = 0; j < n_report; j++) {
		maxScores[j] = -100;
		maxPos[j] = -1;
		maxSpacer1[j] = 0;
		maxSpacer2[j] = 0;
		maxStrand[j] = 1;
	}
	PyObject*  return_list = PyList_New(0);
	
	// Combine scores for different spacer lengths
	int spacer1, spacer2;
	int p,q;
	for (spacer1 = spacer_min; spacer1 <= spacer_max; spacer1++) {
		for (spacer2 = spacer_min; spacer2 <= spacer_max; spacer2++) {
			for (j = 0; j < j_max - spacer1 - spacer2; j++) {
				score = score_matrix_left[j] + score_matrix_middle[j + pwm_left_len + spacer1] + score_matrix_right[j + pwm_left_len + spacer1 + pwm_middle_len + spacer2];
				if (n_report > 0) {
					if (score >= cutoff) {
						p = n_report - 1;
						while ((p >= 0) && (score > maxScores[p])) {
							p--;
						}
						if (p < (n_report-1)) {
							for (q = n_report - 1; q > (p + 1); q--) {
								maxScores[q] = maxScores[q - 1];
								maxPos[q] = maxPos[q - 1];
								maxStrand[q] = maxStrand[q - 1];
								maxSpacer1[q] = maxSpacer1[q - 1];
								maxSpacer2[q] = maxSpacer2[q - 1];
							}
							maxScores[p + 1] = score;
							maxPos[p + 1] = j;
							maxStrand[p + 1] = 1;
							maxSpacer1[p + 1] = spacer1;
							maxSpacer2[p + 1] = spacer2;
						}
					}
				}
				else {
					if (score >= cutoff) {
						PyObject* row = PyList_New(0);
						PyList_Append(row, PyFloat_FromDouble(score));
						PyList_Append(row, PyInt_FromLong((long) j));
						PyList_Append(row, PyInt_FromLong((long) spacer1));
						PyList_Append(row, PyInt_FromLong((long) spacer2));
						PyList_Append(row, PyInt_FromLong((long) 1));
						PyList_Append(return_list, row);
					}
				}
			}
		}
	}	

	for (i = 0; i < n_report; i++) {
		if (maxPos[i] > - 1) {
			PyObject* row = PyList_New(0);
			PyList_Append(row, PyFloat_FromDouble(maxScores[i]));
			PyList_Append(row, PyInt_FromLong((long)maxPos[i]));
			PyList_Append(row, PyInt_FromLong((long)maxSpacer1[i]));
			PyList_Append(row, PyInt_FromLong((long)maxSpacer2[i]));
			PyList_Append(row, PyInt_FromLong((long)maxStrand[i]));
			PyList_Append(return_list, row);
		}
	}
	return return_list;
}


static PyObject * core_scan(PyObject *self, PyObject * args)
{
	PyObject *seq_o;
	PyObject *pwm_left_o;
	PyObject *pwm_right_o;
	PyObject *mult_o;
	PyObject *cutoff_o;
	char *seq;
	int seq_len;
	int n_report;
	int pwm_left_len, pwm_right_len;
	int spacer_min, spacer_max;
	int i, j;
	double zeta = 0.01;

	if (!PyArg_ParseTuple(args, "OOOiiOi", &seq_o, &pwm_left_o, &pwm_right_o, &spacer_min, &spacer_max, &cutoff_o, &n_report))
		return NULL;

	// Sequence and length
	if (!PyString_Check(seq_o))
		return NULL;
	seq = PyString_AsString(seq_o);
	seq_len = PyString_Size(seq_o);

	// Retrieve frequency matrix
	if (!PyList_Check(pwm_left_o) || !PyList_Check(pwm_right_o))
		return NULL;

	// Weight matrices
	pwm_left_len = PyList_Size(pwm_left_o);
	pwm_right_len = PyList_Size(pwm_right_o);
	double pwm_left[pwm_left_len][4];
	double pwm_right[pwm_right_len][4];
	fill_pwm(pwm_left, pwm_left_o);
	fill_pwm(pwm_right, pwm_right_o);

	// Cutoff for every spacer length
	int cutoff_size = PyList_Size(cutoff_o);
	double cutoff[cutoff_size];
	for (i = 0; i < cutoff_size; i++) {
		cutoff[i] = PyFloat_AsDouble(PyList_GetItem(cutoff_o, i));
	}
	
	// Scan sequence
	int j_max = seq_len - pwm_right_len - pwm_left_len;
	double score_matrix_left[j_max];
	double score_matrix_right[j_max];
	double score, rc_score;
	
	if (j_max < 0) { j_max = 0;}
	int m;
	for (j = 0; j < j_max; j++) {
		score = 0;
		for (m = 0; m < pwm_left_len; m++) {
			score += calc_score(seq[j + m], m, pwm_left);
		}
		score_matrix_left[j] = score;
		score = 0;
		for (m = 0; m < pwm_right_len; m++) {
			score += calc_score(seq[j + m + pwm_left_len], m, pwm_right);
		}
		score_matrix_right[j] = score;
	}

	// Initialize matrices of n_report highest scores and corresponding positions + strands
	double maxScores[n_report];
	double maxPos[n_report];
	int maxStrand[n_report];
	int maxSpacer[n_report];
	for (j = 0; j < n_report; j++) {
		maxScores[j] = -100;
		maxPos[j] = -1;
		maxSpacer[j] = 0;
		maxStrand[j] = 1;
	}
	PyObject*  return_list = PyList_New(0);
	
	// Combine scores for different spacer lengths
	int spacer;
	int p,q;
	for (spacer = spacer_min; spacer <= spacer_max; spacer++) {
		for (j = 0; j < j_max - spacer; j++) {
			score = score_matrix_left[j] + score_matrix_right[j + spacer];
			if (n_report > 0) {
				if (cutoff_size == 0 || score >= cutoff[spacer - spacer_min]) {
					p = n_report - 1;
					while ((p >= 0) && (score > maxScores[p])) {
						p--;
					}
					if (p < (n_report-1)) {
						for (q = n_report - 1; q > (p + 1); q--) {
							maxScores[q] = maxScores[q - 1];
							maxPos[q] = maxPos[q - 1];
							maxStrand[q] = maxStrand[q - 1];
							maxSpacer[q] = maxSpacer[q - 1];
						}
						maxScores[p + 1] = score;
						maxPos[p + 1] = j;
						maxStrand[p + 1] = 1;
						maxSpacer[p + 1] = spacer;
					}
				}
			}
			else {
				if (score >= cutoff[spacer - spacer_min]) {
					PyObject* row = PyList_New(0);
					PyList_Append(row, PyFloat_FromDouble(score));
					PyList_Append(row, PyInt_FromLong((long) j));
					PyList_Append(row, PyInt_FromLong((long) spacer));
					PyList_Append(row, PyInt_FromLong((long) 1));
					PyList_Append(return_list, row);
				}
			}
		}
	}	

	for (i = 0; i < n_report; i++) {
		if (maxPos[i] > - 1) {
			PyObject* row = PyList_New(0);
			PyList_Append(row, PyFloat_FromDouble(maxScores[i]));
			PyList_Append(row, PyInt_FromLong((long)maxPos[i]));
			PyList_Append(row, PyInt_FromLong((long)maxSpacer[i]));
			PyList_Append(row, PyInt_FromLong((long)maxStrand[i]));
			PyList_Append(return_list, row);
		}
	}
	return return_list;
}

static PyObject * core_match(PyObject *self, PyObject * args)
/* Scan sequence with weight matrix and return n best matches */ 
{
	PyObject *seq_o;
	PyObject *pwm_o;
	PyObject *mult_o;
	char *seq;
	int seq_len;
	int n_report;
	int pwm_len;
	int i,j;
	double zeta = 0.01;

	if (!PyArg_ParseTuple(args, "OOi", &seq_o, &pwm_o, &n_report))
		return NULL;

	// Retrieve sequence and length
	if (!PyString_Check(seq_o))
		return NULL;
	seq = PyString_AsString(seq_o);
	seq_len = PyString_Size(seq_o);

	// Retrieve frequency matrix
	if (!PyList_Check(pwm_o))
		return NULL;

	pwm_len = PyList_Size(pwm_o);
	double pwm[pwm_len][4];
	fill_pwm(pwm, pwm_o);

	
	double inf_vector[pwm_len];
	double inf;
	double min;
	double max;
	double tmp_min;
	double tmp_max;
	int pos_min;
	int pos_max;
	min = 0;
	max = 0;
	for (i = 0;i < pwm_len; i++) {
		//printf("i: %d\n", i);
		inf = 0;
		tmp_max = 0;
		tmp_min = 1;
		for (j = 0; j < 4; j++) {
			//printf("j: %d\n", j);
			inf += pwm[i][j] * log(4 * pwm[i][j]);
			//printf("inf: %0.5f\n", inf);
			if (pwm[i][j] > tmp_max) {
				//printf("j: %d pwm: %0.2f\n", j, pwm[i][j]);
				tmp_max = pwm[i][j];
			}	
			if (pwm[i][j] < tmp_min) {
				//printf("j: %d pwm: %0.2f\n", j, pwm[i][j]);
				tmp_min = pwm[i][j];
			}
		}
		inf_vector[i] = inf;
		//printf("min ebrij: %0.4f pos: %d\n", inf * pos_min, pos_min);
		min += inf_vector[i] * tmp_min;
		//printf("max ebrij: %0.4f\n", inf * pos_max);
		max += inf * tmp_max;
	}
		
	//printf("min: %0.5f max: %0.5f\n", min, max);
	

	// Initialize matrices of n_report highest scores and corresponding positions + strands
	double maxScores[n_report];
	double maxPos[n_report];
	int maxStrand[n_report];
	for (j = 0; j < n_report; j++) {
		maxScores[j] = -100;
		maxPos[j] = -1;
		maxStrand[j] = 1;
	}

	// Scan sequence
	int j_max = seq_len - pwm_len;
	if (j_max < 0) { j_max = 0;}
	//	printf("a:%f c:%f g:%f t:%f\n", g[0], g[1], g[2], g[3]);	

	for (j = 0; j < j_max; j++) {
		double score = 0;
		double rc_score = 0;
		//printf("score: %f\n", score);
		//printf("f: %f, g:%f\n",pwm[0][0] , g[0]);
		int m;
		for (m = 0; m < pwm_len; m++) {
			//printf("score: %f\n", score);
			switch(seq[j + m]) {
				case 'A': 
					score += pwm[m][0] * inf_vector[m]; 
					rc_score += pwm[pwm_len - m - 1][3] * inf_vector[pwm_len - m - 1]; 
					break;
				case 'C': 
					score += pwm[m][1] * inf_vector[m]; 
					rc_score += pwm[pwm_len - m - 1][2] * inf_vector[pwm_len - m - 1]; 
					break;
				case 'G': 
					score += pwm[m][2] * inf_vector[m]; 
					rc_score += pwm[pwm_len - m - 1][1] * inf_vector[pwm_len - m - 1]; 
					break;
				case 'T': 
					score += pwm[m][3] * inf_vector[m]; 
					rc_score += pwm[pwm_len - m - 1][0] * inf_vector[pwm_len - m - 1]; 
					break;
			}
		}
		//printf("final score: %f\n", score);
		score = (score - min)/(max - min);
		rc_score = (rc_score - min)/(max - min);
		int strand = 1;
		if (rc_score > score) {
			score = rc_score;
			strand = -1;
		}
		
		int p;
		p = n_report - 1;
		while ((p >= 0) && (score > maxScores[p])) {
			p--;
		}
		if (p < (n_report-1)) {
			int q;
			for (q = n_report - 1; q > (p + 1); q--) {
				maxScores[q] = maxScores[q - 1];
				maxPos[q] = maxPos[q - 1];
				maxStrand[q] = maxStrand[q - 1];
			}
			maxScores[p + 1] = score;
			maxPos[p + 1] = j;
			maxStrand[p + 1] = strand;
		}
	}
	//printf("%f\t%f\t%f\n", maxScores[0], maxPos[0], maxStrand[0]);	
	PyObject*  return_list = PyList_New(0);
	for (i = 0; i < n_report; i++) {
		if (maxPos[i] > - 1) {
			PyObject* row = PyList_New(0);
			PyList_Append(row, PyFloat_FromDouble(maxScores[i]));
			PyList_Append(row, PyInt_FromLong((long)maxPos[i]));
			PyList_Append(row, PyInt_FromLong((long)maxStrand[i]));
			PyList_Append(return_list, row);
		}
	}
	return return_list;
}

static PyObject * core_scan_seq(PyObject *self, PyObject * args) {
	/* Scan sequence sequentially */
	PyObject *seq_o;
	PyObject *pwm_o;
	char *seq;
	char *seq_name;
	int seq_len;
	double min_score;
	int pwm_len;
	int i,j;
	double zeta = 0.01;
	int win;

	if (!PyArg_ParseTuple(args, "sOOd", &seq_name, &seq_o, &pwm_o, &min_score))
		return NULL;
	//printf("n_report: %d\n", n_report);

	// Retrieve sequence and length
	if (!PyString_Check(seq_o))
		return NULL;
	seq = PyString_AsString(seq_o);
	seq_len = PyString_Size(seq_o);
	//printf("seq: %s\n", seq);

	// Retrieve frequency matrix
	if (!PyList_Check(pwm_o))
		return NULL;

	pwm_len = PyList_Size(pwm_o);
	double pwm[pwm_len][4];
	fill_pwm(pwm, pwm_o);

	int j_max = seq_len - pwm_len;
	if (j_max < 0) { j_max = 0;}

	double g[4];
	get_g(g, "A", 1, 4);
	
	for (j = 0; j < j_max; j++) {
		double score = 0;
		double rc_score = 0;
		int m;
		for (m = 0; m < pwm_len; m++) {
		
			switch(seq[j + m]) {
				case 'A': 
					score += log(pwm[m][0] / g[0] + zeta); 
					rc_score += log(pwm[pwm_len - m - 1][3] / g[0] + zeta); 
					break;
				case 'C': 
					score += log(pwm[m][1] / g[1] + zeta); 
					rc_score += log(pwm[pwm_len - m - 1][2] / g[1] + zeta); 
					break;
				case 'G': 
					score += log(pwm[m][2] / g[2] + zeta); 
					rc_score += log(pwm[pwm_len - m - 1][1] / g[2] + zeta); 
					break;
				case 'T': 
					score += log(pwm[m][3] / g[3] + zeta); 
					rc_score += log(pwm[pwm_len - m - 1][0] / g[3] + zeta); 
					break;
				//default: printf("No!\n");
			}
			//printf("score: %f\n", score);
				
		}
		char *strand = "+";
		if (rc_score > score) {
			score = rc_score;
			strand = "-";
		}
		
		//printf("%0.2f %0.2f\n", score, min_score);
		if (score > min_score) {
			printf("%s\t%d\t%d\t%f\t%s\n", seq_name, j + 1, j + 21, score, strand);
		}
	//printf("end\n");
	}
	return PyInt_FromLong(1);
}

static PyMethodDef CoreMethods[] = {
	{"scan", core_scan, METH_VARARGS,"Test"},
	{"scan_triplet", core_scan_triplet, METH_VARARGS,"Test"},
	{"pwmscan", core_pwmscan, METH_VARARGS,"Test"},
	{"match", core_match, METH_VARARGS,"Test"},
	{"scan_seq", core_scan_seq, METH_VARARGS,"Test"},
	//{NULL, NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcore(void)
{
	(void) Py_InitModule("core", CoreMethods);
};
