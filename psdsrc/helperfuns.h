#include <stdio.h>
#include <math.h>
#include <string.h>

int eye(float *M, int n);
int genrotmat(char axis, float angle, float* R);
int multmat(int M, int N, int P, float* mat1, float* mat2, float* mat3);
int orthonormalize(float* mat, int M, int N);
int printmat(FILE *fID, int M, int N, float* mat);
int conv(float* x, int lenx, float* h, int lenh, float* y);
int diff(float* x, int lenx, float di, float* y);
float fsumarr(float *x, int lenx);
float getmaxabs(float *x, int lenx);
float trap(float t, float t_start, float t_ramp, float t_plateau);
int center_out_idx(int length, int idx);
int reverseArray(float *arr, int n);
int catArray(float *arr1, int n1, float *arr2, int n2, int npad, float *arr3);

int eye(float *M, int n) {
	
	int i;
	for (i = 0; i < pow(n,2); i++) {
		if ((int)fmod(i, n+1) == 0)
			M[i] = 1.0;
		else
			M[i] = 0.0;
	}
	
	return 1;
}

int genrotmat(char axis, float angle, float* R) {

        switch (axis) {
                case 'x' :
                        R[0] = 1.0;
                        R[1] = 0.0;
                        R[2] = 0.0;
                        R[3] = 0.0;
                        R[4] = cos(angle);
                        R[5] = -sin(angle);
                        R[6] = 0.0;
                        R[7] = sin(angle);
                        R[8] = cos(angle);
                        break;
                case 'y' :
                        R[0] = cos(angle);
                        R[1] = 0.0;
                        R[2] = sin(angle);
                        R[3] = 0.0;
                        R[4] = 1.0;
                        R[5] = 0.0;
                        R[6] = -sin(angle);
                        R[7] = 0.0;
                        R[8] = cos(angle);
                        break;
                case 'z' :
                        R[0] = cos(angle);
                        R[1] = -sin(angle);
                        R[2] = 0.0;
                        R[3] = sin(angle);
                        R[4] = cos(angle);
                        R[5] = 0.0;
                        R[6] = 0.0;
                        R[7] = 0.0;
                        R[8] = 1.0;
                        break;
                default :
                        fprintf(stderr,"Error: Invalid axis\n");
                        break;
        }

	return 1;
}

int multmat(int M, int N, int P, float* mat1, float* mat2, float* mat3) {
   
        int row, col, v;
        float sum;
	float mat4[M*P];
        for (row = 0; row < M; row++) {
                for (col = 0; col < P; col++) {
                        sum = 0;
                        for (v = 0; v < N; v++)
                                sum += mat1[row * N + v] * mat2[v * P + col];
                        mat4[row * P + col] = sum;
                }
        }

	for (v = 0; v < M*P; v++)
		mat3[v] = mat4[v];

	return 1;
}

int orthonormalize(float* mat, int M, int N) {

	int row, col, n;
	float norm;
	float tmpmat[M*N];
	for (n = 0; n < M*N; n++) tmpmat[n] = mat[n];
	for (col = 0; col < N; col++) {
		norm = 0;
		for (row = 0; row < M; row++)
			norm += pow(tmpmat[col*M + row], 2);
		norm = sqrt(norm);
		for (row = 0; row < M; row++)
			mat[col*M + row] = tmpmat[col*M + row] / norm;
	}

	return 1;
}

int printmat(FILE *fID, int M, int N, float* mat) {

        int row, col;
        for (row = 0; row < M; row++) {
                for (col = 0; col < N; col++)
                        fprintf(fID, "%f \t", mat[row * N + col]);
                fprintf(fID, "\n");
        }

	return 1;
}

int conv(float* x, int lenx, float* h, int lenh, float* y) {

	int leny = lenx + lenh;
	int i, j;
	int h_start, x_start, x_end;
	float val;

	float normx = 0.0;
	for (i = 0; i < lenx; i++)
		normx += fabs(x[i]);

	float normy = 0.0;
	for (i = 0; i < leny; i++) {
		x_start = fmax(0, i-lenh+1);
		x_end = fmin(i+1, lenx);
		h_start = fmin(i, lenh-1);

		val = 0;
		for (j = x_start; j < x_end; j++)
			val += h[h_start--] * x[j];

		y[i] = val;
		normy += fabs(val);
	}

	for (i = 0; i < leny; i++)
		y[i] *= normx/normy;

	return 1;
};

int diff(float* x, int lenx, float di, float* y) {

	int i;
	for (i = 0; i < lenx; i++) y[i] = 0;
	for (i = 1; i < lenx; i++)
		y[i] = (x[i] - x[i-1]) / di;

	return 1;
};

float fsumarr(float *x, int lenx) {
	int i;
	float total = 0;
	for (i = 0; i < lenx; i++)
		total += x[i];
	return total;
}

float getmaxabs(float *x, int lenx) {

	int i;
	float recordmax = fabs(x[0]);
	for (i = 1; i < lenx; i++)
		recordmax = fmax(recordmax, fabs(x[i]));		

	return recordmax;
};

float trap(float t, float t_start, float t_ramp, float t_plateau) {

	if (t < t_start)
		return 0.0;
	else if (t < t_start + t_ramp)
		return (t - t_start) / t_ramp;
	else if (t < t_start + t_ramp + t_plateau)
		return 1.0;
	else if (t < t_start + 2*t_ramp + t_plateau)
		return 1.0 - (t - t_start - t_ramp - t_plateau) / t_ramp;
	else
		return 0.0;

}

int center_out_idx(int length, int index) {

	int center = length / 2;
	if (index % 2)
		return center + ((index+1) / 2);
	else
		return center - (index / 2);
};

int reverseArray(float *arr, int n, float *arr2) {
	for (int i = 0; i < n/2 + 1; ++i) {
		arr2[i] = arr[n - i - 1];
		arr2[n - i - 1] = arr[i];
	}

	return 0;
}

int catArray(float *arr1, int n1, float *arr2, int n2, int npad, float *arr3) {

	if (n1 > 0)
		memcpy(arr3, arr1, n1 * sizeof(float));
	if (npad > 0)
		memset((char*)arr3 + n1 * sizeof(float), 0, npad * sizeof(float));
	if (n2 > 0)
		memcpy((char*)arr3 + (n1 + npad) * sizeof(float), arr2, n2 * sizeof(float));

	return 0;
}
