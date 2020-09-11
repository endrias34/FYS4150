

void offdiagonal(mat A, int *p, int *q, int n){
	double max = 0;
	for (int = 0 ; i < n ; i++){
		for (int j = i + 1; j < n; j++){
			if (fabs(A[i,j]) > max){
				max = A[i,j]; p = i; q = j;
			}
		}
	}
}
