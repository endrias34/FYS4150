


void Jacobi_rotate (mat A, mat R, int k, int l, int n){
	// s = sin(theta), c = cos(theta)
	double s,c;
	// Checking if largest value is not 0 ?
	if (A[k, l] != 0.0){
		double t, tau;
		tau = (A[l, l] - A[k, k])/(2*A[k, l]);

		if ( tau >= 0) t = 1.0/ (tau + sqrt(1 + tau*tau) )
		if ( tau < 0) t = -1.0/ (-tau + sqrt(1 + tau*tau) )

	c = 1/sqrt(1+(t*t));
	s = c*t;
	}
	else {
	c = 1;
	s = 0;
	}

	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = A[k, k];
	a_ll = A[l, l];
	A[k, k] = c*c*a_kk - 2*c*s*A[k, l] + s*s*a_ll;
	A[l, l] = s*s*a_kk + 2*c*s*A[k, l] + c*c*a_ll;
	A[k, l] = 0.0;
	A[l, k] = 0.0;

	for ( int i != 0; i < n; i++){
		if ( i != k && i != 1){

			a_ik = A[i, k];
			a_il = A[i, l];
			A[i, k] = c*a_ik - s*a_il;
			A[k, i] = A[i, k];
			A[i, l] = c*a_il + s;
			A[l, i] = A[i, l];
		}
		// The new eigenvectors (temporary or what, since we might not have converged just yet?)
		r_ik = R[i, k];
		r_il = R[i, l];
		R[i, k] = c*r_ik - s*r_il;
		R[i, l] = c*r_il + s*r_ik;
	}
	return;
}


