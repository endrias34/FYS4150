
void Gauss seidel(){

  X_new = []
	X_past = []
	b = []
	a = [][]

	while dist(X_new- X_past) > epsilon k++:

		for (i=0; i < row length to a; i++){

			for (j = 0; j < column length to a; j++){
				if (j < i){
					sum_part_1 = a[i][j]*X_new[i]
				}
				if (j > i){
					sum_part_2 = a[i][j]*X_past[i]
				}
			}
			X_new[i] = b[i] - sum_part_1 - sum_part_2


	}
}


