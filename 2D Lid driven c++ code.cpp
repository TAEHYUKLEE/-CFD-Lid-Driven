#include <stdio.h>
#include <math.h>
#define MAX_SIZE 16382

double u_o[130][130], u_p[130][130], u_t[130][130], u_n[130][130], u_c[130][130], v_o[130][130], v_p[130][130], v_t[130][130], v_n[130][130], v_c[130][130], stream[130][130], vorticity[130][130];

double g[130][130], a[130], b[130], c[130], d[130], pi_o[130][130], pi_n[130][130];

//poisson ������ ��
double f_1[130][130];


int main(void)
{
	double factor, tdiff, error, Re, sum_err, sum_stream = 0.0;
	int i = 0, j = 0, m, n, p, k = 0, T = 0;
	double delx, dely, delt, mu, convective, viscosity, med, rotation;
	double right_o, left_o, right_p, left_p, up_o, down_o, up_p, down_p, vis_le, vis_ri;
	double coef_a, coef_b, coef_c, coef_i;
	double epsil = pow(10, -6);
	const char* s_output = "cavity 129 x 129, Re1000.plt";
	const char* s_result = "cavity result Collocated.csv";
	double sum = 0.0;
	double h;

	FILE* pfo;
	FILE* presult;
	fopen_s(&pfo, s_output, "w");
	fopen_s(&presult, s_result, "w");

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	n = 129; /*numbers of x-panel --> ���ڴ� ���ƾ� �Ѵ�*/
	m = 129; /* numbers of y-panel */
	p = 2500;  /* numbers of t-panel */
	Re = 1000;
	delx = (1.0) / double(n);  // ����ȭ�� ������ ��Ÿ���� (X��������) 
	dely = (1.0) / double(m);  // ����ȭ�� ������ ��Ÿ���� (Y��������)	
	delt = (1.0) / double(p);
	mu = delt / (2 * Re * pow(delx, 2));
	///////////////////// Crank Nicholson �Ҷ��� ����� �̸� ����� ���� /////////////////////////////////////
	coef_a = -mu;
	coef_b = 1.0 + 2.0 * mu;
	coef_c = -mu;
	coef_i = 1.0 + 3.0 * mu;
	//////////////////////////////////////////////////////////////////////////////////////////////

	/*Initial condition x-axis & Boundary (t=0)*/
	for (i = 1; i <= n - 1; i++)
		for (j = 1; j <= m - 1; j++) {
			pi_n[i][j] = 0.0, pi_o[i][j] = 0.0;

		}

	for (i = 1; i <= n - 1; i++)
		for (j = 1; j <= m - 1; j++) {
			u_o[n][j] = 0.0, u_p[n][j] = 0.0, u_t[n][j] = 0.0;
			v_o[i][m] = 0.0, v_p[i][m] = 0.0, v_t[i][m] = 0.0;
			u_o[1][j] = 0.0, u_p[1][j] = 0.0, u_t[1][j] = 0.0;
			v_o[i][1] = 0.0, v_p[i][1] = 0.0, v_t[i][1] = 0.0;
		}

	////////////////////////////////////////////////////////// Navier stokes Solver ////////////////////////////////////////////////////////

	while (1) {

		/*if (k == 500) {
			printf("\nNumber of expression: %d\n\n", T);

			printf("\npi_o: \n");

			for (j = m - 1; j >= 1; j--) {
				for (i = 1; i <= n - 1; i++) {
					//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��

					printf("%f ", pi_o[i][j]);

				}
				printf("\n");
			}


			printf("\n");

			printf("\nu_o:\n");

			for (j = m - 1; j >= 1; j--) {
				for (i = 1; i <= n - 1; i++) {
					//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��

					printf("%f ", u_o[i][j]);

				}
				printf("\n");
			}

			printf("\nv_o:\n");

			for (j = m - 1; j >= 1; j--) {
				for (i = 1; i <= n - 1; i++) {
					//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��

					printf("%f ", v_o[i][j]);

				}
				printf("\n");
			}

			printf("type in your decision 'out':0, 'keep going':1");
			scanf_s("%lf", &h);

			if (h == 0) {

				break;
			}


				T = T + 1;
				k = 0;
		}



			for (j = m - 1; j >= 1; j--) {
				for (i = 1; i <= n - 1; i++) {
					//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��


					sum = sum + fabs(u_o[i][j]);

					if (sum > 1000) {

						printf("fail\n");
						break;
					}

			}
				sum = 0;

			//printf("\n");

		}*/
		//������ �ľ�1. ���� ���� u_n[1][j] �� 0�� �ƴ� ����, v_n[i][1]���� 0�� �ƴ�

		//printf("\n");
		/*Boundary condition setting*/
		for (i = 1; i <= n - 1; i++)
			for (j = 1; j <= m - 1; j++) {
				//printf("u_n[%d][%d] %f\n", i, j, u_n[i][j]);
			}

		//printf("\n");
		for (j = 1; j <= m - 1; j++)
			for (i = 1; i <= n - 1; i++) {

				//printf("v_n[%d][%d] %f\n", i, j, v_n[i][j]);
			}

		//���ſ� ���縦 ���⼭ ������Ʈ�Ѵ�.
			//printf("\n");
		for (j = 1; j <= m - 1; j++) {
			for (i = 1; i <= n - 1; i++) {
				//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��
				u_p[i][j] = u_o[i][j]; //�� step�� ���¸� u_p�� (n term) 
				//printf("u_p[%d][%d] %f\n", i, j, u_p[i][j]);
				u_o[i][j] = u_n[i][j]; //���� �ð��� ���ϱ� ���� �� �� step���� ���� ���� old term���� �Է��Ѵ�. (n+1 term)
				//printf("u_o[%d][%d] %f\n", i, j, u_o[i][j]);
			}
		}

		for (j = 1; j <= m - 1; j++) {
			for (i = 1; i <= n - 1; i++) {
				//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��
				v_p[i][j] = v_o[i][j]; // n-1 step
				v_o[i][j] = v_n[i][j]; // n step

			}
		}

		for (j = m - 1; j >= 1; j--) {
			for (i = 1; i <= n - 1; i++) {
				//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��

				//printf("%f ", u_o[i][j]);

			}
			//printf("\n");
		}

		//printf("\nv_o:\n");

		for (j = m - 1; j >= 1; j--) {
			for (i = 1; i <= n - 1; i++) {
				//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��

			//printf("%f ", v_o[i][j]);

			}
			//printf("\n");
		}


		//printf("\n");

		// Boundar calculation point

		for (i = 2; i <= n - 1; i++) {

			u_o[i][m] = 2 - u_o[i][m - 1], u_p[i][m] = 2 - u_p[i][m - 1];
			u_o[i][0] = -u_o[i][1], u_p[i][0] = -u_p[i][1];//3�� 10�� u_p[i][1] = -u_p[i][1] <--�̷��� �Ǿ� �־���
		}
		for (j = 2; j <= m - 1; j++) {
			v_o[n][j] = -v_o[n - 1][j], v_p[n][j] = -v_p[n - 1][j];
			v_o[0][j] = -v_o[1][j], v_p[0][j] = -v_p[1][j];

		}




		/////////////////////////////////////////////////// First TDMA for (x-direction) /////////////////////////////////////////////////////
		for (j = 1; j <= m - 1; j++) {

			for (i = 2; i <= n - 1; i++) { //x���� �� �� ������ y�� �ö󰡸鼭 ����ؼ� udate

				a[i] = coef_a;
				b[i] = coef_b;
				c[i] = coef_c;

				right_o = (pow(u_o[i + 1][j], 2) + pow(u_o[i][j], 2)) / (2 * delx);
				left_o = (pow(u_o[i][j], 2) + pow(u_o[i - 1][j], 2)) / (2 * delx);
				right_p = (pow(u_p[i + 1][j], 2) + pow(u_p[i][j], 2)) / (2 * delx);
				left_p = (pow(u_p[i][j], 2) + pow(u_p[i - 1][j], 2)) / (2 * delx);
				/*right_o = (u_o[i+1][j]+u_o[i][j])* (u_o[i + 1][j] + u_o[i][j]) / (4 * delx);
				left_o = (u_o[i][j] + u_o[i - 1][j]) * (u_o[i][j] + u_o[i - 1][j]) / (4 * delx);
				right_p = (u_p[i + 1][j] + u_p[i][j]) * (u_p[i + 1][j] + u_p[i][j]) / (4 * delx);
				left_p = (u_p[i][j] + u_p[i - 1][j]) * (u_p[i][j] + u_p[i - 1][j]) / (4 * delx);*/
				up_o = (u_o[i][j + 1] + u_o[i][j]) * (v_o[i - 1][j + 1] + v_o[i][j + 1]) / (4 * dely);
				down_o = (u_o[i][j] + u_o[i][j - 1]) * (v_o[i - 1][j] + v_o[i][j]) / (4 * dely);
				up_p = (u_p[i][j + 1] + u_p[i][j]) * (v_p[i - 1][j + 1] + v_p[i][j + 1]) / (4 * dely);
				down_p = (u_p[i][j] + u_p[i][j - 1]) * (v_p[i - 1][j] + v_p[i][j]) / (4 * dely);
				vis_ri = delt * (u_o[i + 1][j] - 2 * u_o[i][j] + u_o[i - 1][j]) / (pow(delx, 2) * Re);
				vis_le = delt * (u_o[i][j + 1] - 2 * u_o[i][j] + u_o[i][j - 1]) / (pow(dely, 2) * Re);


				convective = -delt / 2 * (3 * (right_o - left_o + up_o - down_o) - (right_p - left_p + up_p - down_p));
				/*-delt / (8 * delx) * (3 * (2 * (pow(u_o[i + 1][j], 2) - pow(u_o[i - 1][j], 2)) + (u_o[i][j + 1] + u_o[i][j]) * (v_o[i - 1][j + 1] + v_o[i][j + 1]) - (u_o[i][j] + u_o[i][j - 1]) * (v_o[i - 1][j] + v_o[i][j]))
				- (2 * (pow(u_p[i + 1][j], 2) - pow(u_p[i - 1][j], 2)) + (u_p[i][j + 1] + u_p[i][j]) * (v_p[i - 1][j + 1] + v_p[i][j + 1]) - (u_p[i][j] + u_p[i][j - 1]) * (v_p[i - 1][j] + v_p[i][j])));*/

				viscosity = vis_ri + vis_le;
				/*delt / (Re * pow(delx, 2)) * (u_o[i+1][j] + u_o[i-1][j] + u_o[i][j+1] + u_o[i][j-1] - 4 * u_o[i][j]);*/

				d[i] = convective + viscosity;

				//printf("numberof j: %d, a[%d]: %f, b[%d]: %f, c[%d]: %f, d[%d]: %f\n", j, i, a[i], i, b[i], i, c[i], i, d[i]);


			}

			//printf("\n");

			for (i = 3; i <= n - 1; i++) { //���⼭�� i=2�����ϱ� i=3�����ϱ�?

				factor = a[i] / b[i - 1]; //a,b,c (b/a)
				b[i] = b[i] - factor * c[i - 1]; //b term
				d[i] = d[i] - factor * d[i - 1]; // d term

				//printf("factor: %f, b[%d]: %f, d[%d]: %f\n", factor, i, b[i], i, d[i]);
			}

			g[n - 1][j] = d[n - 1] / b[n - 1];

			//printf("g[n-1][%d]:%f\n", j, g[n-1][j]);

			for (i = n - 2; i >= 2; i--) {
				g[i][j] = (d[i] - c[i] * g[i + 1][j]) / b[i];
				//printf("g[%d][%d]:%f\n", i, j, g[i][j]);
			}
		}
		//printf("\n");
		/////////////////////////////////////////////////// Second TDMA /////////////////////////////////////////////////////
		for (i = 2; i <= n - 1; i++) {

			for (j = 1; j <= m - 1; j++) { //x���� �� �� ������ y�� �ö󰡸鼭 ����ؼ� udate

				a[j] = coef_a;
				b[j] = coef_b;
				c[j] = coef_c;

				b[1] = coef_i;
				b[m - 1] = coef_i;

				d[j] = g[i][j];

				//printf("numberof i: %d, a[%d]: %f, b[%d]: %f, c[%d]: %f, d[%d]: %f\n", i, j, a[j], j, b[j], j, c[j], i, d[j]);

			}

			//printf("\n");
			for (j = 2; j <= m - 1; j++) {
				factor = a[j] / b[j - 1]; //a,b,c (b/a)
				b[j] = b[j] - factor * c[j - 1]; //b term
				d[j] = d[j] - factor * d[j - 1]; // d term
				//printf("factor: %f, b[%d]: %f, d[%d]: %f\n", factor, i, b[i], i, d[i]);
			}

			u_n[i][m - 1] = d[m - 1] / b[m - 1];

			for (j = m - 2; j >= 1; j--) {
				u_n[i][j] = (d[j] - c[j] * u_n[i][j + 1]) / b[j];
			}

		}

		for (i = 2; i <= n - 1; i++) {
			for (j = 1; j <= m - 1; j++) {

				//u tilda�� ���ϴ� �κ�
				u_t[i][j] = u_n[i][j] + u_o[i][j];

				//printf("check the u_n[%d][%d]: %f\n", i, j, u_n[i][j]);
			}
		}

		//printf("\n");


		/////////////////////////////////////////////////// First TDMA for (y-direction) /////////////////////////////////////////////////////
		for (j = 2; j <= m - 1; j++) {

			for (i = 1; i <= n - 1; i++) { //x���� �� �� ������ y�� �ö󰡸鼭 ����ؼ� udate

				a[i] = coef_a;
				b[i] = coef_b;
				c[i] = coef_c;

				right_o = (u_o[i + 1][j] + u_o[i + 1][j - 1]) * (v_o[i + 1][j] + v_o[i][j]) / (4 * delx);
				left_o = (u_o[i][j] + u_o[i][j - 1]) * (v_o[i - 1][j] + v_o[i][j]) / (4 * delx);
				right_p = (u_p[i + 1][j] + u_p[i + 1][j - 1]) * (v_p[i + 1][j] + v_p[i][j]) / (4 * delx);
				left_p = (u_p[i][j] + u_p[i][j - 1]) * (v_p[i - 1][j] + v_p[i][j]) / (4 * delx); //3�� 11�� ����� left_o�� �Ǿ� �־��� 9�� 46 �����Բ��� / (4 * delx)�̰� ���������� �߰���
				/*up_o = (v_o[i][j+1] + v_o[i][j])* (v_o[i][j + 1] + v_o[i][j]) / (4 * dely);
				down_o = (v_o[i][j] + v_o[i][j-1]) * (v_o[i][j] + v_o[i][j-1]) / (4 * dely);
				up_p = (v_p[i][j + 1] + v_p[i][j]) * (v_p[i][j + 1] + v_p[i][j]) / (4 * dely);
				down_p = (v_p[i][j] + v_p[i][j - 1]) * (v_p[i][j] + v_p[i][j - 1]) / (4 * dely);*/

				up_o = (pow(v_o[i][j + 1], 2) + pow(v_o[i][j], 2)) / (2 * dely);
				down_o = (pow(v_o[i][j], 2) + pow(v_o[i][j - 1], 2)) / (2 * dely);
				up_p = (pow(v_p[i][j + 1], 2) + pow(v_p[i][j], 2)) / (2 * dely);
				down_p = (pow(v_p[i][j], 2) + pow(v_p[i][j - 1], 2)) / (2 * dely);
				vis_ri = delt * (v_o[i + 1][j] - 2 * v_o[i][j] + v_o[i - 1][j]) / (pow(delx, 2) * Re);
				vis_le = delt * (v_o[i][j + 1] - 2 * v_o[i][j] + v_o[i][j - 1]) / (pow(dely, 2) * Re);



				convective = -delt / 2 * (3 * (right_o - left_o + up_o - down_o) - (right_p - left_p + up_p - down_p));

				/*-delt / (8 * delx) * (3 * ((u_o[i+1][j]+u_o[i+1][j-1])*(v_o[i+1][j]+v_o[i][j])-(u_o[i][j]+u_o[i][j-1])*(v_o[i-1][j]+v_o[i][j])+2*(pow(v_o[i][j+1],2) - pow(v_o[i][j-1],2)))
				- ((u_p[i + 1][j] + u_p[i + 1][j - 1]) * (v_p[i + 1][j] + v_p[i][j]) - (u_p[i][j] + u_p[i][j - 1]) * (v_p[i - 1][j] + v_p[i][j]) + 2 * (pow(v_p[i][j + 1], 2) - pow(v_p[i][j - 1], 2))));*/

				viscosity = vis_ri + vis_le;

				d[i] = convective + viscosity;

				//printf("a[%d]: %f, b[%d]: %f, c[%d]: %f, d[%d]: %f\n", i, a[i], i, b[i], i, c[i], i, d[i]);

				//printf("numberof j: %d, a[%d]: %f, b[%d]: %f, c[%d]: %f, d[%d]: %f\n", j, i, a[i], i, b[i], i, c[i], i, d[i]);
			}

			b[1] = coef_i;
			b[n - 1] = coef_i;

			for (i = 2; i <= n - 1; i++) {

				factor = a[i] / b[i - 1]; //a,b,c (b/a)
				b[i] = b[i] - factor * c[i - 1]; //b term
				d[i] = d[i] - factor * d[i - 1]; // d term

				//printf("factor: %f, b[%d]: %f, d[%d]: %f\n", factor, i, b[i], i, d[i]);
			}

			g[n - 1][j] = d[n - 1] / b[n - 1];

			for (i = n - 2; i >= 1; i--) {
				g[i][j] = (d[i] - c[i] * g[i + 1][j]) / b[i];
			}
		}

		/////////////////////////////////////////////////// Second TDMA /////////////////////////////////////////////////////
		for (i = 1; i <= n - 1; i++) {

			for (j = 2; j <= m - 1; j++) { //x���� �� �� ������ y�� �ö󰡸鼭 ����ؼ� udate

				a[j] = coef_a;
				b[j] = coef_b;
				c[j] = coef_c;

				d[j] = g[i][j];

			}

			for (j = 3; j <= m - 1; j++) {//j=3�����ϱ�? 2�����ϱ�?
				factor = a[j] / b[j - 1]; //a,b,c (b/a)
				b[j] = b[j] - factor * c[j - 1]; //b term
				d[j] = d[j] - factor * d[j - 1]; // d term
			}

			v_n[i][m - 1] = d[m - 1] / b[m - 1];

			for (j = m - 2; j >= 2; j--) {
				v_n[i][j] = (d[j] - c[j] * v_n[i][j + 1]) / b[j];
			}

		}

		for (j = 2; j <= m - 1; j++) {
			for (i = 1; i <= n - 1; i++) { // v�� ���� update�� j=2���� �ϴ� �Ŵ�.

				//u tilda�� ���ϴ� �κ�
				v_t[i][j] = v_n[i][j] + v_o[i][j];
				//printf("check the v_t[%d][%d]: %f\n", i, j, v_t[i][j-1]);
				//printf("v_t[%d][%d]: %f\n", i, j, v_t[i][j]);
			}
		}




		///////////////////////////////////////////// Poisson Equation /////////////////////////////////////////////////////

		/////////// Source term /////////////////
		for (j = 1; j <= m - 1; j++) { //v�� ���� update�� j=2���� �ϴ°�???
			for (i = 1; i <= n - 1; i++) {

				f_1[i][j] = ((u_t[i + 1][j] - u_t[i][j]) + (v_t[i][j + 1] - v_t[i][j])) / (delt * delx);

			}
		}

		/////////////////  Gauss siedel method  ///////////////////

		while (1) {

			sum_err = 0.0;

			for (j = 1; j <= m - 1; j++) { //v�� ���� update�� j=2���� �ϴ°�
				pi_o[0][j] = pi_o[1][j];
				pi_o[n][j] = pi_o[n - 1][j];
			}
			for (i = 1; i <= n - 1; i++) {
				pi_o[i][0] = pi_o[i][1];
				pi_o[i][m] = pi_o[i][m - 1];
			}




			for (j = 1; j <= m - 1; j++) { //v�� ���� update�� j=2���� �ϴ°�

				for (i = 1; i <= n - 1; i++) {


						pi_n[i][j] = (pi_o[i + 1][j] + pi_o[i - 1][j] + pi_o[i][j + 1] + pi_o[i][j - 1] - pow(delx, 2) * (f_1[i][j])) / 4;
					

					sum_err = sum_err + fabs(pi_n[i][j] - pi_o[i][j]) / pow(double(n) - 1, 2);

					/*// Neumann boundary condition check point////
					for (int k = 1; k <= n - 1; k++) {

						if ( pi_o[k][0] == pi_o[k][1]	&&	pi_o[k][m] == pi_o[k][m - 1]) {

							printf("Ok\n");
						}
					}

					for (int k = 1; k <= n - 1; k++) {

						if (pi_o[0][k] == pi_o[1][k] && pi_o[n][k] == pi_o[n - 1][k]) {

							printf("Ok\n");
						}
					}
					/////////////////////////////////////////////*/

					pi_o[i][j] = pi_n[i][j]; //�ٷ� �ٷ� ������Ʈ

				}

			}

			if (sum_err < epsil)
			{
				break;
			}
		}
		//printf("pi_o: \n");

		for (j = m - 1; j >= 1; j--) {
			for (i = 1; i <= n - 1; i++) {
				//���� �ð��� �� ������ �ð��� �ʿ��� u_o�� ���� �޾ƿ� ��

				//printf("%f ", pi_o[i][j]);

			}
			//printf("\n");
		}


		//printf("\n");



		for (j = 2; j <= m - 1; j++) { //v�� ���� update�� j=2���� �ϴ°� staggered grid �������� 0�κп����� �з� ����� �ȵǴµ� j=1�϶��� �� ���� ���� �ǹǷ�
			for (i = 1; i <= n - 1; i++) {

				v_n[i][j] = v_t[i][j] - delt * (pi_o[i][j] - pi_o[i][j - 1]) / (dely); //Update old velocity into new velocity
			}
		}


		for (j = 1; j <= m - 1; j++) {
			for (i = 2; i <= n - 1; i++) {

				u_n[i][j] = u_t[i][j] - delt * (pi_o[i][j] - pi_o[i - 1][j]) / (delx); //Update old velocity into new velocity

			}
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




				/////////////////////////////// time difference Calcuation //////////////////////////////////////////////

		tdiff = 0;

		for (i = 2; i <= n - 1; i++)
			for (j = 1; j <= m - 1; j++) {

				tdiff = tdiff + fabs(u_n[i][j] - u_o[i][j]);

				//+ fabs(v_n[i][j] - v_p[i][j])
				//printf("exact[%d][%d]: %f u[%d][%d]: %f, time step %f\n", i, j, exact[i][j], i, j, u_o[i][j], timer);
				//printf("u[%d][%d]: %f\n", i,j,u_n[i][j]);
			}

		double tdiffv = 0;

		for (i = 2; i <= n - 1; i++)
			for (j = 1; j <= m - 1; j++) {

				tdiffv = tdiffv + fabs(v_n[i][j] - v_o[i][j]);

				//+ fabs(v_n[i][j] - v_p[i][j])
				//printf("exact[%d][%d]: %f u[%d][%d]: %f, time step %f\n", i, j, exact[i][j], i, j, u_o[i][j], timer);
				//printf("u[%d][%d]: %f\n", i,j,u_n[i][j]);
			}

		//////////////t count///////////////

		tdiff = (tdiff + tdiffv) / (pow(double(n) - 1, 2)); //n=m�̹Ƿ� 

		printf("tdiff: %f\n", tdiff);

		if (tdiff < epsil) {
			break;
		}




		//check point for boundary (Debugging)

		for (j = 1; j <= m - 1; j++) {
			for (i = 1; i <= n - 1; i++) {


				/*printf("u_o[n][%d]: %f, u_p[n][%d]: %f, u_t[n][%d]: %f\n", j, u_o[n][j], j, u_p[n][j], j, u_t[n][j]);
				printf("v_o[%d][m]: %f, v_p[%d][m]: %f, v_t[%d][m]: %f\n", i, v_o[i][m], i, v_p[i][m], i, v_t[i][m]);
				printf("u_o[1][%d]: %f, u_p[1][%d]: %f, u_t[1][%d]: %f\n", j, u_o[1][j], j, u_p[1][j], j, u_t[1][j]);
				printf("v_o[%d][1]: %f, v_p[%d][1]: %f, v_t[%d][1]: %f\n", i, v_o[i][1], i, v_p[i][1], i, v_t[i][1]);

				if (u_o[i][m] == 2 - u_o[i][m - 1] && u_p[i][m] == 2 - u_p[i][m - 1] && u_o[i][0] == -u_o[i][1] && u_p[i][0] == -u_p[i][1]) {

					printf("Ok\n");
				}

				if (v_o[n][j] == -v_o[n - 1][j] && v_p[n][j] == -v_p[n - 1][j] && v_o[0][j] == -v_o[1][j] && v_p[0][j] == -v_p[1][j]) {

					printf("Ok\n");
				}*/



			}

		}



		k = k + 1;



	}


	/////////////////////////////// last update //////////////////////////////////////////////

	for (i = 1; i <= n - 1; i++)
		for (j = 1; j <= m - 1; j++) {
			u_n[n][j] = 0.0;
			v_n[i][m] = 0.0;
			u_n[1][j] = 0.0;
			v_n[i][1] = 0.0;
		}

	for (i = 1; i <= n - 1; i++) {

		u_n[i][m] = 2 - u_n[i][m - 1];
		u_n[i][0] = -u_n[i][1];//3�� 10�� u_p[i][1] = -u_p[i][1] <--�̷��� �Ǿ� �־���
	}
	for (j = 1; j <= m - 1; j++) {
		v_n[n][j] = -v_n[n - 1][j];
		v_n[0][j] = -v_n[1][j];
	}



	//////////////////////////////////////////////////////////////////////////////////////////////

	// Move points based on staggered grid in collocated grid//

	for (i = 1; i <= n; i++) {

		for (j = 1; j <= m; j++) {

			u_c[i - 1][j - 1] = (u_n[i][j] + u_n[i][j - 1]) / 2;
			v_c[i - 1][j - 1] = (v_n[i][j] + v_n[i - 1][j]) / 2;
			//Collocated grid�� (0,0)���� �����ϱ⶧���� i-1, j-1�� ���� ���� --> Collocated grid �������δ� 0,0���� (n-1, m-1)���� �����Ѵ�.

		}
	}

	u_c[0][m - 1] = 1.0;
	u_c[n - 1][m - 1] = 1.0;

	/// stream line calculation

	for (i = 0; i <= n - 1; i++) {

		for (j = 0; j <= m - 1; j++) {

			sum_stream = sum_stream + u_c[i][j] * dely;
			stream[i][j] = sum_stream;

		}
		sum_stream = 0.0;
	}



	/*	for (i = 0; i <= n - 1; i++) {

			med = stream[i][0];

			for (j = 1; j <= m-1; j++) {



				stream[i][j] = stream[i][j] - med;


			}
			//stream[i][m - 1] = stream[i][m - 1] - st


			stream[i][m - 1];
		}*/


		//Vorticity //

	for (i = 1; i <= n; i++) {

		for (j = 1; j <= m; j++) {

			vorticity[i - 1][j - 1] = ((v_n[i][j] - v_n[i - 1][j]) / delx - (u_n[i][j] - u_n[i][j - 1]) / dely);

		}
	}



	/////// Output /////////
	fprintf(pfo, "VARIABLES = \"x\",\"y\",\"u\",\"v\",\"stream\",\"vorticity\"\n");

	fprintf(pfo, "Zone T=\"%d\"\n", (1));
	fprintf(pfo, "I=%d J=%d\n", (n), (m));
	for (j = 0; j <= m - 1; j++) //time ��
	{
		for (i = 0; i <= n - 1; i++) //y�� (x)
		{

			fprintf(pfo, "%d %d %f %f %f %f\n", i, j, u_c[i][j], v_c[i][j], stream[i][j], vorticity[i][j]);


		}
	}


	for (j = 0; j <= m - 1; j++) //time ��
	{
		for (i = 0; i <= n - 1; i++) //y�� (x)
		{

			fprintf(presult, "%d, %d, %f, %d, %d, %f\n", i, j, u_o[i][j], i, j, v_o[i][j]);


		}
	}

	fclose(pfo);

	return 0;
}