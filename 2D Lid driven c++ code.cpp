#include <stdio.h>
#include <math.h>
#define MAX_SIZE 16382

double u_o[130][130], u_p[130][130], u_t[130][130], u_n[130][130], u_c[130][130], v_o[130][130], v_p[130][130], v_t[130][130], v_n[130][130], v_c[130][130], stream[130][130], vorticity[130][130];

double g[130][130], a[130], b[130], c[130], d[130], pi_o[130][130], pi_n[130][130];

//poisson 정리할 것
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
	n = 129; /*numbers of x-panel --> 격자는 같아야 한다*/
	m = 129; /* numbers of y-panel */
	p = 2500;  /* numbers of t-panel */
	Re = 1000;
	delx = (1.0) / double(n);  // 차분화한 간격을 나타낸다 (X방향으로) 
	dely = (1.0) / double(m);  // 차분화한 간격을 나타낸다 (Y방향으로)	
	delt = (1.0) / double(p);
	mu = delt / (2 * Re * pow(delx, 2));
	///////////////////// Crank Nicholson 할때의 계수들 미리 계산해 놓기 /////////////////////////////////////
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
					//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것

					printf("%f ", pi_o[i][j]);

				}
				printf("\n");
			}


			printf("\n");

			printf("\nu_o:\n");

			for (j = m - 1; j >= 1; j--) {
				for (i = 1; i <= n - 1; i++) {
					//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것

					printf("%f ", u_o[i][j]);

				}
				printf("\n");
			}

			printf("\nv_o:\n");

			for (j = m - 1; j >= 1; j--) {
				for (i = 1; i <= n - 1; i++) {
					//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것

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
					//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것


					sum = sum + fabs(u_o[i][j]);

					if (sum > 1000) {

						printf("fail\n");
						break;
					}

			}
				sum = 0;

			//printf("\n");

		}*/
		//문제점 파악1. 문제 생김 u_n[1][j] 가 0이 아님 현재, v_n[i][1]또한 0이 아님

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

		//과거와 현재를 여기서 업데이트한다.
			//printf("\n");
		for (j = 1; j <= m - 1; j++) {
			for (i = 1; i <= n - 1; i++) {
				//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것
				u_p[i][j] = u_o[i][j]; //전 step의 상태를 u_p로 (n term) 
				//printf("u_p[%d][%d] %f\n", i, j, u_p[i][j]);
				u_o[i][j] = u_n[i][j]; //다음 시간을 구하기 위해 그 전 step에서 구한 것을 old term으로 입력한다. (n+1 term)
				//printf("u_o[%d][%d] %f\n", i, j, u_o[i][j]);
			}
		}

		for (j = 1; j <= m - 1; j++) {
			for (i = 1; i <= n - 1; i++) {
				//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것
				v_p[i][j] = v_o[i][j]; // n-1 step
				v_o[i][j] = v_n[i][j]; // n step

			}
		}

		for (j = m - 1; j >= 1; j--) {
			for (i = 1; i <= n - 1; i++) {
				//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것

				//printf("%f ", u_o[i][j]);

			}
			//printf("\n");
		}

		//printf("\nv_o:\n");

		for (j = m - 1; j >= 1; j--) {
			for (i = 1; i <= n - 1; i++) {
				//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것

			//printf("%f ", v_o[i][j]);

			}
			//printf("\n");
		}


		//printf("\n");

		// Boundar calculation point

		for (i = 2; i <= n - 1; i++) {

			u_o[i][m] = 2 - u_o[i][m - 1], u_p[i][m] = 2 - u_p[i][m - 1];
			u_o[i][0] = -u_o[i][1], u_p[i][0] = -u_p[i][1];//3월 10일 u_p[i][1] = -u_p[i][1] <--이렇게 되어 있었음
		}
		for (j = 2; j <= m - 1; j++) {
			v_o[n][j] = -v_o[n - 1][j], v_p[n][j] = -v_p[n - 1][j];
			v_o[0][j] = -v_o[1][j], v_p[0][j] = -v_p[1][j];

		}




		/////////////////////////////////////////////////// First TDMA for (x-direction) /////////////////////////////////////////////////////
		for (j = 1; j <= m - 1; j++) {

			for (i = 2; i <= n - 1; i++) { //x부터 한 번 돌리고 y로 올라가면서 계속해서 udate

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

			for (i = 3; i <= n - 1; i++) { //여기서는 i=2부터일까 i=3부터일까?

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

			for (j = 1; j <= m - 1; j++) { //x부터 한 번 돌리고 y로 올라가면서 계속해서 udate

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

				//u tilda를 구하는 부분
				u_t[i][j] = u_n[i][j] + u_o[i][j];

				//printf("check the u_n[%d][%d]: %f\n", i, j, u_n[i][j]);
			}
		}

		//printf("\n");


		/////////////////////////////////////////////////// First TDMA for (y-direction) /////////////////////////////////////////////////////
		for (j = 2; j <= m - 1; j++) {

			for (i = 1; i <= n - 1; i++) { //x부터 한 번 돌리고 y로 올라가면서 계속해서 udate

				a[i] = coef_a;
				b[i] = coef_b;
				c[i] = coef_c;

				right_o = (u_o[i + 1][j] + u_o[i + 1][j - 1]) * (v_o[i + 1][j] + v_o[i][j]) / (4 * delx);
				left_o = (u_o[i][j] + u_o[i][j - 1]) * (v_o[i - 1][j] + v_o[i][j]) / (4 * delx);
				right_p = (u_p[i + 1][j] + u_p[i + 1][j - 1]) * (v_p[i + 1][j] + v_p[i][j]) / (4 * delx);
				left_p = (u_p[i][j] + u_p[i][j - 1]) * (v_p[i - 1][j] + v_p[i][j]) / (4 * delx); //3월 11일 디버깅 left_o로 되어 있었음 9시 46 준혁님께서 / (4 * delx)이게 빠져있음을 발견함
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

			for (j = 2; j <= m - 1; j++) { //x부터 한 번 돌리고 y로 올라가면서 계속해서 udate

				a[j] = coef_a;
				b[j] = coef_b;
				c[j] = coef_c;

				d[j] = g[i][j];

			}

			for (j = 3; j <= m - 1; j++) {//j=3부터일까? 2부터일까?
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
			for (i = 1; i <= n - 1; i++) { // v에 대한 update는 j=2부터 하는 거다.

				//u tilda를 구하는 부분
				v_t[i][j] = v_n[i][j] + v_o[i][j];
				//printf("check the v_t[%d][%d]: %f\n", i, j, v_t[i][j-1]);
				//printf("v_t[%d][%d]: %f\n", i, j, v_t[i][j]);
			}
		}




		///////////////////////////////////////////// Poisson Equation /////////////////////////////////////////////////////

		/////////// Source term /////////////////
		for (j = 1; j <= m - 1; j++) { //v에 대한 update는 j=2부터 하는거???
			for (i = 1; i <= n - 1; i++) {

				f_1[i][j] = ((u_t[i + 1][j] - u_t[i][j]) + (v_t[i][j + 1] - v_t[i][j])) / (delt * delx);

			}
		}

		/////////////////  Gauss siedel method  ///////////////////

		while (1) {

			sum_err = 0.0;

			for (j = 1; j <= m - 1; j++) { //v에 대한 update는 j=2부터 하는거
				pi_o[0][j] = pi_o[1][j];
				pi_o[n][j] = pi_o[n - 1][j];
			}
			for (i = 1; i <= n - 1; i++) {
				pi_o[i][0] = pi_o[i][1];
				pi_o[i][m] = pi_o[i][m - 1];
			}




			for (j = 1; j <= m - 1; j++) { //v에 대한 update는 j=2부터 하는거

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

					pi_o[i][j] = pi_n[i][j]; //바로 바로 업데이트

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
				//과거 시간과 그 다음의 시간이 필요함 u_o는 새로 받아온 것

				//printf("%f ", pi_o[i][j]);

			}
			//printf("\n");
		}


		//printf("\n");



		for (j = 2; j <= m - 1; j++) { //v에 대한 update는 j=2부터 하는거 staggered grid 기준으로 0부분에서는 압력 계산이 안되는데 j=1일때는 그 값을 쓰게 되므로
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

		tdiff = (tdiff + tdiffv) / (pow(double(n) - 1, 2)); //n=m이므로 

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
		u_n[i][0] = -u_n[i][1];//3월 10일 u_p[i][1] = -u_p[i][1] <--이렇게 되어 있었음
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
			//Collocated grid는 (0,0)부터 시작하기때문에 i-1, j-1로 해준 것임 --> Collocated grid 기준으로는 0,0부터 (n-1, m-1)까지 존재한다.

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
	for (j = 0; j <= m - 1; j++) //time 축
	{
		for (i = 0; i <= n - 1; i++) //y축 (x)
		{

			fprintf(pfo, "%d %d %f %f %f %f\n", i, j, u_c[i][j], v_c[i][j], stream[i][j], vorticity[i][j]);


		}
	}


	for (j = 0; j <= m - 1; j++) //time 축
	{
		for (i = 0; i <= n - 1; i++) //y축 (x)
		{

			fprintf(presult, "%d, %d, %f, %d, %d, %f\n", i, j, u_o[i][j], i, j, v_o[i][j]);


		}
	}

	fclose(pfo);

	return 0;
}