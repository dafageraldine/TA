#include "stdafx.h"
#include "fftw3.h"
#include "pbPlots.hpp"
#include "supportLib.hpp"
#include <Windows.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <time.h>
#include <fstream>
#include "svm.h"
#include <ctype.h>
#include <string.h>
#include <errno.h>


extern "C" {
#include "extApi.h"
#include "extApiPlatform.h"
}

using namespace std;
using namespace cv;
#define PI 3.14
#define imgArr 100
#define major_axis 50
#define minor_axis 50

fftw_complex x_[100];
fftw_complex y_[100];
//fftw_plan plan = fftw_plan_dft_1d(32, x_, y_, FFTW_FORWARD, FFTW_MEASURE);
//ofstream MyFile("DATASETS.txt");
//ofstream MyFile("dummy.txt");
//ofstream Mygraph("dum.txt");
ofstream Mygraph("Datagraph.txt");
int clientID = 0;

simxUChar* laserdepan;
simxUChar* laserbelakang;

int roda1 = 0;
int roda2 = 0;
//int robot;

simxInt datasize;
simxInt datasize_b;

Mat imgSave;

char file_name[255];

float tetha = (float)180 / (float)1024;

float xlaser;
float ylaser;

float xlaser_b;
float ylaser_b;



// Class for handling robot's properties



struct laser
{
	float x, y, x_b, y_b;
};

struct lsr{
	float x, y, z, dist_from_urg, dist_from_robot_center, angle;
	cv::Point urg, from_robot_center[700], from_urg, obj[200];
	float dist[700], alpha[700], min_dist1, min_dist2, beta1, beta2;
};

struct lsr lrrf[5];

int robothandle;
struct laser lrf[5];
vector<double> xhasil;
vector<double> yhasil;

int container; int container_b; float alldatavalue[5000]; float alldatavalue_b[5000]; float resadv[5000];
float resadv_b[5000]; float resadv2[2000]; float resadv_b2[2000]; float resadv3[2000]; float resadv_b3[2000];
float resadv4[2000]; float resadv_b4[2000]; float resadv5[2000]; float resadv_b5[2000]; float resadv6[2000];
float resadv_b6[2000]; float resadv7[2000]; float resadv_b7[2000]; float hasilfft[5000][32];
int REAL = 0; int Imaginer = 1; int pros;int start_loop; int finish_loop; int pengurangan;

vector<double> xll;
vector<double> yll;
vector<double> x_1{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, };
vector<double> y_1{
	5.5 ,5.5, 5.5, 1.74492, 1.72362, 1.72238, 1.72818, 1.75987, 5.5, 5.5, 5.5, 5.5, 5.5, 1.73633, 1.72435, 5.5, 1.72847, 1.77033, 5.5, 5.5 };
RGBABitmapImageReference *imageRes = CreateRGBABitmapImageReference();
RGBABitmapImageReference *imageRezz = CreateRGBABitmapImageReference();
RGBABitmapImageReference *imageOrik = CreateRGBABitmapImageReference();
RGBABitmapImageReference *imageori = CreateRGBABitmapImageReference();

double datares[100];
double dataori[100];


Mat rotate(Mat src, double angle)
{
	Mat dst;
	Point2f pt(src.cols / 2., src.rows / 2.);
	Mat r = getRotationMatrix2D(pt, angle, 1.0);
	warpAffine(src, dst, r, Size(src.cols, src.rows));
	return dst;
}

void slidingwindow(int resized,int finish, int mode){
	//cout << "mode : " << mode << "\n";
	int tes = 0;
	start_loop = 0; finish_loop = 20; pengurangan = 0; pros = 1;
	for (int k = start_loop; k < finish_loop; k++){
		if (k < (resized / 2)){
			//tes = tes + k - pengurangan;
			//k - pengurangan > 19 ? cout << "ERROR": cout <<"";
			//mode == 1 ? cout << k - pengurangan << "\n" : cout << "";
			x_[k - pengurangan][REAL] =  mode == 0 ? alldatavalue[k] : mode == 1 ? resadv[k] : mode == 2 ? resadv2[k] : mode == 3 ? resadv3[k] : mode == 4 ? resadv4[k] : mode == 5 ? resadv5[k] : mode == 6 ? resadv6[k]: resadv7[k];
			x_[k - pengurangan][Imaginer] = 0;
		}

		if (k >= (resized / 2)){
			//tes = tes + k - pengurangan;
			//k - pengurangan > 19 ? cout << "ERROR" : cout << "";
			//mode == 1 ? cout << k - pengurangan << "\n" : cout << "";
			x_[k - pengurangan][REAL] = mode == 0 ? alldatavalue_b[k - (resized / 2)] : mode == 1 ? resadv_b[k - (resized / 2)] : mode == 2 ? resadv_b2[k - (resized / 2)] : mode == 3 ? resadv_b3[k - (resized / 2)] : mode == 4 ? resadv_b4[k - (resized / 2)] : mode == 5 ? resadv_b5[k - (resized / 2)] : mode == 6 ? resadv_b6[k - (resized / 2)] : resadv_b7[k - (resized / 2)];
			x_[k - pengurangan][Imaginer] = 0;
		}

		if (k == finish_loop - 1){

			
			for (int l = 0; l < 20; l++){
				if (l == 0){
					Mygraph << "+1 ";
				}
				l == 19 ? Mygraph << l+1 << ":" << x_[l][REAL] << "\n" : Mygraph << l+1 << ":" << x_[l][REAL] << " ";
			}

			for (int a = 20; a < 32; a++){
				x_[a][REAL] = 0;
				x_[a][Imaginer] = 0;
			}
			/*fftw_execute(plan);

			for (int i = 1; i <= 31; i++){
				
				hasilfft[pros-1][i-1] =  sqrt(y_[i][REAL] * y_[i][REAL] + y_[i][Imaginer] * y_[i][Imaginer]);
				if (i == 1){
					MyFile << "+1 ";
				}
				i == 31 ? MyFile << i << ":" << hasilfft[pros - 1][i - 1] << "\n" : MyFile << i << ":" << hasilfft[pros - 1][i - 1] << " ";
				
			}*/
			
			
			pros = pros + 1;
			
			if (pros == finish){
				
				break;

			}
			start_loop = start_loop + 2;
			finish_loop = finish_loop + 2;
			pengurangan = pengurangan + 2;
			k = start_loop - 1;
		}

	}
}


void pyramidscanning(int resized,int looplength,int mode){
	//cout << "mode : " << mode << "\n";
	for (int k = 0; k < resized / 2; k++){
		
		double doubleindex1 = (double)k * (looplength/2) / (resized / 2);
		int index = (int)floor(doubleindex1);
		int hasil = doubleindex1 - index;
		mode == 0 ? resadv[k] = ((1.0 - hasil) * alldatavalue[index] + hasil * alldatavalue[index + 1]) : mode == 1 ? resadv2[k] = ((1.0 - hasil) * resadv[index] + hasil * resadv[index + 1]) : mode == 2 ? resadv3[k] = ((1.0 - hasil) * resadv2[index] + hasil * resadv2[index + 1]) : mode == 3 ? resadv4[k] = ((1.0 - hasil) * resadv3[index] + hasil * resadv3[index + 1]) : mode == 4 ? resadv5[k] = ((1.0 - hasil) * resadv4[index] + hasil * resadv4[index + 1]) : mode == 5 ? resadv6[k] = ((1.0 - hasil) * resadv5[index] + hasil * resadv5[index + 1]) : resadv7[k] = ((1.0 - hasil) * resadv6[index] + hasil * resadv6[index + 1]);
		
	}
	for (int k = 0; k < resized / 2; k++){

		double doubleindex1 = (double)k * (looplength / 2) / (resized / 2);
		int index = (int)floor(doubleindex1);
		int hasil = doubleindex1 - index;
		mode == 0 ? resadv_b[k] = ((1.0 - hasil) * alldatavalue_b[index] + hasil * alldatavalue_b[index + 1]) : mode == 1 ? resadv_b2[k] = ((1.0 - hasil) * resadv_b[index] + hasil * resadv_b[index + 1]) : mode == 2 ? resadv_b3[k] = ((1.0 - hasil) * resadv_b2[index] + hasil * resadv_b2[index + 1]) : mode == 3 ? resadv_b4[k] = ((1.0 - hasil) * resadv_b3[index] + hasil * resadv_b3[index + 1]) : mode == 4 ? resadv_b5[k] = ((1.0 - hasil) * resadv_b4[index] + hasil * resadv_b4[index + 1]) : mode == 5 ? resadv_b6[k] = ((1.0 - hasil) * resadv_b5[index] + hasil * resadv_b5[index + 1]) : resadv_b7[k] = ((1.0 - hasil) * resadv_b6[index] + hasil * resadv_b6[index + 1]);

	}
}


////svm
int print_null(const char *s, ...) { return 0; }

static int(*info)(const char *fmt, ...) = &printf;
int classified = 0;
struct svm_node *K;
int max_nr_attr = 64;

struct svm_model* model;
int predict_probability = 0;
unsigned int klasifikasi;
static char *lined = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
	int len;

	if (fgets(lined, max_line_len, input) == NULL)
		return NULL;

	while (strrchr(lined, '\n') == NULL)
	{
		max_line_len *= 2;
		lined = (char *)realloc(lined, max_line_len);
		len = (int)strlen(lined);
		if (fgets(lined + len, max_line_len - len, input) == NULL)
			break;
	}
	return lined;
}

void exit_input_error(int line_num)
{
	fprintf(stderr, "Wrong input format at line %d\n", line_num);
	exit(1);
}

void predict(FILE *input, FILE *output)
{
	int correct = 0;
	int total = 0;
	double error = 0;
	double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;

	int svm_type = svm_get_svm_type(model);
	int nr_class = svm_get_nr_class(model);
	double *prob_estimates = NULL;
	int j;

	if (predict_probability)
	{
		if (svm_type == NU_SVR || svm_type == EPSILON_SVR)
			info("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n", svm_get_svr_probability(model));
		else
		{
			int *labels = (int *)malloc(nr_class*sizeof(int));
			svm_get_labels(model, labels);
			prob_estimates = (double *)malloc(nr_class*sizeof(double));
			fprintf(output, "labels");
			for (j = 0; j<nr_class; j++)
				fprintf(output, " %d", labels[j]);
			fprintf(output, "\n");
			free(labels);
		}
	}

	max_line_len = 1024;
	lined = (char *)malloc(max_line_len*sizeof(char));
	while (readline(input) != NULL)
	{
		int i = 0;
		double target_label, predict_label;
		char *idx, *val, *label, *endptr;
		int inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0

		label = strtok(lined, " \t\n");
		if (label == NULL) // empty line
			exit_input_error(total + 1);

		target_label = strtod(label, &endptr);
		if (endptr == label || *endptr != '\0')
			exit_input_error(total + 1);

		while (1)
		{
			if (i >= max_nr_attr - 1)	// need one more for index = -1
			{
				max_nr_attr *= 2;
				K = (struct svm_node *) realloc(K, max_nr_attr*sizeof(struct svm_node));
			}

			idx = strtok(NULL, ":");
			val = strtok(NULL, " \t");

			if (val == NULL)
				break;
			errno = 0;
			K[i].index = (int)strtol(idx, &endptr, 10);
			if (endptr == idx || errno != 0 || *endptr != '\0' || K[i].index <= inst_max_index)
				exit_input_error(total + 1);
			else
				inst_max_index = K[i].index;

			errno = 0;
			K[i].value = strtod(val, &endptr);
			if (endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(total + 1);

			++i;
		}
		K[i].index = -1;

		if (predict_probability && (svm_type == C_SVC || svm_type == NU_SVC))
		{
			predict_label = svm_predict_probability(model, K, prob_estimates);
			fprintf(output, "%g", predict_label);
			for (j = 0; j<nr_class; j++)
				fprintf(output, " %g", prob_estimates[j]);
			fprintf(output, "\n");
		}
		else
		{
			predict_label = svm_predict(model, K);
			klasifikasi = predict_label;
			fprintf(output, "%.17g\n", predict_label);
		}

		if (predict_label == target_label)
			++correct;
		error += (predict_label - target_label)*(predict_label - target_label);
		sump += predict_label;
		sumt += target_label;
		sumpp += predict_label*predict_label;
		sumtt += target_label*target_label;
		sumpt += predict_label*target_label;
		++total;
	}
	if (svm_type == NU_SVR || svm_type == EPSILON_SVR)
	{
		info("Mean squared error = %g (regression)\n", error / total);
		info("Squared correlation coefficient = %g (regression)\n",
			((total*sumpt - sump*sumt)*(total*sumpt - sump*sumt)) /
			((total*sumpp - sump*sump)*(total*sumtt - sumt*sumt))
			);
	}
	else
		info("Accuracy = %g%% (%d/%d) (classification)\n",
		(double)correct / total * 100, correct, total);
	classified = correct;
	if (predict_probability)
		free(prob_estimates);
}

void exit_with_help()
{
	printf(
		"Usage: svm-predict [options] test_file model_file output_file\n"
		"options:\n"
		"-b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0); for one-class SVM only 0 is supported\n"
		"-q : quiet mode (no outputs)\n"
		);
	exit(1);
}
////


int main()
{

	/*DrawScatterPlot(imageRes, 1366, 768, &x_1, &y_1);
	vector<double> *pngData = ConvertToPNG(imageRes->image);
	WriteToFile(pngData, "grafiktess.png");
	DeleteImage(imageRes->image);
	cout << "finish";
	*/
	/*int a = 1998;
		int val = a - 1702;
		int value = val * 1.2 *1.2*1.2;
		cout << value;
	getchar();*/
	
	//close previous window and make connection to vrep
	//simxFinish(-1);//! Close any previously unfinished business
	clientID = simxStart((simxChar*)"127.0.0.1", 19000, true, true, 5000, 5);  //koneksi ke vrep
	simxStartSimulation(clientID, simx_opmode_oneshot);
	
	
	if (clientID != -1)// connection success
	{
		
		cout << " Connection status to VREP: SUCCESS" << endl;

		simxGetObjectHandle(clientID, "Pioneer_p3dx", &robothandle, simx_opmode_oneshot_wait);
		////int n = 1000-684; int r = 1060-684;
		//simxInt syncho = simxSynchronous(clientID, 1);
		//int start = simxStartSimulation(clientID, simx_opmode_oneshot_wait);

		//object handle
		simxGetObjectHandle(clientID, "Pioneer_p3dx_leftMotor", &roda1, simx_opmode_oneshot_wait);
		simxGetObjectHandle(clientID, "Pioneer_p3dx_rightMotor", &roda2, simx_opmode_oneshot_wait);
		/*simxSetJointTargetVelocity(clientID, roda1, 2, simx_opmode_oneshot_wait);
		simxSetJointTargetVelocity(clientID, roda2, 2, simx_opmode_oneshot_wait);
*/
		//get robot position
		simxFloat* _pos = new float[3];
		simxGetObjectPosition(clientID, robothandle, -1, _pos, simx_opmode_streaming);

		//get lidar data
		simxGetStringSignal(clientID, "LidarDepan", &laserdepan, &datasize, simx_opmode_streaming);
		simxGetStringSignal(clientID, "LidarBelakang", &laserbelakang, &datasize_b, simx_opmode_streaming);

		
		float dx; float dy;

		while (clientID != -1){
			clock_t tStart = clock();
			imgSave = Mat(640, 480, CV_8UC3, Scalar(0));
			//draw image size + rectangle background
			rectangle(imgSave, Point(0, 0), Point(480, 640), Scalar(255, 255, 255), -1);

			//draw robot position
			cv::circle(imgSave, Point(_pos[0]+80, 300 - _pos[1]), 8, Scalar(255, 0, 0), -1);

			//processing lidardepan data
			if (simxGetStringSignal(clientID, "LidarDepan", &laserdepan, &datasize, simx_opmode_streaming) == simx_error_noerror){

				container = datasize / (4 * 3);

				for (int i = 0; i < container; i++){
					lrf[0].x = ((simxFloat*)(laserdepan + 4 * 3 * i))[0];
					lrf[0].y = ((simxFloat*)(laserdepan + 4 * 3 * i))[1];
					

					//memasukkan semua data lidardepan ke array alldatavalue
					alldatavalue[i] = lrf[0].x <= 0 ? 5.5 : lrf[0].x;

					//yhasil.push_back(alldatavalue[i]);
					//rumus xlaser = r . cos tetha
					//rumus ylaser = r . sin tetha

					//data laser dikali 100 untuk memperbesar gambar
					//xlaser = (lrf[0].x * 100) * cos(((tetha*(i + 1)) - 49) * PI / 180);
					//ylaser = (lrf[0].y * 100) * sin(((tetha*(i + 1)) + 61) *  PI / 180);

					//xlaser = 200 + xlaser;
					//ylaser = 300 - ylaser;

					//rumus jarak real akar kuadrat dari (xlaser - xrobot)2 + (ylaser-yrobot)2
					//float calculate = ((xlaser - _pos[0]) * (xlaser - _pos[0])) + ((ylaser - _pos[1]) * (ylaser - _pos[1]));
					//float realdist = sqrt(calculate);

					//if (xlaser != 200 || ylaser != 300){
					//cv::circle(imgSave, Point(xlaser, ylaser), 2, Scalar(0, 0, 255), -1);
					//}

				}
				/*
				Mat rotated;
				rotated = rotate(imgSave, 180);*/

				//processing lidarbelakang data
				if (simxGetStringSignal(clientID, "LidarBelakang", &laserbelakang, &datasize_b, simx_opmode_streaming) == simx_error_noerror){
					container_b = datasize_b / (4 * 3);
					for (int i = 0; i < container_b; i++){
						lrf[1].x_b = ((simxFloat*)(laserbelakang + 4 * 3 * i))[0];
						lrf[1].y_b = ((simxFloat*)(laserbelakang + 4 * 3 * i))[1];

						//memasukkan semua data lidarbelakang ke array alldatavalue_b
						alldatavalue_b[i] = lrf[1].x_b <= 0 ? 5.5 : lrf[1].x_b;
						yll.push_back(alldatavalue_b[i]);
						xll.push_back(i);

						xlaser_b = (lrf[1].x_b * 100) * cos(((tetha*(i + 1)) - 30) *  PI / 180);
						ylaser_b = (lrf[1].y_b * 100) * sin(((tetha*(i + 1)) + 32) *  PI / 180);

						xlaser_b = 80 + xlaser_b;
						ylaser_b = 300 - ylaser_b;

						if (xlaser_b != 80 || ylaser_b != 300){

							cv::circle(imgSave, Point(xlaser_b, ylaser_b), 2, Scalar(0, 0, 255), -1);
						}
					}
				}

				
				//sprintf(file_name, "hasil%d.jpg",10);
				

				//Mat gethasil = imread("E://TA and OPENCV//vrepapi//vrepapi//Capture.JPG");
				//imshow("hasil", gethasil);

				/*DrawScatterPlot(imageori,1366,768,&xll,&yll);
				vector<double> *pngDataMm = ConvertToPNG(imageori->image);
				WriteToFile(pngDataMm, "grafikscanLRF.png");
				DeleteImage(imageori->image);
				cout << "finish";
				getchar();*/

				//DrawScatterPlot(imageRes, 1366, 768, &x, &y);
				//vector<double> *pngData = ConvertToPNG(imageRes->image);
				//WriteToFile(pngData, "tes2res.png");
				//DeleteImage(imageRes->image);
				//Mat result;
				//result = rotate(rotated, - 90);
				//Mat cropped;
				//cv::Rect crop_region(0, 80, 480, 480);
				//cropped = result(crop_region);

				//menjumlah banyak data lidardepan dan belakang
				int datalength = container + container_b;
				//new method
				int looplength = datalength;

				slidingwindow(looplength, 676, 0);//slidingwindow pertama 675

				int resized = looplength / 1.2;

				pyramidscanning(resized, looplength, 0);//pyramid scan pertama 1236
				slidingwindow(resized, 562, 1);//slidingwindow kedua

				looplength = resized;
				resized = looplength / 1.2;

				pyramidscanning(resized, looplength, 1);//pyramid scan kedua 1702
				slidingwindow(resized, 467, 2);//slidingwindow ketiga

				looplength = resized;
				resized = looplength / 1.2;

				pyramidscanning(resized, looplength, 2);//pyramid scan ketiga 2088
				slidingwindow(resized, 387, 3);//slidingwindow keempat

				looplength = resized;
				resized = looplength / 1.2;

				pyramidscanning(resized, looplength, 3);//pyramid scan keempat 2408
				slidingwindow(resized, 321, 4);//slidingwindow kelima

				looplength = resized;
				resized = looplength / 1.2;

				pyramidscanning(resized, looplength, 4);//pyramid scan kelima 2673
				slidingwindow(resized, 266, 5);//slidingwindow keenam

				looplength = resized;
				resized = looplength / 1.2;

				pyramidscanning(resized, looplength, 5);//pyramid scan keenam 2892
				slidingwindow(resized, 220, 6);//slidingwindow ketujuh

				looplength = resized;
				resized = looplength / 1.2;

				pyramidscanning(resized, looplength, 6);//pyramid scan ketujuh 3073
				slidingwindow(resized, 182, 7);//slidingwindow kedelapan

				//MyFile.close();
				Mygraph.close();
				

				FILE *input, *output;
				input = fopen("Datagraph.txt", "r"); //data input
				if (input == NULL)
				{
					cout << "1" << endl;
					exit(1);
				}
				output = fopen("hasil.txt", "w");
				if (output == NULL)
				{
					cout << "2" << endl;
					exit(2);
				}
				if ((model = svm_load_model("datasetsnew.txt.model")) == 0)
				{
					cout << "3" << endl;
					exit(3);
				}
				K = (struct svm_node *) malloc(max_nr_attr*sizeof(struct svm_node));
				if (predict_probability)
				{
					if (svm_check_probability_model(model) == 0)
					{
						fprintf(stderr, "Model does not support probabiliy estimates\n");
						exit(4);
					}
				}
				else
				{
					if (svm_check_probability_model(model) != 0)
						info("Model supports probability estimates, but disabled in prediction.\n");
				}

				predict(input, output);
				//cout << klasifikasi << endl;
				//cout << "\n";
				//svm_free_and_destroy_model(&model);
				free(K);
				free(lined);
				fclose(input);
				fclose(output);
				int run=0;
				int indeks = -1;
				int arr[150];
				string myText;
				ifstream MyReadFile("hasil.txt");
				while (getline(MyReadFile, myText)) {
					// Output the text from the file
					run = run + 1;
					if (myText == "1"){
						indeks = indeks + 1;
						if (run <= 675){
							arr[indeks] = run;
						}
						else if (run > 675 && run <= 1236){
							int value = run - 675;
							int lvalue = value * 1.2;
							arr[indeks] = lvalue;
						}
						else if (run > 1236 && run <= 1702){
							int value = run - 1236;
							int lvalue = value * 1.2 * 1.2;
							arr[indeks] = lvalue;
						}
						else if (run > 1702 && run <= 2088){
							int value = run - 1702;
							int lvalue = value * 1.2 * 1.2 * 1.2;
							arr[indeks] = lvalue;
						}
						else if (run > 2088 && run <= 2408){
							int value = run - 2088;
							int lvalue = value * 1.2 * 1.2 * 1.2 * 1.2;
							arr[indeks] = lvalue;
						}
						else if (run > 2408 && run <= 2673){
							int value = run - 2408;
							int lvalue = value * 1.2 * 1.2 * 1.2 * 1.2 * 1.2;
							arr[indeks] = lvalue;
						}
						else if (run > 2673 && run <= 2892){
							int value = run - 2673;
							int lvalue = value * 1.2 * 1.2 * 1.2 * 1.2 * 1.2 * 1.2;
							arr[indeks] = lvalue;
						}
						else if (run > 2892 && run <= 3073){
							int value = run - 2892;
							int lvalue = value * 1.2 * 1.2 * 1.2 * 1.2 * 1.2 * 1.2 *1.2;
							arr[indeks] = lvalue;
						}
						//cout << "value: " <<myText << "baris ke : " << run << "\n";
						
					}
				}
				MyReadFile.close();
				int arrcek[50];
				int index = -1;
				for (int k = 0; k < classified; k++){
					//cout << arr[k] << " ";
					if (count(std::begin(arr), std::end(arr), arr[k]) > 1 && count(std::begin(arrcek), std::end(arrcek), arr[k]) < 1){
						index = index + 1;
						arrcek[index] = arr[k];
					};
				}
				
				//draw result detected human foot
				int koor[5];
				int detectedman = 0;
				/*for (int i = 0; i < index; i++){
					for (int k = i + 1;k < index; k++){
						if (arrcek[i] > arrcek[k]){
							int tmp = arrcek[k];
							arrcek[k] = arrcek[i];
							arrcek[i] = tmp;
						}
					}
				}*/

				int param = 0;
				int counting = 0;
				int compare = 0;
				int iterateat = 0;
				for (int i = 0; i < index-1; i++){
					if (arrcek[i + 1]-arrcek[compare]  < 10){
						//cout << arrcek[i] << "\n"<<"\n";
						counting = counting + arrcek[i];
						//cout << counting << "\n";
						param = param + 1;
						//cout << param << "\n";
					}
					else if (arrcek[i + 1] - arrcek[compare] >= 10){
						/*cout << arrcek[compare] << "\n";*/
						//cout << "\n";/*
						//cout << "hm";
						//cout << arrcek[i+1] << "\n";
						if (param != 0){
							int data = counting / param;
							koor[0] = data;
							//cout << koor[0];
							param = 0;
							counting = 0;
						}
						compare = compare + i;
						detectedman = detectedman + 1;
						iterateat = iterateat + 1;
						
						//compare = compare + i;
					}
				}
				if (param != 0){
					detectedman = detectedman + 1;
					//cout <<"\n"<< counting << "\n";
					int data = counting / param;
					koor[1] = data;
					for (int i = 0; i < 2; i++){
						cout << koor[i]<<"\n";
						//cout << "oi";
					}
				}
					//for (int k = 0; k < detectedman; k++){
						int iterasi = ((arrcek[0] - 1) * 2) - 684;

						while (1){
							if (lrf[1].x_b = ((simxFloat*)(laserbelakang + 4 * 3 * iterasi))[0] == 0){
								iterasi = iterasi + 1;
							}
							else{
								break;
							}
						}

						lrf[1].x_b = ((simxFloat*)(laserbelakang + 4 * 3 * iterasi))[0];
						lrf[1].y_b = ((simxFloat*)(laserbelakang + 4 * 3 * iterasi))[1];

						/*xlaser_b = (lrf[1].x_b * 100) * cos(((tetha*(iterasi + 1)) - 30) *  PI / 180);
						ylaser_b = (lrf[1].y_b * 100) * sin(((tetha*(iterasi + 1)) + 32) *  PI / 180);
						xlaser_b = 80 + xlaser_b;
						ylaser_b = 300 - ylaser_b;
						rectangle(imgSave, Point(xlaser_b, ylaser_b), Point(xlaser_b - 20, ylaser_b - 40), Scalar(0, 255, 0));*/
					//}
				printf("Time taken in ms : %.2f ms\n FPS : %g fps\n", (double)(clock() - tStart) / CLOCKS_PER_SEC * 1000, 1000 / ((double)(clock() - tStart) / CLOCKS_PER_SEC * 1000));
				
				//menampilkan hasil
				//imshow("ialahasil", imgSave);
				//imwrite("foto.jpg", imgSave);
				//waitKey(1);
				//getchar();
			}

		}

	}
	else
	{
		cout << " Connection status to VREP: FAILED" << endl;
	}
	simxFinish(clientID);
	return clientID;
}
