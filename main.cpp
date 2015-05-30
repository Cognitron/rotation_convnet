#include <stdio.h>
#include <vector>
#include <complex>
#include <math.h>
#include <time.h>

#include "neuron.h"

#define M_PI       3.14159265358979323846



struct pixel
{
	float x, y, a;
	int n;

	std::vector< std::pair< float, std::pair<int, int> > > substrate;
	std::vector< std::vector< std::pair<float, float>  > > sub_pixel;

	~pixel(){
		for (int i = 0; i < sub_pixel.size(); i++){
			sub_pixel.at(i).clear();
		}
		sub_pixel.clear();
		substrate.clear();
	}

	pixel(float cx, float cy, int num){
		x = cx;
		y = cy;
		a = 0.0;
		n = num;

		float h = 1.0 / (2 * n);
		for (int i = -n; i <= n; i++){
			std::vector< std::pair<float, float> > v;
			for (int j = -n; j <= n; j++){
				std::pair<float, float> c;
				c.first = j*h + x;
				c.second = i*h + y;

				if (j == -n) c.first += 0.00001f;
				if (j == n) c.first -= 0.00001f;
				if (i == -n) c.second += 0.00001f;
				if (i == n) c.second -= 0.00001f;
				v.push_back(c);
			}
			sub_pixel.push_back(v);
		}
	}

	void Rotate(float rot){
		a = rot * (M_PI / 180);
		float cs = cos(a);
		float sn = sin(a);
		for (int i = 0; i < sub_pixel.size(); i++){
			for (int j = 0; j < sub_pixel.at(i).size(); j++){
				float xx = sub_pixel.at(i).at(j).first;
				float yy = sub_pixel.at(i).at(j).second;

				sub_pixel.at(i).at(j).first = xx*cs + yy*sn;
				sub_pixel.at(i).at(j).second = xx*(-sn) + yy*cs;
			}
		}
	}

	void Substrate(){
		float dis = 0.0f;
		if ((x - floor(x))<0.1 && (y - floor(y))<0.1)
		{
			dis = 0.5;
		}

		std::vector< std::pair< float, std::pair<int, int> > > v;
		for (int i = 0; i < sub_pixel.size(); i++){
			for (int j = 0; j < sub_pixel.at(i).size(); j++){
				std::pair< float, std::pair<int, int> > c;
				c.first = 1.0 / ((2 * n + 1)*(2 * n + 1));
				c.second = std::pair<int, int>(floor(sub_pixel.at(i).at(j).first + dis), floor(sub_pixel.at(i).at(j).second + dis));
				v.push_back(c);
			}
		}

		for (int i = 0; i < v.size(); i++){
			bool flag = true;
			float val = v.at(i).first;
			std::pair<int, int> c = v.at(i).second;

			for (int j = 0; j < substrate.size(); j++){
				if (substrate.at(j).second == c){
					substrate.at(j).first += val;
					flag = false;
					break;
				}
			}
			if (flag){
				substrate.push_back(std::pair< float, std::pair<int, int> >(val, c));
			}
		}
	}

	void Print(){
		for (int i = 0; i < sub_pixel.size(); i++){
			for (int j = 0; j < sub_pixel.at(i).size(); j++){
				if (n<5)printf(" % .5lf % .5lf |", sub_pixel.at(i).at(j).first, sub_pixel.at(i).at(j).second);
			}
			if (n<5)printf("\n");
		}
		printf("\n---------------------------\n");
		for (int i = 0; i < substrate.size(); i++){
			printf(" (%d;%d) = % .3lf \n", substrate.at(i).second.first, substrate.at(i).second.second, substrate.at(i).first);
		}
	}

	void Draw()
	{
		int hw = 200;
		unsigned char *image = new unsigned char[hw*hw];

		float dis = 0.0f;
		if ((x - floor(x))<0.1 && (y - floor(y))<0.1)
		{
			dis = 0.5;
		}
		int d = dis*(hw / 10);

		for (int i = 0; i < hw; i++){
			for (int j = 0; j < hw; j++){
				image[j + i*hw] = 0;
				if (i>hw / 2) image[j + i*hw] += 32;
				if (j>hw / 2) image[j + i*hw] += 32;

				if (((i + d) % (hw / 10)) == 0 || ((j + d) % (hw / 10)) == 0)	image[j + i*hw] = 128;
			}
		}

		for (int i = 0; i < sub_pixel.size(); i++){
			for (int j = 0; j < sub_pixel.at(i).size(); j++){
				float xx = sub_pixel.at(i).at(j).first;
				float yy = sub_pixel.at(i).at(j).second;

				int ix = hw / 2 + (hw / 10)*xx;
				int iy = hw / 2 + (hw / 10)*yy;
				image[ix + (hw - iy - 1)*hw] = 255;
			}
		}
		/*/
		for( int k = 0 ; k < substrate.size() ; k++ ){
		int ax = 5 + substrate.at(k).second.first;
		int ay = 5 + substrate.at(k).second.second;
		float val = substrate.at(k).first;
		for( int i = 0 ; i < hw/10 ; i++ ){
		for( int j = 0 ; j < hw/10 ; j++ ){
		int ix = ax*(hw/10) + i;
		int iy = ay*(hw/10) + j;
		image[ix + (hw-iy-1)*hw] = 255*val;
		}
		}
		}
		//*/
		FILE *f = fopen("image.raw", "wb");
		fwrite(image, sizeof(unsigned char), hw*hw, f);
		fclose(f);

		delete(image);
	}
};


template <int fx, int fy, int fz, int rot >
class plane
{
	wave output, energy, error;
public:
	neuron<fx, fy, fz, rot> nrn;
};

int main()
{
	cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE);

	//cv::Mat img = cv::imread("F:\\_EC1M_v2\\images_00\\image_00037.jpg");

	FILE *f = fopen("C:\\Users\\Opex\\OneDrive\\database\\MNIST\\train-images.idx3-ubyte", "rb");
	char mnist[1024];
	fread(mnist, 16, sizeof(unsigned char), f);

	for (int i = 0; i < 1341; i++){
		fread(mnist, 28 * 28, sizeof(unsigned char), f);
	}

	cv::Mat img(256, 256, CV_8UC1);
	memset(img.data, 0, 256 * 256);

	for (int y = 0; y < 28; y++){
		for (int x = 0; x < 28; x++){
			img.data[(x + 114) + (y + 114) * 256] = mnist[x + 28 * y];
		}
	}


	//cv::Point start( 0, 256/2 );
	//cv::Point end( 256-1, 256/2 );
	//cv::line( img, start, end, cv::Scalar( 255, 255, 255 ) );
	cv::imwrite("src.bmp", img);

	int rows = img.rows;
	int cols = img.cols;


	cv::Mat dst;


	for (int i = 0; i < 361; i++){
		cv::Mat M = cv::getRotationMatrix2D(cv::Point2f(cols / 2.0, rows / 2.0), i, 1);
		cv::warpAffine(img, dst, M, cv::Point(cols, rows), CV_INTER_AREA);
		cv::imshow("Display window", dst);
		cv::waitKey(1);
	}

	cv::Mat M = cv::getRotationMatrix2D(cv::Point2f(cols / 2.0, rows / 2.0), 10, 1);
	cv::warpAffine(img, dst, M, cv::Point(cols, rows), CV_INTER_AREA);
	cv::imwrite("dst.bmp", dst);

	neuron<7, 7, 1, 32> nrn;
	nrn.CreateMask();
	int timer = -clock();
	for (int i = 0; i < 1; i++){
		nrn.Forward();
		nrn.Backward();
		nrn.Draw(i);
	}

	timer += clock();

	nrn.output.draw();


	printf("\n TIMER: %d ms\n", timer);

	return 0;
}
