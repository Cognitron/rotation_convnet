#pragma once

#include <vector>
#include <complex>
#include <math.h>
#include <opencv2/opencv.hpp>


#define M_PI       3.14159265358979323846

#define SAMPLING 32
#define FREQUENCY 5

std::complex<double> KF[SAMPLING][SAMPLING];
std::complex<double> iKF[SAMPLING][SAMPLING];

void Fourier(int n)
{
	static int r = 0;
	if (r) return;
	printf("\n Calculation fourier\n");
	for (int i = 0; i < SAMPLING; i++){
		for (int j = 0; j < SAMPLING; j++){
			double re = 0.0;
			double im = 2.0*M_PI*double(i*j) / double(n);
			KF[i][j] = exp(std::complex<double>(re, -im));
			iKF[i][j] = exp(std::complex<double>(re, im));
		}
	}
	r++;
}

void ifft(std::complex<double> val[SAMPLING], std::complex<double> params[SAMPLING])
{
	for (int i = 0; i < SAMPLING; i++){
		std::complex<double> sum(0.0, 0.0);
		for (int j = 0; j < SAMPLING; j++){
			sum += params[j] * iKF[i][j];
		}
		val[i] = sum;
	}
}

void fft(std::complex<double> val[SAMPLING], std::complex<double> params[SAMPLING])
{
	static int num = 0;
	if (!num) Fourier(SAMPLING);
	for (int i = 0; i < SAMPLING; i++){
		std::complex<double> sum(0.0, 0.0);
		for (int j = 0; j < SAMPLING; j++){
			sum += val[j].real() * KF[i][j];
		}
		params[i] = sum / double(SAMPLING);
	}
	num = 1;
}


class wave
{
public:
	double h;
	double ampl[FREQUENCY], phase[FREQUENCY];
	std::complex<double> val[SAMPLING], params[SAMPLING];

	friend const wave operator+(const wave& left, const wave& right)
	{
		wave res;
		for (int i = 0; i < SAMPLING; i++){
			res.val[i] = left.val[i] + right.val[i];
		}
		return res;
	}

	friend wave& operator+=(wave& left, const wave& right)
	{
		for (int i = 0; i < SAMPLING; i++){
			left.val[i] = left.val[i] + right.val[i];
		}
		left.update();
		return left;
	}

	friend const wave operator-(const wave& left, const wave& right)
	{
		wave res;
		for (int i = 0; i < SAMPLING; i++){
			res.val[i] = left.val[i] - right.val[i];
		}
		res.update();
		return res;
	}

	friend wave operator>>(const wave& left, const int phase)
	{
		wave res;
		for (int i = 0; i < SAMPLING; i++){
			res.val[i] = left.val[(i + phase) % SAMPLING];
		}
		res.update();
		return res;
	}

	friend wave operator<<(const wave& left, const int phase)
	{
		wave res;
		for (int i = 0; i < SAMPLING; i++){
			res.val[i] = left.val[(i - phase) % SAMPLING];
		}
		res.update();
		return res;
	}

	friend const double operator*(const wave& left, const wave& right)
	{
		double sum = 0.0;
		for (int i = 0; i < SAMPLING; i++){
			sum += left.val[i].real() * right.val[i].real();
		}
		return sum / double(SAMPLING);
	}

	friend wave operator*(const wave& left, double k)
	{
		wave res;
		for (int i = 0; i < SAMPLING; i++){
			res.val[i] = k*left.val[i];
		}
		res.update();
		return res;
	}

	friend wave operator*(double k, const wave& right)
	{
		wave res;
		for (int i = 0; i < SAMPLING; i++){
			res.val[i] = k*right.val[i];
		}
		res.update();
		return res;
	}

	wave()
	{
		h = 1.0 / double(SAMPLING);
		for (int i = 0; i < FREQUENCY; i++){
			ampl[i] = 0.0;
			phase[i] = 0.0;
		}
		for (int i = 0; i < SAMPLING; i++){
			val[i] = std::complex<double>(0.0, 0.0);
			params[i] = std::complex<double>(0.0, 0.0);
		}
	};

	void update(bool f = false)
	{
		fft(val, params);
		int ind = 1;
		double mxm = abs(ampl[1]);

		for (int i = 0; i < FREQUENCY; i++){
			ampl[i] = std::abs(params[i]);
			phase[i] = std::arg(params[i]);
			if (mxm < abs(ampl[i]) && i){
				mxm = abs(ampl[i]);
				ind = i;
			}
		}
		for (int i = 1; i < SAMPLING; i++){
			if (i >= FREQUENCY)params[i] = std::complex<double>(0.0, 0.0);
			if (ind != i && 0){
				if (i<FREQUENCY){
					ampl[i] = 0.0;
					phase[i] = 0.0;
				}
				params[i] = std::complex<double>(0.0, 0.0);
			}
		}
		ifft(val, params);
	}

	cv::Mat draw(int sz = 256, int background = 0, char *fname = "wave")
	{
		cv::Mat img(sz, sz, CV_32FC3, cv::Scalar(background, background, background));

		cv::Scalar colors[6];
		colors[0] = cv::Scalar(255, 32, 32);
		colors[1] = cv::Scalar(0, 255, 0);
		colors[2] = cv::Scalar(0, 0, 255);
		colors[3] = cv::Scalar(255, 0, 255);
		colors[4] = cv::Scalar(127, 127, 255);
		colors[5] = cv::Scalar(185, 200, 60);

		for (int c = 0; c < FREQUENCY; c++){
			for (int i = 0; i <= sz; i++){
				float x0 = 2.0*M_PI*float(i) / float(sz);
				float x1 = 2.0*M_PI*float(i + 1) / float(sz);
				float y0 = ampl[c] * cos(c*x0 + phase[c]);
				float y1 = ampl[c] * cos(c*x1 + phase[c]);

				cv::Point start(1.0f*i, sz / 2 - (sz / 4)*y0);
				cv::Point end(1.0f*(i + 1), sz / 2 - (sz / 4)*y1);

				cv::line(img, start, end, colors[c]);
			}
		}


		for (int i = 0; i <= sz; i++){
			float y0 = 0.0f;
			float y1 = 0.0f;
			for (int c = 0; c < FREQUENCY; c++){
				float x0 = 2.0*M_PI*float(i) / float(sz);
				float x1 = 2.0*M_PI*float(i + 1) / float(sz);
				y0 += ampl[c] * cos(c*x0 + phase[c]);
				y1 += ampl[c] * cos(c*x1 + phase[c]);
			}

			cv::Point start(1.0f*i, sz / 2 - (sz / 4)*y0);
			cv::Point end(1.0f*(i + 1), sz / 2 - (sz / 4)*y1);

			cv::line(img, start, end, colors[5]);
		}

		cv::Point start(0, sz / 2);
		cv::Point end(sz - 1, sz / 2);
		cv::line(img, start, end, cv::Scalar(255, 255, 255));

		if (fname != "wave"){
			char filename[1024];
			sprintf(filename, "%s.jpg", fname);
			cv::imwrite(filename, img);
		}

		return img;
	}
};
