#pragma once

#include <vector>
#include <complex>
#include <math.h>
#include <opencv2/opencv.hpp>

#include "wave.h"

template <int fx, int fy, int fz, int rot >
class neuron
{
public:
	float alfa;
	float cx, cy;
	wave output, energy, error;
	wave weight[fz][fy][fx];
	wave deltaw[fz][fy][fx];
	wave input[fz][3 * fy][3 * fx];

	std::vector< std::pair< float, std::pair<int, int> > > mask[rot][fy][fx];


	neuron()
	{
		alfa = 0.0;
		cx = 0.5f*(fx - 1);
		cy = 0.5f*(fy - 1);

		for (int z = 0; z < fz; z++){
			for (int y = 0; y < 3 * fy; y++){
				for (int x = 0; x < 3 * fx; x++){
					for (int i = 0; i < 32; i++){
						input[z][y][x].val[i] = 0.0;

						if (x == int(1.5*fx))		input[z][y][x].val[i] = 1.0;
						if (y == int(1.5*fy))		input[z][y][x].val[i] = 1.0;

					}
					input[z][y][x].update();
				}
			}
		}

		double k = 50.0 / double(fx*fy*fz);
		for (int z = 0; z < fz; z++){
			for (int y = 0; y < fy; y++){
				for (int x = 0; x < fx; x++){
					double ampl = 1.0 - (rand() & 0xfff) / 4096.0;
					double phase = 1.0 - (rand() & 0xfff) / 4096.0;
					int period = rand() % 4;
					for (int i = 0; i < 32; i++){
						weight[z][y][x].val[i] = -0.10;

						if (y == int(fy / 2))	weight[z][y][x].val[i] = 1.0;

						weight[z][y][x].val[i] *= k;

						//weight[z][y][x].val[i] = ampl*cos( period*2.0*M_PI*( double(i)/double(32) + phase ) );

						deltaw[z][y][x].val[i] = 0.0;
					}
					weight[z][y][x].update();
					deltaw[z][y][x].update();
				}
			}
		}
	};


	void CreateMask()
	{
		float dx = cx - floor(cx);
		float dy = cy - floor(cy);

		for (int y = -fy / 2, i = 0; y <= fy / 2; y++, i++){
			for (int x = -fx / 2, j = 0; x <= fx / 2; x++, j++){
				printf("\n ( % .2lf, % .2lf )", x + dx, y + dy);
				for (int r = 0; r < rot; r++){

					float rt = 360.0*float(r) / float(rot);
					pixel px = pixel(x + dx, y + dy, 100);
					px.Rotate(rt);
					px.Substrate();

					mask[r][i][j] = px.substrate;
				}
			}
		}
	}

	void Forward()
	{
#define printf //
		printf("\n");
		for (int z = 0; z < fz; z++){
			for (int y = 0; y < 3 * fy; y++){
				for (int x = 0; x < 3 * fx; x++){
					printf(" % .1lf", input[z][y][x].val[0].real());
				}
				printf("\n");
			}
		}
		printf("\n");
		printf("\n");
		wave sum;
		for (int z = 0; z < fz; z++){
			for (int y = 0; y < fy; y++){
				for (int x = 0; x < fx; x++){

					for (int r = 0; r < rot; r++){
						for (int i = 0; i < mask[r][y][x].size(); i++){
							float kf = mask[r][y][x].at(i).first;
							int xx = 1.5*fx + mask[r][y][x].at(i).second.first;
							int yy = 1.5*fy + mask[r][y][x].at(i).second.second;
							sum.val[r] += kf * (weight[z][y][x] << r) * input[z][yy][xx];
						}
					}
					printf(" % .1lf", weight[z][y][x].val[0].real());

				}
				printf("\n");
			}
		}

		for (int r = 0; r < rot; r++){
			printf("\nrot % .5lf", sum.val[r].real());
			//sum.val[r].real( std::max( sum.val[r].real(), 0.0 ) );
		}

		sum.update(true);
		printf("\n");
		for (int i = 0; i < 5; i++){
			printf("\n params: %.5lf | %.5lf", sum.ampl[i], sum.phase[i]);
		}
		printf("\n");
		for (int r = 0; r < rot; r++){
			printf("\nrot<< % .5lf", sum.val[r].real());
		}

		output = sum;
#undef printf
	}
	//*/
	void Backward()
	{
		for (int z = 0; z < fz; z++){
			for (int y = 0; y < fy; y++){
				for (int x = 0; x < fx; x++){
					//alfa * error.val[r].real() * 
					for (int r = 0; r < rot; r++){
						for (int i = 0; i < mask[r][y][x].size(); i++){
							float kf = mask[r][y][x].at(i).first;
							int xx = 1.5*fx + mask[r][y][x].at(i).second.first;
							int yy = 1.5*fy + mask[r][y][x].at(i).second.second;
							//alfa * error.val[r].real() * input[z][yy][xx];
							deltaw[z][y][x] += alfa * error.val[r].real() * ((input[z][yy][xx]) << r);
						}
					}
					for (int i = 0; i < SAMPLING; i++){
						//deltaw[z][y][x].val[i] = sin(2.0*M_PI*float(i)/float(SAMPLING));
					}
					deltaw[z][y][x].update();

					weight[z][y][x] = deltaw[z][y][x];
				}
			}
		}
	}
	//*/
	void Print()
	{
		for (int i = fy - 1; i < fy; i++){
			for (int j = fx - 1; j < fx; j++){
				for (int r = 0; r < rot; r++){
					for (int m = 0; m < mask[r][i][j].size(); m++){
						printf("\n ( %d , %d ) = % .5lf", mask[r][i][j].at(m).second.first, mask[r][i][j].at(m).second.second, mask[r][i][j].at(m).first);
					}
					printf("\n------------------------------------------\n");
					getchar();
				}
			}
		}
	}

	void Draw(int num)
	{
		int sz = 256;
		cv::Mat image(fx*sz, fy*sz, CV_32FC3);

		for (int y = 0; y < fy; y++){
			for (int x = 0; x < fx; x++){
				//cv::Mat src = output.draw(sz);
				cv::Mat src = weight[0][y][x].draw(sz, ((x % 2) + (y % 2) + 1) * 16);

				src.copyTo(image(cv::Range(x*sz, (x + 1)*sz), cv::Range(y*sz, (y + 1)*sz)));
			}
		}

		for (int y = 0; y < fy; y++){
			for (int x = 0; x < fx; x++){
				int thickness = 3;
				cv::Point start0(x*sz, y*sz);
				cv::Point end0(x*sz, fy*sz - 1);
				cv::line(image, start0, end0, cv::Scalar(127, 127, 127), thickness);

				cv::Point start1(x*sz, y*sz);
				cv::Point end1(fx*sz - 1, y*sz);
				cv::line(image, start1, end1, cv::Scalar(127, 127, 127), thickness);
			}
		}
		char filename[1024];
		sprintf(filename, "images\\neuron_%.3d.jpg", num);
		cv::imwrite(filename, image);
	}

};
