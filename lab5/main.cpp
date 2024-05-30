#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <set>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>
#include <chrono>

double I(std::vector<double> &image, int i, int j, int H, int W)
{
	if (i < 0 || i > H - 1 || j < 0 || j > W - 1)
		return 0;
	return image[(i * W + j) * 3] + image[(i * W + j) * 3 + 1] + image[(i * W + j) * 3 + 2];
}

void compute_energy_map(std::vector<double> &E, std::vector<double> &image, int H, int W)
{
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			E[i * W + j] = abs(I(image, i, j + 1, H, W) - I(image, i, j - 1, H, W)) + abs(I(image, i + 1, j, H, W) - I(image, i - 1, j, H, W));
		}
	}
}

void compute_cumulated_emap(std::vector<double> &C, std::vector<double> &E, int H, int W)
{
#pragma omp parallel for schedule(dynamic, 1)
	for (int j = 0; j < W; j++)
	{
		C[j] = E[j];
	}
	// #pragma omp parallel for schedule(dynamic, 1)
	for (int i = 1; i < H; i++)
	{
		C[i * W] = std::min(C[(i - 1) * W], C[(i - 1) * W + 1]) + E[i * W];
		for (int j = 1; j < W - 1; j++)
		{
			C[i * W + j] = std::min(std::min(C[(i - 1) * W + j - 1], C[(i - 1) * W + j]), C[(i - 1) * W + j + 1]) + E[i * W + j];
		}
		C[i * W + W - 1] = std::min(C[(i - 1) * W + W - 2], C[(i - 1) * W + W - 1]) + E[i * W + W - 1];
	}
}

void backtracking(std::vector<double> &C, std::set<int> &path, int H, int W)
{
	std::vector<double>::iterator min_id = std::min_element(C.end() - W, C.end());
	int cur_i, cur_j;
	cur_i = H - 1, cur_j = min_id - (C.end() - W);

	std::cout << cur_j << "\n";
	std::cout << C.end() - min_id << "\n";
	std::cout << *(min_id - 1) << " " << *(min_id) << " " << *(min_id + 1) << "\n";
	for (auto itr = C.end() - W; itr != C.end(); itr++)
	{
		std::cout << *(itr) << " ";
	}
	path.insert(cur_i * W + cur_j);
	while (cur_i > 0)
	{
		cur_i--;
		if (C[(cur_i)*W + cur_j - 1] < C[(cur_i)*W + cur_j] && C[(cur_i)*W + cur_j - 1] < C[(cur_i)*W + cur_j + 1])
		{
			if (cur_j > 0)
			{
				cur_j--;
			}
			path.insert(cur_i * W + cur_j);
		}
		else if (C[(cur_i)*W + cur_j + 1] < C[(cur_i)*W + cur_j - 1] && C[(cur_i)*W + cur_j + 1] < C[(cur_i)*W + cur_j])
		{
			if (cur_j < W - 1)
			{
				cur_j++;
			}
			path.insert(cur_i * W + cur_j);
		}
		else
		{
			path.insert(cur_i * W + cur_j);
		}
	}
}

void shift_image(std::vector<double> &image, std::vector<double> &new_image, std::set<int> path, int H, int W)
{
	int total_length = W * H * 3;
	std::vector<double>::iterator start_src, end_src, start_dest;
	start_src = image.begin();
	start_dest = new_image.begin();
	for (auto itr = path.begin(); itr != path.end(); itr++)
	{
		int x = *itr;
		end_src = image.begin() + x * 3;
		std::copy(start_src, end_src, start_dest);
		start_dest += (end_src - start_src);
		start_src = end_src + 3;
	}
}

void draw_cmap(std::vector<double> &CM, std::set<int> &path, int W, int H, unsigned int idx)
{
	std::vector<unsigned char> cmap(W * H, 0);
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			cmap[i * W + j] = std::min(CM[i * W + j] / 500., 255.);
			if (path.find(i * W + j) != path.end())
			{
				cmap[i * W + j] = 255;
			}
		}
	}

	std::string filename = "./cmap/cmap_redim" + std::to_string(idx) + ".png";

	stbi_write_png(&filename[0], W, H, 1, &cmap[0], 0);
}

int main()
{
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	int W, H, C;
	const int reduce_w = 900;

	// stbi_set_intertically_on_load(true);
	unsigned char *image = stbi_load("redim.jpg",
									 &W,
									 &H,
									 &C,
									 STBI_rgb);

	std::vector<double> image_double(W * H * 3);
	for (int i = 0; i < W * H * 3; i++)
		image_double[i] = image[i];

	std::vector<double> I(W * H), E(W * H), CM(W * H);

	std::cout << "Input Size " << W << "x" << H << std::endl;

	std::vector<double> new_image(image_double.size());
	std::set<int> path;
	for (int i = 0; i < reduce_w; i++)
	{
		path = {};
		compute_energy_map(E, image_double, H, W);
		compute_cumulated_emap(CM, E, H, W);
		backtracking(CM, path, H, W);
		shift_image(image_double, new_image, path, H, W);
		image_double = new_image;
		if (i % 50 == 0)
		{
			std::cout << "i count: " << i << "\n";
			draw_cmap(CM, path, W, H, i);
		}
		W--;
	}

	std::cout << "Output Size " << W << "x" << H << std::endl;
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	// std::cout << "Size " << path.size() << std::endl;

	// W++;
	// std::vector<unsigned char> emap(W * H, 0);
	// for (int i = 0; i < H; i++)
	// {
	// 	for (int j = 0; j < W; j++)
	// 	{
	// 		E[i * W + j] /= 5;
	// 		emap[i * W + j] = std::min(E[i * W + j], 255.);
	// 	}
	// }

	// stbi_write_png("emap_redim.png", W, H, 1, &emap[0], 0);

	// std::vector<unsigned char> cmap(W * H, 0);
	// for (int i = 0; i < H; i++)
	// {
	// 	for (int j = 0; j < W; j++)
	// 	{
	// 		cmap[i * W + j] = std::min(CM[i * W + j] / 200, 255.);
	// 		if (path.find(i * W + j) != path.end())
	// 		{
	// 			cmap[i * W + j] = 255;
	// 		}
	// 	}
	// }

	// stbi_write_png("cmap_redim.png", W, H, 1, &cmap[0], 0);
	// W--;

	std::vector<unsigned char> image_result(W * H * 3, 0);
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{

			image_result[(i * W + j) * 3 + 0] = image_double[(i * W + j) * 3 + 0];
			image_result[(i * W + j) * 3 + 1] = image_double[(i * W + j) * 3 + 1];
			image_result[(i * W + j) * 3 + 2] = image_double[(i * W + j) * 3 + 2];
		}
	}
	stbi_write_png("redim_result.png", W, H, 3, &image_result[0], 0);

	return 0;
}